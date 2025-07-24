#!/bin/bash

# Based on: Targeted capture of complete coding regions across divergent species
# Ryan K Schott, Bhawandeep Panesar, Daren C Card, Matthew Preston, Todd A Castoe, Belinda SW Chang

# Updates: Ryan K Schott and Taiye Estwick
# Dependancies: bwa-mem, samtools, bcftools, bedtools, mosdepth, perl, bio-perl, python 3, biopython
# Recommend using conda/mamba to manage packages
# Blast searches currently disabled. Add any dummy variable for the bast database.
# Usage example:
# bash /mnt/c/scripts/BT2_SeqCap_Assembly_v6.sh -s /mnt/c/scripts/ /mnt/c/target_capture/ref_genome.fas blastdb 12 /mnt/c/target_capture/trimmed/*.fq

clean_up() {
	#Deletes broken files if they were created and ends program
	printf "Creating $1" >&2
	rm -f $1
	shift
	for file in $@
	do
		printf ", $file" >&2
		rm -f $file
	done
	printf " failed\n" >&2
	exit 1
}

USAGE="Usage: bash BWA_SeqCap_Assembly_v6.sh [-s scriptdir] ref_genome.fas blastdb threads [fwd_read.fastq rvrse_read.fastq]..."
scriptdir="$(pwd)"
seqs=1

while getopts ':hs:' flag
do
	case $flag in
		h) # Help
			echo "$USAGE"
			exit 1
			;;
		s) # Change script directory
			scriptdir=$OPTARG
			;;
		\?)# Unexpected flag
			echo "Unexpected option -$OPTARG" >&2
			exit 1
			;;
		:) # Missing argument
			echo "Option -$OPTARG requires an argument" >&2
			exit 1
			;;
		*) # Unexpected error
		    echo "Unexpected error in getopts" >&2
			exit 1
			;;
	esac
done
shift $( expr $OPTIND - 1 )

if [ $# -le 1 ]
then
	echo $USAGE
	exit 1
fi

date +"Time started: %F %T"
echo "Using Reference Genome: $1"
ref="$1"
shift

echo "Using blast db $1"
db="$1"
shift

echo "Using $1 threads"
threads=$1
shift

if [ -f $ref.bwt ]
then
   echo "Reference index exists, skipping reference index build"
else
   echo "Making bwa index"
   bwa index $ref || clean_up $ref.bwt
fi

if [ -f $ref.fai ]
then
   echo "Samtools index exists, skipping"
else
   echo "Making Samtools Index"
   samtools faidx $ref || clean_up $ref.fai
fi

if [ $# -gt 0 ]
then
	tmp=$( dirname $ref | cut -d / -f 1 )
	dirs=("." "..")
	dirs+=( $(ls -F | grep ".*/" | sed 's|/||') )
	for d in ${dirs[@]}
	do
		if [ $tmp = $d ]
		then
			ref="../../$ref"
			break
		fi
	done

	tmp=$(basename $ref)
	refname=${tmp%.*}
	echo "Making ${refname}_BWA directory"
	if [ ! -d ${refname}_BWA ]
	then
		mkdir ${refname}_BWA
	fi
	cd ${refname}_BWA || { echo "Changing Directory failed"; exit 1; }
fi

while test ${#} -gt 0
do
	date +"Pair start: %T"
	echo "Making Directories"

	tmp="$(basename "$1")"
	echo "$tmp"
	seqname="${tmp%_R*}"
	echo "Using Sequence Name: $seqname"
	name=$(basename "${seqname%%_*}")
	echo "Using Name: $name"
	if [ ! -d "$seqname" ]
	then
		mkdir "$seqname"
	fi

	forward="$1"
	echo "Forward read: $forward"
	shift
	reverse="$1"
	echo "Reverse read: $reverse"
	shift
	echo "Changing to directory $seqname"
	cd "$seqname" || { echo "Changing Directory failed"; exit 1; }

	if [ -f $seqname.bam.bai ]
	then
	   echo "BWA SAMTOOLS output exists, skipping alignment"
	else
	   echo "Running BWA SAMTOOLS pipe: bwa mem -B 2 -M -t $threads $ref $forward $reverse"
	   bwa mem -B 2 -M -t $threads $ref "$forward" "$reverse" | samtools fixmate -O bam -@ $threads - - | samtools sort -@ $threads -o $seqname.bam - 
	fi

	if [ -f $seqname.bam.bai ]
	then
	    echo "Sorted BAM Index exists, skipping"
	else
	    echo "Indexing Sorted BAM file"
	    samtools index -@ $threads $seqname.bam || clean_up $seqname.bam.bai
	fi

	if [ -f stats_$seqname.txt ]
	then
	   echo "Stats outputs exist, skipping"
	else
	   echo "Generating Stats Files"
	   samtools flagstat -@ $threads $seqname.bam > flagstat_$seqname.txt || clean_up flagstat_$seqname.txt
	   samtools stats -@ $threads $seqname.bam > stats_$seqname.txt || clean_up stats_$seqname.txt
	   samtools idxstats $seqname.bam > idxstats_$seqname.txt || clean_up idxstats_$seqname.txt
	fi
	
	if [ -f depth_$seqname.txt ]
	then
	   echo "Stats outputs exist, skipping"
	else
	   echo "Generating Stats Files"
	   samtools depth -d 1000000 -aa $seqname.bam > depth_$seqname.txt || clean_up depth_$seqname.txt
	   mosdepth -t $threads $seqname $seqname.bam
	fi

	if [ -f mapped_$seqname.bam ]
	then
	    echo "Mapped output exists, skipping"
	else
	    echo "Filtering out only the mapped seqs"
	    samtools view -b -F 4 -@ $threads $seqname.bam > mapped_$seqname.bam || clean_up mapped_$seqname.bam
		echo "Removing full bam file"
		rm $seqname.bam
	fi

	if [ -f bcftools_consensus_masked_fixed_${seqname}.fas ]
	then
		echo "bcftools masked consensus exists, skipping"
	else
		echo "Running bcftools consensus pipeline with low-coverage masked reference"
		cat depth_$seqname.txt | awk '$3 < 5 {print}' | awk -v OFS='\t' '{print $1,int($2)-1,$2}' > low_coverage_$seqname.bed
		bedtools maskfasta -fi $ref -bed low_coverage_$seqname.bed -fo masked_ref_${seqname}.fas
		bcftools mpileup --threads $threads -Ou -f masked_ref_${seqname}.fas mapped_$seqname.bam | bcftools call --threads $threads -Ou -mv | bcftools norm --threads $threads -Oz -f masked_ref_${seqname}.fas > normalized_calls_masked_${seqname}.vcf.gz
		bcftools index --threads $threads normalized_calls_masked_${seqname}.vcf.gz
		bcftools consensus -f masked_ref_${seqname}.fas normalized_calls_masked_${seqname}.vcf.gz > bcftools_consensus_masked_${seqname}.fas
		echo "Fixing fasta headers"
		cat bcftools_consensus_masked_${seqname}.fas | sed "s/>.*_/>${seqname}_/" > bcftools_consensus_masked_fixed_${seqname}.fas
		rm bcftools_consensus_masked_${seqname}.fas
	fi

	if [ -f masked_coverage_$seqname.tsv ]
	then
		echo "Masked Coverage file exists, skipping"
	else
		echo "Calculating Masked Coverage"
		perl "$scriptdir"/Coverage2.pl bcftools_consensus_masked_fixed_$seqname.fas > masked_coverage_$seqname.tsv || clean_up masked_coverage_$seqname.tsv
	fi


	date +"Pair end: %T"
	echo "_________________________________________________"
	echo "FINISHED. Moving on to next pair"
	echo "_________________________________________________"
	cd ..  || { echo "Changing Directory failed"; exit 1; }
	seqs=$((seqs+1))
done


if [ $seqs -gt 0 ]
then
	echo "Creating Summary Files"
    #cd ${refname}_BWA
	python "$scriptdir"/Masked_Coverage.py */masked_coverage*.tsv
	python "$scriptdir"/Mosdepth.py */*mosdepth.summary.txt
	
fi

echo "All Done here!"
date +"Time ended: %F %T"
