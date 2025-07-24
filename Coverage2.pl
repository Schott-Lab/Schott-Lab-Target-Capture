#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $usage = "$0 infile\n";
my $in = shift or die $usage;

my $seq_in = Bio::SeqIO->new(-file => $in, -format => 'fasta');

my ($n_count, $total, $seq_ct, $cov_count, $n_perc, $cov_perc, $cov_perc_t) = (0, 0, 0, 0, 0, 0, 0);

print "Seq_ID\tLength\tN_count\tcov_count\tN_perc\tCov_perc\n";

while(my $seq = $seq_in->next_seq) {
    if ($seq->length > 0) {
        $seq_ct++;
        $total += $seq->length;
        $n_count = ($seq->seq =~ tr/Nn//);
        $cov_count = sprintf("%.3f",$seq->length-$n_count);
        $n_perc  = sprintf("%.3f",$n_count/$seq->length * 100);
        $cov_perc = sprintf("%.3f",$cov_count/$seq->length * 100);
        $cov_perc_t += $cov_perc;
        print join("\t",($seq->id,$seq->length,$n_count,$cov_count,$n_perc,$cov_perc)),"\n";
    }
}

my $len_ave     = sprintf("%.0f",$total/$seq_ct);
my $cov_perc_ave  = sprintf("%.3f",$cov_perc_t/$seq_ct);

print "\nAverage Coverage for\t$in\t$cov_perc_ave\n";