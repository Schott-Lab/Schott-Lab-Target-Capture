from Bio import SeqIO
import sys
import getopt


def main(argv):
    printout = False
    try:
        opts, args = getopt.getopt(
            argv, "phr:i:o:", ["--print","ref=", "ifile=", "ofile="])
    except getopt.GetoptError:
        print ('CoverageMap.py -r <reference> -i <inputfile> -o <outputfile.csv>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-r", "--ref"):
            reference = arg
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-p", "--print"):
            printout = True
        
    referenceFile = SeqIO.to_dict(SeqIO.parse(open(reference, "r"), "fasta"))
    consensusFile = SeqIO.to_dict(SeqIO.parse(open(inputfile, "r"), "fasta"))

    outputlist = []
    for key in referenceFile.keys():
        thisSeq = consensusFile.get(key, None)
        if (len(key.split('_')) > 1):
            tup = (key.split('_')[1],)
        else:
            tup = (key,)
        if (not thisSeq):
            tup = tup + ("0",)
        else:
            percentCoverage = 100 * \
                float(len(str(thisSeq.seq).replace("n","").replace("N",""))) / \
                len(referenceFile[key].seq)
            tup = tup + (str(percentCoverage),)
        outputlist.append(tup)
    
    outputlist.sort(key=lambda k: k[0].lower())
    with open(outputfile, 'w') as output:
        output.write("Gene,Coverage\n")
        for i in outputlist:
            if (printout):
                print (i)
            output.write("%s,%s\r\n" % i)

    # for refSeq, conSeq in itertools.izip(SeqIO.parse(referenceFile, "fasta"), SeqIO.parse(consensusFile, "fasta")):
    #     print(refSeq.id)
    #     print(conSeq.id)
    # referenceFile.close()
    # consensusFile.close()

if __name__ == "__main__":
    main(sys.argv[1:])
