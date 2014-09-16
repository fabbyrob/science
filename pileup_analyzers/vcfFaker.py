import sys
import getopt

_L = 12
_N = 13
afs = [1]+[0]*2*_N

def __main__():
    #check aruguments
    processArgs(1)
    
    if len(afs) != 2*_N+1:
        sys.stderr.write("AFS length does not match sample size. Need %s AFS categories, had %s. Exiting.\n" % (2*_N+1, len(afs)))
        sys.exit(0)
    
    if sum(afs) != 1:
        sys.stderr.write("AFS does not sum to 1, putting excess in invariant category.\n")
        afs[0] += 1-sum(afs)
        sys.stderr.write("-a %s\n" % ",".join(map(str, afs)))
    
    afs_counts = []
    
    for a in afs:
        afs_counts.append(int(round(a*_L)))
    
    curr = 0
    print("##fileformat=VCFv4.0")
    header = "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=PL,Number=3,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">\n\
##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n\
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n\
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">\n\
##INFO=<ID=DS,Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">\n\
##INFO=<ID=Dels,Number=1,Type=Float,Description=\"Fraction of Reads Containing Spanning Deletions\">\n\
##INFO=<ID=HRun,Number=1,Type=Integer,Description=\"Largest Contiguous Homopolymer Run of Variant Allele In Either Direction\">\n\
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=\"Consistency of the site with at most two segregating haplotypes\">\n\
##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n\
##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">\n\
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">\n\
##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">\n\
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">\n\
##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">"
    
    print(header)
    header2 = "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  "
    for i in range(_N):
        header2 += "sample_"+str(i)+" "
    print(header2)
    
    for i, a in enumerate(afs_counts):
        for l in range(a):
            curr += 1
            #build int individual info
            inds = ""
            alt_count = i
            for n in range(_N):
                ind = ""
                if alt_count:
                    ind += "1/"
                    alt_count -= 1
                else:
                    ind += "0/"
                
                if alt_count:
                    ind += "1"
                    alt_count -= 1 
                else:
                    ind += "0"
                
                dp = "60"
                ml = "50"
                if ind == "0/0":
                    pl = [0, 80, 100]
                elif ind == "1/0":
                    pl = [80, 0, 100]
                else:
                    pl = [80, 100, 0]
                pls = ",".join(map(str, pl))
                
                ind += ":%s:%s:%s\t" % (dp, ml, pls)
                inds += ind
                
            #build the site info
            if i:
                alt = "A"
            else:
                alt = "."
            qual = 100
            af = float(i)/(2*_N)
            dels = 0.00
            
            #print the line
            print("scaffold_1 %s  .  C  %s  %s   .  AC=0;AF=%.2f;AN=24;Dels=%s;MQ=19.88;MQ0=7     GT:DP:GQ:PL     %s" % (str(curr), alt, qual, af, dels, inds))
            

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"N:L:a:h")
    except getopt.GetoptError:
        details()
    
    for opt, arg in opts:
        if opt == "-h":
            details()
            sys.exit(0)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-L":
            global _L
            _L = int(arg)
        elif opt == "-a":
            global afs
            afs = arg.split(",")
            afs = map(float, afs)
        else:
            sys.stderr.write ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("N: %s -L %s -a %s\n" % (_N, _L, ",".join(map(str, afs))))
   
use = "python "+__file__.split("/")[-1]+" [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print("Prints a fake VCF with a given AFS.")
    print (use)
    print("Option - type - Default - Description")
    print("h - NONE - NONE - help, prints out details of all options.")
    print("N - int - %s - number of diploid individuals to produce data for." % _N)
    print("L - int - %s - number of lines of a VCF to make." % _L)
    print("a - comma separated floats - %s - the AFS that this VCF should have when complete. The AFS needs to sum to 1. The AFS must have 2*N+1 entries." % (",".join(map(str, afs))))

    
if __name__ == "__main__":   
    __main__()