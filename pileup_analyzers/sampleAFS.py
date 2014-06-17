'''
Takes in a vcfSummary (see vcfSummarizer.py) and outputs a folded AFS for a given site type
X alleles are sampled at each locus, any site with fewer than X alleles is discarded.

usage:

python sampleAFS.py mySummary.txt -X 10 -t 4fold

'''

import sys
import getopt
from random import sample
from math import ceil

types = {'0fold': 3, 'stop': 8, 'intergene': 0, 'downstream': 6, 'upstream': 5, 'exon': 2, 'intron': 1, 'istop': 7, '4fold': 4, 'unknown': 9}

_t = '4fold'
_X = 10

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(2)
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    length = int(ceil(_X/2))+1
    AFS = [0]*(length)    
    #read the infile file...
    for line in infile:
        line = line.rstrip()
        
        if line[0] == "#":
            continue
        
        sline = line.split()
        
        #CHROM    POS    REF    ALT    REF_NUMBER    ALT_NUMBER    TOTAL    SITE_TYPE
        if len(sline) != 8:#bad line
            sys.stderr.write("Bad line:\n\t"+line+"\n")
            continue
        
        scaf = sline[0]
        pos = int(sline[1])
        ref = int(sline[4])
        alt = int(sline[5])
        tot = int(sline[6])
        type = int(sline[7])
        
        if tot < _X:#not enough sites 
            sys.stderr.write("Not enough alleles at:\n\t"+line+"\n")
            continue
        
        af = sampleSites(ref, alt)
        AFS[af] += 1
        
    mystr = ""
    for af in AFS:
        mystr += str(af)+"\t"
        
    print(mystr)
    
def sampleSites(ref, alt):
    samps = [0]*ref+[1]*alt
    alts = sum(sample(samps, _X))
    return min(alts, _X-alts)

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"X:t:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-t":
            global _t 
            _t = arg
        elif opt == "-X":
            global _X 
            _X = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" INFILE [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("Takes in a vcfSummary (see vcfSummarizer.py) and outputs AFS for a given site type\n\nX alleles are sampled at each locus, any site with fewer than X alleles is discarded.\n\nusage:\n\npython sampleAFS.py mySummary.txt -X 10 -t 4fold")
    print()
    print("option - argument type - default - description")
    print("X - INT - "+str(_X)+" - the number of alleles to sample at each site")
    print("t - STR - "+str(_t)+" - the type of sites of interest must be one of: "+str(types.keys()))
    
if __name__ == "__main__":   
    __main__()