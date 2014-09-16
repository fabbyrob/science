'''
This program takes in a vcf summary file (see vcfSummarizer.py) and outputs
data the minor allele frequency and distance to nearest non-intergenic site 
(including UTRs, CNCs, and anything else annotated as not intergenic in the summary file).

example use:

python diversityNearCoding.py mySummary.txt 
#CHROM    POS    MAF    DIST
scaffold_1 100 0 100
scaffold_1 101 3 99
scaffold_1 102 0 98
scaffold_1 103 0 97

the output has 4 columns: scaffold, position, minor allele frequency, and distance to nearest coding region.
'''

import sys
import getopt
types = {'0fold': 3, 'stop': 8, 'intergene': 0, '3utr': 6, '5utr': 5, 'exon': 2, 'intron': 1, 'istop': 7, '4fold': 4, 'unknown': 9, 'cnc':10}
_t = ['intergene']
_c = False
     
def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()

    pCode = -1
    nCode = -1
    sites = []
    types_num = []
    for t in _t:
        types_num.append(types[t])
     
    #read the infile file...
    for line in infile:
        line = line.rstrip()
        if line[0] == "#":#header
            continue
        sline = line.split()
        
        #CHROM    POS    REF    ALT    REF_NUMBER    ALT_NUMBER    TOTAL    SITE_TYPE
        if len(sline) != 8:#bad line
            sys.stderr.write("Bad line:\n\t"+line+"\n")
            continue
        
        scaf = sline[0]
        pos = int(sline[1])
        maf = min(int(sline[4]), int(sline[5]))
        type = int(sline[7])
        
        if type in types_num:
            sites.append((scaf, pos, maf))
            continue
        
        nCode = pos
        if sites:
            processSites(sites, pCode, nCode)
            sites = []
        
        if not _c or (_c and type == types['cnc']):#if we only care about cncs then only make pCode the latest CNC site
            pCode = pos
        nCode = -1
        
    if sites:#get the last few sites
        processSites(sites, pCode, nCode)
        
def processSites(sites, pCode, nCode):
    for site in sites:
        if pCode == -1:
            dist = nCode - site[1]
        elif nCode == -1:
            dist = site[1] - pCode
        else:
            dist = min(site[1]-pCode, nCode-site[1])
            
        print(site[0]+"\t"+str(site[1])+"\t"+str(site[2])+"\t"+str(dist))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"t:c")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-t":
            global _t 
            _t = arg.split(",")
        elif opt == "-c":
            global _c
            _c = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" INFILE [OPTINS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print('This program takes in a vcf summary file (see vcfSummarizer.py) and outputs\ndata the minor allele frequency and distance to nearest non-intergenic site \n(including UTRs, CNCs, and anything else annotated as not intergenic in the summary file).\n\nexample use:\n\npython diversityNearCoding.py mySummary.txt \n#CHROM    POS    MAF    DIST\nscaffold_1 100 0 100\nscaffold_1 101 3 99\nscaffold_1 102 0 98\nscaffold_1 103 0 97\n\nthe output has 4 columns: scaffold, position, minor allele frequency, and distance to nearest coding region.')
    print()
    print("option - argument type - default - description")
    print("t - STR - "+str(_t)+" - the type of sites of interest must be one of (or a comma seperated list of): "+str(types.keys()))
    print("c - flag - "+str(_c)+" - tells the program to only calculate distance from CNCs, not genes")
    
if __name__ == "__main__":   
    __main__()