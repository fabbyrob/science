import sys
import re
import getopt
import cPickle as pickle
import vcf
from math import sqrt

_g = ""
_f = 500
_w = 100
_s = 50
_i = 0
_d = 20
_D = 60
_q = 90
_L = 40
_N = 24

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    infile = open(sys.argv[1], "r")
    outfile = open(sys.argv[2], "w")
    processArgs(3)
    
    
    outfile.write("position\tx\tn\tfolded\n")
    for line in infile:
        line = line.rstrip()
        sline = line.split()

        scaf = sline[0]
        base = int(sline[1])
        genos = sline[3:3+_N]
        alt = int(sline[-2])
        ref = int(sline[-1])

        safe = 0        
        for n in genos:
            if n in ['A','T','C','G']:
                safe += 1
        
        if safe == _N and alt != 0 and alt != _N:
            if alt > _N-alt:
                alt = _N-alt
            outfile.write(str(base)+"\t"+str(alt)+"\t"+str(_N)+"\t1\n")
    outfile.close()

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"f:w:s:i:d:D:q:L:g:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-f":
            global _f 
            _f = int(arg)
        elif opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-s":
            global _s 
            _s = int(arg)
        elif opt == "-i":
            global _i 
            _i = float(arg)
        elif opt == "-d":
            global _d 
            _d = int(arg)
        elif opt == "-D":
            global _D 
            _D = int(arg)
        elif opt == "-q":
            global _q 
            _q = int(arg)
        elif opt == "-L":
            global _L 
            _L = int(arg)
        elif opt == "-g":
            global _g 
            _g = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" pickled_annotation VCF [-g gene_list -f flanking_size -w window_size -s step_size -i indel_cutoff -d min_depth -D max_depth -q min_qual -L min_likelihood]"
def usage():
    print (use)
    sys.exit()
    
def details():
    descrip = "This program takes an annotation, and a VCF and analyzes diversity around all genes present in that annotation. \nIt prints out \
Watterson's theta, pi, and Tajima's D for sliding windows around and within genes. \nIf a list of genes is provided it does this analysis only for a subset of genes.\n\
The annotation should be pre-processed using pickleAnnotation.py.\n\
The program outputs data for each window and each gene formatted like:\n\n\
    gene_name,midpoint,pi,theta,TajD\n\n\
    Midpoint is between -1*flanking size and 0 for windows upstream of the gene, \n\
    0 to 100 for windows within the gene (representing a % of the length of the gene). So midpoint of 50 would be halfway across the gene.\n\
    Or between 101 and 100+flanking size for windows downstream of the gene.\n"
    print(descrip)
    print (use)
    print("\nOption details below are printed in the form:\n\t OPTION_NAME ARGUMENT_TYPE (DEFAULT OPTION) - DESCRIPTION\n")
    print("-g STR ("+str(_g)+") - Used to provide a list containing a subset of the genes to analyze.")
    print("-f INT ("+str(_f)+") - The number of bases on either side of a gene to consider. If two genes are closer together than this number, then the full distance is used.")
    print("-w INT ("+str(_w)+") - Size of the sliding window to analyze sequence in.")
    print("-s INT ("+str(_s)+") - Step size for the sliding window.")
    print("-i FLOAT ("+str(_i)+") - Maximum 'Dels' value a site can have and be included in analysis.")
    print("-d INT ("+str(_d)+") - Minimum depth to include a site in the analysis, this is the individual, not total site, depth.")
    print("-D INT ("+str(_D)+") - Maximum depth to include a site in the analysis, this is the individual, not total site, depth.")
    print("-q INT ("+str(_q)+") - Minimum quality to include a site in the analysis.")
    print("-L INT ("+str(_L)+") - Minimum likelihood of a particular individual's genotype call to include a site in the analysis.")

    
if __name__ == "__main__":   
    __main__()
