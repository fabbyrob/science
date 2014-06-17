import sys
import re
import getopt
from collections import defaultdict

_N = 13

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[2:],"N:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-N":
            global _N
            _N = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    vcf = open(sys.argv[1],"r")
    
    if (vcf == None):
        print("Bad vcf name: "+sys.argv[1])
        sys.exit()
             
    #reg ex
    vcf_pat = re.compile("\w+_(\d+)\s+(\d+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(.+)")
    
    depths = []
    for i in range(0, _N):
        depths.append(defaultdict(int))
    
    #read the vcf file...
    for line in vcf:
        line = line.rstrip()
        m = vcf_pat.match(line)
        
        if (m == None):# if we find a weird line...
            sys.stderr.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            scaf = int(m.group(1))
            base = int(m.group(2))
            alt = m.group(3)
            qual = m.group(4)
            if qual == ".":
                qual = 0
            else:
                qual = float(qual)
            depth = m.group(5).split(":")
            genos = m.group(6)
            
            if 'DP' in depth:
                depthIndex = depth.index('DP')
            else:
                continue
            
            i = 0
            for ind in genos.split():
                if ind == "./.":
                    continue
                ind = ind.split(":")
                iDepth = int(ind[depthIndex])
                #print(iDepth)
                depths[i][iDepth] += 1
                i += 1

    print("Depth\t#samps")
    i = 0
    for ind in depths:
        indStr = "Ind"+str(i)
        indDPs = "Depth\t"
        indCounts = "Num_Sites\t"
        for dp in ind.keys(): 
            indDPs += str(dp)+"\t"
            indCounts += str(ind[dp])+"\t"
        print(indStr)
        print(indDPs)
        print(indCounts)
        print()
        i += 1
   
use = "python "+__file__.split("/")[-1]+" VCF [-N <Samp size>]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()