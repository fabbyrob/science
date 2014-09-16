import sys
import getopt
import random

_b = 200
_w = 20

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    infile = open(sys.argv[1],"r")
    
    #read sites
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    #read the infile file...
    size = 200
    my_min = 0
    my_max = 3000
    bins = []
    for i in range(int((my_max/size)+1)):
        bins.append([])
    start = None
    bin = None
    end = None
    window = _w
    snps = 0
    sites = 0
    
    for line in infile:
        line = line.rstrip()
        sline = line.split()
        if len(sline) != 4:
            continue
        
        scaf, pos, af, dist = sline
        pos = int(pos)
        af = int(af)
        dist = int(dist)
        
        if start == None:
            start = pos
            bin = min(len(bins)-1, int(start/size))
            if bin == len(bins)-1:
                end = start + window
            else:
                end = min((bin+1)*size, start+window)
            
        
        if pos > end:
            bins[bin].append((snps, sites))
            start = pos
            bin = min(len(bins)-1, int(start/size))
            if bin == len(bins)-1:
                end = start + window
            else:
                end = min((bin+1)*size, start+window)
                
            snps = 0
            sites = 0
            
        if af > 0:
            snps += 1
        sites += 1
    
    bins[bin].append((snps, sites))

    #calc the point estimate totals
    point_estimates = []
    for ls in bins:
        point_estimates.append(calcDiversity(ls))
    

    #do the bootstraps
    boots = []
    for i in range(len(bins)):
        boots.append([])
        
    for i in range(_b):
        this_boot = []
        j = 0
        for ls in bins:
            this_bin = []
            while len(this_bin) < len(ls):
                num = random.randint(0,len(ls)-1)
                this_bin.append(ls[num])
            boots[j].append(calcDiversity(this_bin))
            j += 1
        
    cis = calcCIs(boots)
    
    print("starting_dist\tpoint_estimate\tlowCI\thighCI")
    #sys.stderr.write(str(point_estimates)+"\n")
    #sys.stderr.write(str(cis)+"\n")
    for i in range(0, len(point_estimates)-1):
        print(str(my_min+size*i)+"\t"+str(point_estimates[i])+"\t"+str(cis[i][0])+"\t"+str(cis[i][1]))

def calcDiversity(ls):
    snps = 0
    tots = 0
    for snp, tot in ls:
        snps += snp
        tots += tot
    
    if tots == 0:
        return 0
    
    return float(snps)/tots

def calcCIs(boots):
    cis = []
    
    for ls in boots:
        ls.sort()
        minimum = ls[int(len(ls)*0.025)]
        maximum = ls[int(len(ls)*0.975-1)]
        cis.append((minimum, maximum))
    
    return cis
        

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:b:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-b":
            global _b 
            _b = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" infile [-w window size -b num_boots]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()