#!/usr/bin/python
import sys
import re
import getopt
from Queue import PriorityQueue

_q = 40#min quality
_d = 20#min depth
_D = 60#max depth
_l = 40#ML cutoff
_G = False#flag for GATK files

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"q:d:D:L:G")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-q":
            global _q
            _q = int(arg)
        elif opt == "-d":
            global _d
            _d = int(arg)
        elif opt == "-D":
            global _D
            _D = int(arg)
        elif opt == "-L":
            global _l 
            _l = int(arg)
        elif opt == "-G":
            global _G 
            _G = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    vcf = open(sys.argv[1],"r")
    
    if (vcf == None):
        print("Bad vcf name: "+sys.argv[1])
        sys.exit()
        
    sites = open(sys.argv[2],"r")  
        
    if (sites == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()
        
    #reg ex
    if not _G:
        vcf_pat = re.compile("\w+_(\S+)\s(\d+)\s\S+\s\S+\s(\S+)\s(\S+)\s\S+.+DP=(\d+);\S+\s\S+\s(.+)")
    else:
        vcf_pat = re.compile("\w+_(\d+)\s+(\d+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(.+)")
    site_pat = re.compile("\w+_(\d+)\s+(\S+)")
         
    #initialize values
    next_site = getNewSite(sites, site_pat)
    
    #read the vcf file...
    for line in vcf:
        line = line.rstrip()
        m = vcf_pat.match(line)
        
        if (m == None):# if we find a weird line...
            sys.stderr.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            #grab the relevant data
            scaf = int(m.group(1))
            base = int(m.group(2))
            alt = m.group(3)
            if not _G:
                qual = float(m.group(4))
                depth = int(m.group(5))
                genos = m.group(6)
            else:
                qual = m.group(4)
                if qual == ".":
                    qual = 0
                else:
                    qual = float(qual)
                depth = m.group(5).split(":")
                genos = m.group(6)
            
            #if we're not at a site we care about...
            if(scaf != next_site[0] or base != next_site[1]):
                while (scaf > next_site[0] or (scaf == next_site[0] and base > next_site[1])):
                    sys.stderr.write("\t FAIL at : "+str(scaf) +", "+str(base)+"\n")
                    sys.stderr.write("site not found: "+str(next_site)+"\n")
                    next_site = getNewSite(sites, site_pat)
                    #sys.stderr.write("site:\t"+str(next_site)+"\n")
                    if (next_site == None):
                        break
                    if (scaf == next_site[0] and base == next_site[1]):
                        safe = callSite(scaf, base, qual, depth, genos)
                        if (safe):
                            print ("scaffold_"+str(scaf) + "\t" + str(base))
                        else:
                            sys.stderr.write("Site did not pass filters:\n\t"+genos+"\n")
                            
                        next_site = getNewSite(sites, site_pat)
                        if (next_site == None):
                            break
            else:
                safe = callSite(scaf, base, qual, depth, genos)
                            
                if (safe):
                    print ("scaffold_"+str(scaf) + "\t" + str(base))
                    
                #grab the next site
                next_site = getNewSite(sites, site_pat)
            if (next_site == None):
                break

def callSite(scaf, base, qual, depth, genos):
    if (qual < _q ):
        safe = 0
        sys.stderr.write("Site did not pass filters (qual):\n\t"+" ("+str(qual)+") "+genos+"\n")
        return safe
    if not _G and (depth > _D or depth < _d):
        safe = 0
        sys.stderr.write("Site did not pass filters (depth):\n\t"+" ("+str(depth)+") "+genos+"\n")
        return safe
    
    if _G:
        if 'DP' in depth:
            depthIndex = depth.index('DP')
        else:
            depthIndex = -1
    
    inds = genos.split()
    safe = 1
    #everybody homozygous reference
    if not _G and (inds[0].split(":") == [inds[0]]):
        safe = 1
    else:
        #other
        c = -1
        for i in inds:
            c += 1
            if depthIndex == -1:
                sys.stderr.write("Site did not pass filters (depthIndex - indel?):\n\t"+str(scaf)+"\t"+str(base)+"\n")
                safe = 0
                break
            
            data = i.split(":")
            if not _G:
                data = data[0].split(",")
            else:
                #in GATK files we check that EVERy individual meets depth requirements
                if data == ["./."]:
                    depth = 0
                else:
                    depth = int(data[depthIndex])
                data = data[-1].split(",")
                if depth < _d or depth > _D:
                    sys.stderr.write("Site did not pass filters (depth - at least one too low or high):\n\t"+" ("+"ind = "+str(c)+" depth = "+str(depth)+") "+str(scaf)+"\t"+str(base)+"\n")
                    safe = 0
                    break
                
            makeInts(data)
            
            (min, max, mid) = getMinMaxIndex(data)
            
            if (data[min] != 0):
                #site is ambiguous, no max likelihood
                safe = 0
                sys.stderr.write("Site did not pass filters (likelihood - no max):\n\t"+" ("+"ind = "+str(c)+") "+str(scaf)+"\t"+str(base)+"\n")
                break
            
            if (data[mid] < _l):
                #dont pass quality cut off missing site
                safe = 0
                sys.stderr.write("Site did not pass filters (likelihood - to low):\n\t"+" ("+"ind = "+str(c) + " likelihood = "+ str(data[mid])+") "+str(scaf)+"\t"+str(base)+"\t"+"\n")
                break
    return safe
        
def makeInts(ls):
    for i in range(0, len(ls)):
        ls[i] = int(ls[i])
   
def getMinMaxIndex(ls):
    min = 0
    minV = ls[0]
    max = 0
    maxV = ls[0]
    for i in range(0, len(ls)):
        if ls[i] < minV:
            min = i
            minV = ls[i]
        if ls[i] > maxV:
            max = i
            maxV = ls[i]
            
    return (min, max, list(set([1,2,0])-set([min, max]))[0])
   
#grabs a new site too look for from the sites input file
def getNewSite(sites, site_pat):
    next_site = None
    while(next_site == None):
        new_site = sites.readline()
        
        if (new_site == ""):
            return None
        
        m = site_pat.match(new_site)
        
        if (m == None):
            sys.stderr.write("Non-standard line (SITE FILE):\n\t"+new_site+"\n")
            continue
        next_site = (int(m.group(1)), int(m.group(2)))
        
    print(next_site)
    return next_site
   
def usage():
    print ("usage: filter.py <VCF> <SITE_LIST> [-q <QUALITY CUTOFF> -d <DEPTH MIN> -D <DEPTH MAX> -L <LIKELIHOOD CUTOFF> -G]]\n")
    sys.exit()
    
def details():
    print ("usage: filter.py <VCF> <SITE_LIST> [-q <QUALITY CUTOFF> -d <DEPTH MIN> -D <DEPTH MAX> -G\n")
    print ("<VCF> - a variant call format file")
    print ("<SITE_LIST> - a file containing a list of sites of interest. Contains 2 columns, Scaffold Name and Position, each site on its own line.")
    print ("q - the minimum quality to include a site in the analysis (40)")
    print ("d - the minimum depth to include a site in the analysis (20)")
    print ("D - the maximum depth to include a site in the analysis (60)")
    print ("L - the cut off for calling a site ambiguous based on heterozygote ML scores (40)")
    print ("G - flag indicating our VCF came from GATK (False)")
    
if __name__ == "__main__":   
    __main__()
