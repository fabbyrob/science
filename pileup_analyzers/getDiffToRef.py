#!/usr/bin/python
import sys
import re
import getopt

_q = 40#min quality
_d = 20#min depth
_D = 60#max depth
_l = 60
_log = "getDiffToRef.log"#log file name
_p = False
_s = 30

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"q:d:D:l:L:s:p")
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
        elif opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-p":
            global _p
            _p = True
        elif opt == "-s":
            global _s
            _s = float(arg)
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
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    #reg ex
    if (_p):
        vcf_pat = re.compile("\w+_(\S+)\s(\d+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s\S+\s(\S+)")
    else:   
        vcf_pat = re.compile("\w+_(\S+)\s(\d+)\s\S+\s\S+\s(\S+)\s(\S+)\s\S+.+DP=(\d+);\S+\s\S+\s(\S+)")
    site_pat = re.compile("\w+_(\d+)\s+(\S+)")
         
    #initialize values
    next_site = getNewSite(sites, site_pat, logfile)
    
    num_sites = 0
    num_diff = 0
         
    if (next_site == None):
        logfile.write("No sites in site file\n")
        return
         
    countLast = 0
    #read the vcf file...
    for line in vcf:
        if (next_site == (-1,-1)):
            break
        line = line.rstrip()
        m = vcf_pat.match(line)
        
        if (m == None):
            logfile.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            if(_p):
                scaf = int(m.group(1))
                base = int(m.group(2))
                ref = m.group(3)
                alt = m.group(4)
                qual = float(m.group(5))
                snpqual = float(m.group(6))
                depth = int(m.group(7))           
            else:
                alt = ref = ""
                scaf = int(m.group(1))
                base = int(m.group(2))
                qual = float(m.group(4))
                depth = int(m.group(5))
                genos = m.group(6)
            
            #incase we pass a site
            while(scaf > next_site[0] or (scaf == next_site[0] and base > next_site[1])):
                logfile.write("Site not fount:\n\t"+str(next_site)+"\n")
                next_site = getNewSite(sites, site_pat, logfile)
                if (next_site == (-1,-1)):
                    break
            
            if(scaf == next_site[0] and base == next_site[1]):
                #print(scaf, " ", base, " ", next_site)
                next_site = getNewSite(sites, site_pat, logfile)
                #if we dont pass filters
                if (depth < _d or depth > _D or qual < _q and alt == ref):
                    countLast = 0
                    logfile.write("Site did not pass filters.\t\n"+str(next_site)+str(depth)+" "+str(qual)+"\n")
                    continue
                
                #print(line)
                
                if(not _p):
                    data = genos.split(":")
                    if([genos] == data):
                        num_sites += 1
                        logfile.write("Site accepted - same:\n\t"+line+"\n")
                        #logfile.write("homozygote, yay!")
                        continue
                    
                    data = data[0].split(",")
                    makeInts(data)
                    (min, max, mid) = getMinMaxIndex(data)
                    
                    if (min == 1):#heterozygotes should not exist
                        logfile.write("Possible heterozygote at:\n\t"+line+"\n")
                        continue
                    
                    if (data[mid] < _l):#site likelihoods are below quality thresholds
                        logfile.write("Likelihood quality too low:\n\t"+line+"\n")
                        continue
                    
                    num_sites += 1
                    if (max == 2):#homozygote non-ref
                        logfile.write("Site accepted - different:\n\t"+line+"\n")
                        num_diff += 1
                else:
                    if (alt == "*" or ref == "*"):
                        if countLast == 1:
                            num_sites -= 1
                            logfile.write("INDEL detected, uncounted site\n")
                        elif countLast == 2:
                            num_sites -= 1
                            num_diff -= 1
                            logfile.write("INDEL detected, uncounted site and divergence\n")
                        countLast = 0
                        continue
                                
                    countLast = 0        
                    if (alt != ref):
                        #non-ref
                        if(snpqual < _s):
                            #did not pass filter
                            logfile.write("Site did not pass filters, snpqual.\t\n"+str(next_site)+" "+str(snpqual)+"\n")
                        else:
                            #pass filter
                            num_sites += 1
                            logfile.write("Site accepted - different:\n\t"+line+"\n")
                            num_diff += 1
                            countLast = 2
                    else:
                        #reference
                        num_sites += 1
                        logfile.write("Site accepted - same:\n\t"+line+"\n")
                        countLast = 1
                        
       
    print ("No. Sites\tNo. Differences")
    print(str(num_sites)+"\t"+ str(num_diff))

#grabs a new site too look for from the sites input file
def getNewSite(sites, site_pat, logfile):
    next_site = None
    while(next_site == None):
        new_site = sites.readline()
        
        if (new_site == ""):
            return (-1,-1)
        
        m = site_pat.match(new_site)
        
        if (m == None):
            logfile.write("Non-standard line (SITE FILE):\n\t"+new_site+"\n")
            continue
        next_site = (int(m.group(1)), int(m.group(2)))
        
    return next_site
            
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
            
use = "usage: getDiffToRef.py <VCF> <SITES> [-q <QUALITY CUTOFF> -d <DEPTH MIN> -D <DEPTH MAX> -L <LIKELIHOOD> -p -s <SNPQUAL> -l <LOG FILE NAME>]\n"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("<VCF> - a variant call format file")
    print ("<SITES> - a file containing a list of sites of interest. Contains 2 columns, Scaffold Name and Position, each site on its own line.")
    print ("q - the minimum quality to include a site in the analysis (40)")
    print ("d - the minimum depth to include a site in the analysis (20)")
    print ("D - the maximum depth to include a site in the analysis (60)")
    print ("L - the minimum likelihood to call a site as homozygous, or heterozygous rather than N (60)")
    print ("p - indicates input is a pileup rather than a vcf (False)")
    print ("s - snp quality cut off (for sue with pileups only) (30)")
    print ("l - the name of the log file (getDiffToRef.log)")
    
if __name__ == "__main__":   
    __main__()
