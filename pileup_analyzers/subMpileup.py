#!/usr/bin/python
import sys
import re
import getopt

_q = 40#min quality
_d = 20#min depth
_D = 60#max depth
_log = "subMpileup.log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"l:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
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
    vcf_pat = re.compile("\w+_(\d+)\s(\d+).+")
    site_pat = re.compile("\w+_(\d+)\s+(\S+)")
        
    #initialize values
    next_site = getNewSite(sites, site_pat, logfile)
    
    if (next_site == None):
        logfile.write("Site file empty\n")
        return
    
    #read the vcf file...
    for line in vcf:
        line = line.rstrip()
        m = vcf_pat.match(line)
        
        if (m == None):# if we find a weird line...
            logfile.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            #grab the relevant data
            scaf = int(m.group(1))
            base = int(m.group(2))
            
            #if we're not at a site we care about...
            if(scaf != next_site[0] or base != next_site[1]):
                while (scaf > next_site[0] or (scaf == next_site[0] and base > next_site[1])):
                    logfile.write("\t FAIL at : "+str(scaf) +", "+str(base)+"\n")
                    logfile.write("site not found: "+str(next_site)+"\n")
                    next_site = getNewSite(sites, site_pat, logfile)
                    #logfile.write("site:\t"+str(next_site)+"\n")
                    if (next_site == None):
                        break
                    if (scaf == next_site[0] and base == next_site[1]):
                        print(line)
                        next_site = getNewSite(sites, site_pat, logfile)
                        if (next_site == None):
                            break
            else: 
                #grab the next site
                while next_site[0] == scaf and next_site[1] == base:
                    next_site = getNewSite(sites, site_pat, logfile)
                    print(line)
                    if (next_site == None):
                        break
                    
                #logfile.write("site:\t"+str(next_site)+"\n")
            if (next_site == None):
                break
   
#grabs a new site too look for from the sites input file
def getNewSite(sites, site_pat, logfile):
    next_site = None
    while(next_site == None):
        new_site = sites.readline()
        
        if (new_site == ""):
            return None
        
        m = site_pat.match(new_site)
        
        if (m == None):
            logfile.write("Non-standard line (SITE FILE):\n\t"+next_site+"\n")
            continue
        next_site = (int(m.group(1)), int(m.group(2)))
        
    return next_site
   
def usage():
    print ("usage: filter.py <VCF> <SITE_LIST> [-l <LOG FILE NAME>]\n")
    sys.exit()
    
def details():
    print ("usage: filter.py <VCF> <SITE_LIST> [-l <LOG FILE NAME>]\n")
    print ("<VCF> - a variant call format file")
    print ("<SITE_LIST> - a file containing a list of sites of interest. Contains 2 columns, Scaffold Name and Position, each site on its own line.")
    #print ("q - the minimum quality to include a site in the analysis (40)")
    #print ("d - the minimum depth to include a site in the analysis (20)")
    #print ("D - the maximum depth to include a site in the analysis (60)")
    print ("l - the name of the log file (subMpileup.log)")
    
if __name__ == "__main__":   
    __main__()
