import sys
import re
import getopt

_N = 24
_log = __file__.split("/")[-1]+".log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 4:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[4:],"l:N:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    input = open(sys.argv[1],"r")
    
    if (input == None):
        print("Bad input name: "+sys.argv[1])
        sys.exit()
        
    sites = open(sys.argv[2],"r")  
        
    if (sites == None):
        print("Bad sites file name: "+sys.argv[2])
        sys.exit()
       
    afsfile = open(sys.argv[3],"w")  
        
    if (afsfile == None):
        print("Bad afs file name: "+sys.argv[3])
        sys.exit()
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    #reg ex
    input_pat = re.compile("\w+_(\S+)\s+(\d+)\s+\S+\s+\S+\s+(.+)\s+(\d+)\s+(\d+)")
    site_pat = re.compile("\w+_(\d+)\s+(\S+)")
    
    #load the list of sites of interest into memory
#    mysites = []
#    for line in sites:
#        line = line.rstrip()
#        m = site_pat.match(line)
#        
#        if (m == None):# if we find a weird line...
#            sys.stderr.write("Non-standard line (sites:\n\t"+line+"\n")
#            continue
#        else:
#            scaf = int(m.group(1))
#            base = int(m.group(2))
#            mysites.append((scaf, base))
    
    next_site = getNextSite(sites, site_pat)
    #initialize AFS
    afs = [0]*(int(_N/2.0+.5)+1)
    
    #read the input file...
    for line in input:
        if next_site == None:
            break
#        if not mysites:#out of sites
#            break
        line = line.rstrip()
        sline = line.split()
        if len(sline) != _N + 6:
            sys.stderr.write("Bad line: \n\t"+line)
            continue
        
        
        scaf = int(sline[0].split("_")[1])
        base = int(sline[1])
        genos = sline[3:3+_N]
        alt = int(sline[-2])
        ref = int(sline[-1])
        
        while next_site != None and scaf == next_site[0] and base == next_site[1]:
            safe = 0
            #print(genos.split())
            for ind in genos:
                if ind in ['A','T','C','G']:
                    safe += 1
                
            if safe == _N:
                print("scaffold_"+str(scaf)+"\t"+str(base)+"\t"+str(alt))
                if alt < _N-alt:
                    count = alt
                else:
                    count = _N-alt
                afs[count] += 1
            else:
                sys.stderr.write("Weird number of alleles "+str(safe) +":\n\t"+line+"\n")
                
            next_site = getNextSite(sites, site_pat)
        if next_site == None or (scaf == next_site[0] and base < next_site[1]) or scaf < next_site[0]:
            continue
        else:# scaf == next_site[0] and base > next_site[0] or scaf > next_site[0]:
            while next_site != None and (next_site[0] < scaf or (next_site[0] == scaf and next_site[1] < base+1)):
                sys.stderr.write("Site not found: \n\t"+str(next_site)+"\n")
                next_site = getNextSite(sites, site_pat)
            
    
    for i in afs:
        afsfile.write(str(i)+"\t")
    afsfile.write("\n")
    
#    for site in mysites:#print sites missing to log
#        sys.stderr.write("Missing line:\n\t"+str(site)+"\n")

def getNextSite(file, pat):
    dat = None
    while(dat == None):
        line = file.readline()
        
        if (line == ""):
            return
        
        line = line.split()
        if len(line) < 2:
            sys.stderr.write("Bad line:\n\t"+str(line)+"\n")
            continue
        
        scaf = int(line[0].split("_")[1])
        site = int(line[1])
        dat = (scaf, site)
    return dat
   
use = "python "+__file__.split("/")[-1]+" inputFile siteFile afsFile -N SampSize"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()