import sys
import re
import getopt

_N = 26
_n = 24
_log = __file__.split("/")[-1]+".log"#log file name
_siteFile = "sites.txt"

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"l:N:n:s:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-n":
            global _n 
            _n = int(arg)
        elif opt == "-s":
            global _siteFile 
            _siteFile = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    grand = open(sys.argv[1],"r")
    
    if (grand == None):
        print("Bad input name: "+sys.argv[1])
        sys.exit()
        
    rub = open(sys.argv[2],"r")  
        
    if (rub == None):
        print("Bad sites file name: "+sys.argv[2])
        sys.exit()
        
    siteFile = open(_siteFile, "w")
        
    #reg ex
    grand_pat = re.compile("\w+_(\S+)\s+(\d+)\s+(\d+)")
    rub_pat = re.compile("\w+_(\S+)\s+(\d+)\s+(\d+)")
    
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
    
    next_site = getNextSite(rub, rub_pat)
    
    counts = {}
    for i in range(0,_N+1):
        counts[i] = [0]*(_n+1)
    
    
    #print("scaffold\tbase\tCgAF\tCrAF")
    #read the input file...
    for line in grand:
        if next_site == None:
            sys.stderr.write("Ran out of rubella sites at this location in grandiflora: \n\t"+line+"\n")
            break
#        if not mysites:#out of sites
#            break
        line = line.rstrip()
        m = grand_pat.match(line)
        
        if (m == None):# if we find a weird line...
            sys.stderr.write("Non-standard line (input):\n\t"+line+"\n")
            continue
        else:
            scaf = int(m.group(1))
            base = int(m.group(2))
            afs = int(m.group(3))
            
            if scaf == next_site[0] and base == next_site[1]:
                siteFile.write(line+"\t"+str(next_site[2])+"\n")
                counts[afs][next_site[2]] += 1
                next_site = getNextSite(rub, rub_pat)
            elif (scaf == next_site[0] and base < next_site[1]) or scaf < next_site[0]:
                sys.stderr.write("Site not found (grandiflora only): \n\t"+line+"\n\t"+str(next_site)+"\n")
                continue
            else:# scaf == next_site[0] and base > next_site[0] or scaf > next_site[0]:
                while next_site != None and (next_site[0] < scaf or (next_site[0] == scaf and next_site[1] < base)):
                    sys.stderr.write("Site not found (rubella only): \n\t"+str(next_site)+"\n")
                    next_site = getNextSite(rub, rub_pat)
                if (scaf == next_site[0] and base == next_site[1]):
                    counts[afs][next_site[2]] += 1
                    next_site = getNextSite(rub, rub_pat)
                else:
                    sys.stderr.write("Site not found (grandiflora only): \n\t"+line+"\n\t"+str(next_site)+"\n")
                    
            
    siteFile.close()
    
    print("Cg\tCr\tNum_Sites")
    for i in range(0,_N+1):
        for j in range(0,_n+1):
            print(str(i)+"\t"+str(j)+"\t"+str(counts[i][j]))
            
            
        
#    for site in mysites:#print sites missing to log
#        sys.stderr.write("Missing line:\n\t"+str(site)+"\n")

def getNextSite(file, pat):
    dat = None
    while(dat == None):
        aLine = file.readline()
        
        if (aLine == ""):
            return
        
        m = pat.match(aLine)
        
        if(m == None):
            continue
        scaf = int(m.group(1))
        site = int(m.group(2))
        afs = int(m.group(3))
        
        dat = (scaf, site, afs)
    return dat
   
use = "python "+__file__.split("/")[-1]+" grandifloraFile rubellaFile"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()