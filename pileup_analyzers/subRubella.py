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
    input_pat = re.compile("\w+_(\S+)\s+(\d+).+\s+(\d+)\s+(\d+)")
    site_pat = re.compile("\w+_(\d+)\s+(\S+)")
    
    #load the list of sites of interest into memory
    mysites = []
    for line in sites:
        line = line.rstrip()
        m = site_pat.match(line)
        
        if (m == None):# if we find a weird line...
            sys.stderr.write("Non-standard line (sites:\n\t"+line+"\n")
            continue
        else:
            scaf = int(m.group(1))
            base = int(m.group(2))
            mysites.append((scaf, base))
    
    #initialize AFS
    afs = [0]*(int(_N/2.0+.5)+1)
    
    #read the input file...
    for line in input:
        if not mysites:#out of sites
            break
        line = line.rstrip()
        m = input_pat.match(line)
        
        if (m == None):# if we find a weird line...
            sys.stderr.write("Non-standard line (input):\n\t"+line+"\n")
            continue
        else:
            scaf = int(m.group(1))
            base = int(m.group(2))
            alt = int(m.group(3))
            ref = int(m.group(4))
            if (scaf, base) in mysites:
                mysites.pop(mysites.index((scaf,base)))
                if alt+ref == _N:
                    print(line)
                    if alt < ref:
                        count = alt
                    else:
                        count = ref
                    afs[count] += 1
                else:
                    sys.stderr.write("Weird number of alleles:\n\t"+line+"\n")
    
    for i in afs:
        afsfile.write(str(i)+"\t")
    afsfile.write("\n")
    
    for site in mysites:#print sites missing to log
        sys.stderr.write("Missing line:\n\t"+str(site)+"\n")

   
use = "python "+__file__.split("/")[-1]+" inputFile siteFile afsFile -N SampSize"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()