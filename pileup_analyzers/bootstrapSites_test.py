import sys
import getopt
from Queue import PriorityQueue
from random import sample

_log = __file__.split("/")[-1]+".log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
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
    
    sites = open(sys.argv[1],"r")
    
    if (sites == None):
        print("Bad vcf name: "+sys.argv[1])
        sys.exit()
        
    
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    outfile = open(sys.argv[2], "w")

    #reg ex
    siteslist = []
    
    #read the vcf file...
    for line in sites:
        line = line.rstrip()
        lineSplit = line.split()
        if len(lineSplit) >= 2:
            base = int(lineSplit[1])
            siteslist.append((base, line))
            continue
        else:
            sys.stderr.write("Weird line: "+str(len(lineSplit))+"..."+line+"\n")

    newSiteList = []
    while len(newSiteList) < len(siteslist):
        newSiteList.append(sample(siteslist, 1)[0])

    newSiteList.sort()
    printList = []
    while len(newSiteList) > 0:
        printList.append(newSiteList[0][1])
        newSiteList.pop(0)
    
    out = "\n".join(printList)
    outfile.write(out)
    outfile.close()

   
use = "python "+__file__.split("/")[-1]+" SITES"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("This program takes a list of sites, and bootstraps it. The same number of sites in the input file will be printed.")
    print(">python bootsrapSites.py mySites.txt > myNewSites.txt")

    
if __name__ == "__main__":   
    __main__()
