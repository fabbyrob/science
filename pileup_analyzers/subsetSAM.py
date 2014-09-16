doc = """
Takes in a SAM file and a list of sites, will then output only lines from the SAM file that are represented in the site list.

Sites list MUST be sorted.
"""

import sys
import getopt

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    sites = open(sys.argv[2],"r")  
        
    if (sites == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()

    mySites = {}
    for line in sites:
        line = line.rstrip()
        sline = line.split()
        if len(sline) < 2:
            sys.stderr.write("Malformed sites line:\n\t%s\n" % line)
            continue
        
        if line.startswith("#"):
            continue
        
        scaf = sline[0]
        pos = int(sline[1])
        
        if scaf in mySites.keys():
            mySites[scaf][pos] = 1#make dictionaries because hash searching is much faster than searching through a list.
        else:
            mySites[scaf] = {pos:1}

    #read the infile file...
    for line in infile:
        line = line.rstrip()
        sline = line.split()

        if line.startswith("@"):
            print(line)
            continue
        
        scaf = sline[2]
        pos = int(sline[3])
        
        if search(mySites, scaf, pos):
            print(line)
        
def search(sites, myScaf, myPos):
    if not sites.has_key(myScaf):
        return False
    return sites[myScaf].has_key(myPos)
        

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"")
    except getopt.GetoptError:
        usage()
#    
#    for opt, arg in opts:
#        if opt == "-d":
#            global _d 
#            _d = True
#        else:
#            print ("Unrecognized option: "+opt+"\n")
#            usage()

    sys.stderr.write("samFile: %s sitesFile: %s\n" % (sys.argv[1], sys.argv[2]))
   
use = "python "+__file__.split("/")[-1]+" myFile.sam mySites.txt"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")

import time
if __name__ == "__main__": 
    start = time.time()
    __main__()
    sys.stderr.write("run time: %.2f\n" % ((time.time()-start)/60))  