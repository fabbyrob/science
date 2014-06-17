doc = """
Takes a list of sites, and randomly splits it into N groups. The sites are grouped by regions, the size of which users can specify, then these regions are randomly assigend to groups. 
This program assumes that sites are sorted in the input file by chromosome and position.
Note that the output site list is NOT sorted. 
"""

import sys
import getopt
import random

_n = 2#number of groups
_w = 10000#region size
_p = "group"#prefix

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(2)
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    regions = []

    currentRegion = []
    currentScaf = ""
    start = 0
    end = start + _w
    #read the infile file...
    for line in infile:
        line = line.rstrip()
        sline = line.split()
        if line[0] == "#" or len(sline) < 2:
            continue
        
        scaf = sline[0]
        pos = int(sline[1])
        
        if currentScaf == "":
            currentScaf = scaf
            start = pos
            end = start + _w
            
        if scaf != currentScaf or (scaf == currentScaf and pos >= end):#end of prev region, start a new one
            regions.append((currentScaf, currentRegion))
            currentScaf = scaf
            start = pos
            end = start + _w
            currentRegion = []
    
        currentRegion.append(pos)
            
    
    if currentRegion:#grab the last one
        regions.append((currentScaf, currentRegion))
        
    #make all the output group files
    files = []
    for i in range(_n):
        files.append(open("%s_%s.txt" % (_p, i), "w"))
        
    #randomly assign the regions to files
    for scaf, pos in regions:
        myFile = random.choice(files)
        for p in pos:
            myFile.write("%s\t%s\n" % (scaf, p))
            
    #close the files
    for f in files:
        f.close()

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"n:w:p:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-n":
            global _n 
            _n = int(arg)
        elif opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-p":
            global _p 
            _p = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("site file: %s -p %s -n %s -w %s\n" % (sys.argv[1], _p, _n, _w))
   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("p - STR - %s - Prefix for the output filenames." % _p)
    print("n - INT - %s - The number of groups to split sites into" % _n)
    print("w - INT - %s - The size of regions in which to group sites (in base pairs)" % _w)

    
if __name__ == "__main__":   
    __main__()