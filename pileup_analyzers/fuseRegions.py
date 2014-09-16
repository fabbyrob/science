doc = """
Fuses regions in a file that are close together. 

Made for use with IBD regions identified from checkIBD, but can be used on any list of regions.
"""

import sys
import getopt

_d = 1000

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()

    regions = {}#regions sorted by ind
    #read the infile file...
    for line in infile:
        line = line.rstrip()
        sline = line.split()
        
        if sline[0] == "CHROM":
            continue
        
        if len(sline) != 4:
            continue
        
        scaf = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        ind = sline[3]
        
        if ind in regions.keys():
            if scaf in regions[ind].keys():
                regions[ind][scaf].append((start,end))
            else:
                regions[ind][scaf] = [(start,end)]
        else:
            regions[ind] = {}
            regions[ind][scaf] = [(start,end)]
            
    #fuse the regions...
    for ind in regions.keys():
        for scaf in regions[ind].keys():
            myRegions = regions[ind][scaf]
            newRegions = []
            prev = myRegions[0]
            
            for start, end in myRegions[1:]:
                if start - prev[1] < _d:
                    sys.stderr.write("Fusing regions on %s in sample %s:\n\t%s\t%s\n\t%s\t%s\n" % (scaf, ind, prev[0], prev[1], start, end))
                    prev = (prev[0], end)
                else:
                    newRegions.append(prev)
                    prev = (start, end)
            newRegions.append(prev)
            
            regions[ind][scaf] = newRegions
            
    #print regions
    print("CHROM\tSTART\tEND\tSAMPLE")
    for ind in regions.keys():
        for scaf in regions[ind].keys():
            for start, end in regions[ind][scaf]:
                print("%s\t%s\t%s\t%s\n" % (scaf, start, end, ind))
        

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"d:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-d":
            global _d 
            _d = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -arg %s\n" % (sys.argv[1], _d))
   
use = "python "+__file__.split("/")[-1]+" regionFile [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("d - INT - %s - the maximum distance between two regions that would allow them to be fused." % _d)

    
if __name__ == "__main__":   
    __main__()