import sys
import re
import getopt

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
    
    input = open(sys.argv[1],"r")
    
    if (input == None):
        print("Bad input name: "+sys.argv[1])
        sys.exit()
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    numInds = -1
    bins = []
    #read the input file...
    for line in input:
        line = line.split()
        if numInds == -1:
            numInds = len(line)-6
            bins = [0]*(numInds+1)
        goodSites = int(line[-1])
        altCount = int(line[-2])
        
        if goodSites != numInds:
            sys.stderr.write("Missing data at site: "+line[0]+" "+line[1]+"\n")
            continue
        
        bins[altCount] += 1
        
    #print(bins)
    #fold the afs
    foldedBins = []
    for i in range(0,int(len(bins)/2+len(bins)%2)):
        if len(bins)%2 == 1 and i == int(len(bins)/2-.5):#if we have an odd-length AFS then the middle one should not be added to anything
            foldedBins.append(bins[i])
        else:
            foldedBins.append(bins[i]+bins[-1*(i+1)])
        
    myStr = ""
    for f in foldedBins:
        myStr += str(f)+"\t"
        
    print(myStr)
        

   
use = "python "+__file__.split("/")[-1]+" rubInput"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()