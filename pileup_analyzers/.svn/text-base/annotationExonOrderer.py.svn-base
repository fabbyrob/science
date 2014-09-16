import sys
import re
import getopt
from Queue import PriorityQueue

_log = __file__.split("/")[-1]+".log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[2:],"l:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    annot = open(sys.argv[1],"r")
    
    if (annot == None):
        print ("Bad annotation name: "+sys.argv[1])
        sys.exit()
        
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print ("Bad logfile name: "+_log)
        sys.exit()
        
    #reg ex
    annotpat = re.compile("\w+_(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(.*)")
    
    negQueue = PriorityQueue()
    
    #read the vcf file...
    for line in annot:
        line = line.rstrip()
        m = annotpat.match(line)
        
        if (m == None):# if we find a weird line...
            sys.stderr.write("Non-standard line (annot):\n\t"+line+"\n")
            continue
        else:
            scaf = int(m.group(1))
            start = int(m.group(2))
            direction = m.group(5)
            
            if direction == "-":
                negQueue.put_nowait((start, line))
                continue
            else:
                while not negQueue.empty():
                    print(negQueue.get()[1])
                print(line)

   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()