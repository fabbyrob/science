import sys
import getopt
import random

_w = 20000
_N = 26

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
        

    prevPolys = []
    prevDivs = []
    polys = 0
    divs = 0
    start = -1
    end = -1
    #read the infile file...
    for line in infile:
        if line[0] == "#":
            continue
        
        line = line.rstrip()
        sline = line.split()
        
        base = int(sline[1])
        alt = int(sline[4])
        ref = int(sline[5])
        tot = int(sline[6])
        div = int(sline[-1])
        
        if start == -1:
            start = base
            end = start + _w
            
        if base >= end:
            prevPolys.append(polys)
            prevDivs.append(divs)
            
            polys = 0
            divs = 0
            
            start = end
            end = start + _w
        
        
        if tot == _N and alt + ref == tot:
            if alt and alt < _N:
                polys += 1
            
            if div == 1:
                divs += 1
        
    #grab the last window
    prevPolys.append(polys)
    prevDivs.append(divs)
        
    ratio = [float(a)/float(b) for a, b in zip(prevPolys, prevDivs) if b > 0]
    
    random.shuffle(prevPolys)
    
    randomRatio = [float(a)/float(b) for a, b in zip(prevPolys, prevDivs) if b > 0]
    
    print("Ratio\tRandomized")
    for vals in zip(ratio, randomRatio):
        print("%s\t%s" % vals)
     
        

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:N:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -w %s -N %s\n" % (sys.argv[1], _w, _N))
   
use = "python "+__file__.split("/")[-1]+" summaryFile [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print("Takes in a summary of a VCF file, with divergence info. Calculates polymorphism and divergence in windows, then takes the ratio of the two in each window. Also randomizes the polymorphism/divergence pairs to approximate a null distribution.")
    print (use)
    print("")
    print("OPTION - TYPE - DEFAULT - Description")
    print("w - INT - %s - the window size to analyze" % _w)
    print("N - INT - %s - sample size" % _N)

    
if __name__ == "__main__":   
    __main__()