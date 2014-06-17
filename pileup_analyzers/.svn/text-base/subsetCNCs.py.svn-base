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
    
    cncfile = open(sys.argv[1],"r")
    
    if (cncfile == None):
        print("Bad cncfile name: "+sys.argv[1])
        sys.exit()
        
    genefile = open(sys.argv[2],"r")  
        
    if (genefile == None):
        print("Bad genefile name: "+sys.argv[2])
        sys.exit()

    myGenes = []
    #read the genes file...
    for line in genefile:
        line = line.rstrip()
        sline = line.split()
        myGenes.append(sline[0])
        
    myCNCs = []
    for line in cncfile:
        line = line.rstrip()
        sline = line.split()
        
        if len(sline) < 2:
            continue
        
        name = sline[0]
        cnc = float(sline[1])
        
        if name in myGenes:
            print(line)
            myCNCs.append(cnc)
            
    mean = sum(myCNCs)/float(len(myCNCs))
    sq_diff = 0
    for c in myCNCs:
        sq_diff += (c - mean)**2
        
    stdiv = (sq_diff/len(myCNCs))**0.5
    print("mean %s\nstdiv %s" % (mean, stdiv))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        continue
#        if opt == "-o":
#            global _o 
#            _o = arg
#        else:
#            print ("Unrecognized option: "+opt+"\n")
#            usage()

    sys.stderr.write("infile: %s -arg %s\n" % (sys.argv[1], 1))
   
use = "python "+__file__.split("/")[-1]+" CNC_FILE GENE_FILE"
def usage():
    print (use)
    sys.exit()
    
def details():
    print("takes a file with the number of CNCs per gene, and a list of genes, and only outputs the CNC info for that subset of genes. Calculates the mean and standard deviation.")
    print("")
    print (use)
    print("\npython subsetCNCs.py myCNCs.txt myGenes.txt")
    print("PAC:202 10\nPAC:204 0\n...\naverage 12.7\nstdiv 32.1")

    
if __name__ == "__main__":   
    __main__()