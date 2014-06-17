import sys
import getopt

_l = 0

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    genes = open(sys.argv[1],"r")
    
    if (genes == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    cncs = open(sys.argv[2],"r")  
        
    if (cncs == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()

    #(name, scaf, start, stop)
    myGenes = []

    #read the infile file...
    for line in genes:
        line = line.rstrip()
        sline = line.split()
        if len(sline) < 4:
            continue
        
        myGenes.append((sline[0], int(sline[1].split("_")[1]), int(sline[2]), int(sline[3])))
        
    myCNCs = [0]*len(myGenes)
    
    current = 0
    #scaf start stop name dir
    for line in cncs:
        line = line.rstrip()
        sline = line.split()
        
        scaf = int(sline[0].split("_")[1])
        start = int(sline[1])
        end = int(sline[2])
        
        while scaf > myGenes[current][1]:
            current += 1
        
        if start > myGenes[current][3] + _l:
            current += 1
        
        myCNCs[current] += calcOverlap(scaf, start, end, myGenes[current])
        
        if current+1 < len(myGenes):
            myCNCs[current+1] += calcOverlap(scaf, start, end, myGenes[current+1])
#            
#        if scaf == myGenes[current][1] and ((end >= myGenes[current][2] - _l and end <= myGenes[current][3] + _l) or (start <= myGenes[current][3] + _l and start >= myGenes[current][2] - _l)):
#            myCNCs[current] += 1
#            
#        if current+1 < len(myGenes) and myGenes[current+1][1] == scaf and (end >= myGenes[current+1][2] - _l and end <= myGenes[current+1][3] + _l) or (start <= myGenes[current+1][3] + _l and start >= myGenes[current+1][2] - _l):
#            myCNCs[current+1] += 1
        
    tot = sum(myCNCs)
    for i in range(len(myGenes)):
        print("%s %s" % (myGenes[i][0], myCNCs[i]))
        
    mean = tot/float(len(myGenes))
    sq_diff = 0
    for c in myCNCs:
        sq_diff += (c - mean)**2
        
    stdiv = (sq_diff/len(myGenes))**0.5
    print("mean %s\nstdiv %s" % (mean, stdiv))

def calcOverlap(cnc_scaf, cnc_start, cnc_end, gene):
    overlap = 0
    if cnc_scaf == gene[1]:
        #whole CNC overlaps gene
        if cnc_start >= gene[2] - _l and cnc_end <= gene[3] + _l:
            overlap += cnc_end - cnc_start + 1
            
        #right of CNC overlaps
        elif cnc_start < gene[2] - _l and cnc_end >= gene[2] - _l:
            overlap += cnc_end - max(1, (gene[2] - _l)) + 1#the max makes sure that if its at the begining of the scaf then we dont count negative overlap
        
        #left of CNC overlaps
        elif cnc_start < gene[3] + _l and cnc_end > gene[3] + _l:
            overlap += (gene[3] + _l) - cnc_start + 1
    
    return overlap

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"l:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _l 
            _l = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -l %s\n" % (sys.argv[1], _l))
   
use = "python "+__file__.split("/")[-1]+" gene_file cnc_file [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("")
    print("Takes in a list of genes, and a list of CNCs. Counts the number of CNCs that overlap each gene, and the regions flanking a gene.")
    print("python cnc_counter.py myGenes.txt myCNCs.txt -l 1000")
    print("gene_1 24\ngene_2 34\ngene_3 105\naverage 88.6\nstd_dev 24.2")
    
if __name__ == "__main__":   
    __main__()