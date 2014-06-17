import sys
import getopt
from subprocess import call, Popen

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
#    try: 
#        opts, args = getopt.getopt(sys.argv[3:],"")
#    except getopt.GetoptError:
#        usage()
#    
#    for opt, arg in opts:
#        if opt == "-l":
#            global _log 
#            _log = str(arg)
#        else:
#            print ("Unrecognized option: "+opt+"\n")
#            usage()
    
#    vcf = open(sys.argv[1],"r")
#    
#    if (vcf == None):
#        print("Bad vcf name: "+sys.argv[1])
#        sys.exit()
        
    annot = open(sys.argv[1],"r")  
        
    if (annot == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()
        
    intron = getIntron(annot)
    print("Gene\tIntron_Ord\tLength\tNum_Sites\tDiv_Sites")#output header
    while intron:
        scaf = intron[0]
        start = intron[1]
        end = intron[2]
        gene = intron[3]
        dir = intron[4]
        ord = intron[5]
        
        sites = genSites(scaf, start, end)
        
        div = getDiv(sites, scaf)
        
        print(gene+"\t"+str(ord)+"\t"+str(end-start)+"\t"+str(div[0])+"\t"+str(div[1]))
        
        intron = getIntron(annot, gene, ord)

def getIntron(annot, prevGene = None, ord = 0):
    line = annot.readline()
    if line == "":
        return None
    
    line = line.rstrip()
    line = line.split()
    
    scaf = line[0]
    start = int(line[1])
    end = int(line[2])
    gene = line[3]
    dir = line[4]
    
    if (prevGene == None or prevGene != gene):
        ord = 0
    else:
        ord += 1
        
    return ((scaf, start, end, gene, dir, ord))
 
def genSites(scaf, start, end, file = "temp.sites"):
    output = open(file, "w")
    for i in range(start, end+1):
         output.write(scaf+"\t"+str(i)+"\n")
    output.close()
    return file

def getDiv(file, scaf, outfile = "temp.div"):
    scafNum = int(scaf.split("_")[1])
    cmd = "python ~/bin/divergence2.py "+file+" /data/mathieu.blanchette/alignmentsStephen/CR_v_NP.chain.ortho.Scaffold"+str(scafNum)+".sorted.maf > "+outfile 
    try:
        p = Popen(cmd, shell=True)
        p.wait()
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)
        
    data = open(outfile, "r")
    data.readline()
    div = data.readline().split()
    data.close()
    return(int(div[0]), int(div[1]))
   
use = "python "+__file__.split("/")[-1]+" <Intron Annotation>"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()