import sys
import getopt
from math import log

_w = 10
_N = 26
_g = ""

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
        
    #read in and store all the SNP info into a list
    snps = []
    afs = []
    end = 0
    infile.readline()#get rid of header
    for line in infile:
        line = line.rstrip()
        sline = line.split()
        
        if len(sline) < 2:
            sys.stderr.write("Weird line:\n\t"+line+"\n")
            continue
        
        pos = int(sline[0])
        af = int(sline[1])
        snps.append(pos)
        afs.append(af)
        end = pos
    
    #calculate the global AFS
    num_snps = len(afs)
    total_afs = [None]
    for i in range(1, _N):
        total_afs.append(float(afs.count(i))/num_snps)

    if _g == "":
        #for each window calculate the AFS and LR
        print("midpoint\tLR")
        for start in range(0, len(snps), _w):
            #grab the window
            window = afs[start:start+_w]
            
            #calc the afs in this window
            window_afs = [None]
            window_snps = len(window)
            
            L_genome = 0
            L_window = 0
            
            for i in range(1, _N):
                curr_count = float(window.count(i))
                window_afs.append(curr_count/window_snps)
                
                #L_genome += log(total_afs[i]**curr_count)#because I need to log the product later, and Log(a*b) = Log(a) + Log(b), keeps the numbers from becoming too small and rounding to 0
                for j in range(int(curr_count)):
                    L_genome += log(total_afs[i])
                #L_window += log(window_afs[i]**curr_count)
                for j in range(int(curr_count)):
                    L_window += log(window_afs[i])
                
            endpos = min(len(afs)-1, start+_w-1)#position of last snp in window
            midpoint = float(snps[endpos]+snps[start])/2#midpoint of window
            LR = 2*(L_window-L_genome)#Likelihood ratio test
            print(str(midpoint)+ "\t"+str(LR))
    else:#we want specific genes
        print("midpoint\tLR\tgene\tnum_snps")
        genes = open(_g, "r")
        myGenes = []
        for line in genes:
            line = line.rstrip()
            sline = line.split()
            if len(sline) < 4:
                sys.stderr.write("Weird line:\n\t"+line+"\n")
                continue
            
            start = int(sline[1])
            end = int(sline[2])
            name = sline[3]
            myGenes.append((start, end, name))
            
        for gene in myGenes:
            start, end, name = gene
            
            starti = -1
            endi = -1
            for i, pos in enumerate(snps):
                if pos > end:
                    break
                if starti == -1 and pos >= start:
                    starti = i
                endi = i

            window = afs[starti:endi+1]
            
            #calc the afs in this window
            window_afs = [None]
            window_snps = len(window)
            
            if window_snps == 0:
                print(str(0)+ "\t"+str(0)+"\t"+name+"\t"+str(0))
                continue
            
            L_genome = 0
            L_window = 0
            
            for i in range(1, _N):
                curr_count = float(window.count(i))
                window_afs.append(curr_count/window_snps)
                
                #L_genome += log(total_afs[i]**curr_count)#because I need to log the product later, and Log(a*b) = Log(a) + Log(b), keeps the numbers from becoming too small and rounding to 0
                for j in range(int(curr_count)):
                    L_genome += log(total_afs[i])
                #L_window += log(window_afs[i]**curr_count)
                for j in range(int(curr_count)):
                    L_window += log(window_afs[i])
                
            endpos = endi#position of last snp in window
            midpoint = float(snps[endpos]+snps[starti])/2#midpoint of window
            LR = 2*(L_window-L_genome)#Likelihood ratio test
            print(str(midpoint)+ "\t"+str(LR)+"\t"+name+"\t"+str(len(window)))
                        
            
            
            
        
def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:N:g:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-N":
            global _N 
            _N = arg
        elif opt == "-g":
            global _g 
            _g = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" INFILE [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()