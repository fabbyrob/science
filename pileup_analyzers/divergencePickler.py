import sys
import re
import getopt
import cPickle as pickle

_w = 5
_pickle = "pickle.p"
_makePickle = True

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()

    try: 
        opts, args = getopt.getopt(sys.argv[2:],"w:o:d:")
    except getopt.GetoptError:
        usage()
    
    global _pickle
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-o": 
            _pickle = arg
        elif opt == "-d":
            global _makePickle 
            _makePickle = False
            _pickle = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
	files = sys.argv[1].split(",")
    
    if _makePickle:
        infile = open(sys.argv[1],"r")
        
        if (infile == None):
            print("Bad vcf name: "+sys.argv[1])
            sys.exit()
     
    if _makePickle:
        pickleDivergence(infile)
    else:
        unpickleDivergence(files)
        
    #infile.close()
        
def unpickleDivergence(infiles):
    divergence = pickle.load(open(_pickle, "rb"))
    #print(divergence)
    
    for file in infiles:
        infile = open(file, "r")
        outfile = open(file+".div", "w")
        sites = 0
        diverged = 0
        for line in infile:
            line = line.rstrip()
            line = line.split()
            if len(line) < 2:
                sys.stderr.write("Weird input line:\n"+file+":\t"+str(line)+"\n")
                continue
        	
            base = int(line[1])
        	#print(base,divergence[base])
            if base >= len(divergence):
                break
        	
            if divergence[base][2] != None:
        		#print(str(base) + ": "+str(divergence[base]))
                sites += 1
                if not divergence[base][2]:
                    diverged += 1
           
        outfile.write("#sites\t#diverged\n")         
        outfile.write(str(sites)+"\t"+str(diverged)+"\n")
        outfile.close()
  
def pickleDivergence(infile):        
    #read the vcf file...
    divergence = [None]
    rubella = ""
    neslia = ""
    start = -1
    for line in infile:
        line = line.rstrip()
        if len(line) == 0:
            continue
        
        if line[0] == "#":
            continue
        elif line[0] == 'a':
            process(rubella, neslia, start, divergence)

            #process
            #reset variables
            rubella = ""
            neslia = ""
            start = -1
        elif line[0] == 's' and rubella == "":
            rubella = line.split()[-1]
            start = int(line.split()[2])+1#add one because this file is zero indexed, but the site file isn't 
        else:
            neslia = line.split()[-1]

    process(rubella, neslia, start, divergence)
    #pickle the scaffold
    #print(len(divergence))
    pickle.dump(divergence, open(_pickle, "wb"))
    #print(divergence)
    #print(len(divergence))

def process(rubella, neslia, start, divergence):
    if rubella == "" or neslia == "":
        return
    
    while len(divergence) < start:
       divergence.append((None, None, None))
    
    window_left = 0
    #print(len(neslia))
    #print(len(rubella))
    for i in range(0,len(rubella)):
        if neslia[i] == "-" or rubella[i] == "-":
            window = min(len(divergence)-start, _w)  
            #divergence[len(divergence)-window:]=[None]*window 
            fixedLoci = []
            for locus in divergence[len(divergence)-window:]:
                fixedLoci.append((locus[0], locus[1], None))
            divergence[len(divergence)-window:]=fixedLoci
                
            current = None
            window_left = _w
        elif window_left > 0:
            current = None
            window_left -= 1
        elif neslia[i] == rubella[i]:
            current = True
        else:
            current = False
        
        if rubella[i] !="-":#if its  a deletion in rubella we don't want to include it in the list
            divergence.append((rubella[i], neslia[i], current))
    
use = "python "+__file__.split("/")[-1]+" <INFILE> [-w <window size> -o <pickle filename> -d <pickle filename>]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("INFILE - either a global alignment to get pickled, or a list of sites to check for divergence in the file provided with -d")
    print("w (5) - the window around indels to ignore")
    print("o (pickle.p) - if in pickle mode, the file to save the pickled genome")
    print("d - indicates that we are checking divergence from a list of sites, should provide a pickle filename from which to check divergence")

    
if __name__ == "__main__":   
    __main__()