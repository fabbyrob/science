import sys
import getopt

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    files = []
    lines = []
    for f in sys.argv[1:]:
        files.append(open(f, "r"))
        lines.append(files[-1].readline())
        
    
    if len(files) > 1:
        while lines[0] != "":
            likelihoods = []
            for l in lines:
                if l.split()[0] == "name":
                    continue
                likelihood = l.split()[1]
                if likelihood == "NA":
                    likelihoods.append(likelihood)
                else:
                    likelihoods.append(float(likelihood))
                
            ML = likelihoods[0]
            ML_i = 0
            for i, ml in enumerate(likelihoods[1:]):
                if (ML == "NA"  and ml != "NA") or (ml != "NA" and ml > ML):
                    ML = ml
                    ML_i = i+1
                    
            print(lines[ML_i].rstrip())
            
            #get the next line
            lines = []
            for f in files:
                lines.append(f.readline())
    else:#one file, with numbers of smaples in order
        last_name = ""
        likelihoods = []
        lines = []
        for line in files[0]:
            line = line.rstrip()
            sline = line.split()
            if sline[0] == "name":
                continue
            
            name = sline[0]
            if last_name != name:
                printBest(lines, likelihoods)
                last_name = name
                lines = []
                likelihoods = []
                
            lines.append(line)
            if sline[1] == "NA":
                likelihoods.append(sline[1])
            else:
                likelihoods.append(float(sline[1]))
                
        #grab the last one
        printBest(lines, likelihoods)
    
def printBest(lines, likelihoods):  
    if len(lines) == 0 or len(likelihoods) == 0:
        return 
    
    ML = likelihoods[0]
    ML_i = 0
    for i,ml in enumerate(likelihoods[1:]):
        if (ML == "NA"  and ml != "NA") or (ml != "NA" and ml > ML):
            ML = ml
            ML_i = i
    
    print(lines[ML_i])  

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
   
use = "python "+__file__.split("/")[-1]+" dfe_out_1 dfe_out_2 .... dfe_out_i"
def usage():
    print (use)
    sys.exit()
    
def details():
    print("Takes several output files from dfealpha.pl and compares them line by line. The run with the highest ML across the files in each line will be printed to stdout.")
    print (use)

    
if __name__ == "__main__":   
    __main__()
