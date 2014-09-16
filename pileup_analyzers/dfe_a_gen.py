import sys
import re
import getopt

_log = __file__.split("/")[-1]+".log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
#    try: 
#        opts, args = getopt.getopt(sys.argv[3:],"l:")
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
    
    afsfile = open(sys.argv[1],"r")
    
    if (afsfile == None):
        print("Bad afs name: "+sys.argv[1])
        sys.exit()
        
    div = open(sys.argv[2],"r")  
        
    if (div == None):
        print("Bad div file name: "+sys.argv[1])
        sys.exit()
       
    namesfile = open(sys.argv[3],"r")  
        
    if (namesfile == None):
        print("Bad names file name: "+sys.argv[1])
        sys.exit()
        
    samps = int(sys.argv[4])
    ref = map(int, sys.argv[5].split(","))
    #print(ref)
    prefix = sys.argv[6]
    
    afs = []
    for line in afsfile:
        #print(line)
        line = line.rstrip().lstrip()
        afs.append(line)
        
    divs = []
    for line in div:
        items = line.split()
        if len(items) == 2:
            divs.append((items[0], items[1]))
         
    names = []
    for line in namesfile:
        #print(line)
        line = line.rstrip()
        names.append(line)
        
    zeros = ""
    for i in range(0, int(samps/2)):
        zeros += " 0"
        
    i = -1
    outfile = open(prefix+".txt", "w")
    outfile.write(str(len(afs)-1)+"\n1\n"+str(samps)+"\n")#write sample size
    written = 1
    lastRef = ref[0]
    for a in afs:
        i += 1
        if i in ref:
            lastRef = i
            continue
        
        outfile.write(str(written)+"\n")#write sample size
        outfile.write(str(divs[i][0])+"\n"+str(divs[i][1])+"\n")#write sample divergence
        outfile.write(str(divs[lastRef][0])+"\n"+str(divs[lastRef][1])+"\n")#write reference divergence
        outfile.write(a+zeros+"\n")#write sample afs
        outfile.write(afs[lastRef]+zeros+"\n")#write ref afs
        written += 1

   
use = "python "+__file__.split("/")[-1]+" AFS DIV NAMES samp_size reference_row(0 indexed) prefix_filename"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()