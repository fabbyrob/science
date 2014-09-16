doc = """

"""

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
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    sites = open(sys.argv[2],"r")  
        
    if (sites == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()

    #read the infile file...
    for line in infile:
        line = line.rstrip()
        sline = line.split()

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
   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")

    
if __name__ == "__main__":   
    __main__()