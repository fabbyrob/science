doc = """

"""

import sys
import getopt
from Bio.AlignIO import MafIO
from Bio import AlignIO

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    for alignment in AlignIO.parse(sys.argv[1], "maf"):
        print("alignment")
        for seqrec in alignment:
            print ("starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % (seqrec.annotations["start"], seqrec.annotations["strand"], seqrec.annotations["srcSize"], seqrec.annotations["size"]))

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