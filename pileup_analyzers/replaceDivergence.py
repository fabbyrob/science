doc = """
Replaces the divergence in a summary file with new divergence estimates.
"""

import sys
import getopt
import summary
import cPickle as pickle

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    summaryReader = summary.Reader(open(sys.argv[1],"rb"))
    div = pickle.load(open(sys.argv[2], "rb"))

    print(summaryReader.Header)

    #read the infile file...
    for site in summaryReader:
        if div.has_key(site.CHROM):
            try:
                site.DIVERGENCE = div[site.CHROM][site.POS]
            except IndexError:
                sys.stderr.write("No divergence for site %s %s in divergence file. Assigning a -1.\n" % (site.CHROM, site.POS))
                site.DIVERGENCE = -1
        else:
            sys.stderr.write("No divergence for site %s %s in divergence file (Missing CHROM in divergence). Assigning a -1.\n" % (site.CHROM, site.POS))
            site.DIVERGENCE = -1
        print(site)

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
   
use = "python "+__file__.split("/")[-1]+" summary.txt divergence.pickle"
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