doc = """
Takes a summary file and pulls out divergence and AFS for each intron.
"""

import sys
import getopt
sys.path.append("/data/robert.williamson/bin")

import summary


_v = False
_N = 26

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    reader = summary.Reader(open(sys.argv[1], "rb"))
    
    afs = [0]*(_N+1)
    div = [0,0]#checked, diverged
    start = 0
    pSite = 0
    
    header = "chrom\tstart\tstop\tlen\tdivergenceSites\tdivergence\tafs_"
    for i in range(_N+1):
        header += str(i)+"\tafs_"
    header = header[:-4]
    print(header)
    
    for record in reader:
        if _v: sys.stderr.write("On line:\t%s\n" % record.prettyStr())
        
        if 'intron' in record.typeNames():
            if start == 0:
                if _v: sys.stderr.write("Start intron: %s.\n" % record.prettyStr())
                start = record.POS
                
            pSite = record.POS
            
            if record.TOTAL != _N:
                continue
            
            afs[record.ALT_NUM] += 1
            
            if record.DIVERGENCE != -1:
                div[0] += 1
                if record.DIVERGENCE:
                    div[1] += 1
        elif pSite != 0:
            if _v: sys.stderr.write("End of intron at %s.\n"% (record.prettyStr()))
            #scaf, start, stop, length, divergenceChecked, divergence, afs
            afsStr = "\t".join(map(str, afs))
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (record.CHROM, start, pSite, pSite-start, div[0], div[1], afsStr))
            pSite = 0
            start = 0
            afs = [0]*(_N+1)
            div = [0,0]
    

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"N:v")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-v":
            global _v 
            _v = True
        elif opt == "-N":
            global _N 
            _N = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -N %s -v %s\n" % (sys.argv[1], _N, _v))
   
use = "python "+__file__.split("/")[-1]+" SummaryFile"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("N - INT - %s - Sample size" % _N)

    
if __name__ == "__main__":   
    __main__()