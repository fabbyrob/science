doc = """
takes a summary and splits it into windows, reports the % of sites in window that are of a certain type.

python codingDensity.py mySummary.txt -t exon -t 4fold -t 3utr > myPercentages.txt
"""

import sys
import getopt
import summary

_w = 1000
_t = []

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    summaryReader = summary.Reader(open(sys.argv[1],"rb"))
    
    global _t
    safeT = []

    myCounts = {}

    for t in _t:
        if t not in summaryReader.summary.typeToCode.keys():
            sys.stderr.write("Type %s not found in summary. It will be ignored. The valid options are: %s.\n" % (t, summaryReader.summary.typeToCode.keys()))
        else:
            safeT.append(summaryReader.summary.typeToCode[t])
            myCounts[summaryReader.summary.typeToCode[t]] = []
    
    _t = safeT
    
    if len(_t) == 0:
        sys.stderr.write("No valid types to count were provided. Use the -t option to specify sites of interest.\n")
        sys.exit(0)

    start = None
    pScaf = None
    pSite = None

    print("TYPE\tCHROM\tMIDPOINT\tCOUNT")

    #read the infile file...
    for site in summaryReader:
        if start == None:
            start = site.POS
            pScaf = site.CHROM
            pSite = site.POS
            
        if pScaf != site.CHROM:
            printWindow(myCounts, summaryReader.summary, (pSite+start)/2.0, pScaf)
            for t in _t:
                myCounts[t] = []
            start = site.POS
            pScaf = site.CHROM
            
        if site.POS - start >= _w:
            printWindow(myCounts, summaryReader.summary, (site.POS+start)/2.0, pScaf)
            for t in _t:
                myCounts[t] = []
            start = site.POS
                
        for t in _t:
            if t in site.Types:
                myCounts[t].append(1)
            else:
                myCounts[t].append(0)
                
        pSite = site.POS
                
    printWindow(myCounts, summaryReader.summary, (pSite-start)/2.0, pScaf)

def printWindow(counts, summary, midpoint, scaf):
    tot = 0
    for t in _t:
        print("%s\t%s\t%s\t%s" % (summary.Types[t], scaf, midpoint, sum(counts[t])))
        tot += sum(counts[t])
        
    print("ALL\t%s\t%s\t%s" % (scaf, midpoint,tot))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:t:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-t":
            global _t 
            _t.append(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("summary: %s -w %s -t %s \n" % (sys.argv[1], _w, _t))
   
use = "python "+__file__.split("/")[-1]+" summaryFile.txt [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("w - INT - %s - window size" % _w)
    print("t - STR - %s - a type name that you want counted. This option can be used several times." % _w)

    
if __name__ == "__main__":   
    __main__()