doc = """
Takes in a summary, for each continuous stretch of CNCs that are also intergenic, the distance up and downstream to the nearest gene is caclulated.

Also, for each region divergence and pi are calculated.

output:

CHROM CNC_START CNC_END CNC_DIV CNC_PI UP_DIST DOWN_DIS
scaffold_1 100 110 0.1 0.15 10000 23141
"""

import sys
import summary
import getopt

_N = 26
_t = []

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    if _N % 2 == 0:
        folded_size = int(_N/2+1)
    else:
        folded_size = int(ceil(float(_N)/2))
    
    summaryReader = summary.Reader(open(sys.argv[1],"rb"))
    
    #read the infile file...
    lastGeneSite = None# only count exon, 4, and 0 fold sites
    currentCNCs = []# list of tuples (start, stop, div, pi)
    cncStart = None#start of current stretch of CNCs
    cncDiv = 0
    cncDivTot = 0
    cncWindow = []
    pScaf = None
    pSite = None
    
    geneTypes = [summaryReader.summary.typeToCode['exon'],summaryReader.summary.typeToCode['4fold'],summaryReader.summary.typeToCode['0fold']]
    cncTypes = [summaryReader.summary.typeToCode['cnc']]
    
    for t in _t:
        cncTypes.append[summaryReader.summary.typeToCode[t]]
    
    print("CHROM\tCNC_START\tCNC_END\tCNC_DIV\tCNC_PI\tUP_DIST\tDOWN_DIST")
    
    for site in summaryReader:
        if pScaf == None:
            pScaf = site.CHROM
            
        if pScaf != site.CHROM:#on next scaf, clean up
            if cncStart != None:#end CNC at last site
                if cncDivTot:
                    div = float(cncDiv)/cncDivTot
                else:
                    div = "NA"
                pi = calcPi(cncWindow)
                currentCNCs.append((cncStart, pSite, div, pi))
                outputCNCs(pScaf, currentCNCs, lastGeneSite, None)
            cncStart = None
            cncWindow = []
            cncDiv = 0
            cncDivTot = 0
            currentCNCs = []
            pScaf = site.CHROM
            
        if isCNC(cncTypes, site.Types):
            if cncStart == None:
                cncStart = site.POS
                
            if site.REF_NUM + site.ALT_NUM != _N or site.TOTAL != _N:
                sys.stderr.write("Bad data at site %s %s, ignoring it.\n" % (site.CHROM, site.POS))
            
            cncWindow.append(site.REF_NUM)
            
            if site.DIVERGENCE == 1:
                cncDiv += 1
                cncDivTot += 1
            elif site.DIVERGENCE == 0:
                cncDivTot += 1
        elif isGene(geneTypes, site.Types):
            if cncStart != None:
                if cncDivTot:
                    div = float(cncDiv)/cncDivTot
                else:
                    div = "NA"
                pi = calcPi(cncWindow)
                currentCNCs.append((cncStart, pSite, div, pi))
                cncStart = None
                cncWindow = []
                cncDiv = 0
                cncDivTot = 0
                
            if currentCNCs:
                outputCNCs(pScaf, currentCNCs, lastGeneSite, site.POS)
                currentCNCs = []
                
            lastGeneSite = site.POS
        else:
            if cncStart != None:
                if cncDivTot:
                    div = float(cncDiv)/cncDivTot
                else:
                    div = "NA"
                    
                pi = calcPi(cncWindow)
                currentCNCs.append((cncStart, pSite, div, pi))
                cncStart = None
                cncWindow = []
                cncDiv = 0
                cncDivTot = 0
            
        pSite = site.POS

def isGene(geneTypes, siteTypes):
    for gt in geneTypes:
        for st in siteTypes:
            if gt == st:
                return True
    return False

def isCNC(cncTypes, siteTypes):
    for ct in cncTypes:
        if ct not in siteTypes:
            return False
    return True
            
def calcPi(window):
    total = 0
    numSites = 0.0
    for frequency in window:
        if frequency == -1:#no data
            continue
        
        minorfreq = float(frequency)/_N
        majorfreq = 1-minorfreq
        
        total += 2*minorfreq*majorfreq
        numSites += 1.0
        
    total *= _N/(_N-1.0)
    if numSites == 0:
        return "NA"
    return (total/numSites)

def outputCNCs(scaf, cncs, prevGene, nextGene):
    for start, stop, div, pi in cncs:
        if prevGene:
            upDist = start-prevGene
        else:
            upDist = "NA"
            
        if nextGene:
            downDist = nextGene-stop
        else:
            downDist = "NA"
            
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (scaf, start, stop, div, pi, upDist, downDist))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"N:t:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-t":
            global _t 
            _t.append(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("summary: %s -N %s -t %s\n" % (sys.argv[1], _N, _t))
   
use = "python "+__file__.split("/")[-1]+" summaryFile.txt"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("N - INT - %s - number of chromosomes in sample" % _N)
    print("t - str - %s - a type that must overlap with CNCs for the CNC to be counted. Can be used multiple times. (e.g. if you wanted intergenic CNCs only add -t intergene)" % _N)

    
if __name__ == "__main__":   
    __main__()