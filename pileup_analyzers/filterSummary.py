import sys
import getopt
import summary

_c = 0.3#minimum fraction of sites that must have good data in a window
_w = 20000#size of window to test coverage
_N = 26 #sample size
_o = "" #output regions mode


_coding = ['4fold', '0fold', 'exon', 'stop']
_ncoding = ['intergene', 'intron', '3utr', '5utr', 'cnc']

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    args = {"fraction coding":_c, "window size":_w, "sample size":_N, "region output file":_o}
    
    argStr = ""
    for arg in args.keys():
        argStr += ("#%s\t%s\n" % (arg, args[arg]))
    sys.stderr.write(argStr)
    
    myReader = summary.Reader(open(sys.argv[1],"rb"))

    codingN = []
    for t in _coding:
        if t in myReader.TypeToCode.keys():
            codingN.append(myReader.TypeToCode[t])
    codingN = set(codingN)
    
    ncodingN = []
    for t in _ncoding:
        if t in myReader.TypeToCode.keys():
            ncodingN.append(myReader.TypeToCode[t])
    ncodingN = set(ncodingN)

    if _o: 
        filteredRegions = open(_o, "w")
    else:
        filteredRegions = None

    window = []
    safe_sites = 0.0

    print(myReader.Header)
    pScaf = ""

    #read the infile file...
    for site in myReader:
        if pScaf == "":
            site.CHROM
        
        if (len(window) == _w or pScaf != site.CHROM) and window:
            processWindow(window, safe_sites, filteredRegions)
            window = []
            safe_sites = 0.0
        
        #check if all individuals are represented
        if site.TOTAL == _N:
            safe_sites += 1
           
        #check if a site is both a coding site and a noncoding site
        #if it is, make the site not safe, set type to unknown
        safe = True
        if len(set(site.Types) & codingN) and len(set(site.Types) & ncodingN):
            safe = False
        
        if not safe:
            site.Types = [myReader.TypeToCode['unknown']]
                        
        window.append(site)
    
        pScaf = site.CHROM
        
    processWindow(window, safe_sites, filteredRegions)#get the last window
    if _o: filteredRegions.close()
        
def processWindow(window, safe_sites, filteredRegions):
    safe = True
    if safe_sites/len(window) < _c:
        safe = False
    
    scaf = None
    start = None
    end = None
        
    for site in window:
        if not safe:
            site.TOTAL = 0
        
        if not scaf:
            scaf = site.CHROM
            start = site.POS
        end = site.POS
        
        print(site)
        
    if _o and not safe:
        filteredRegions.write("%s\t%s\t%s\n" % (scaf, start, end))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"c:w:N:o:")
    except getopt.GetoptError:
        usage()
    print(sys.argv[num:])
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-c":
            global _c 
            _c = float(arg)
        elif opt == "-N":
            global _N 
            _N= int(arg)
        elif opt == "-o":
            global _o 
            _o = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" INFILE [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("Takes in a summary file and filters out sites in windows where the coverage is too low. Or sites where coding and noncoding site types overlap.")
    print("option - type - default - description\n****************************")
    print("c - FLOAT - %0.2f - the minimum fraction of sites in a window that must have good data." % _c)
    print("w - INT - %d - the size of the window to analyze coverage in." % _w)
    print("N - INT - %d - sample size of your data (# of alleles)." % _N)
    print("o - STR - %s - optional output file to put the list of regions that are filtered out." % _o)

    
if __name__ == "__main__":   
    __main__()
