doc = """
Reads in a vcf and calculates heterozygosity in windows from raw data. No filtering is applied.

OUTPUT:
CHROM\tPOS\tNum_Hets\tNum_Genos
scaffold_1    50000    576    28999
scaffold_1    50000    123    1000
scaffold_2    50000    543    2333

columns are the midpoint of the window, the number of heterozygous sites, and the total number of sites with called genotypes.
"""

import sys
import getopt
import vcf

_w = 100000

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    reader = vcf.Reader(open(sys.argv[1], 'rb'))

    #read the infile file...
    print("CHROM\tPOS\tNum_Hets\tNum_Genos")
    window = []
    start = 0
    end = 0
    scaf = ""
    for record in reader:
        if not scaf:
            scaf = record.CHROM
            start = record.POS
        
        if scaf != record.CHROM:
            processWindow(scaf, (start+end)/2.0, window)
            window = []
            start = record.POS
            scaf = record.CHROM
        elif len(window) == _w:
            processWindow(scaf, (start+end)/2.0, window)
            window = []
            start = record.POS
            
        hets = 0
        samps = 0
        for s in record.samples:
            if 'GT' in s.keys():
                if '0/1' == s['GT']:
                    hets += 1
                    samps += 1
                else:
                    samps += 1
        
        window.append((hets, samps))
        end = record.POS
            
    processWindow(scaf, (start+end)/2.0, window)

def processWindow(scaf, pos, window):
    numhets = sum([x for x,y in window])
    numsamps = sum([y for x,y in window])
    print("%s\t%s\t%s\t%s" % (scaf, pos, numhets, numsamps))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -w %s\n" % (sys.argv[1], _w))
   
use = "python "+__file__.split("/")[-1]+" VCF [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("w - INT - %s - the window size to use" % _w)

    
if __name__ == "__main__":   
    __main__()