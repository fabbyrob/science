doc = """
Takes in a list of regions in each individual that are IBD and a VCF. 
Down samples alleles at each locus to a given sample size. 
Individuals with IBD can only be sampled once.
"""

import sys
import getopt
import summary
import random

_N = 23

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    summaryReader = summary.Reader(open(sys.argv[1],"r"))
        
    regions = open(sys.argv[2],"r") 
    
    #read the regions
    #regions{IND{SCAF{(start, end)}}}
    regionsList = {}
    for ind in summaryReader.summary.Samples:
        regionsList[ind] = {}

    for line in regions:
        line = line.rstrip()
        sline = line.split()
        
        if line.startswith("CHROM"):
            continue
        
        if len(sline) != 4:
            sys.stderr.write("Weird line:\n\t%s\n" % line)
            continue
        
        chrom = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        ind = sline[3]
        
        if chrom in regionsList[ind].keys():
            regionsList[ind][chrom].append((start,end))
        else:
            regionsList[ind][chrom] = [(start, end)]

    header = summaryReader.summary.prettyHeader()
    header = header[0:header.index("DIVERGENCE")+10]
    print(header)

    genoReverse = {v:k for k, v in summaryReader.summary.Genotypes.items()}
    #read the infile file...
    for site in summaryReader:
        ibdInds = getUnsafeInds(site, regionsList)
        #print("%s %s" % (site.POS, ibdInds))
        alleles = []
        for ind in site.Genotypes.keys():
            if ind in ibdInds:
                if site.Genotypes[ind] == genoReverse["homozygote reference"]:
                    alleles.append(site.REF)
                elif site.Genotypes[ind] == genoReverse["homozygote alternate"]:
                    alleles.append(site.ALT)
                elif site.Genotypes[ind] == genoReverse["heterozygote"]:#assume ref for all IBD individuals with data
                    alleles.append(site.REF)
                elif site.Genotypes[ind] == genoReverse["unknown"]:
                    alleles.append(genoReverse["unknown"])
            else:
                if site.Genotypes[ind] == genoReverse["homozygote reference"]:
                    alleles.append(site.REF)
                    alleles.append(site.REF)
                elif site.Genotypes[ind] == genoReverse["homozygote alternate"]:
                    alleles.append(site.ALT)
                    alleles.append(site.ALT)
                elif site.Genotypes[ind] == genoReverse["heterozygote"]:
                    alleles.append(site.REF)
                    alleles.append(site.ALT)
                elif site.Genotypes[ind] == genoReverse["unknown"]:
                    alleles.append(genoReverse["unknown"])
                    alleles.append(genoReverse["unknown"])
                    
        newAlleles = random.sample(alleles, _N)
        ref = newAlleles.count(site.REF)
        alt = newAlleles.count(site.ALT)
        site.REF_NUM = ref
        site.ALT_NUM = alt
        
        if site.TOTAL != 0:#if total was 0 for some reason, like being filtered, keep it that way
            site.TOTAL = ref+alt
        
        site.Genos = []
        
        print(site)        

def getUnsafeInds(site, regions):
    notSafe = []
    for ind, myDict in regions.items():
        if site.CHROM in myDict.keys():
            myRegions = myDict[site.CHROM]
            for start, end in myRegions:
                if site.POS >= start and site.POS <= end:
                    notSafe.append(ind)
                    break
    return notSafe

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"N:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-N":
            global _N 
            _N = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("summary: %s regions: %s -N %s\n" % (sys.argv[1], sys.argv[2], _N))
   
use = "python "+__file__.split("/")[-1]+" SummaryFile.txt Regions.txt [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("N - INT - %s - number of alleles to down sample to." % _N)

    
if __name__ == "__main__":   
    __main__()
