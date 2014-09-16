doc = """
Takes a VCFSummary and an annotated reference and outputs a new summary with the GENE and DIR columns added.
"""

import sys
import getopt
import summary

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    summaryReader = summary.Reader(open(sys.argv[1],"r"))
    sites = open(sys.argv[2],"r")  
    
    summaryReader.addGenes()
    print(summaryReader.Header)
    
    next_site = getNextSite(sites)
    #read the infile file...
    for site in summaryReader:
        if next_site == (None, None, None, None, None):#print the remaining sites
            site.GENE = "N"
            site.DIR = "0"
            print(site)
            continue
    
        while next_site[0] != None and site.CHROM > next_site[0]:#we missed this site
            sys.stderr.write("Missed a site (CHROM): %s %s\n" % (next_site[0], next_site[1]))
            next_site = getNextSite(sites)
            
        if next_site == (None, None, None, None, None):
            site.GENE = "N"
            site.DIR = "0"
            print(site)
            continue
        
        if site.CHROM < next_site[0]:#we havent found the right scaf yet, go to the next site. 
            site.GENE = "N"
            site.DIR = "0"
            print(site)
            continue
        
        #on the right scaf, now deal with pos info
        
        while next_site[1] != None and site.POS > next_site[1]:#we missed a site WARNING, this while loop could bring us to the next chromosome....
            sys.stderr.write("Missed a site (POS): "+str(next_site)+"\n")
            next_site = getNextSite(sites)
            
        if next_site[1] > site.POS:#haven't reached the site yet
            site.GENE = "N"
            site.DIR = "0"
            print(site)
            continue
        
        if next_site == (None, None, None, None, None):
            site.GENE = "N"
            site.DIR = "0"
            print(site)
            continue
        
        if next_site[1] == site.POS:#at the site
            site.GENE = next_site[3]
            site.DIR = next_site[4]
            print(site)
            next_site = getNextSite(sites)
            
def getNextSite(sitefile):
    line = sitefile.readline()
    if line == "":
        return (None, None, None, None, None)
    
    line = line.rstrip()
    
    while line[0] == "{":#dictionary line
        return getNextSite(sitefile)
    
    sline = line.split()
    if len(sline) < 6:
        if _v:
            sys.stderr.write("Line not formatted correctly\n\t"+str(line)+"\n")
        return getNextSite(sitefile)
    
    scaf = sline[0]
    pos = int(sline[1])
    gene = sline[3]
    dir = sline[4]
    type = sline[5]
    return (scaf, pos, type, gene, dir)

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

    sys.stderr.write("summary: %s annotation: %s \n" % (sys.argv[1], sys.argv[2]))
   
use = "python "+__file__.split("/")[-1]+" Summary Annotation [OPTIONS]"
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