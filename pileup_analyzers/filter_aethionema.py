#!/usr/bin/python
import sys
import re
import getopt
import vcf
#from Queue import PriorityQueue

_q = 40#min quality
_d = 20#min depth
_D = 60#max depth
_l = 60#genotype quality cutoff

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"q:d:D:l:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-q":
            global _q
            _q = int(arg)
        elif opt == "-d":
            global _d
            _d = int(arg)
        elif opt == "-D":
            global _D
            _D = int(arg)
        elif opt == "-l":
            global _l 
            _l = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    reader = vcf.Reader(open(sys.argv[1],"rb"))
    sites = open(sys.argv[2], "r")
    
    next_site = getNewSite(sites)
    
    prev = None
    want = False
    for record in reader:
        if next_site == None:
            sys.stderr.write("Ran out of sites.")
            break
        
        site = (record.CHROM, record.POS)

        #check indel
        if prev != None and prev.CHROM == record.CHROM and prev.POS == record.POS:#INDEL
            #do not print this site ,or the previous one
            sys.stderr.write("Indel site: %s %s\n" % site)
            want = False
            continue
            
        elif prev != None and want:
            print("%s\t%s\t%s" % (prev.CHROM, prev.POS, prev.AF))#previous site was safe, use it
            want = False
          
        #check position  
        if record.CHROM < next_site[0]:#aren't on the right scaf yet
            continue
        
        if record.POS < next_site[1]:#aren't at the right pos yet
            continue
        
                  
        #make sure we haven't passed up the last site we were looking for
        while next_site != None and record.CHROM > next_site[0]:
            sys.stderr.write("Missed site: %s %s\n" % next_site)
            next_site = getNewSite(sites)#we missed a chromosome, grab a new site
            
        while next_site != None and record.POS > next_site[1]:
            sys.stderr.write("Missed site: %s %s\n" % next_site)
            next_site = getNewSite(sites)#we missed a position, grab a new site
        
        if site == next_site:
            #check each filter
            safe = True
            if record.QUAL < _q:
                sys.stderr.write("Quality too low: %s %s \n" % site)
                safe = False
            
            if "DP" not in record.INFO.keys() or record.INFO["DP"] < _d or record.INFO["DP"] > _D:
                sys.stderr.write("Depth missing or outside range: %s %s \n" % site)
                safe = False
        
            if record.ALT[0] == ".":
                minorFreq = 0
            else:
                minorFreq = 0
                for sample in record.samples:
                    if "GQ" in sample.keys() and sample["GQ"] >= _l:
                        minorFreq += sample["GT"].count("1")
                    else:
                        safe = False
                        sys.stderr.write("At least one sample's genotype is bad: %s %s\n" % site)
                        break
                
            if safe:
                want = True
                record.AF = minorFreq
                prev = record
                
            next_site = getNewSite(sites)
            continue

    if prev != None and want:
        print("%s\t%s\t%s" % (prev.CHROM, prev.POS, pref.AF))#previous site was safe, use it
        
    
#grabs a new site too look for from the sites input file
def getNewSite(sites):
    next_site = None
    while(next_site == None):
        new_site = sites.readline()
        
        if (new_site == ""):
            return None
        
        sline = new_site.split()
        
        if len(sline) < 2:
            sys.stderr.write("Non-standard line (SITE FILE):\n\t"+new_site+"\n")
            continue
        next_site = (sline[0], int(sline[1]))

    return next_site
   
def usage():
    print ("usage: filter.py <VCF> <SITE_LIST> [-q <QUALITY CUTOFF> -d <DEPTH MIN> -D <DEPTH MAX> -L <LIKELIHOOD CUTOFF>]\n")
    sys.exit()
    
def details():
    print ("usage: filter.py <VCF> <SITE_LIST> [-q <QUALITY CUTOFF> -d <DEPTH MIN> -D <DEPTH MAX> -G -l <LOG FILE NAME>]\n")
    print ("<VCF> - a variant call format file")
    print ("<SITE_LIST> - a file containing a list of sites of interest. Contains 2 columns, Scaffold Name and Position, each site on its own line.")
    print ("q - the minimum quality to include a site in the analysis (40)")
    print ("d - the minimum depth to include a site in the analysis (20)")
    print ("D - the maximum depth to include a site in the analysis (60)")
    print ("L - the minimum genotype likelihood score (60)")
    
if __name__ == "__main__":   
    __main__()
