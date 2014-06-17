#!/usr/bin/python
import sys
import re
import getopt


_log = "vcfToFastq.log"
_w = 5

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"l:w:")
    except getopt.GetError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-w":
            global _w 
            _w = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    vcf = open(sys.argv[2],"r")
    
    if (vcf == None):
        print("Bad vcf name: "+sys.argv[2])
        sys.exit()
        
    sites = open(sys.argv[1],"r")  
        
    if (sites == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    #reg ex
    #1 = scaf, 2 = pos
    indel = re.compile("\w+_(\d+)\s(\S+).+INDEL.+")
    #1 = scaf, 2 = pos, 3 = ref, 4 = alt, 5 = qual, 6 = depth, 7 = inds
    vcf_pat = re.compile("\w+_(\d+)\s(\d+)\s\S+\s(\S+)\s(\S+)\s(\S+)\s\S+\s.*DP=(\d+);\S+\s\S+\s(.+)")
    #1 = scaffold, 2 = start, 3 = end, 4 = model, 5 = +/-
    sitepat = re.compile(">\S+\s(\S+)\s\S+\s\S+\s(\S+)\s(\S+)\s(\S+)\s(\S+).+")
    
    #grab our first model to search for
    (scaf, region) = getNextLocus(sites, sitepat)
    
    #read the vcf file...
    for line in vcf: 
        indelCount += 1
  
def getNextLocus(file, pat):
    mod = None
    while (mod == None):
        aLine = file.readline()
        
        if (aLine == ""):
            return
        
        m = pat.match(aLine)
        
        if(m == None):
            continue
        modP = mod
        
        scaf = int(m.group(1))
        sites = [int(m.group(2)),int(m.group(3)),int(m.group(4)),int(m.group(5))]
            
    return (scaf, sites)
        
use = "usage: vcfToFastq.py <VCF> <SITE LIST> [-w <WINDOW SIZE> -l <LOG FILE NAME>]\n"
        
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("<VCF> - a vcf file")
    print ("<SITE LIST> - a list of sites")
    print ("w - number of bases to ignore around INDELs (5)")
    
if __name__ == "__main__":   
    __main__()
