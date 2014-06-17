#!/usr/bin/python
import sys
import re
import getopt

from subprocess import *

_log = "afsForRanges.log"
_d = 20
_D = 60
_q = 40
_u = -1

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[2:],"l:d:D:q:u:")
    except getopt.GetError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-d":
            global _d 
            _d = int(arg)
        elif opt == "-D":
            global _D 
            _D = int(arg)
        elif opt == "-q":
            global _q 
            _q = int(arg)
        elif opt == "-u":
            global _u 
            _u = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    sites = open(sys.argv[1],"r")  
        
    if (sites == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    #reg ex
    #1 = scaffold, 2 = start, 3 = end, 4 = model, 5 = +/-
    sitepat = re.compile("(\S+)\s+(\S+)\s+(\S+).+")
    
    #grab our first model to search for
    myRange = getNextRange(sites, sitepat)
    siteName = "afsForRanges.sites"
    vcfName = "afsForRanges.vcf"
    afsName = "afsForRanges.afs"
    
    outfile = open("afsForRanges.csv", "w")
    
    #read the vcf file...
    c = 0
    while(myRange != None):
        #make the site file
        siteFile = open(siteName, "w")
        for i in range(myRange[1][0],myRange[1][1]+1):
            siteFile.write(myRange[0]+"\t"+str(i)+"\n") 
        
        siteFile.close()
        #filter sites
        command = "python /home/robert.williamson/repos/robert.williamson/trunk/pileup_analyzers/filter.py /data/stephen.wright/genomic_SE_runs/final_Cg_fastq/variantcalls_allsampspileup_cgf /data/robert.williamson/temp/"+siteName+" -q 40 -d 260 -D 500 > "+siteName+".flt"
        call(command, shell=True)
        #make the vcf
        command = "python /home/robert.williamson/repos/robert.williamson/trunk/pileup_analyzers/subMpileup.py /data/stephen.wright/genomic_SE_runs/final_Cg_fastq/variantcalls_allsampspileup_cgf /data/robert.williamson/temp/"+siteName+".flt > /data/robert.williamson/temp/"+vcfName
        call(command, shell=True)
        
        #make the afs 
        command = "perl /home/robert.williamson/repos/robert.williamson/trunk/pileup_analyzers/getDepthQual.pl /data/robert.williamson/temp/"+vcfName+" /data/robert.williamson/temp/trash1 /data/robert.williamson/temp/trash2 /data/robert.williamson/temp/trash3 /data/robert.williamson/temp/"+afsName
        call(command, shell=True)

        #compile output
        afs = open("/data/robert.williamson/temp/"+afsName, "r")
        afs.readline()
        outfile.write(myRange[0]+","+str(myRange[1][0])+","+ str(myRange[1][1])+",")
        t = 0
        for line in afs:
            t =1
            data = line.split()
            outfile.write(str(data[1])+",")
            
        if t == 0:
            for i in range(0,27):
                outfile.write("0,")
            
        if(_u > -1 and c > _u):
            break
        
        c += 1
            
        outfile.write("\n")
        afs.close()
        
        myRange = getNextRange(sites, sitepat)
  
def getNextRange(file, pat):
    sites = None
    while (sites == None):
        aLine = file.readline()
        
        if (aLine == ""):
            return
        
        m = pat.match(aLine)
        
        if(m == None):
            continue
        
        scaf = m.group(1)
        sites = [int(m.group(2)), int(m.group(3))]
            
    return (scaf, sites)
        
use = "usage: afsForRanges.py <sites> [-d <DEPTH> -D <DEPTH> -q <QUALITY> -l <LOG FILE NAME>]\n"
        
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("<SITEs> - a list of sites")
    print ("q - the minimum quality to include a site in the analysis (40)")
    print ("d - the minimum depth to include a site in the analysis (20)")
    print ("D - the maximum depth to include a site in the analysis (60)")
    print ("l - the name of the log file (vcfToFastq.log)")
    
if __name__ == "__main__":   
    __main__()
