import sys
import re
import getopt
import vcf
from math import floor 

_log = __file__.split("/")[-1]+".log"#log file name
_N = 26
_v = False

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[5:],"l:N:v")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-v":
            global _v 
            _v = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    if _v: sys.stderr.write("-N = " +str(_N)+"\n")
    
    vcf_reader = vcf.Reader(open(sys.argv[1],"rb"))
    
    depthFile = open(sys.argv[2],"w")  
    afsFile = open(sys.argv[3],"w") 
    qualFile = open(sys.argv[4],"w")     
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    #reg ex
    #vcf_pat = re.compile("\w+\s+\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+AF=([^;]+);\S+\s+(\S+)\s+(.+)")
    
    quals = {}
    afs = {}
    depths = {}
    
    for i in range(0,_N+1):
        afs[i] = 0
    
    if _v: sys.stderr.write("afs: %s\n" % afs)
    
    #read the vcf file...
    ctr = 0
    for record in vcf_reader:
        if ctr % 5000 == 0:
            sys.stderr.write("Processed up in here yo to %s lines...\n" % ctr)
        ctr += 1
        
        qual = int(record.QUAL)#sets the max qual to 200, to limit the file size
        if qual > 200:
            qual = 200
        af = 0
        
        if not quals.has_key(qual):
            quals[qual] = 1
        else:
            quals[qual] += 1
            
        
        for samp in record.samples:
            if 'DP' in samp.keys():
                depth = samp['DP'][0]
                
                if not depths.has_key(depth):
                    depths[depth] = 1
                else:
                    depths[depth] += 1
            
            if 'GT' in samp.keys():
                af += samp['GT'].count("1")
                    
        afs[af] += 1
    
    if _v: sys.stderr.write("Finished reading in all data.\nSorting depth.\n")
    
    for k in sorted(depths.keys()):
        depthFile.write(str(k)+"\t"+str(depths[k])+"\n")
        
    if _v: sys.stderr.write("Finished outputting depth.\nSorting afs.\n")
    for k in sorted(afs.keys()):
        afsFile.write(str(k)+"\t"+str(afs[k])+"\n")
        
    if _v: sys.stderr.write("Finished outputting afs.\nSorting quals.\n")
    for k in sorted(quals.keys()):
        qualFile.write(str(k)+"\t"+str(quals[k])+"\n")
        
    if _v: sys.stderr.write("Finished outputting quals.\n")

   
use = "python "+__file__.split("/")[-1]+" <VCF> <DEPTH FILE> <AFS FILE> <QUAL FILE> [-N <Sample size>]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()