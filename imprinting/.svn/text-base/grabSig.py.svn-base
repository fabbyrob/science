#!/usr/bin/python
import sys
import re
import getopt

_log = __file__.split("/")[-1]+".log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[2:],"l:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()

        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    heliSig = []
    mudrSig = []
    
    #read the vcf file...
    for line in infile:
        line=line.rstrip()
        line = line.split(",")

        if len(line) > 7:
            if line[-1] == "1":
                mudrSig.append((line[0],line[-4]))
            elif line[-3] == "1":
                heliSig.append((line[0],line[2]))

    hfile = open("helis.csv","w")
    mfile = open("mudrs.csv","w")
    
    for i in heliSig:
        hfile.write(i[0]+","+i[1]+"\n")
        
    for i in mudrSig:
        mfile.write(i[0]+","+i[1]+"\n")

    print(heliSig)
    print(mudrSig)
   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()
