#!/usr/bin/python2.7
import sys
import re
import getopt
from collections import defaultdict
_window = 100
_log = "CG_v_depth.log"
_hist = ""

AD = []
CG = []
scafs = []
bins = defaultdict(list)

def __main__():
    #check aruguments
    if len(sys.argv) < 2:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[2:],"w:h:l")
    except getopt.GetError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _window 
            _window = int(arg)
        elif opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-h":
            global _hist 
            _hist = str(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    #open reference file
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad infile name: "+_log)
        sys.exit()
        
    if (_hist != ""):
        histfile = open(_hist, "w")
        if (histfile == None):
            print("Bad infile name: "+_hist)
            sys.exit()
        
    #print("FILES OPEN")
    print("#Window size - "+str(_window)+"\n#CHROM\tAVERAGE_DEPTH\tCG_CONTENT")
        
    pat = re.compile("(\S+)\s\d+\s\S+\s(\S+)\s(\S+)\s.+DP=(\d+);.+")

    cgCount = 0
    depths = []
    prevscaf = ""
        
    for line in infile:
        line = line.rstrip()
        m = pat.match(line)
        
        if (m != None):
            scaf = m.group(1)
            base = m.group(2)
            depth = int(m.group(4))
            #print("ref = "+str(base)+" depth = "+str(depth))
                
            if(scaf == prevscaf or prevscaf == ""):
                if (len(depths) < _window):
                    depths.append(depth)
                    if (base == "C" or base == "G"):
                        cgCount += 1
                else:
                    calc(depths, cgCount, scaf)
                    depths = []
                    cgCount = 0
            else:
                calc(depths, cgCount, prevscaf)
                depths = []
                cgCount = 0
            prevscaf = scaf
        else:
            logfile.write("Non-standard line:\n\t"+line+"\n")
            
    calc(depths, cgCount, prevscaf) 
    depths = []
    cgCount = 0
    
    if (_hist != ""):
        for k in bins.keys():
            histfile.write(str(k) + ": " + str(len(bins[k]))+"\n")
        
        for x in bins[.1]:
            histfile.write("0.1\t"+x[0]+"\t"+x[1]+"\n")
        for x in bins[.2]:
            histfile.write("0.2\t"+x[0]+"\t"+x[1]+"\n")
        for x in bins[.9]:
            histfile.write("0.9\t"+x[0]+"\t"+x[1]+"\n")
        for x in bins[1]:
            histfile.write("1\t"+x[0]+"\t"+x[1]+"\n")
    
    infile.close()
    logfile.close()
    
def setbin(depth, cgcont):
    if (cgcont < .1):
        bins[0.1].append((depth, cgcont))
    elif (cgcont <.2):
        bins[0.2].append((depth, cgcont))
    elif (cgcont <.3):
        bins[0.3].append((depth, cgcont))
    elif (cgcont <.4):
        bins[0.4].append((depth, cgcont))
    elif (cgcont <.5):
        bins[0.5].append((depth, cgcont))
    elif (cgcont <.6):
        bins[0.6].append((depth, cgcont))
    elif (cgcont <.7):
        bins[0.7].append((depth, cgcont))
    elif (cgcont <.8):
        bins[0.8].append((depth, cgcont))
    elif (cgcont <.9):
        bins[0.9].append((depth, cgcont))
    else:
        bins[1].append((depth, cgcont))
    
def calc(depths, cgCount, scaf):
    if (len(depths) == 0):
        return
    avgDepth = float(sum(depths))/len(depths)
    AD.append(avgDepth)
    cgContent = float(cgCount)/len(depths) 
    CG.append(cgContent)
    scafs.append(scaf)
    print(scaf+"\t"+str(avgDepth)+"\t"+str(cgContent))
    setbin(avgDepth, cgContent)
    
def usage():
    print ("usage: CG_v_depth.py <INFILE> [-w <WINDOW SIZE> -l <LOG FILE NAME>]\n")
    sys.exit()
    
if __name__ == "__main__":   
    __main__()
