#!/usr/bin/python
import sys
import re
import getopt

_log = "filter.log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    input = open(sys.argv[1],"r")
    
    if (input == None):
        print("Bad input name: "+sys.argv[1])
        sys.exit()
    
    #reg ex
    pat = re.compile("(\S+);(\S+);(\S+)")
    bads = ['Low_complexity','Simple_repeat', 'rRNA', 'snRNA', 'Satellite', 'Other/Composite', 'Unknown', 'buffer']
    
    
    #read the vcf file...
    for line in input:
        line = line.rstrip()
        m = pat.match(line)
        
        if (m == None):# if we find a weird line...
            continue
        else:
            fam = m.group(1)
            start = m.group(2)
            end = m.group(3)
            
            skip = False
            for b in bads:
                if b in fam:
                    skip = True
                    
            if skip:
                continue
            
            print (fam+","+start+","+end)
            
use = "python "+__file__.split("/")[-1]+"<TE locations>"

def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("<TE locations> - file listing TE locations and types")

    
if __name__ == "__main__":   
    __main__()
