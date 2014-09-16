#!/usr/bin/python
#Robert Williamson
#U of T
#2011

#This program takes in one text file with a series of command line commands and executes each of them using the OS.
#There are 2 options for this program.
#-d : Debug mode. This prints out information at the beginning of each execution to let you know what command is being executed.

#example input file included in /data/robert.williamson/examples/wrappertest.txt
#Program will ignore lines in the input file starting with a pound '#' symbol.
#All paths used in the input file must be ABSOLUTE, no matter where you call the program from.

import sys
import getopt
import re

_debug = False
_exit = True

def __main__():
    if len(sys.argv) < 2:
        usage()
        
    infile = open(sys.argv[1],"r")
    
    try: 
        opts, args = getopt.getopt(sys.argv[2:],"d")
    except getopt.GetError:
        usage()
    
    for opt, arg in opts:
        if opt == "-d":
            global _debug 
            _debug = True
        else:
            print ("***Unrecognized option: "+opt+"\n")
            usage()
    
    if(infile == None):
        print("***Bad infile name: "+sys.argv[1])
        sys.exit()
        
    pat = re.compile("#")
        
    if _debug:
        print("***File open, beginning execution.")
        
    for line in infile:
        
        m = pat.match(line)
        
        if (m != None):
            if _debug:
                print("***Comment encountered:\n\t" + line + "\n***Skipping to next line.")
            continue
        
        if _debug:
            print("***Now executing line:\n\t"+ line)
        try:
            retcode = call(line, shell=True)
            if retcode < 0 and _debug:
                print ("***Child was terminated by signal " + str(retcode))
            elif _debug:
                print ("***Child returned " + str(retcode))
        except OSError, e:
            print ("***Execution failed: " + e)
            if _exit:
                sys.exit()
            
    if _debug:
        print("***Commands completed.")
        
def usage():
    print("use: wrapper.py <INFILE> [-d]")
    sys.exit()
    
if __name__ == "__main__":
    __main__()
