import sys
from subprocess import Popen
import getopt
import time

_p = False
_e = True
_w = 0.0

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
            
    #read in the loop file
    loopfile = open(sys.argv[1])
    
    program = ""
    folder = ""
    maxProcs = 10
    replicates = 1
    
    lineEnding = ".txt"
    
    params = []
    
    runs = 1
    
    for line in loopfile:
        line = line.split("#")[0]#removes everything after a #
        line = line.rstrip()
        
        sline = line.split()
        if len(sline) < 1:#blank line
            continue
        
        if sline[0] == "program":
            if len(sline) < 2:
                sys.stderr.write("No program path provided on \"program\" line.\n")
                sys.exit()
            program = " ".join(sline[1:])
        elif sline[0] == "output":
            if len(sline) < 2:
                sys.stderr.write("No output folder path provided on \"output\" line.\n")
                sys.exit()
            folder = sline[1]
            
            if folder[-1] != "/":#make sure folder path ends in a /
                folder += "/"
        elif sline[0] == "processes":
            if len(sline) < 2:
                sys.stderr.write("No max processes provided on \"processes\" line.\n")
                sys.exit()
                
            maxProcs = int(sline[1])
        elif sline[0] == "replicate":
            if len(sline) < 2:
                sys.stderr.write("No number of replicates provided on \"replicate\" line.\n")
                sys.exit()
                
            replicates = int(sline[1])
            runs *= replicates
        else:#it must be a parameter
            param = sline[0]
            vals = sline[1:]
            
            sys.stderr.write("For parameter %s (parameter number %s) will use input values: %s\n" % (param, len(params), vals))
            
            if len(vals) > 1 and vals[0] == "OUTPUT":
                lineEnding = vals[1]
                vals = ["OUTPUT"]
            else:
                runs *= len(vals)
                
            params.append((param, vals))
    
    sys.stderr.write("Will run a total of %s instances of %s.\n" % (runs, program))
    
    cmds = buildCmd(program, params, "", folder, lineEnding, replicates)
    
    if _p:#only print out the commands, dont run them
        print("\n".join(cmds))
    else:
        procs = []
        for cmd in cmds:
            time.sleep(_w)
            print(cmd)
            try:
                p = Popen(cmd, shell=True)
                procs.append(p)
            except OSError, e:
                print ("***Execution failed: " + e + "\t"+cmd)
            
            while len(procs) >= maxProcs:
                procs[0].wait()
                #print("one done")
                procs.pop(0)
                
        while len(procs) > 0:
                procs[0].wait()
                #print("one done")
                procs.pop(0)

def buildCmd(cmd, params, fileName, outFolder, lineEnding, replicates):
    if not params:#if we're out of parameters
        #then make a command for each replicate
        final_cmds = []
        fileName = outFolder+fileName
        for i in range(replicates):
            mycmd = cmd.replace("OUTPUT", fileName+"replicate"+str(i)+lineEnding)
            
            if _e:
                mycmd += " 2> "+fileName+"replicate"+str(i)+".err"
            
            #print(mycmd)
            final_cmds.append(mycmd)
            
        return final_cmds
    else:  
        paramName = params[0][0].replace("-","")#removes -'s from param names, for readability
        final_cmds = []
        for val in params[0][1]:
            if paramName == "?" or paramName == "":
                param = ""
                paramName = ""
            else:
                param = params[0][0]
            
            
            mycmd = cmd + " %s %s" % (param, val)
            if val == "OUTPUT":#dont modify the filename for the OUTPUT parameter...
                myfileName = fileName
            else:
                myfileName = fileName + "%s%s_" % (paramName, val)
            final_cmds += buildCmd(mycmd, params[1:], myfileName, outFolder, lineEnding, replicates)
        return final_cmds

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:pe")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-p":
            global _p 
            _p = True
        elif opt == "-e":
            global _e 
            _e = False
        elif opt == "-w":
            global _w 
            _w = float(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -p %s\n" % (sys.argv[1], _p))

use = "python "+__file__.split("/")[-1]+" LOOP_FILE.txt [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print("Reads in a loop file, which specifies a program to run, and what combinations of parameters to run it for. Then executes each instance of that program.")
    print (use)
    
    print("OPTION - ARG - DEFAULT - Description")
    print("p - NONE - %s - If this flag is set, then commands will be printed out but not executed. Use this to test your input file and check its formatting." % _p)
    print("e - NONE - %s - If this flag is set, then error files will not automatically be made for every instance." % _e)
    print("w - FLOAT - %s - The amount of time specified by this flag is waited between each run (in seconds)." % _w)
    
    print("\n-----------------\nAn example loop file:")
    file = "program /opt/python27/bin/python27 ~/bin/TEsim.py #program to run\noutput /data/robert.williamson/temp/loopTest/ #output folder\nprocesses 2#maximum number of processes\nreplicate 2#number of replicates\n#will use parameters in order\n#if \"?\" is provided as param name, it will not use a flag on the parameter, jut put it directly in\n#i.e if the first line was \"? 2\"\n#the program here would always execute: TEsim.py 2 -T ... -N ... -G ... > OUTPUT\n#if the first line was: ? OUTPUT .txt\n#then the first thing in the line would be the aggregated output filename (see below)\n#i.e. TEsim.py T2_N10_G5.txt -T 2 -N 10 -G 5 (assuming the last line in this file did not exist)\n-T 2 5#list of parameters followed by all options\n-N 10 100\n-G 5 10\n> OUTPUT .txt#if the argument is OUTPUT it will make an output filename by concatenating all the options for this run together, it will give the file ending given after OUTPUT\n"

    print(file)
if __name__ == "__main__":   
    __main__()