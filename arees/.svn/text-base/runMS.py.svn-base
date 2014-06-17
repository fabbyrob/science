import sys
from subprocess import call, Popen

procs = []
procLimit = 10

infile = open(sys.argv[1],"r")
outname = sys.argv[2]

def waitAll(procs):
    while len(procs) > 0:
        procs[0].wait()
        #print("one done")
        procs.pop(0)
        
ctr = 0
for line in infile:
    sline = line.split()
    theta = sline[0]
    length = sline[1]
    cmd = "/data/apps/msdir/ms 27 1 -t "+theta+" -I 2 1 26 -ej 4.3 1 2 "
    cmd2 = "/bin/echo \"length: "+length+"\" >> "+outname
    if ctr == 0:
        cmd += "> "+outname
    else:
        cmd += ">> "+outname
    try:
        p = Popen(cmd, shell=True)
        procs.append(p)
        p.wait()
        p = Popen(cmd2, shell=True)
        p.wait()
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)
    
        
    ctr += 1
            
while len(procs) > procLimit:
    procs[0].wait()
    #print("one done")
    procs.pop(0)