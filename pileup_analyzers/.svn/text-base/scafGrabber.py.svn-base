import sys
from subprocess import Popen

samps = [11, 12, 13, 15, 16, 168]
samps2 = [172,173,171,182,20,22,30,36,38,40,52,53,59,63,64,71,73,77,79,81,82,83]
scafs = [8,7,6,5,4,3,2,1]

#subset sites
procs = []
procLimit = 4

def waitAll(procs):
    while len(procs) > 0:
        procs[0].wait()
        #print("one done")
        procs.pop(0)
        
#for s in samps2:
#    cmd = "/data/apps/align/samtools-0.1.18/samtools view -b ~/rubella_samps/aligns/"+str(s)+".bam scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 > "+str(s)+".sub.bam"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#    while len(procs) > procLimit:
#        procs[0].wait()
#        #print("one done")
#        procs.pop(0)
#waitAll(procs)

#for s in samps2:
#    cmd = "/data/apps/align/samtools-0.1.18/samtools view -Sb ~/rubella_samps/aligns/dan.sample."+str(s)+".sam > "+str(s)+".bam"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#    while len(procs) > procLimit:
#        procs[0].wait()
#        #print("one done")
#        procs.pop(0)
#waitAll(procs)

for s in samps2:
    cmd = "java -Xmx2g -jar /data/apps/conversion_tools/picard-tools-1.67/BuildBamIndex.jar \
    INPUT="+str(s)+".bam  \
    2> buildIndex."+str(s)+".err"
    try:
        p = Popen(cmd, shell=True)
        procs.append(p)
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)
            
    while len(procs) > procLimit:
        procs[0].wait()
        #print("one done")
        procs.pop(0)
waitAll(procs)