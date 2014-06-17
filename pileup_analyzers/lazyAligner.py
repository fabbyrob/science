import sys
from subprocess import Popen

samps = [11, 12, 13, 15, 16, 168,172,173,171,182,20,22,30,36,38,40,52,53,59,63,64,71,73,77,79,81,82,83]
#REDO 171
#REDO 171
#REDO 171
#REDO 171
scafs = [8,7,6,5,4,3,2,1]

#subset sites
procs = []
procLimit = 10

def waitAll(procs):
    while len(procs) > 0:
        procs[0].wait()
        #print("one done")
        procs.pop(0)
        
#for s in samps:
#    cmd = "/data/robert.williamson/python27/Python-2.7.3/python /data/wei.wang/apps/stampy-1.0.19/stampy.py -g /data/robert.williamson/rubella_samps/aligns/rubella \
#    -h /data/robert.williamson/rubella_samps/aligns/rubella \
#    -t 8 -M /data/arvid.agren/dan.rubella/dan.sample."+str(s)+".1.txt /data/arvid.agren/dan.rubella/dan.sample."+str(s)+".2.txt \
#    > "+str(s)+".sam 2> stampy."+str(s)+".err"
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

for s in samps:
    cmd = "java -Xmx2g -jar /data/apps/conversion_tools/picard-tools-1.67/ReorderSam.jar \
    INPUT="+str(s)+".sam \
    OUTPUT="+str(s)+".reordered.bam  \
    REFERENCE=/data/robert.williamson/Capsella_rubella_v1.0_combined.fasta 2> reorder."+str(s)+".err"
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
#
for s in samps:
    cmd = "java -Xmx2g -jar /data/apps/conversion_tools/picard-tools-1.67/AddOrReplaceReadGroups.jar \
    INPUT="+str(s)+".reordered.bam \
    OUTPUT="+str(s)+".readgroup.bam  \
    RGLB=rubella \
    RGPL=illumina \
    RGPU=1234 \
    RGSM=sample_"+str(s)+" \
    2> readgroups."+str(s)+".err"
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
##
for s in samps:
    cmd = "java -Xmx2g -jar /data/apps/conversion_tools/picard-tools-1.67/SortSam.jar \
    INPUT="+str(s)+".readgroup.bam  \
    OUTPUT="+str(s)+".sorted.bam  \
    SORT_ORDER=coordinate 2> sort."+str(s)+".err"
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
#
for s in samps:
    cmd = "java -Xmx2g -jar /data/apps/conversion_tools/picard-tools-1.67/BuildBamIndex.jar \
    INPUT="+str(s)+".sorted.bam  \
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
#
for s in samps:
    cmd = "java -Xmx2g -jar /data/apps/align/GenomeAnalysisTK-1.6-5-g557da77/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /data/robert.williamson/Capsella_rubella_v1.0_combined.fasta \
    -o "+str(s)+".intervals \
    -I "+str(s)+".sorted.bam \
    2> target."+str(s)+".err"
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
#
for s in samps:
    cmd = "java -Xmx2g -jar /data/apps/align/GenomeAnalysisTK-1.6-5-g557da77/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /data/robert.williamson/Capsella_rubella_v1.0_combined.fasta \
    -targetIntervals "+str(s)+".intervals \
    -I "+str(s)+".sorted.bam \
    -o "+str(s)+".realigned.bam \
    2> realigner."+str(s)+".err"
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