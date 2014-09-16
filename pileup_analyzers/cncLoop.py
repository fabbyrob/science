import sys
from subprocess import call, Popen

sites = ['3utr','5utr','ambig','down','up','intergene','intron','snc']

scafs = [8,7,6,5,4,3,2,1]
procs = []
procLimit = 10

def waitAll(procs):
    while len(procs) > 0:
        procs[0].wait()
        #print("one done")
        procs.pop(0)
            

#for s in sites:
#    cmd = "perl ~/bin/rangeToSites.pl /data/robert.williamson/CNC/"+s+".CNS.range /data/robert.williamson/CNC/"+s+".CNS.sites"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#while len(procs) > procLimit:
#    procs[0].wait()
#    #print("one done")
#    procs.pop(0)
#        
#for s in sites:
#    for i in scafs:
#        cmd = "grep -P \"scaffold_"+str(i)+"\s\" /data/robert.williamson/CNC/"+s+".CNS.sites > /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".CNS.sites"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
#                
#    while len(procs) > procLimit:
#        procs[0].wait()
#        #print("one done")
#        procs.pop(0)
   
#sites = ["allSites"]
       
#for s in sites:
#    for i in scafs:
#        cmd = "sort -gk 2 /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".CNS.sites > /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".sorted.sites"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
#                
#    while len(procs) > procLimit:
#        procs[0].wait()
#        #print("one done")
#        procs.pop(0)
#  
#waitAll(procs)
#         
##filter sites
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/filter.py /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".sorted.sites -q 40 -d 20 -D 60 -G -l scaf"+str(i)+"."+s+".flt.log > /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".gatk.flt.sites 2>"+str(i)+s+".filter.err"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
#            
#    while len(procs) > procLimit:
#        procs[0].wait()
#        #print("one done")
#        procs.pop(0)
#
#waitAll(procs)
##for range(num bootstraps)
##bootsrap that beetch!
#          
##sub vcf
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/subMpileup.py /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".gatk.flt.sites "+s+".sub.log >/data/robert.williamson/CNC/scaf"+str(i)+"."+s+".gatk.vcf 2>"+str(i)+s+".sub.err"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
#            
#    while len(procs) > procLimit:
#        procs[0].wait()
#        #print("one done")
#        procs.pop(0)
##  
#waitAll(procs)        
###get hists
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/getHistsGATK.py /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".gatk.vcf scaf"+str(i)+"."+s+".gatk.depth scaf"+str(i)+"."+s+".gatk.af1s scaf"+str(i)+"."+s+".gatk.quals"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
#            
#    while len(procs) > procLimit:
#        procs[0].wait()
#        #print("one done")
#        procs.pop(0)
#
#waitAll(procs)
###ADD IN DIVERGENCE
for s in sites:
    for i in scafs:
        cmd = "python ~/bin/divergence2.py /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".gatk.flt.sites -s /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".div.flt.sites /data/mathieu.blanchette/alignmentsStephen/CR_v_NP.chain.ortho.Scaffold"+str(i)+".sorted.maf > /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".div" 
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
#generate DIV totals
cmd = "tail -n 1 /data/robert.williamson/finalGenome_analysis/div/4fold.div > all.div"
try:
    p = Popen(cmd, shell=True)
    procs.append(p)
except OSError, e:
    print ("***Execution failed: " + e + "\t"+cmd)
 
waitAll(procs)
    
for s in sites:
    cmd = "python /data/robert.williamson/finalGenome_analysis/div/grabber.py /data/robert.williamson/CNC/scaf*."+s+".div >> all.div"
    try:
        p = Popen(cmd, shell=True)
        procs.append(p)
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)
    waitAll(procs)

waitAll(procs)

#generate AFS totals
for s in sites:
    cmd = "perl ~/bin/afsCleaner2.pl "+s+".gatk.af1s scaf*."+s+".gatk.af1s"
    try:
        p = Popen(cmd, shell=True)
        procs.append(p)
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)

waitAll(procs)

#sum the totals and glue together
cmd = "python /data/robert.williamson/finalGenome_analysis/wraps/summer.py /data/robert.williamson/finalGenome_analysis/wraps/GATK/total.4fold.af1s > all.gatk.af1s"
try:
    p = Popen(cmd, shell=True)
    procs.append(p)
except OSError, e:
    print ("***Execution failed: " + e + "\t"+cmd)
waitAll(procs)
   
for s in sites:
    cmd = "python /data/robert.williamson/finalGenome_analysis/wraps/summer.py "+s+".gatk.af1s >> all.gatk.af1s"
    try:
        p = Popen(cmd, shell=True)
        procs.append(p)
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)
    waitAll(procs)
waitAll(procs)


#generate the DFE inputs
cmd = "python ~/bin/dfe_a_gen.py all.gatk.af1s all.div names 26 0 gatk_CNC"
try:
    p = Popen(cmd, shell=True)
    procs.append(p)
except OSError, e:
    print ("***Execution failed: " + e + "\t"+cmd)
waitAll(procs)

##REPEAT FOR DIV FILTERED SITES
#sub vcf
for s in sites:
    for i in scafs:
        cmd = "python ~/bin/subMpileup.py /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".div.flt.sites "+s+".sub.log >/data/robert.williamson/CNC/scaf"+str(i)+"."+s+".div.vcf 2>"+str(i)+s+".sub.err"
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
##get hists
for s in sites:
    for i in scafs:
        cmd = "python ~/bin/getHistsGATK.py /data/robert.williamson/CNC/scaf"+str(i)+"."+s+".div.vcf scaf"+str(i)+"."+s+".div.depth scaf"+str(i)+"."+s+".div.af1s scaf"+str(i)+"."+s+".div.quals"
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
#generate AFS totals
for s in sites:
    cmd = "perl ~/bin/afsCleaner2.pl "+s+".div.af1s scaf*."+s+".div.af1s"
    try:
        p = Popen(cmd, shell=True)
        procs.append(p)
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)

waitAll(procs)

#sum the totals and glue together
cmd = "python /data/robert.williamson/finalGenome_analysis/wraps/summer.py /data/robert.williamson/finalGenome_analysis/wraps/GATK/total.4fold.af1s > all.div.af1s"
try:
    p = Popen(cmd, shell=True)
    procs.append(p)
except OSError, e:
    print ("***Execution failed: " + e + "\t"+cmd)
waitAll(procs)

for s in sites:
    cmd = "python /data/robert.williamson/finalGenome_analysis/wraps/summer.py "+s+".div.af1s >> all.div.af1s"
    try:
        p = Popen(cmd, shell=True)
        procs.append(p)
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)
    waitAll(procs)

waitAll(procs)

#generate the DFE inputs
cmd = "python ~/bin/dfe_a_gen.py all.div.af1s all.div names 26 0 div_CNC"
try:
    p = Popen(cmd, shell=True)
    procs.append(p)
except OSError, e:
    print ("***Execution failed: " + e + "\t"+cmd)

#for i in scafs:
#    cmd = "perl Polymorphorama_newPopTags.txt ../finalGenome_analysis/fastas/GATK/scaf"+str(i)+"/ 0.0 26 > scaf"+str(i)+".poly.out 2> scaf"+str(i)+".err"
#    try:
#        p = Popen(cmd, shell=True)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#    p.wait()
#            
#    cmd = "cp frequencies.xls /data/robert.williamson/finalGenome_analysis/polyruns/frequencies.scaf"+str(i)+".xls 2> cp"+str(i)+".err"
#    try:
#        p = Popen(cmd, shell=True)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#    p.wait()
#    
#    cmd = "cp summarystats.xls /data/robert.williamson/finalGenome_analysis/polyruns/summarystats.scaf"+str(i)+".xls 2> cp"+str(i)+".err"
#    try:
#        p = Popen(cmd, shell=True)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#    p.wait()