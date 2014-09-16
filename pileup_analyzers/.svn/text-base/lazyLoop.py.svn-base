import sys
from subprocess import call, Popen
#
sites = ["4fold","0fold","intron","intergene","3utr","5utr"]#,"CNC"]
#sites = ["4fold", "first","non"]
#sites = ['ambig', 'down', 'intergene', 'intron', 'snc', 'up', '3utr', '5utr']
scafs = [8,7,6,5,4,3,2,1]

#subset sites
procs = []
procLimit = 10

def waitAll(procs):
    while len(procs) > 0:
        procs[0].wait()
        #print("one done")
        procs.pop(0)

#for s in sites:
#    for i in scafs:
#        cmd = "/bin/grep -P \"scaffold_"+str(i)+"\\s\" /data/robert.williamson/finalGenome_analysis/rubella/"+s+".rub2.txt > /data/robert.williamson/finalGenome_analysis/rubella/scaf"+str(i)+"."+s+".sites"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
#            
#while len(procs) > procLimit:
#    procs[0].wait()
#    #print("one done")
#    procs.pop(0)
#
##print("all DONE!")
#filter sites
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/filter.py /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".sites -q 90 -d 20 -D 60 -G -l scaf"+str(i)+"."+s+".flt.log > /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.flt.sites 2>"+str(i)+s+".filter.err"
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

#filter ALL data
#for i in scafs:
#    cmd = "python ~/bin/filter.py /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf -q 90 -d 20 -D 60 -G -l scaf"+str(i)+".flt.log > /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+".gatk.flt.sites 2>"+str(i)+".filter.err"
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
##
#waitAll(procs)

###
####for range(num bootstraps)
####bootsrap that beetch!
###          
####sub vcf
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/subMpileup.py /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.flt.sites "+s+".sub.log >/data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.vcf 2>"+str(i)+s+".sub.err"
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
####          
#####get hists
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/getHistsGATK.py /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.vcf scaf"+str(i)+"."+s+".gatk.depth scaf"+str(i)+"."+s+".gatk.af1s scaf"+str(i)+"."+s+".gatk.quals -N 26"
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
###
####sites = ["intergene","intron","0fold","4fold","3utr","5utr"]
#######ADD IN DIVERGENCE
##for s in sites:
##    for i in scafs:
##        cmd = "python ~/bin/divergence2.py /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.flt.sites /data/mathieu.blanchette/alignmentsStephen/CR_v_NP.chain.ortho.Scaffold"+str(i)+".sorted.maf > /data/robert.williamson/finalGenome_analysis/div/scaf"+str(i)+"."+s+".gatk.div" 
##        try:
##            p = Popen(cmd, shell=True)
##            procs.append(p)
##        except OSError, e:
##            print ("***Execution failed: " + e + "\t"+cmd)
##            
##    while len(procs) > procLimit:
##        procs[0].wait()
##        #print("one done")
##        procs.pop(0)
## 
##waitAll(procs)
# 
###DIVERGENCE PILEUP       
##for s in sites:
##    for i in scafs:
##        cmd = "python ~/bin/getDiffToRef.py /data/robert.williamson/finalGenome_analysis/neslia/scaf"+str(i)+".pile /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.flt.sites -q 40 -d 10 -D 70 -s 30 -p -l scaf"+str(i)+"."+s+".log > /data/robert.williamson/finalGenome_analysis/div/scaf"+str(i)+"."+s+".gatk.sam.div" 
##        try:
##            p = Popen(cmd, shell=True)
##            procs.append(p)
##        except OSError, e:
##            print ("***Execution failed: " + e + "\t"+cmd)
##            
##    while len(procs) > procLimit:
##        procs[0].wait()
##        #print("one done")
##        procs.pop(0)
##
##waitAll(procs)
#
###DIVERGENCE PICKLE       
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/divergencePickler.py /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.flt.sites -d /data/robert.williamson/finalGenome_analysis/div/pickle/scaf"+str(i)+"_div.p > /data/robert.williamson/finalGenome_analysis/div/scaf"+str(i)+"."+s+".pickle.div 2> err"+str(s)+str(i)+".log" 
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
##   
#for s in sites:
#    cmd = "python /data/robert.williamson/finalGenome_analysis/div/grabber.py /data/robert.williamson/finalGenome_analysis/div/scaf*."+s+".pickle.div "
#    if s == "4fold":
#        cmd += "> all.div"
#    else:
#        cmd += ">> all.div"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#    waitAll(procs)
#
#waitAll(procs)
##
###generate AFS totals
#for s in sites:
#    cmd = "perl ~/bin/afsCleaner2.pl "+s+".gatk.af1s scaf*."+s+".gatk.af1s"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#
#waitAll(procs)
##
###sum the totals and glue together  
#for s in sites:
#    cmd = "python /data/robert.williamson/finalGenome_analysis/wraps/summer.py "+s+".gatk.af1s "
#    if s == "4fold":
#        cmd += "> all.gatk.af1s"
#    else:
#        cmd += ">> all.gatk.af1s"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#    waitAll(procs)
#    #
#waitAll(procs)
##
##
###generate the DFE inputs
#cmd = "python ~/bin/dfe_a_gen.py all.gatk.af1s all.div names 26 0 grandiflora_introns_3_6"
#try:
#    p = Popen(cmd, shell=True)
#    procs.append(p)
#except OSError, e:
#    print ("***Execution failed: " + e + "\t"+cmd)
#waitAll(procs)

#for i in scafs:
#    cmd = "perl Polymorphorama_newPopTags.txt ../finalGenome_analysis/fastas/GATK_new/scaf"+str(i)+"/ 0.0 26 > scaf"+str(i)+".poly.out 2> scaf"+str(i)+".err"
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

#for s in sites:
#    for i in scafs:
#        cmd = "perl ~/bin/getAF.pl /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.vcf /data/robert.williamson/finalGenome_analysis/sites/per_site_afs/scaf"+str(i)+"."+s+".sites.afs 26"
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

#for s in sites:
#    cmd = "python ~/bin/subRubella2.py /data/daniel.koenig/finished_mat.txt /data/robert.williamson/finalGenome_analysis/sites/"+s+".sites "+s+".rubella.afs -N 24 > "+s+".rub.flt.sites"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#        

#sliding window stats
#for i in scafs:
#        cmd = "python ~/bin/slidingWindowStats.py /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf -w 100 -s 100 -q 90 -d 20 -D 60 -l 40 -N 13 > scaf"+str(i)+".100100SNP.csv"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
#
#
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/slidingWindowStats.py /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".gatk.vcf -w 100 -s 100 -q 90 -d 20 -D 60 -l 40 -N 13 > scaf"+str(i)+"."+s+".100100SNP.csv"
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
#filter ALL data

#cmd = "python ~/bin/subRubella2.py /data/daniel.koenig/finished_mat.txt /data/daniel.koenig/finished_mat.txt rubella.afs -N 24 > all.rub.flt.sites"
#try:
#    p = Popen(cmd, shell=True)
#    procs.append(p)
#except OSError, e:
#    print ("***Execution failed: " + e + "\t"+cmd)

#for s in sites:
#        cmd = "python ~/bin/alignAFS.py /data/robert.williamson/finalGenome_analysis/sites/per_site_afs/"+s+".sites.afs /data/robert.williamson/finalGenome_analysis/rubella/"+s+".rub.flt.sites -s /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"+s+".cg.cr.sites > /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"+s+".aligned 2> "+s+".err"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)

#for i in scafs:
#    cmd = "perl ~/bin/getAF.pl /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf /data/robert.williamson/finalGenome_analysis/sites/per_site_afs/scaf"+str(i)+".sites.afs 26"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#
#waitAll(procs)
#            
#for i in scafs:
#        cmd = "python ~/bin/alignAFS.py /data/robert.williamson/finalGenome_analysis/sites/per_site_afs/scaf"+str(i)+".sites.afs /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+".rub.flt.sites -s /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/scaf"+str(i)+".cg.cr.sites > /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/scaf"+str(i)+".aligned 2> "+str(i)+".err"
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)            
#
#waitAll(procs)
            
##
#waitAll(procs)
##split up sites for divergence
#for s in sites:
#    for i in scafs:
#        cmd = "/bin/grep -P \"scaffold_"+str(i)+"\\s\" "+s+".rub.flt.sites > scaf"+str(i)+"."+s+".rub.flt.sites"
#        print(cmd)
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
# 
#waitAll(procs)   
##DIVERGENCE
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/divergencePickler.py scaf"+str(i)+"."+s+".rub.flt.sites -d /data/robert.williamson/finalGenome_analysis/div/pickle/scaf"+str(i)+"_div.p > /data/robert.williamson/finalGenome_analysis/rubella/scaf"+str(i)+"."+s+".pickle.div 2> err"+str(s)+str(i)+".log" 
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

#glue all the divs together
#for s in sites:
#    cmd = "python /data/robert.williamson/finalGenome_analysis/div/grabber.py /data/robert.williamson/finalGenome_analysis/rubella/scaf*."+s+".pickle.div "
#    if s == "4fold":
#        cmd += "> all.div"
#    else:
#        cmd += ">> all.div"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#    waitAll(procs)
    
#put the AFS all into one file
#for s in sites:
#    cmd = "tail -n 1 "+s+".rubella.afs "
#    if s == "4fold":
#        cmd += "> all.afs"
#    else:
#        cmd += ">> all.afs"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#    waitAll(procs)
#        
##generate the DFE inputs
#cmd = "python ~/bin/dfe_a_gen.py all.afs all.div names 24 0 rubella_shared_2_27"
#try:
#    p = Popen(cmd, shell=True)
#    procs.append(p)
#except OSError, e:
#    print ("***Execution failed: " + e + "\t"+cmd)
#waitAll(procs)
        
        
#get only sites that pass filters in both grandiflora AND rubella
#for s in sites:
#    for i in scafs:
#        cmd = "python ~/bin/subRubella2.py /data/daniel.koenig/finished_mat.txt /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".sites /data/robert.williamson/finalGenome_analysis/sites/scaf"+str(i)+"."+s+".rubella.afs -N 24 > /data/robert.williamson/finalGenome_analysis/sharedSites/scaf"+str(i)+"."+s+".shared.sites"
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
        
#get out unique sites in each species
#for s in sites:
#    cmd = "/bin/grep -P \"^0\" /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"+s+".aligned > /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"+s+".cr.unique"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#    p.wait()
#    cmd = "/bin/grep -P \"^\d+\s+0\s+\d+\" /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"+s+".aligned > /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"+s+".cg.unique"
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#        
#    p.wait()
#    
## 
#waitAll(procs)

#clean the AFS
#for s in sites:
#    cmd = "python ~/bin/uniqueAFSCleaner.py /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"+s+".cr.unique"
#    if s == "4fold":
#        cmd += " > cr.afs "
#    else:
#        cmd += " >> cr.afs "
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#    p.wait()
#    cmd = "python ~/bin/uniqueAFSCleaner.py /data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"+s+".cg.unique"
#    if s == "4fold":
#        cmd += " > cg.afs "
#    else:
#        cmd += " >> cg.afs "
#    try:
#        p = Popen(cmd, shell=True)
#        procs.append(p)
#    except OSError, e:
#        print ("***Execution failed: " + e + "\t"+cmd)
#            
#    p.wait()
#
#    
## 
#waitAll(procs)

#bootstrapping grandiflora
#make names
#for i in range(4):
#    for s in sites:
#        cmd = "echo \""+str(s)+"_"+str(i)+"\""
#        if i == 0 and s == "4fold":
#            cmd += " > names"
#        else:
#            cmd += " >> names"
#            
#        try:
#            p = Popen(cmd, shell=True)
#            procs.append(p)
#        except OSError, e:
#            print ("***Execution failed: " + e + "\t"+cmd)
#bootstrap
boots = 100
files = -1
for j in range(int(boots/4)):
    files += 1
    for num in range(4):
        #bootstrap
        for s in sites:
            for i in scafs:
                cmd = "/opt/python27/bin/python27 ~/bin/bootstrapSites_rubella.py /data/robert.williamson/finalGenome_analysis/rubella/sites/scaf"+str(i)+"."+s+".flt.sites /storage/robert.williamson/finalGenome_analysis/scaf"+str(i)+"."+s+".boot.sites > scaf"+str(i)+"."+s+".boot.afs"
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

        ###sub vcf
        #GRAND
#        for s in sites:
#            for i in scafs:
#                cmd = "python ~/bin/subMpileup.py /data/robert.williamson/finalGenome_analysis/GATK/scaf"+str(i)+".gatk.vcf /data/robert.williamson/finalGenome_analysis/boot/scaf"+str(i)+"."+s+".boot.sites "+s+".sub.log > scaf"+str(i)+"."+s+".gatk.vcf 2>"+str(i)+s+".sub.err"
#                try:
#                    p = Popen(cmd, shell=True)
#                    procs.append(p)
#                except OSError, e:
#                    print ("***Execution failed: " + e + "\t"+cmd)
#                    
#            while len(procs) > procLimit:
#                procs[0].wait()
#                #print("one done")
#                procs.pop(0)
#        
#        waitAll(procs)
        ###          
        
        #RUB
#        for s in sites:
#            for i in scafs:
#                cmd = "python ~/bin/subRubella2.py /data/daniel.koenig/finished_mat.txt /storage/robert.williamson/finalGenome_analysis/scaf"+str(i)+"."+s+".boot.sites /storage/robert.williamson/finalGenome_analysis/scaf"+str(i)+"."+s+".boot.afs -N 24 > scaf"+str(i)+"."+s+".boot.sites.flt 2>"+str(i)+s+".sub.err"
#                try:
#                    p = Popen(cmd, shell=True)
#                    procs.append(p)
#                except OSError, e:
#                    print ("***Execution failed: " + e + "\t"+cmd)
#                    
#            while len(procs) > procLimit:
#                procs[0].wait()
#                #print("one done")
#                procs.pop(0)
#        
#        waitAll(procs)        
        
        ####get hists
        #GRAND
#        for s in sites:
#            for i in scafs:
#                cmd = "python ~/bin/getHistsGATK.py scaf"+str(i)+"."+s+".gatk.vcf scaf"+str(i)+"."+s+".gatk.depth scaf"+str(i)+"."+s+".gatk.af1s scaf"+str(i)+"."+s+".gatk.quals -N 26"
#                try:
#                    p = Popen(cmd, shell=True)
#                    procs.append(p)
#                except OSError, e:
#                    print ("***Execution failed: " + e + "\t"+cmd)
#                    
#            while len(procs) > procLimit:
#                procs[0].wait()
#                #print("one done")
#                procs.pop(0)
#        
#        waitAll(procs)
        
        ##DIVERGENCE PICKLE   
		
        for i in scafs:
            input = ""
            for s in sites:
                    input += ",/storage/robert.williamson/finalGenome_analysis/scaf"+str(i)+"."+s+".boot.sites"
                
            input = input[1:]
            cmd = "python ~/bin/divergencePickler.py "+input+" -d /data/robert.williamson/finalGenome_analysis/div/pickle/scaf"+str(i)+"_div.p 2> err"+str(s)+str(i)+".log" 
            try:
                p = Popen(cmd, shell=True)
                procs.append(p)
            except OSError, e:
                print ("***Execution failed: " + e + "\t"+cmd)
			

        waitAll(procs)
        #   
        for s in sites:
            cmd = "python /data/robert.williamson/finalGenome_analysis/div/grabber.py scaf*."+s+".boot.sites.div "
            if s == "4fold" and num == 0:
                cmd += "> all.div"
            else:
                cmd += ">> all.div"
            try:
                p = Popen(cmd, shell=True)
                procs.append(p)
            except OSError, e:
                print ("***Execution failed: " + e + "\t"+cmd)
            waitAll(procs)
        
        waitAll(procs)
        
        #generate AFS totals
        #GRAND
#        for s in sites:
#            cmd = "perl ~/bin/afsCleaner2.pl "+s+".gatk.af1s scaf*."+s+".gatk.af1s"
#            try:
#                p = Popen(cmd, shell=True)
#                procs.append(p)
#            except OSError, e:
#                print ("***Execution failed: " + e + "\t"+cmd)
#        
#        waitAll(procs)
        
        #RUB
        for s in sites:
            cmd = "python /data/robert.williamson/finalGenome_analysis/rubella/boot/afs_glue.py scaf*."+s+".boot.afs > "+s+".afs"
            try:
                p = Popen(cmd, shell=True)
                procs.append(p)
            except OSError, e:
                print ("***Execution failed: " + e + "\t"+cmd)
            waitAll(procs)
        
        #sum the totals and glue together  
        for s in sites:
            cmd = "python /data/robert.williamson/finalGenome_analysis/wraps/summer.py "+s+".afs "
            if s == "4fold" and num == 0:
                cmd += "> all.gatk.af1s"
            else:
                cmd += ">> all.gatk.af1s"
            try:
                p = Popen(cmd, shell=True)
                procs.append(p)
            except OSError, e:
                print ("***Execution failed: " + e + "\t"+cmd)
            waitAll(procs)
            
        waitAll(procs)
        
        
        #generate the DFE inputs
    cmd = "python ~/bin/dfe_a_gen.py all.gatk.af1s all.div names 24 0,6,12,18 rubella_boot_5_20_"+str(files)
    try:
        p = Popen(cmd, shell=True)
        procs.append(p)
    except OSError, e:
        print ("***Execution failed: " + e + "\t"+cmd)
    waitAll(procs)