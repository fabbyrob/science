import sys
from subprocess import call, Popen

folder = 'confirm_wright_test'
selection_tes = [0.001]
selection_silencing = [0]
replicates = 5
N = [1000]
loci = [200]
selfing = [0.1, 0.5, 0.9]
silence = [0]
silenceMutation = [0]
silenceEpsilon = [0]
teNum = [10]#% of number of loci to start as TEs
generations = [1000]
crossovers = [2]#% of number of loci
calc = ['ALL','HET','HOMO']#g for fitness
transrate = [0.005]
exponent = [1.5]
seltype = ['R']
maxClasses = [1]
classMutation = [0]
excision = [0]

l = [selection_tes, selection_silencing, N, loci, selfing, silence, silenceMutation, silenceEpsilon, teNum, generations, crossovers, calc, transrate, exponent, seltype, maxClasses, classMutation, excision]

tot = replicates
for thingy in l:
    tot *= len(thingy)
    
print(tot)
#sys.exit()
procs = []
maxProcs = 15

for n in N:
    for l in loci:
        pathStart = "/data/robert.williamson/TEsim/results/"+folder+"/N"+str(n)+"_loci"+str(l)+"_"
        for G in generations:
            path1 = pathStart + "_generations"+str(G)
            for r in selfing:
                path2 = path1 + "_selfing"+str(r)
                for s in silence:
                    path3 = path2 + "_silence"+str(s)
                    for sd in silenceMutation:
                        path4 = path3 +"_silMut"+str(sd)
                        for t in teNum:
                            num = t#int(l*t)
                            path5 = path4 + "_TEs"+str(num)
                            for c in crossovers:
                                crosses = c#int(c*l)
                                path6 = path5 + "_cross"+str(crosses)
                                for g in calc:
                                    path7 = path6 + "_g"+g
                                    for u in transrate:
                                        path8 = path7 + "_u"+str(u)
                                        for e in exponent:
                                            path9 = path8 + "_t"+str(e)
                                            for b in maxClasses:
                                                path10 = path9 + "_b"+str(b)
                                                for c in classMutation:
                                                    path11 = path10 + "_c"+str(c)
                                                    for st in seltype:
                                                        path12 = path11 + "_S"+st
                                                        for z in selection_silencing:
                                                            path13 = path12 + "_z"+str(z)
                                                            for x in selection_tes:
                                                                path14 = path13 + "_x"+str(x)
                                                                for epsilon in silenceEpsilon:
                                                                    path15 = path14 + "_e"+str(epsilon)
                                                                    for v in excision:
                                                                        path16 = path15 + "_v"+str(v)
                                                                        for R in range(0, replicates):
                                                                            path17 = path16 + "_R"+str(R)#may need python 2.7...
                                                                            cmd = "/opt/python27/bin/python27 ~/bin/TEsim.py "+" -T "+str(num)+" -l "+str(s)+" -N "+str(n)+" -G "+str(l)+" -g "+str(g)+" -X "+str(crosses)+" -s "+str(sd)+" -r "+str(r)+" -K "+str(G)+" -t "+str(e)+" -u "+str(u)+" -S "+st+" -b "+str(b)+" -c "+str(c)+" -x "+str(x)+" -z "+str(z)+" -e "+str(epsilon)+" -v "+str(v)+" > "+path17
                                                                            #print(cmd)
                                                                            try:
                                                                                p = Popen(cmd, shell=True)
                                                                                procs.append(p)
                                                                            except OSError, e:
                                                                                print ("***Execution failed: " + e + "\t"+cmd)
                                                        
                                                                            #only run up to 10 jobs at once
                                                                            while len(procs) > maxProcs:
                                                                                procs[0].wait()
                                                                                #print("one done")
                                                                                procs.pop(0)
