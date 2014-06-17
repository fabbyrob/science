import sys
from subprocess import call, Popen

folder = 'rewrite_small_pop'
selection_tes = [0.001]
replicates = 100
N = [100]
loci = [200]
selfing = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
teNum = [2]#% of number of loci to start as TEs
generations = [5000]
crossovers = [2]#% of number of loci
calc = ['HOMO','ALL']#g for fitness
transrate = [0.005]
excisionrate = [0]
exponent = [1.5]

params = [selection_tes, N, loci, selfing, teNum, generations, crossovers, calc, transrate, exponent, excisionrate]
tot = 1
for x in params:
    tot *= len(x)
tot *= replicates
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
                for t in teNum:
                    num = t#int(l*t)
                    path3 = path2 + "_TEs"+str(num)
                    for c in crossovers:
                        crosses = c#int(c*l)
                        path4 = path3 + "_cross"+str(crosses)
                        for g in calc:
                            path5 = path4 + "_g"+g
                            for u in transrate:
                                path6 = path5 + "_u"+str(u)
                                for e in exponent:
                                    path7 = path6 + "_t"+str(e)
                                    for x in selection_tes:
                                        path8 = path7 + "_x"+str(x)
                                        for v in excisionrate:
                                            path9 = path8 + "_v"+str(v)
                                            for R in range(0, replicates):
                                                path10 = path9 + "_R"+str(R)#may need python 2.7...
                                                cmd = "/opt/python27/bin/python27 ~/repos/robert.williamson/TE_sim_simple/TEsim_simp.py "+" -T "+str(num)+" -N "+str(n)+" -G "+str(l)+" -g "+str(g)+" -X "+str(crosses)+" -r "+str(r)+" -K "+str(G)+" -t "+str(e)+" -u "+str(u)+" -x "+str(x)+" > "+path10
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
