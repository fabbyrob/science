import sys
import random
import matplotlib.pyplot as plt
import numpy
from multiprocessing import Pool
import time

class individual():
    def __init__(self, mito, nuc):
        self.mito = mito
        self.nuc = nuc
        self.fitness = None
        return
    
    def calcFitness(self, nucCost = 0.00, mitoCost = 0.00):
        self.fitness = 1-(nucCost*sum(self.nuc)+self.mito*mitoCost)
        return self.fitness
        
    def mate(self, father, mu, mu_mito):
        #if father.mito and 1 not in father.nuc:#mito in father has CMS gene and nuc has no restorer
        #    return None
        
        child = individual(self.mito, [random.choice(self.nuc), random.choice(father.nuc)])
        
        #mutate kids genome
        if random.random() <= mu_mito:#mutate mito
            child.mito = -1*(child.mito-1)#switch mito allele
            
        for i in range(0,2):#try and mutate each allele
            if random.random() <= mu:
                child.nuc[i] = -1*(child.nuc[i]-1)
        
        return child
        
    def __lt__(self, other):
        return self.fitness < other.fitness
    
    def __gt__(self, other):
        return self.fitness > other.fitness
    
    def __eq__(self, other):
        return self.fitness == other.fitness
        
    def __str__(self):
        return ("Individual M:%s N:%s F:%.2f" % (self.mito, self.nuc, self.fitness))

class population():
    def __init__(self, N, nucCost, mitoCost, mu, mu_mito, plots):
        self.N = N
        self.nucCost = nucCost
        self.mitoCost = mitoCost
        self.mu = mu
        self.mu_mito = mu_mito
        self.fitnesses = []
        self.mitoFreqs = []
        self.nucFreqs = []
        self.numFar = []
        self.resFars = []
        
        self.pop = []
        while len(self.pop) < self.N:
            self.pop.append(individual(random.randint(0,1), [random.randint(0,1), random.randint(0,1)])) 
            
        self.meanFitness()
        self.mitoFreq()
        self.restorerFreq() 
        
        if plots:
            self.fitPlot, self.mitPlot, self.nucPlot, self.farPlot, self.farResPlot = plots
            self.update_plot(self.fitPlot, [0, self.fitnesses[-1]])
            self.update_plot(self.mitPlot, [0, self.mitoFreqs[-1]])
            self.update_plot(self.nucPlot, [0, self.nucFreqs[-1]])
            self.plots = 1
        else:
            self.plots = 0
            
        #print("%s,%s,%s,%s" % ("generation","fitness","CMS","Restorer"))
        #print("%s,%s,%s,%s" % (0, self.fitnesses[-1], self.mitoFreqs[-1], self.nucFreqs[-1]))
            
    def meanFitness(self):
        fits = sum([ind.calcFitness(self.nucCost, self.mitoCost) for ind in self.pop])
        self.fitnesses.append(float(fits)/self.N)
        return self.fitnesses[-1]
     
    def mitoFreq(self):
        mitos = sum([ind.mito for ind in self.pop])
        self.mitoFreqs.append(float(mitos)/self.N)
        return self.mitoFreqs[-1]
    
    def restorerFreq(self):
        restorers = sum([sum(ind.nuc) for ind in self.pop])
        self.nucFreqs.append(float(restorers)/(2*self.N))
        return self.nucFreqs[-1]
    
    def pickParent(self, parents, other = None):
        if len(parents) == 1:
            return parents[0]
        
        tries = 0
        p1 = random.choice(parents)
        while (p1 is other):
            tries += 1
            p1 = random.choice(parents)
            if tries > 200:
                sys.stderr.write("Could not find another parent. Trying again.\n")
                return None
         
        tries = 0
        p2 = random.choice(parents)
        while (p2 is p1 or p2 is other):
            tries += 1
            p2 = random.choice(parents)
            if tries > 200:
                sys.stderr.write("Could not find another parent. Trying again.\n")
                return None
            
        if p1 < p2:
            return p2
        else:
            return p1
        
    
    def generation(self, T):
        newPop = []
        
        #pull out all individuals that can produce pollen
        fathers = []
        resFar = 0.0
        
        for ind in self.pop:
            if not ind.mito:
                fathers.append(ind)
            elif (ind.mito and 1 in ind.nuc):
                fathers.append(ind)
                resFar += 1
        self.numFar.append(len(fathers)/float(self.N))
        self.resFars.append(resFar/float(self.N))#len(fathers))
                
        if len(fathers) == 0:
            sys.stderr.write("No fathers in the population. Exiting.\n")
            return 0
        
        fails = 0
        while len(newPop) < self.N:
            tries = 0
            father = self.pickParent(fathers)
            mother = self.pickParent(self.pop, father) 
            while not father or not mother:
                tries += 1
                father = self.pickParent(fathers)
                mother = self.pickParent(self.pop, father)  
                if tries > 200:
                    sys.stderr.write("Could not pick two parents. Exiting.\n")
                    return 0
            
            offspring = mother.mate(father, self.mu, self.mu_mito)
            if offspring:
                newPop.append(offspring)
            else:
                fails += 1
                
            if fails > 4*self.N:
                sys.stderr.write("over 4N failed offspring, population not viable. Exiting.\n")
                return 0
                
        self.pop = newPop
        
        
        self.meanFitness()
        self.mitoFreq()
        self.restorerFreq() 
        
        if self.plots and T % max(1, int(T*0.1)) == 0:
            self.update_plot(self.fitPlot, [T, self.fitnesses[-1]])
            self.update_plot(self.mitPlot, [T, self.mitoFreqs[-1]])
            self.update_plot(self.nucPlot, [T, self.nucFreqs[-1]])
            self.update_plot(self.farPlot, [T-1, self.numFar[-1]])
            self.update_plot(self.farResPlot, [T-1, self.resFars[-1]])
            
        #print("%s,%s,%s,%s" % (T, self.fitnesses[-1], self.mitoFreqs[-1], self.nucFreqs[-1]))
        return 1

    def update_plot(self, hl, new_data):
         hl.set_xdata(numpy.append(hl.get_xdata(), new_data[0]))
         hl.set_ydata(numpy.append(hl.get_ydata(), new_data[1]))
         plt.draw()

def run_sim(args):
    start = time.time()
        #plotting stuff
    K, p = args
    if p:
        plt.ion()
    
        plt.subplot(221)
        fitPlot, = plt.plot([],[],"r")
        plt.ylim([.95,1.01])
        plt.xlim([0, K])
        plt.ylabel("fitness")
        
        plt.subplot(222)
        mitPlot, = plt.plot([],[],"g")
        plt.ylim([0,1.11])
        plt.xlim([0, K])
        plt.ylabel("CMS")
        
        plt.subplot(224)
        nucPlot, = plt.plot([],[],"g")
        plt.ylim([0,1.11])
        plt.xlim([0, K])
        plt.ylabel("Restorer")
        plt.xlabel("generation")
        
        plt.subplot(223)
        farPlot,farResPlot, = plt.plot([],[],"b",[],[],"r")
        plt.ylim([0,1.11])
        plt.xlim([0, K])
        plt.ylabel("Max Fathers")
        plt.xlabel("generation")
        
        p = (fitPlot, mitPlot, nucPlot, farPlot, farResPlot)
    
    pop = population(N=500, nucCost=0.000, mitoCost=0.0, mu=0.00, mu_mito=0.000, plots=p)
        
    for i in range(K):
        status = pop.generation(i+1)
        if not status:
            break
        
    sys.stderr.write("Run time: %.2f min\n" % ((time.time()-start)/60.0))
    return pop

def __main__():
    K =100#number of generations to run
    p = 0#turns on/off plotting while simulation is running (its really damn slow if it's on)
    multi = 1#turns on/off multiprocessing
    if multi:#cant plot while simulating if multiprocessing is being used
        p = 0
    R = 30#number of replicates
    pops = []
    
    if multi:
        pool = Pool(processes=2)
        pops = pool.map(run_sim, [(K, p)]*R)
    else:
        for i in range(R):
            #simulation start
            pop = run_sim((K, p))
            
            pops.append(pop)
        
    if not p:
        for pop in pops:
            plt.subplot(221)
            fitPlot, = plt.plot(pop.fitnesses,"r")
            #plt.ylim([.999,1.001])
            plt.xlim([0, K])
            plt.ylabel("fitness")
            
            plt.subplot(222)
            mitPlot, = plt.plot(pop.mitoFreqs,"g")
            plt.ylim([0,1.11])
            plt.xlim([0, K])
            plt.ylabel("CMS")
            
            plt.subplot(224)
            nucPlot, = plt.plot(pop.nucFreqs,"g")
            plt.ylim([0,1.11])
            plt.xlim([0, K])
            plt.ylabel("Restorer")
            plt.xlabel("generation")
            
            plt.subplot(223)
            farPlot,farResPlot, = plt.plot(pop.numFar,"b", pop.resFars, "r")
            plt.ylim([0,1.11])
            plt.xlim([0, K])
            plt.ylabel("Max Fathers")
            plt.xlabel("generation")
        
    plt.show()

if __name__ == "__main__":
    __main__()