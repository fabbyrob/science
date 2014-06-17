doc = '''
Simulates a population with a quantitative trait under stabilizing selection, controlled by many loci. 

Lets the population reach equilibrium, then the population is split into two new populations. 

Evolution continues, and the proportion of adaptive fixations (alpha) is calculated each generation.
'''

import sys
import random
import time
import getopt

_N = 100
_L = 10
_u = 0.01
_E = 100
_K = 100
_o = 0
_t = 1

def __main__():
    #check aruguments
    processArgs(1)
    
    #initialize pop with same params
    sys.stderr.write("Initializing population...\n")
    initialPop = Pop(N = _N, loci = _L, mu = _u, optimum = _o)
    initialPop.generate()
    initialPop.printGeneration(0, "I", "NA")
    
    #reach equilibrium
    for i in range(1, _E):
        initialPop.generation()
        initialPop.printGeneration(i, "I", "NA")
    sys.stderr.write("Finished letting population equilibrate...\n")
    
    #sample from pop with replacement, to make 2 new populations
    popA = initialPop.makeNewPop()
    popB = initialPop.makeNewPop()
    popA.calcMeans()
    popB.calcMeans()
    alpha = popA.alpha(popB)
    popA.printGeneration(_E, "A", alpha)
    
    sys.stderr.write("Populations split, divergence begining...\n")
    #run for K generations
    #every generation calculate alpha by comparing 1 individual from each population to one from the other
    for i in range(_E, _E+_K):
        popA.generation()
        popB.generation()
        alpha = popA.alpha(popB)
        
        popA.printGeneration(i, "A", alpha)

class Pop:
    def __init__(self, N = 100, alleles = [-0.02,-0.01,0.01,0.02], loci = 10, mu = 0.01, optimum = 0):
        self.alleles = alleles
        self.N = N
        self.numLoci = loci
        self.mu = mu
        self.optimum = optimum
        
        self.max = loci*max(alleles)
        self.min = loci*min(alleles)
        self.maxDiff = max([abs(self.optimum-self.max), abs(self.optimum-self.min)])
        
        self.meanFit = None
        self.meanNeut = None
        self.meanPheno = None
        
        self.inds = []
        
    def generate(self):
        totPheno = 0
        totNeut = 0
        totFitness = 0
        
        while len(self.inds) < self.N:
            newInd = Ind(genome = [], pop = self)
            totFitness += newInd.fitness()
            totPheno += newInd.pheno
            totNeut += newInd.neut
            self.inds.append(newInd)
            
        self.meanFit = totFitness/self.N
        self.meanNeut = totNeut/self.N
        self.meanPheno = totPheno/self.N
    
    def generation(self):
        newInds = []
        totPheno = 0
        totNeut = 0
        totFitness = 0
        
        tries = 0
        while len(newInds) < self.N:
            mom = random.randrange(0,self.N-1)
            dad = random.randrange(0,self.N-1)
            while mom == dad:#no selfing
                dad = random.randrange(0,self.N-1)
            
            mom = self.inds[mom]
            dad = self.inds[dad]
                
            off = mom.mate(dad)
            if off.fit > random.random():#selection!
                newInds.append(off)
                totFitness += off.fit
                totPheno += off.pheno
                totNeut += off.neut
                tries = 0
            else:
                tries += 1
                if tries > _N:#failed at making offspring N times, population is extinct
                    sys.stderr.write("Failed to make a viable offspring in N tries. Population extinct at generation %s.\n" % K)
    
        self.meanFit = totFitness/self.N
        self.meanNeut = totNeut/self.N
        self.meanPheno = totPheno/self.N
      
    def calcMeans(self):
        totPheno = 0
        totNeut = 0
        totFitness = 0
        
        for ind in self.inds:
            ind.fitness()
            totFitness += ind.fit
            totPheno += ind.pheno
            totNeut += ind.neut
        
        self.meanFit = totFitness/self.N
        self.meanNeut = totNeut/self.N
        self.meanPheno = totPheno/self.N
        
    def makeNewPop(self):
        newPop = Pop(N = self.N, alleles = self.alleles, loci = self.numLoci, mu = self.mu, optimum = self.optimum)
        while len(newPop.inds) < newPop.N:
            ind = random.choice(self.inds)
            newInd = Ind(genome = ind.genome[:], pop = newPop)
            newPop.inds.append(newInd)
        return newPop
        
    def alpha(self, other):
        myInd = random.choice(self.inds)
        otherInd = random.choice(other.inds)
        
        neutPoly, selPoly, polyList = self.polymorphism()
        neutDiv, selDiv = myInd.divergence(otherInd, polyList)
        
        if neutDiv == 0:
            sys.stderr.write("No neutral divergence.\n")
            a = "NA"
        elif selPoly == 0:
            sys.stderr.write("No selected polymorphism.\n")
            a = "NA"
        else:
            a = 1.0 - (float(selDiv)*neutPoly)/(float(neutDiv)*selPoly)
        return a
        
    def polymorphism(self):
        sites = [0]*(2*self.numLoci)
        
        ref = self.inds[0].genome
        for i in range(1, self.N):
            if sum(sites) == 2*self.numLoci:#all polymorphic
                break
            
            diffs = [0 if x == y else 1 for x, y in zip(ref, self.inds[i].genome)]
            sites = [1 if x+y != 0 else 0 for z, y in zip(sites, diffs)]
            
        sel = sum([sites[i] for i in range(0, 2*self.numLoci) if i % 2 == 0])
        neut = sum(sites)-sel
        
        return neut, sel, sites
        
    def printGeneration(self, K, pop, alpha):
        if K == 0 and pop == "I":
            print("Population, Generation, meanFitness, meanPhenotype, meanNeutral, alpha")
        print("%s,%s,%s,%s,%s,%s" % (pop, K, self.meanFit, self.meanPheno, self.meanNeut, alpha))
    
    def __str__(self):
        return "Pop(meanFit = %s, maxPheno = %s, minPheno = %s, optimum = %s)\n" % (self.meanFit, self.max, self.min, self.optimum)

class Ind:
    def __init__(self, genome = [], pop = None):
        #have N loci, that affect a quantitative trait
        self.pop = pop
        
        if genome:
            self.genome = genome
        else:
            self.genome = self.generate()
            
        self.fit = None
    
    def generate(self):
        newGenome = []
        for i in range(2*self.pop.numLoci):
            newGenome.append(random.choice(self.pop.alleles))
            
        return newGenome
    
    def fitness(self):
        if self.fit != None:
            return self.fit
        
        pheno = sum([self.genome[i] for i in range(0, 2*self.pop.numLoci) if i % 2 == 0])#only even positions count for fitness
        
        self.pheno = pheno
        self.neut = sum(self.genome)-pheno
        
        self.fit = (1 - abs(self.pop.optimum - pheno)/self.pop.maxDiff)**_t
        return (self.fit)
    
    def mutate(self):
        for i in range(2*self.pop.numLoci):
            if random.random() > self.pop.mu:
                self.genome[i] = random.choice(self.pop.alleles)
                
        if self.fit != None:
            self.fit = None
    
    def mate(self, other):
        #free recombination
#        for i in 2*range(self.pop.numLoci):
#            parent = random.randrange(0,2)
#            if parent:
#                newGenome.append(self.genome[i])
#            else:
#                newGenome.append(other.genome[i])
        #exactly one recombination point
        point = random.randrange(0,self.pop.numLoci*2-1)
        newGenome = self.genome[0:point]+other.genome[point:]
                
        off = Ind(genome=newGenome, pop = self.pop)
        off.mutate()
        off.fitness()
        return off
    
    def divergence(self, other, polyList):
        neutDiv = 0
        selDiv = 0
        for i in range(0, 2*self.pop.numLoci):
            if self.genome[i] != other.genome[i] and not polyList[i]:#dont count it if it is polymorphic in my population
                if i % 2 == 0:
                    selDiv += 1
                else:
                    neutDiv += 1
                    
        return neutDiv, selDiv
    
    def __str__(self):
        return "Ind(%s, fit=%s)\n"%(self.genome, self.fit)  

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"N:L:u:K:E:o:t:h")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-L":
            global _L 
            _L = int(arg)
        elif opt == "-K":
            global _K 
            _K = int(arg)
        elif opt == "-E":
            global _E 
            _E = int(arg)
        elif opt == "-u":
            global _u 
            _u = float(arg)
        elif opt == "-h":
            details()
            sys.exit()
        elif opt == "-o":
            global _o 
            _o = float(arg)
        elif opt == "-t":
            global _t 
            _t = float(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("-N %s -L %s -u %s -K %s -E %s -o %s -t %s\n" % (_N, _L, _u, _K, _E, _o, _t))
   
use = "python "+__file__.split("/")[-1]+" [options]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("h - None - NA - prints help")
    print("N - INT - %s - population size" % _N)
    print("L - INT - %s - number of loci affecting phenotype" % _L)
    print("E - INT - %s - number of generations population is allowed to evolve and reach equilibrium before it is split" % _E)
    print("K - INT - %s - number of generations of evolution after population split" % _K)
    print("u - FLOAT - %s - per site mutation probability" % _u)
    print("o - FLOAT - %s - optimum phenotype" % _o)
    print("t - FLOAT - %s - exponent in fitness function. Bigger number makes selection stronger." % _t)

if __name__ == "__main__":
    start = time.time()
    __main__()
    sys.stderr.write("Run time (min): "+str((time.time()-start)/60.0)+"\n")