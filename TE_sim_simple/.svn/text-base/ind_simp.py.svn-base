import sys
import random
sys.path.append("/opt/python27/lib")
sys.path.append("/usr/lib64/python2.6/site-packages/") 
from math import pow
from numpy.random import poisson
import pop_simp

class Individual:
    def __init__(self, population, genome = [], silence = []):#, loci = 10, u = 0.05, tes = 2, v = 0.05, silence = 0):
        self.population = population
        self.verboseOutputter("__init__")
        if genome:
            self.genome = genome
        else:
            self.generate()
            
        if silence:
            self.silence = silence
        else:
            self.silence = self.population.silStart
            
        self.fitness = None
        self.k = 0
    
    def verboseOutputter(self, function = "verboseOutputter", other = ""):
        if self.population.verbose:
            sys.stderr.write("Individual.%s()\t%s\n" % (function, other))
        
    def generate(self):
        self.verboseOutputter("generate")
        loci = self.population.loci
        tes = self.population.startTEs
        self.genome = [[0]*loci, [0]*loci]
        myTEs = [poisson(t) for t in tes]#for each class of TEs pick the starting number

        for i, t in enumerate(myTEs):#put each class on the genome
            myClass = i + 1
            for j in range(t):
                chr, locus = self.pickSite(self.genome)
                    
                if chr == None:#no empty sites left
                    break
                
                self.genome[chr][locus] = myClass
        
    def pickSite(self, genome, type = 0):
        #self.verboseOutputter("pickSite")
        chr = random.randint(0,1)
        sites = [k for k,x in enumerate(genome[chr]) if x == type]
        if len(sites) == 0:
            chr = (chr-1)*-1
            sites = [k for k,x in enumerate(genome[chr]) if x == type]
            
            if len(sites) == 0:
                #sys.stderr.write("No sites of type %s remain in the genome.\n" % type)
                return (None, None)
            
        locus = random.choice(sites)
        return (chr, locus)
        
        
    def transpose(self):
        self.verboseOutputter("transpose")
        tes = [self.genome[0].count(i)+self.genome[1].count(i) for i in range(1, len(self.silence)+1)]
        total = sum(tes)
        if total == 2*self.population.loci or total == 0:#all TEs all the time or no TEs, no need to transpose
            return
        
        newTEs = [poisson(t*(1.0-sum(self.silence[i]))*self.population.u) for i, t in enumerate(tes)]#get number of new TEs for each class
        
        if sum(newTEs) == 0:
            return
        
        newGenome = [self.genome[0][:], self.genome[1][:]]
        for i, t in enumerate(newTEs):
            myClass = i+1
            for j in range(t):
                chr, locus = self.pickSite(newGenome)
                
                if chr == None:#no empty sites left
                    break
                
                newGenome[chr][locus] = myClass
            
        self.genome = newGenome
    
    def excise(self):
        self.verboseOutputter("excise")
        tes = [self.genome[0].count(i)+self.genome[1].count(i) for i in range(1, self.population.classes+1)]
        total = sum(tes)
        if total == 0:#no TEs, no need to excise
            return
        
        removals = [poisson(t*self.population.v) for t in tes]#get number of removed TEs for each class
        
        if sum(removals) == 0:
            return
        
        newGenome = [self.genome[0][:], self.genome[1][:]]
        for i, t in enumerate(removals):
            myClass = i + 1
            for j in range(t):
                chr, locus = self.pickSite(newGenome, type = myClass)
                if chr == None:
                    break
                
                newGenome[chr][locus] = 0
            
        self.genome = newGenome
     
    def mutateSilence(self):
        self.verboseOutputter("mutateSilence")
        if self.population.silIntro:#don't mutate in -I mode!
            return
        
        newSilence = []
        for allele1, allele2 in self.silence:
            if random.random() < self.population.silenceU:
                allele1 = min(0.5, max(0, allele1+self.population.silenceMag*random.choice([1,-1])))#keep it between 0 and 0.5
                
            if random.random() < self.population.silenceU:
                allele2 = min(0.5, max(0, allele2+self.population.silenceMag*random.choice([1,-1])))#keep it between 0 and 0.5
                
            newSilence.append([allele1, allele2])
            
        self.silence = newSilence
                
        
    def calcFitness(self):
        if self.fitness != None:
            return self.fitness
        
        self.verboseOutputter("calcFitness")
        x = self.population.selection
        t = self.population.t
        mode = self.population.mode
        
        x_sil = self.population.selectionSilence
        
        if mode == "ALL":
            g = 2*self.population.loci - (self.genome[0].count(0)+self.genome[1].count(0))
        elif mode == "HET":
            g = 0
            for i in range(self.population.loci):
                if self.genome[0][i] != self.genome[1][i]:
                    g += 1
        elif mode == "HOMO":
            g = 0
            for i in range(self.population.loci):
                if self.genome[0][i] != 0 and self.genome[0][i] == self.genome[1][i]:
                    g += 2
        
        fitness = 1
        fitness -= x*pow(g,t)
        fitness -= x_sil*sum([allele1+allele2 for allele1, allele2 in self.silence])
        
        self.fitness = max(0, fitness)
        return self.fitness
        
#    def makeGamete(self, crosses = 2):
#        numCrosses = poisson(crosses)
#        
#        chr1 = self.genome[0][:]
#        chr2 = self.genome[1][:]
#        
#        for i in range(numCrosses):
#            point = random.randint(0,self.loci-1)
#            temp = chr1[0:point]+chr2[point:]
#            chr2 = chr2[0:point]+chr1[point:]
#            chr1 = temp
#        
#        if random.randint(0,1):
#            return chr1
#        else:
#            return chr2
        
    #pick N sites, then swap as appropriate
    def makeGamete(self):
        self.verboseOutputter("makeGamete")
        #if self.genome[0] == self.genome[1]:#both the same, crossover doesn't matter
        #    return self.genome[0][:]
        
        numCrosses = poisson(self.population.crosses)
        points = []
        numSil = self.population.classes
        
        for i in range(numCrosses):
            p = random.randint(0, self.population.loci-1+numSil)#bigger to include the silencing locus
            if p in points:
                points.remove(p)
            else:
                points.append(p)
              
        points.sort()
        
        sites = []
        last = 0
        chr = 0
        for p in points:
            sites += [chr]*(p-last)
            chr = (chr-1)*-1
            last = p
            
        sites += [chr]*(self.population.loci-last+numSil)#bigger for silencing locus

        chr1 = [0]*self.population.loci
        chr2 = [0]*self.population.loci
        
        silencingPoint = int((self.population.loci)/2)
        silencing1 = []
        silencing2 = []
        
        for i, chr in enumerate(sites):
             if i < silencingPoint:
                chr1[i] = self.genome[chr][i]
                chr2[i] = self.genome[(chr-1)*-1][i]
             elif i >= silencingPoint and i < silencingPoint + numSil:
                 myClass = i - silencingPoint
                 silencing1.append(self.silence[myClass][chr])
                 silencing2.append(self.silence[myClass][(chr-1)*-1])
             else:
                 i = i - numSil#- to account for the silencing locus displacement
                 chr1[i] = self.genome[chr][i]
                 chr2[i] = self.genome[(chr-1)*-1][i]
        
        if random.randint(0,1):
            return (chr1, silencing1)
        else:
            return (chr2, silencing2)
        
    def mate(self, other):
        self.verboseOutputter("mate")
        mine = self.makeGamete()
        others = other.makeGamete()
        
        return Individual(genome = [mine[0], others[0]], silence = zip(mine[1], others[1]), population = self.population)
        
    def numTEs(self):
        self.verboseOutputter("numTEs")
        return [self.genome[0].count(i)+self.genome[1].count(i) for i in range(1, self.population.classes+1)]
    
    def totalSilencing(self):
        self.verboseOutputter("totalSilencing")
        return [allele1+allele2 for allele1, allele2 in self.silence]
        
    def silenceFrequencies(self):
        freqs = []
        for allele1, allele2 in self.silence:
            freq = 0
            if allele1:
                freq += 1
            if allele2:
                freq += 1
            freqs.append(freq)
        return freqs
        
    def __str__(self):
        tes = [self.genome[0].count(i)+self.genome[1].count(i) for i in range(0, self.population.classes+1)]
        return "Individual::genome= %s\n\t tes = %s\n\t silence= %s\n\t fitness= %s" % (self.genome, tes, self.silence, self.fitness)
        
    def __lt__(self, other):
        return self.fitness < other.fitness
    
    def __gt__(self, other):
        return self.fitness > other.fitness
        
    
if __name__ == "__main__":
    p = pop_simp.Population(loci=20,u = 0.5, v= 0.25, startTEs = [10], mode = "HOMO", crosses = 1, silStart = [0.1], silenceMag = 0.05, silenceU = 0.5)
    ind = Individual(population = p)
    ind.calcFitness()
    print(ind)
    
    t = 0
    for i in range(100):
        ind = Individual(population = p)
        print(ind)
        nb4 = ind.numTEs()[0]
        ind.excise()
        naf = ind.numTEs()[0]
        ind.mutateSilence()
        print(ind)
        t += float(float(nb4-naf)/nb4)
        print("")
        
    print(t)
    t = t/100.0
    print("excision rate: %s\n" % t)
    
    ind.excise()
    print(ind)
    ind.transpose()
    print(ind)
    ind.calcFitness()
    print(ind)
    
    ind.mutateSilence()
    print(ind)
    print("")
    
    ind2 = Individual(population = p)
    ind2.calcFitness()
    print(ind)
    print(ind2)
    offspring = ind.mate(ind2)
    offspring.calcFitness()
    
    print("lt")
    print(ind < ind2)
    print("gt")
    print(ind > ind2)
    
    print(offspring)
    print("")
    offspring.excise()
    print(offspring)