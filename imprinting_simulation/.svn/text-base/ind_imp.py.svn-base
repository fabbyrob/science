import random

class Individual:
    def __init__(self, genome = []):
        if genome:
            self.genome = genome
            #2 tuples, one for each locus
            #first in tuple is mothers second is fathers
        else:
            self.genome = self.generate()
          
        self.mutate()
                
        self.growth = self.genome[0][0]+self.genome[1][1]#first locus from mom, second from dad
            
        self.offspring = []
        self.resources = 0
        return
    
    def generate(self):
        genes = []
        for i in range(2):
            genes.append((random.random(), random.random()))
            
        return genes
    
    def mutate(self, m = 0.1):
        newGenome = []
        for a1, a2 in self.genome:
            coin = random.random()
            if coin <= m:
                a1 = self.mutateAllele(a1)
            coin = random.random()
            
            if coin <= m:
                a2 = self.mutateAllele(a2)
            newGenome.append((a1, a2))
            
        self.genome = newGenome
    
    def mutateAllele(self, allele, e = 0.05):
        allele = allele + e*random.choice([-1,1])
        allele = min(1, max(0, allele))
        return allele
    
    def mate(self, other):
        genome = []
        for i, g in enumerate(self.genome):
            g = (random.choice(g), random.choice(other.genome[i]))
            genome.append(g)

        self.offspring.append(Individual(genome))
        return
    
    def birth(self):
        resources = 5.0
        born = []
        tgrowth = 0
        for b in self.offspring:
            tgrowth += max(0, b.growth)
            
        if tgrowth == 0:
            return []
            
        biggest = None
        for n in self.offspring:
            n.resources = max(0, resources*(n.growth/tgrowth))
            if n.resources > 1.0:
                born.append(n)
            if biggest == None or biggest.resources < n.resources:
                biggest = n
        
        if len(born) == 0 and biggest != None and biggest.resources > 0:
            born.append(biggest)
        
        return born
            
    def heterozygosity(self):
        h = []
        for a1, a2 in self.genome:
            h.append(abs(a1-a2))
        return h