import random
class Individual:
    def __init__(self, genome = None):
        if genome:
            self.genome = genome
        else:
            self.genome = []
            bases = ["A","T","G","C"]
            for i in range(10):
                self.genome.append(random.choice(bases))
        
        self.fitness = self.calcFit()
    
    def calcFit(self):
        cgs = 0
        for b in self.genome:
            if b == "C" or b == "G":
                cgs += 1
        return float(cgs)/len(self.genome)
    
    def mutate(self):
        bases = ["A","T","G","C"]
        for i in range(len(self.genome)):
            p = random.random()
            if p < 0.1:
                self.genome[i] = random.choice(bases)
        self.fitness = self.calcFit()
        
    def crossover(self, other):
        offspring1 = self.genome[0:(len(self.genome)/2)]+other.genome[(len(other.genome)/2):]
        offspring2 = other.genome[0:(len(other.genome)/2)]+self.genome[(len(self.genome)/2):]
        
        p = random.random()
        if p < 0.5:
            return Individual(offspring1)
        else:
            return Individual(offspring2)
        
    
    def __str__(self):
        return "Individual(genome="+str(self.genome)+", fitness="+str(self.fitness)+")"