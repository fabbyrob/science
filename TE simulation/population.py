from ind import Individual
import random
import sys
from math import floor, ceil
from collections import Counter

class Population:
    def __init__(self, inds = [], size = 50, selfing = 0, silence = [0.1], numCrosses = 2, selection = (0.001,0.001), t = 1, mode = 'ALL', transRate = 0.01, maxClasses = 1, classMutation = 0, silenceMutation = 0.01, epsilon = 0.01, excision = 0.01):
        self.size = size#population size
        self.inds = inds#number of inds in the population
        self.selfing = selfing#population selfing rate
        
        self.silence = silence#initial silencing rates for each class
        self.numCrosses = numCrosses#mean number of crossovers per individual
        self.selection = selection#selection coefficients for both TEs and silencing
        self.t = t#shape parameter for fitness function
        self.mode = mode#selection mode (i.e. roulette vs. random vs. tournament)
        self.transRate = transRate#maximum transposition rate
        self.maxClasses = maxClasses#maximum number of classes
        self.classMutation = classMutation#mutation rate between classes of TEs
        self.silenceMutation = silenceMutation#mutation rate in silencing
        self.epsilon = epsilon#mutation magnitude in silencing
        self.excision = excision#excision rate
        
        while len(self.silence) < self.maxClasses:
            sys.stderr.write("Not enough initial silences provided. Assuming 0 for remainder.\n")
            self.silence.append(0.0)
  
        self.fitnesses = []#population's mean fitness each generation
        self.TEnums = []#population's mean # of TEs each generation
        self.silencenums = []#population's mean silencing rate each generation
  
        print("Generation\tClass\tFitness\tTEnums\tSilence\tNe")
    
    #initializes the first generation of the population
    def firstGen(self, loci = 10, numTEs = 5):
        #create all individuals
        while(len(self.inds) < self.size):
            newind = Individual(loci = loci, TEnum = numTEs, initSilence = self.silence, transRate = self.transRate, maxClasses = self.maxClasses, classMutation = self.classMutation, silenceMutation = self.silenceMutation, epsilon = self.epsilon, excision = self.excision)
            self.inds.append(newind)
            
        #calculate and store the mean population parameters
        fit = self.meanFitness()
        tenums = self.meanTEnum()
        silencenums = self.meanSilence()
        silencediffs = self.meanSilenceDiff()
        
        self.fitnesses.append(fit)
        self.TEnums.append(tenums)
        self.silencenums.append(silencenums)
        
        self.printGenerationInfo(0, fit, tenums, silencenums, silencediffs, -1)
        
    #function that performs all the steps of a single generation
    def generation(self, sel = 'O', gen = 0):
        self.parents = None#reset the parents used in the last generation
        newPop = []#this will contain our new population
        
        #do this while our new population is not yet full
        while(len(newPop) < self.size):
            if (sel != 'R'):#selection acting on parents, i.e. roulette or tournament selection
                p1 = self.pickParent(type = sel)
                if(random.random() < self.selfing):
                    p2 = p1
                else:
                    p2 = self.pickParent(other = p1, type = sel)
                    
                child = p1.mate(p2, self.numCrosses)#make parents mate
                child.mutateSilencing()#mutate the silencing rate
                child.transpose()#transpose in offspring
                child.excise()#excise in offspring
                newPop.append(child)#store offspring
            else:#selection acting on zygotes
                p1 = self.pickParent(type = "R")
                if(random.random() < self.selfing):
                    p2 = p1
                else:
                    p2 = self.pickParent(other = p1, type = "R")
                    
                child = p1.mate(p2, self.numCrosses)#parents mate
                fit = child.calcFitness(selection = self.selection, t = self.t, mode = self.mode)#calculate offspring fitness
                if (random.random() < fit):#flip a coin to see if we keep this offspring
                    #keep the offspring, transpose and excise, mutate silencing, then store it
                    child.mutateSilencing()
                    child.transpose()
                    child.excise()
                    newPop.append(child)
                else:
                    #offspring did not survive
                    #decrease the parent's offspring count 'k'
                    p1.k -= 1
                    p2.k -= 1
                    
        Ne = self.calcNe()#calculate Ne of parental population, based on variation in family size
        self.inds = newPop#replace the population with the new offspring we just created
        
        #calculate and store mean population parameters
        fit = self.meanFitness()
        tenums = self.meanTEnum()
        silencenums = self.meanSilence()
        silencediffs = self.meanSilenceDiff()
        self.fitnesses.append(fit)
        self.TEnums.append(tenums)
        self.silencenums.append(silencenums)
        
        self.printGenerationInfo(gen, fit, tenums, silencenums, silencediffs, Ne)
            
    #prints the summary info for one generation
    def printGenerationInfo(self, gen, fitness, tenums, silencenums, silencediffs, Ne):
        for i in range(0, len(silencenums)):
            print(str(gen)+","+str(i)+","+str(fitness)+","+str(tenums[i])+","+str(silencenums[i])+","+str(silencediffs[i])+","+str(Ne))
    
    #picks a parent, makes sre the parent chosen is not the same as 'other' (to avoid accidental selfing)
    def pickParent(self, other = None, type = 'O'):
        if (type == 'R'):#randomly choose parents
            p = random.randint(0, len(self.inds)-1)  
            p = self.inds[p]
            while (p is other):
                p = random.randint(0, len(self.inds)-1)
                p = self.inds[p]    
            return p
        elif (type == 'T'):#pick by tournament selection
            p = self.tournament()   
            while (p is other):
                p = self.tournament()  
            return p
        elif (type == 'O'):#pick by roulette
            p = self.roulette()   
            while (p is other):
                p = self.roulette()  
            return p
        elif (type == 'P'):#pick by truncation selection
            p = self.truncate()   
            while (p is other):
                p = self.truncate() 
            return p
    
    #tournament selection. Pick two random potential parents, the one with the higher fitness becomes the parent.
    def tournament(self):
        p1 = random.randint(0,len(self.inds)-1)
        p2 = random.randint(0,len(self.inds)-1)
        ind1 = self.inds[p1]
        ind2 = self.inds[p2]
        if ind1.fitness < ind2.fitness:
            return ind2
        else:
            return ind1

    #pick parents randomly, weighted by their fitness, such that higher fitness parents are more likely to be chosen
    def roulette(self):
        prob = random.random()
        total = 0
        self.inds.sort()
        for ind in self.inds:
            total += ind.fitness/self.totalFit
            if (total>=prob):
                return ind
    
    #pick a parent randomly from the parents with the 10% highest fitness
    def truncate(self, p=0.1):
        if not self.parents:
            sortedInds = sorted(self.inds)
            topPercent = int(len(self.inds)-ceil(p*len(self.inds)))
            self.parents = sortedInds[topPercent:]
        parent = random.randint(0,len(self.parents)-1)
        return self.parents[parent]
    
    #calculate the mean fitness of the population
    def meanFitness(self):
        total = 0
        maxFit = 0
        for i in self.inds:#for every individual
            fit = i.calcFitness(selection = self.selection, t = self.t, mode = self.mode)#calculate its fitness, based on the given selection parameters, shape parameter, and selection mode
            total += fit
            if fit > maxFit:
                maxFit = fit
        self.totalFit = total
        
        mean = float(total)/self.size#get the mean
            
        return mean
    
    #calculate the mean number of TEs per individual for each class of TEs
    def meanTEnum(self):
        total = [0]*self.maxClasses#we will store the total number of TEs of each class in the population in this list
        for i in self.inds:#for every individual
            nums = Counter(i.teGenome[0]+i.teGenome[1])#count the number of TEs of each class in this individual's genome

            for k in nums.keys():
                if k == 0:
                    continue

                total[k-1] += nums[k]#add the number of TEs to the total
        return map(lambda t: float(t)/self.size, total)#divide each total by the population size, return the list
    
    #calculate the mean silencing rate of each class of TEs 
    def meanSilence(self):
        total = [0]*self.maxClasses#make a list to store the population 'total' silencing rate
        for i in self.inds:#for each individual
            for j in range(0, len(total)):#for each class of TEs
                total[j] += i.silence[j]#add that class's silencing rate to total
            #total = list(sum(t) for t in zip(total, i.silence))# doesnt work beacuse min isnt 0
        return map(lambda s: float(s)/self.size, total)#divide the totals by the population size
    
    #calculate the mean difference in the silencing rate of an individual's two silencing alleles (~heterozygosity)
    def meanSilenceDiff(self):
        total = [0]*self.maxClasses
        for i in self.inds:
            diffs = i.silenceDiff()
            total = [sum(t) for t in zip(total, i.silence)]
        return map(lambda s: float(s)/self.size, total)
    
    #calculate the population's Ne, based on variation in family size
    def calcNe(self):
        meanK = 2#mean number of gametes contributed by each individual, will always be 2 because population size is constant
        
        #calculate the variance in gametic contribution
        totSqDiffK = 0.0
        
        for i in self.inds:
            totSqDiffK += abs(i.k-meanK)**2
            
        varK = float(totSqDiffK)/len(self.inds)
        
        Ne = float(len(self.inds)*4)/(varK+2)
        return Ne

        
#    def splitPop(self, sizes = [50,50], number = 2):
#        pops = []
#        parents = int(self.size/number)
#        random.shuffle(self.inds)
#        start = 0
#        end = 0
#        for i in range(0, number):
#            end = start + parents
#            subPop = self.inds[start:end]
#            start = end
#            pop = Population(subPop, sizes[i], self.selfing, self.mutation)
#            pop.generation()
#            pops.append(pop)
#        return pops
         
#    def merge(self, other, size = 50):
#        per = int(size/2)
#        random.shuffle(self.inds)
#        random.shuffle(other.inds)
#        
#        new = self.inds[0:per]
#        new += other.inds[per:size]
#        
#        pop = Population(inds = new, size = size, selfing = self.selfing, mutation = self.mutation)
#        pop.fitnesses.append(pop.meanFitness())
#        return pop
            
    def __str__(self):
        message = "__________________________________\n"
        message += str(self.size)+":"+str(self.fitnesses)+"\n"
        for i in self.inds:
            message += str(i)+"\n"
            
        message += "__________________________________"   
        
        return message
    
    def __len__(self):
        return len(self.inds)
 
def printPop(pop):
    stri = ""
    for ind in pop:
        stri += str(ind) + ", "
    print(stri)
    
if __name__=="__main__":
    pop = Population(size=5, maxClasses = 5, classMutation = 0.8, transRate = 0.9)
    pop.firstGen(loci = 4, numTEs = 2)
    
    print(pop)
    
    pop.generation()
    pop.inds.sort()
    printPop(pop.inds)
    print()
    pop.generation(gen = 1)
    pop.inds.sort()
    printPop(pop.inds)
    print()
    pop.generation(gen = 2)
    pop.inds.sort()
    printPop(pop.inds)
    print()
    pop.generation(gen = 3)
    pop.inds.sort()
    printPop(pop.inds)
    print()
    
#    print(pop)
#    print(pop.meanTEnum())
#    print(pop.meanSilence())
#    print("sorted:")
#    pop.truncate()
    

