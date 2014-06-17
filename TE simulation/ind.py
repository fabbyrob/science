import random
import sys
sys.path.append("/opt/python27/lib")
sys.path.append("/usr/lib64/python2.6/site-packages/") 
from math import pow, floor
from numpy.random import poisson
from collections import Counter

class Individual:
    def __init__(self, teGenome = [], silencingGenome = [], loci = 10, TEnum = 2,  initSilence = [0.1], transRate = 0.01, maxClasses = 1, classMutation = 0, silenceMutation = 0.01, epsilon = 0.01, excision = 0.01):
        if(teGenome):
            self.teGenome = teGenome
            self.silencingGenome = silencingGenome
        else:#if no genome is provided, generate a random one
            self.generate(loci, TEnum, initSilence)
            
        #mutate our silencing rates
        self.silenceMutation = silenceMutation
        self.silenceEpsilon = epsilon
            
        self.loci = loci
            
        self.maxClasses = maxClasses
        self.classMutation = classMutation

        #caclulate the transposition rate of each category
        self.maxRate = transRate

        self.excision = excision
        
        self.calculateTransposition()#calculate transposition rates and silencing rates BEFORE mutation has occurred
        
        self.fitness = -1
        self.k = 0#number of gametes contributed to the next generation

    #generates a random genome
    def generate(self, loci, numTEs, silencing):
        genome = [[0]*loci, [0]*loci]#makes two chromosomes each of length 'loci'
        numTEs = poisson(numTEs)#samples the number of TEs from a poisson, with mean 'numTEs'
        
        while sum(genome[0]+genome[1]) < numTEs:#picks random locations in the genome to insert those TEs
            chr = random.randint(0,1)
            locus = random.randint(0,loci-1)
            genome[chr][locus] = 1
            
        #generate genome of silencing
        silencingLoci = [[],[]]
        for n in silencing:#for each solencing rate
            a1 = min(0.5, max(0, random.gauss(n/2.0, 0.05*n)))#sample a random number for allele 1 from a gauss
            a2 = min(0.5, max(0, random.gauss(n/02.0, 0.05*n)))#sampe a random number for allele 2
            silencingLoci[0].append(a1)
            silencingLoci[1].append(a2)
            
        #store both the TE numbers and teh silencing numbers
        self.teGenome = genome
        self.silencingGenome = silencingLoci
    
    #calculates the silencing rate and transposition rate of each class of TEs
    def calculateTransposition(self):
        #caclulate our total silencing for each category
        self.silence = []
        for i in range(len(self.silencingGenome[0])):
            self.silence.append(self.silencingGenome[0][i]+self.silencingGenome[1][i])
            
        #calculate transposition for each class based on total silencing
        self.transposeRate = []
        for i in range(0, len(self.silence)):
            self.transposeRate.append(self.maxRate*(1-self.silence[i]))
    
    #mutate silencing rate
    def mutateSilencing(self):
        for i, r in enumerate(self.silencingGenome[0]):#for each silencing locus
            prob = random.random()#filp a coin to mutate allele 1
            if prob < self.silenceMutation:#if we do mutate
                dir = random.choice([1,-1])#decide to eitehr increase or decrease silencing
                self.silencingGenome[0][i] = max(0, min(0.5, (r + dir*self.silenceEpsilon)))#modify silencing by +/- the silencing mutation magnitude.
                #one silencing allele can never give less than 0 silencing, or more than 0.5 silencing
                              
            r = self.silencingGenome[1][i]                
            prob = random.random()#flip a coin to mutate allele 2
            if prob < self.silenceMutation:
                dir = random.choice([1,-1])
                self.silencingGenome[1][i] = max(0, min(0.5, (r + dir*self.silenceEpsilon)))
                
        self.calculateTransposition()#recalculate new transposition and silencing rates POST mutation
       
        
    #transpose TEs within the genome
    def transpose(self):
        Ns = Counter(self.teGenome[0]+self.teGenome[1])#count the number of TEs in each class
        Ls = [0]*len(self.silence)#lambdas for poissons for each class
        newTEs = [0]*len(self.silence)#the number of new TEs we should add
        
        for k in Ns.keys():#for each class of TEs
            if k == 0:
                continue
            Ls[k-1] = Ns[k]*self.transposeRate[k-1]#lambda = numTEs*transposition rate
            newTEs[k-1] = poisson(Ls[k-1])#sample the number of new TEs from a poisson

        newGenome = [self.teGenome[0][:],self.teGenome[1][:]]#copy the old genome
        
        for i, teClass in enumerate(newTEs):#for every class of TEs
            while teClass > 0:#for each new insertion
                #Check to see if the new insertion will mutate to a new class
                if random.random() > self.classMutation and i < self.maxClasses:
                    myclass = i+1#if it does mutate change it to a 'higher' class
                else:
                    myclass = i
                    
                chr = random.randint(0,1)#pick one of the two chromosomes to transpose to
                
                sites = [j for j,x in enumerate(newGenome[chr]) if x == 0]#get a list of all the sites on this chrom that have no TEs
                
                if len(sites) == 0:#if there are no sites on this chr
                    chr = (chr-1)*-1#use the other chromosome
                    sites = [j for j,x in enumerate(newGenome[chr]) if x == 0]
                    
                    if len(sites) == 0:#no sites on the other chromosome too
                        break#give up transposing there are no sites left
                
                locus = random.choice(sites)#pick a random locus to transpose to

                if newGenome[chr][locus] == 0:#if there is no TE at that location
                    newGenome[chr][locus] = myclass#transpose into that site
                    
                teClass -= 1#decrease the number of transpositions remaining
        
        #replace this individual's genome with it's new genome, post transposition
        self.teGenome = newGenome
        
        
        return(newGenome)

    #excise exons
    def excise(self):
        if sum(self.teGenome[0])+sum(self.teGenome[1]) == 0:#there are no TEs, don't bother trying
            return
        #pick the number of TEs that should be excised from a poisson
        #where lambda = excistion rate * number of TEs
        numToExcise = poisson(self.excision*(2*self.loci - (self.teGenome[0].count(0)+self.teGenome[1].count(0))))
        for i in range(numToExcise):
            chr = random.randint(0,1)#pick a random chromosome
            
            sites = [i for i,x in enumerate(self.teGenome[chr]) if x != 0]#get a list of all the sites on this chrom that ARE NOT TEs
            
            if len(sites) == 0:#if there are no sites
                chr = (chr-1)*-1#use other chrom
                sites = [i for i,x in enumerate(self.teGenome[chr]) if x != 0]
                    
                if len(sites) == 0:#no sites on the other chromosome too
                    break#give up on excising, there are no TEs
            
            site = random.choice(sites)
            
            self.teGenome[chr][site] = 0

    #calculates fitness
    def calcFitness(self, selection = (0.001, 0.001), t = 1, mode = 'ALL'):
        #calculate 'g' the number of deleterious TEs
        if mode == 'ALL':#if all TEs are bad
            g = float(self.loci*2-(self.teGenome[0].count(0)+self.teGenome[1].count(0)))#/(2*self.loci)
        elif mode == 'HOMO':#only homozygous TEs are bad
            g = 0
            for i in range(0, self.loci):
                if self.teGenome[0][i] != 0 and self.teGenome[0][i] == self.teGenome[1][i]:
                    g += 1
            g = float(g)#/(self.loci)
        elif mode == 'HET':#only heterozygous TEs are bad
            g = 0
            for i in range(0, self.loci):
                if self.teGenome[0][i] != self.teGenome[1][i]:
                    g += 1
            g = float(g)#/(self.loci)
        else:
            sys.stderr.write("Invalid fitness mode: "+mode+"\n")
            return -1
    
        #fitness = 1 - x*g^t - x_silencing*silencing_rate
        fitness = 1
        fitness -= selection[0]*pow(g, t)
        fitness -= selection[1]*(sum(self.silence))
        fitness = max(0, fitness)#makes sure that fitness is always at least 0
        self.fitness = fitness
        return fitness
    
    #mate this individual with 'other' to produce an offspring
    def mate(self, other, numCrosses = 2):
        myCross = self.cross(numCrosses)#cross over my genome and produce gametes
        otherCross = other.cross(numCrosses)#cross over partner's genome and produce gametes
        
        myChr = random.randint(0,1)#randomly pick one of my two gametes
        otherChr = random.randint(0,1)#randomly pick on of my partners two gametes
        
        self.k += 1#increase my contribution to next gen
        other.k += 1#increase my partners contribution
        
        tes = [myCross[0][myChr], otherCross[0][otherChr]]#take my 'teGenome' and my partners'teGenome' from our gametes and combine them
        silencing = [myCross[1][myChr], otherCross[1][otherChr]]#combine our 'silencingGenomes as well
        
        #make a new individual with these te and silencing genomes and my other parameters
        offspring = Individual(teGenome = tes, silencingGenome = silencing, loci = self.loci, transRate = self.maxRate, maxClasses = self.maxClasses, classMutation = self.classMutation, silenceMutation = self.silenceMutation, epsilon = self.silenceEpsilon)
        return offspring
    
    #crossover my chromosomes and produce gametes, assumes the silencing loci are in the middle of the TE loci
    def cross(self, crosses):
        #fuse the two things together
        midPoint = int(floor(len(self.teGenome[0])/2))#determines the midpoint of the chromosome
        g1 = self.teGenome[0][0:midPoint]+self.silencingGenome[0]+self.teGenome[0][midPoint:]#fuses the silencing loci into the middle of the chromosome
        g2 = self.teGenome[1][0:midPoint]+self.silencingGenome[1]+self.teGenome[1][midPoint:]
        
        #crossover
        loci = len(g1)
        numCrosses = poisson(crosses)#samples a number of crossovers to performs
        for i in range(0, numCrosses):
            point = random.randint(0, loci-1)#pick a random crossover point
            #crossover
            temp = g1[0:point]+g2[point:]
            g1 = g2[0:point]+g1[point:]
            g2 = temp
        
        #yank them apart again
        crossedTEs = [g1[0:midPoint]+g1[midPoint+len(self.silencingGenome[0]):], g2[0:midPoint]+g2[midPoint+len(self.silencingGenome[0]):]]
        crossedSilence = [g1[midPoint:midPoint+len(self.silencingGenome[0])], g2[midPoint:midPoint+len(self.silencingGenome[0])]]
        
        return [crossedTEs, crossedSilence]
    
    
    #DEPRICATED
#    def crossSilence(self, stdev):
#        silence = random.gauss(self.silence, stdev)
##        if silence < 0:
##            silence = 0
##        elif silence > 1:
##            silence = 1
#        self.crossedSilence = silence
        
#    def genomeSize(self, genome):
#        return len(sum(genome, []))
       
#    def mutate(self, genome, rate = 0.001):
#        mutGenome = genome
#        size = self.genomeSize(genome)
#        mean = size*rate
#        stdev = 2*mean
#        numMut = int(random.gauss(mean, stdev))
#        
#        while(numMut > 0):
#            locus = random.randint(0,len(genome)-1)
#            site = random.randint(0,len(genome[locus])-1)
#            if (genome[locus][site]):
#                mutGenome[locus][site] = 0
#            else:
#                mutGenome[locus][site] = 1
#                
#            numMut -= 1
#        return mutGenome
            
    def shorter(self, other):
        if (len(self.teGenome)>len(other.teGenome)):
            return len(other.teGenome)
        return len(self.teGenome)
    
    def normalizeFitness(self, meanFitness):
        self.fitness = self.fitness/meanFitness
        return self.fitness

    #counts the number of TEs in my genome
    def teNum(self):
        return 2*self.loci-(self.teGenome[0].count(0)+self.teGenome[1].count(0))
    
    #calculates the difference between my two silencing alleles at each silencing locus
    def silenceDiff(self):
        diffs = []
        for i in range(0, len(self.silence)):
            diffs.append(abs(self.silencingGenome[0][i]-self.silencingGenome[1][i]))
        return diffs
    
    def __lt__(self, other):
        if isinstance(other, Individual):
            return self.fitness < other.fitness
        return False
    
    def __gt__(self, other):
        if isinstance(other, Individual):
            return self.fitness > other.fitness
        return False
    
    def __eq__(self, other):
        if isinstance(other, Individual):
            return self.fitness == other.fitness
        return False

    def __str__(self):
        return "Individual(TEs="+str(self.teGenome)+", silencing="+str(self.silence)+", fitness="+str(self.fitness)+")"
 
if __name__ == "__main__":
    ind1 = Individual(loci=9, transRate = 1, silenceMutation=0.5, silencingGenome=[[0.2,0.3],[0.1,0.1]], teGenome=[[0,1,1,0,1,1,0,1,1],[0,0,0,0,0,0,0,0,0]])
    ind1.calcFitness(selection = (.1,.1), t=2)
    print(ind1)
    print(ind1.silencingGenome)
    print(ind1.silenceDiff())
    ind1.mutateSilencing()
    print(ind1.silencingGenome)
    print(ind1.silenceDiff())
#    ind1.teGenome = ind1.transpose()
#    ind1.calcFitness(selection = (.1,.1), t=2)
#    print(ind1)
#    print(ind1.calcFitness(selection=(0.1,0.2)))
#    print(ind1.calcFitness(mode = 'HOMO',selection=(0.1,0.2)))
#    print(ind1.calcFitness(mode = 'HET',selection=(0.1,0.2)))
#    inds = []
#    print("")
#    ind1 = Individual()
#    print(ind1)
#    for i in range(0,20):
#        print("")
#        inds.append(Individual())
#        print(inds[i].calcFitness())
#        print(inds[i].calcFitness(mode = 'HOMO'))
#        print(inds[i].calcFitness(mode = 'HET'))
#        print(inds[i])
#        print(inds[i].mate(ind1))

    
    
    
