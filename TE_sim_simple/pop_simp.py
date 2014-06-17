import sys
import random
import ind_simp
import gc

class Population():
    def __init__(self, u = 0.05, v = 0.05, crosses = 2, loci = 10, startTEs = [2], size = 50, selection = 0.001, t = 1.5, r = 0.1, mode = "ALL", classes = 1, silenceU = 0.0, silenceMag = 0.0, silStart =[0], selectionSilence = 0.001, classU = 0.0, verbose = False, selType = "R", silIntro = None, relative = False):
        self.verbose = verbose
        self.verboseOutputter("__init__")
        
        self.crosses = crosses
        self.u = u
        self.v = v
        self.size = size
        self.selection = selection
        self.t = t
        self.r = r
        self.mode = mode
        self.classes = classes
        self.silenceMag = silenceMag
        self.loci = loci
        self.startTEs = startTEs
        self.relative = relative
        
        self.silIntro = silIntro
        if not self.silIntro:
            self.silenceU = silenceU
            self.silStart = [[s/2.0]*2 for s in silStart]
            self.silLater = silStart
        else:
            self.silenceU = 0
            self.silStart = [[0,0]*self.classes]
            self.silLater = silStart
        
        self.selectionSilence = selectionSilence
        self.classU = classU
        
        self.selType = selType
        
        self.firstGeneration()
    
    def verboseOutputter(self, function = "verboseOutputter", other = ""):
        if self.verbose:
            sys.stderr.write("Population.%s()\t%s\n" % (function, other))
    
    def firstGeneration(self):
        self.verboseOutputter("firstGeneration")
        self.inds = []
        
        while len(self.inds) < self.size:
            ind = ind_simp.Individual(population = self)
            ind.transpose()
            ind.excise()
            ind.mutateSilence()
            ind.calcFitness()
            self.inds.append(ind)
            
        self.fit = self.meanFitness()
        self.tes = self.meanTEs()
        self.silencing = self.meanSilencing()
        self.Ne = self.calcNe()
        
        self.printGeneration(g = 0)
      
    def randomParents(self):
        p1 = self.inds[random.randint(0, self.size-1)]
        
        if random.random() < self.r:
            p2 = p1
        else:
            p2 = self.inds[random.randint(0, self.size-1)]
            while p2 is p1:
                p2 = self.inds[random.randint(0, self.size-1)]
        
        off = p1.mate(p2)
        off.calcFitness()
        val = random.random()
        offFit = off.fitness
        if self.relative:
            offFit = offFit/self.maxFitness
        
        if val < offFit:
            off.transpose()
            off.excise()
            off.mutateSilence()
            p1.k += 1
            p2.k += 1
            return off
        elif self.verbose:
            sys.stderr.write("Individual's fitness was too low. \n\tNumber=%s \n\t%s.\n" % (val, off))
            
        return None
      
    def rouletteParents(self):
        r = random.random()
        p1 = self.pickParent()
            
        selfing = random.random()
        if selfing < self.r:
            p2 = p1
        else:
            p2 = self.pickParent(p1)
                    
        p1 = self.inds[p1]
        p2 = self.inds[p2]
        off = p1.mate(p2)
        p1.k += 1
        p2.k += 1
        return off
    
    def pickParent(self, other = None):
        r = random.random()
        p = 0
        tot = self.inds[0].calcFitness()/self.totFitness
        while tot < r:
            tot += self.inds[p].calcFitness()/self.totFitness
            p += 1
            
        if p >= self.size:
            p = self.size-1
            
        while p == other:
            p = self.pickParent()
            
        return p
        
    def generation(self, g):
        self.verboseOutputter("generation", g)
        newPop = []
        fails = 0
        attempts = 0
        if self.selType == "O":
            self.inds.sort()
        
        while len(newPop) < self.size:
            attempts += 1
            if fails > self.size*2:
                sys.stderr.write("%s offspring have been inviable, population meltdown. Exiting.\n" % fails)
                sys.exit(0)
            
            if attempts % 100 == 0:
                gc.collect()
            
            if self.selType == "R":
                off = self.randomParents()
                if not off:
                    fails += 1
                else:
                    newPop.append(off)
            elif self.selType == "O":
                off = self.rouletteParents()
                newPop.append(off)
        
        if g == self.silIntro:#time to introduce silencing
            self.verboseOutputter("generation", "introducing silencing")
            ind = random.choice(newPop)
            ind.silence = [[x,0] for x in self.silLater]
                
        self.inds = newPop
                
        self.fit = self.meanFitness()
        self.tes = self.meanTEs()
        self.silencing = self.meanSilencing()
        self.Ne = self.calcNe()
        
        self.printGeneration(g = g)
        
        self.verboseOutputter("generation", str(self.inds[0]))
        
        if self.silIntro and g > self.silIntro:
            if self.freqs[0] > 0.99 or self.freqs[0] == 0:
                return False
        return True
            
    def meanFitness(self):
        self.verboseOutputter("meanFitness")
        tot = 0
        maxFit = 0
        
        for ind in self.inds:
            fit = ind.calcFitness()
            tot += fit
            maxFit = max(maxFit, fit)
            
        self.totFitness = tot
        self.maxFitness = maxFit
        return float(tot)/self.size
    
    def meanTEs(self):
        self.verboseOutputter("meanTEs")
        tots = [0]*self.classes
        
        for ind in self.inds:
            myTEs = ind.numTEs()
            tots = [tot + mine for tot, mine in zip(tots, myTEs)]
            
        tots = [float(tot)/self.size for tot in tots]
            
        return tots
    
    def meanSilencing(self):
        self.verboseOutputter("meanSilencing")
        tots = [0]*self.classes
        freqs = [0]*self.classes
        for ind in self.inds:
            mySil = ind.totalSilencing()
            tots = [tot + mine for tot, mine in zip(tots, mySil)]
            myFreq = ind.silenceFrequencies()
            freqs = [freq + mine for freq, mine in zip(freqs, myFreq)]
        
#        if self.silIntro:
#            freqs = [float(tots[i])/float(self.silLater[i])/self.size for i in range(self.classes)]
#            self.freqs = freqs
            
        self.freqs = [float(tot)/(2*self.size) for tot in freqs]
        tots = [float(tot)/self.size for tot in tots]
        return tots
        
    def calcNe(self):
        meanK = 2
        
        #variance in contribution
        totSqDiffK = sum([abs(ind.k-meanK)**2 for ind in self.inds])
        
        varK = totSqDiffK/self.size
        
        Ne = float(self.size*4)/(varK+2)
        return Ne
        
        
    def printGeneration(self, g):
        self.verboseOutputter("printGeneration", g)
        for i in range(0, self.classes):
            if self.silIntro:
                print(str(g)+","+str(i+1)+","+str(self.fit)+","+str(self.tes[i])+","+str(self.freqs[i])+",0,"+str(self.Ne))#the zeroes are silencing diff, and Ne respectively
            else:
                print(str(g)+","+str(i+1)+","+str(self.fit)+","+str(self.tes[i])+","+str(self.silencing[i])+",0,"+str(self.Ne))#the zeroes are silencing diff, and Ne respectively
    
    def __str__(self):
        myStr = "meanFitness = "+str(self.fit)+" meanTEnum = "+str(self.tes)+" meanSilence = "+str(self.silencing)
        myStr += "\n__________________\n"
        for ind in self.inds:
            myStr += str(ind)+"\n"
            
        myStr += "-----------------"
        return(myStr)
        
if __name__ == "__main__":
    pop = Population(size = 5, silStart = [1], selType = "O", silIntro = 2)
    print(pop)
    for i in range(1, 5):
        pop.generation(i)
        print(i)
        print(pop)
    
    print(pop)