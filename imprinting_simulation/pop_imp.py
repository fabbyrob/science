from ind_imp import Individual
import random

class Population:
    def __init__(self, N = 100, r = 0):
        self.inds = []
        self.N = N
        self.r = r
        
        while len(self.inds) < N:
            self.inds.append(Individual())
            
    def generation(self, g = 0):
        births = []
        numBirths = 0
        numOffs = 0
        while len(births) < self.N:
            #do 2*N matings
            for i in range(2*self.N):
                mom = random.choice(self.inds)
                coin = random.random()
                if coin <= self.r:#if we self then same parent
                    dad = mom
                else:
                    dad = random.choice(self.inds)
                mom.mate(dad)
            
            #birth offspring
            for mom in self.inds:
                mbirths = mom.birth()
                numOffs += len(mbirths)
                numBirths += 1
                births += mbirths
            
        #randomly pick N to survive
        newInds = random.sample(births, self.N)
        
        meanAllelesMat, meanAllelesPat = self.getMeanAlleles()
        meanHetMat, meanHetPat = self.getMeanHets()
        if g == 0:
            print("Generation,meanMaternal,meanPaternal,meanMhet,meanPhet,meanBirths")
        print(str(g)+","+str(meanAllelesMat)+","+str(meanAllelesPat)+","+str(meanHetMat)+","+str(meanHetPat)+","+str(numOffs/float(numBirths)))
        
        #replace out pop
        self.inds = newInds
        
    def getMeanAlleles(self):
        t1 = 0.0
        t2 = 0.0
        for i in self.inds:
            t1 += i.genome[0][0]
            t2 += i.genome[1][1]
            
        return (t1/self.N, t2/self.N)
    
    def getMeanHets(self):
        h1 = 0.0
        h2 = 0.0
        for i in self.inds:
            h = i.heterozygosity()
            h1 += h[0]
            h2 += h[1]
        
        return (h1/self.N, h2/self.N)