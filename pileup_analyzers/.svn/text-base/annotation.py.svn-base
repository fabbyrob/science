'''
A parser for annotation files. Should make a reader that will return one entire gene, consisting of a list of exons.


    >>> import annotation
    >>> annotation_reader = annotation.Reader(open('example.annot', 'rb'))
    >>> for record in annotation_reader:
    >>>     print record
    
    Gene.name
    Gene.scaf
    Gene.start
    Gene.end
    Gene.exons
    Gene.direction
'''

import sys

class Gene(object):
    def __init__(self, name, scaf, start, end, direction, exons):
        self.name = name
        self.scaf = scaf
        self.start = start
        self.end = end
        self.direction = direction
        self.exons = exons #list of tuples
        
    def sortExons(self):
        self.exons.sort()
        self.start = self.exons[0][0]
        self.end = self.exons[-1][1]
        
    def addExon(self, newExon):
        self.exons.append(newExon)
        
    def makeIntrons(self):
        introns = []
        for i, exon in enumerate(self.exons):
            if i == len(self.exons)-1:
                break
            intronStart = exon[1]+1#one after the end of the exon
            intronEnd = self.exons[i+1][0]-1#one before the start of the next exon
            introns.append((intronStart, intronEnd))
        self.introns = introns
        return self.introns
    
    #return the within-gene index of the given genomic position
    #returns None if it's not within an exon
    def getIndex(self, scaf, pos):
        if scaf != self.scaf:
            return None
        
        if not (pos >= self.start and pos <= self.end):
            return None
        
        i = 0
        
        for start, end in self.exons:
            if pos < start:
                return None
            
            if pos >= start and pos <= end:
                i += pos - start
                return i
            else:
                i += end - start
                
        return None
        
        
    def __lt__(self, other):
        if isinstance(other, Gene):
            return self.start < other.start
        return False
    
    def __gt__(self, other):
        if isinstance(other, Gene):
            return self.start > other.start
        return False
    
    def __eq__(self, other):
        if isinstance(other, Gene):
            return self.start == other.start
        return False
        
    def __str__(self):
        return "Gene(name="+self.name+", scaf="+str(self.scaf)+", start="+str(self.start)+", end="+str(self.end)+", num_exons="+str(len(self.exons))+")"
        
class Reader(object):
    def __init__(self, filename):
        self.reader = filename
        self.myFile = None
        
    def __iter__(self):
        myGene = None
        for line in self.reader:
            line = line.rstrip()
            sline = line.split()
            
            if not (len(sline) == 5 or len(sline) == 6):
                continue
            
            scaf = sline[0]
            start = int(sline[1])
            end = int(sline[2])
            name = sline[3]
            direction = sline[4]       
            
            if myGene == None:
                myGene = Gene(name, scaf, start, end, direction, [(start, end)])
            elif myGene.name == name:
                myGene.addExon((start, end))
            else:
                yield (myGene)
                myGene = Gene(name, scaf, start, end, direction, [(start, end)])   
   

if __name__ == "__main__":
#TEST
    myReader = Reader(open('io/annot.txt','rb'))

#for g in myReader:
#    print(g)
    
    for gene in myReader:
        print(gene)
    
