class Parser:
    def __init__(self, file):
        self.file = open(file, "r")
        
    def __iter__(self):
        for line in self.file:
            yield Exon(line)
        
class Exon:
    def __init__(self, data):
        data = data.split()
        self.CHROM = data[0]
        self.START = int(data[1])
        self.END = int(data[2])
        self.GENE = data[3]
        self.DIRECTION = data[4]
        
    def __str__(self):
        return "Record("+self.GENE+")"
  
class Gene:
    def __init__(self, exons):
        self.exons = exons
        self.name = exons[0].GENE
        self.start = exons[0].START
        self.end = exons[-1].END
        self.direction = exons[0].DIRECTION
    
    def sortExons(self):
        return
    
    def __str__(self):
        return "Gene("+self.name+")"#+", "+str(self.start)+", "+str(self.end)+")"
    
myFile = "annot.txt"
myGene = ""
myExons = []
genes = []
for exon in Parser(myFile):
    if myGene != exon.GENE:
        if myGene != "":
            genes.append(Gene(myExons))
        myGene = exon.GENE
        myExons = [exon]
    else:
        myExons.append(exon)
        
for gene in genes:
    print(gene)