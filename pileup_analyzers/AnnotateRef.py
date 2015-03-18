import sys
import re
import getopt
from Queue import PriorityQueue

#_log = __file__.split("/")[-1]+".log"#log file name

#does not distinguish between 2 and 3 fold, but it can later
codes = {'2fold':2, '3fold':2,'intergene':0 , 'intron':1, 'exon':2, '0fold':3, '4fold':4, '3utr':5, '5utr':6, 'istop':7, 'stop':8, 'unknown':9}

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"l:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    reference = open(sys.argv[1],"r")
    
    if (reference == None):
        print("Bad reference name: "+sys.argv[1])
        sys.exit()

    annotation = open(sys.argv[2],"r")  
        
    if (annotation == None):
        print("Bad annotation file name: "+sys.argv[1])
        sys.exit()
        
 #   logfile = open(_log,"w")
    
  #  if (logfile == None):
   #     print("Bad logfile name: "+_log)
    #    sys.exit()
       
    print(codes)
        
    #patterns to extract data from input files
    #1 = scaffold
    refpat = re.compile(">(\S+)")
    #1 = scaffold, 2 = start, 3 = end, 4 = model, 5 = +/-
    annotpat = re.compile("(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
    
    annot = {}
    #read in the entire annotation
    for line in annotation:
        line = line.rstrip()
        m = annotpat.match(line)
    
        if (m == None):
            sys.stderr.write("Non-standard line (annotation):\n\t"+line+"\n")
            continue
        else:
            scaf = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            name = m.group(5)
            dir = m.group(4)
            
            if scaf not in annot.keys():
                annot[scaf] = []
                gene = Gene(name, dir)
                gene.putExon((start, end))
                annot[scaf].append(gene)
            elif annot[scaf][-1].name == name:#same gene, new exon
                annot[scaf][-1].putExon((start,end))
            else:#new gene
                gene = Gene(name, dir)
                gene.putExon((start, end))
                annot[scaf].append(gene)
    
    #sort the genes in the annotation
    for k in annot.keys():
        annot[k].sort()
    
    #read the reference file...
    fullQueue = PriorityQueue() #items are: (baseNum, scaf, code, gene_name, gene_direction)
    exonQueue = [] #items are: (baseNum, base)
    scaf = ""
    
    for line in reference:
        line = line.rstrip()
        m = refpat.match(line)
        
        if (m == None):# we're at a line of bases
            if scaf not in annot.keys():
                sys.stderr.write("Missing scaffold from annotation: "+scaf+"\n")
                continue
            for base in line:
                #print (base)
                baseNum += 1
                if len(annot[scaf]) > 0:
                    gene = annot[scaf][0]
         #           print(gene.getStart())
          #          print(baseNum)
                    if baseNum < gene.getStart():#we havent reached a gene yet!
                        fullQueue.put((baseNum, scaf, base, codes['intergene'], 'N', '0'))
            #            print('not gene yet')
                    elif baseNum >= gene.getNextStart() and baseNum < gene.getNextEnd():#middle of the exon stick it on the to-process queue!
                        exonQueue.append((baseNum, base))
           #             print('exon')
                    elif baseNum == gene.getNextEnd():#end of the exon
                        sys.stderr.write(str(baseNum) + 'end of exon '+str(gene.name)+"\n")
                        exonQueue.append((baseNum, base))
                        annot[scaf][0].popExon()
                        if annot[scaf][0].empty():#end of the gene
                            sys.stderr.write('end of gene\n')
                            #sys.stderr.write("tiger "+str(baseNum)+"\n")
                            results = process(exonQueue, scaf, gene.name, gene.dir)
                            for i in results:
                                fullQueue.put(i)
                            fullQueue = printQueue(fullQueue)
                            exonQueue = []
                            annot[scaf].pop(0)#get rid of the gene in our list
                            
                            while len(annot[scaf]) > 0 and annot[scaf][0].getStart() <= baseNum:
                                sys.stderr.write("Gene start before the end of previous gene (skipping)! "+str(annot[scaf][0])+"\n")
                                annot[scaf].pop(0)#get rid of the gene in our list
       #                     sys.stderr.write("next: "+str(annot[scaf][0])+"\n")
                            
                    else:#after the start of the gene, but not within the next exon == intron
                        fullQueue.put((baseNum, scaf, base, codes['intron'], gene.name, '0'))
                        #print('intron')
                else:#no genes left on this scaffold, it must be intergenic!
                    fullQueue.put((baseNum, scaf, base, codes['intergene'], 'N', '0'))
                    #print('end of scaffold')
        else:#we're at a title line
            fullQueue = printQueue(fullQueue)
            scaf = m.group(1)
            baseNum = 0
    #print the leftovers
    fullQueue = printQueue(fullQueue)

def process(queue, scaf, name, dir):
    if dir == '-':
        queue.reverse()

    output = []
    
    for i in range(0, len(queue), 3):
        bases = queue[i:i+3]
        
        #make the codon into a string
        codon = ""
        for i in bases:
            codon += str(i[1])
       
        codon = codon.upper()

        #flip the codon if we're on a - strand gene
        if dir == "-":
            codon = sense(codon)
        
        #get the codes for the codons
        myCodes = degeneracy(codon)
        
        #add each base to the output
        c = 0
       # print(bases)
        for b in bases:
       #     print(b)
       #     print(myCodes, c)
            result = (b[0], scaf, b[1], myCodes[c], name, dir)
            output.append(result)
            c += 1
            
    #replace the last values with stop if they were flagged istop
    for i in [-1,-2,-3]:
        base = output[i]
        if base[3] == codes['istop']:
            output[i] = (base[0], base[1], base[2], codes['stop'],base[4], base[5])
            
    #reverse the output if we reversed it to start with
    if dir == '-':
        output.reverse()
            
    return output

#returns the opposite strand of a codon
def sense(codon):
    ncodon = ""
    for b in codon:
        if b == "A":
            ncodon += "T"
        elif b == "T":
            ncodon += "A"
        elif b == "G":
            ncodon += "C"
        elif b == "C":
            ncodon += "G"
        elif b == "N":
            ncodon += "N"

    return ncodon

def printQueue(queue):
    while not queue.empty():
        data = queue.get()
        
        base = str(data[0])
        scaf = data[1]
        ref = data[2]
        code = str(data[3])
        name = data[4]
        dir = data[5]
        
        print (scaf+"\t"+base+"\t"+ref+"\t"+name+"\t"+dir+"\t"+code)
        
    return queue

codon_table = '''
#These are DNA codons
#CODON BASE_1 BASE_2 BASE_3
#Ala
GCT 0fold 0fold 4fold
GCC 0fold 0fold 4fold
GCA 0fold 0fold 4fold
GCG 0fold 0fold 4fold
#Arg
#the first base here is ambiguous, it can be 2 fold (NGA or NGG) or it can be 0fold (CGT and CGC)
CGT 0fold 0fold 4fold *
CGC 0fold 0fold 4fold *
CGA 2fold 0fold 4fold
CGG 2fold 0fold 4fold
AGA 2fold 0fold 2fold
AGG 2fold 0fold 2fold
#Asn
AAT 0fold 0fold 2fold
AAC 0fold 0fold 2fold
#Asp
GAT 0fold 0fold 2fold
GAC 0fold 0fold 2fold
#Cys
TGT 0fold 0fold 2fold
TGC 0fold 0fold 2fold
#Gin
CAA 0fold 0fold 2fold
CAG 0fold 0fold 2fold
#Glu
GAA 0fold 0fold 2fold
GAG 0fold 0fold 2fold
#Gly
GGT 0fold 0fold 4fold
GGC 0fold 0fold 4fold
GGA 0fold 0fold 4fold
GGG 0fold 0fold 4fold
#His
CAT 0fold 0fold 2fold
CAC 0fold 0fold 2fold
#Ile
ATT 0fold 0fold 3fold
ATC 0fold 0fold 3fold
ATA 0fold 0fold 3fold
#Met - Start
#this one was just wrong before (last base was "exon")
ATG 0fold 0fold 0fold *
#Leu
#this first one is ambiguous 2fold (for NTA and NTG) and 0fold (for CTT and CTC)
TTA 2fold 0fold 2fold *
TTG 2fold 0fold 2fold *
CTT 0fold 0fold 4fold *
CTC 0fold 0fold 4fold *
CTA 2fold 0fold 4fold
CTG 2fold 0fold 4fold
#Lys
AAA 0fold 0fold 2fold
AAG 0fold 0fold 2fold
#Phe
TTT 0fold 0fold 2fold
TTC 0fold 0fold 2fold
#Pro
#this one was just wrong before (last base was "exon")
CCT 0fold 0fold 4fold *
CCC 0fold 0fold 4fold *
CCA 0fold 0fold 4fold *
CCG 0fold 0fold 4fold *
#Ser
#These are all 0fold because the two codon types (TCN and AGN) do not share bases 
TCT 0fold 0fold 4fold * 
TCC 0fold 0fold 4fold * 
TCA 0fold 0fold 4fold * 
TCG 0fold 0fold 4fold * 
AGT 0fold 0fold 2fold * 
AGC 0fold 0fold 2fold * 
#Thr
ACT 0fold 0fold 4fold
ACC 0fold 0fold 4fold
ACA 0fold 0fold 4fold
ACG 0fold 0fold 4fold
#Trp
#this one was just wrong before (last base was "exon")
TGG 0fold 0fold 0fold *
#Tyr
TAT 0fold 0fold 2fold
TAC 0fold 0fold 2fold
#Val
GTT 0fold 0fold 4fold
GTC 0fold 0fold 4fold
GTA 0fold 0fold 4fold
GTG 0fold 0fold 4fold
#STOP
TAA stop stop stop
TGA stop stop stop
TAG stop stop stop 
'''
_codonDict = None
def degeneracy(codon):
    #This should be dumped into an external file, not hardcoded
    global _codonDict
    if not _codonDict:
        _codonDict = {}
        for f in codon_table.split("\n"):
            #print(f)
            sf = f.split()
            #print(sf)
            if f.startswith("#") or not f:
                continue
            _codonDict[sf[0]] = [codes[sf[1]], codes[sf[2]], codes[sf[3]]]
    
    #print(_codonDict)

    #print codon
    #if a base is ambiguous then call it all exon, and move on
    if len(codon) < 3:
        res = []
        for i in range (0, len(codon)):
            res.append(codes['exon'])
        # ("Error: Codon must be length 3 - " + codon)
        return res
        #exit(0)

    if ("N" in codon):
        return (codes['unknown'],codes['unknown'],codes['unknown'])
    if codon in _codonDict:
        return _codonDict[codon]
    else:
        sys.stderr.write("Unknown codon %s\n" % codon)
        #sys.exit()
        return(codes['unknown'],codes['unknown'],codes['unknown'])
      
   
use = "python "+__file__.split("/")[-1]+" ReferenceGenome.fa Annotation.txt"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)


class Gene:
    def __init__(self, name, dir):
        self.exons = []
        self.name = name
        self.dir = dir
        self.start = -1
        
    def __lt__(self, other):
        if len(self.exons) > 0 and len(other.exons) > 0:
            return self.exons[0][0] < other.exons[0][0]
        else:
            return 0
        
    def __eq__(self, other):
        if len(self.exons) > 0 and len(other.exons) > 0:
            return self.exons[0][0] == other.exons[0][0]
        else:
            return 0
        
    def __gt__(self, other):
        if len(self.exons) > 0 and len(other.exons) > 0:
            return self.exons[0][0] > other.exons[0][0]
        else:
            return 0
        
    def getStart(self):
        return self.start
        
    def getNextStart(self):
        if len(self.exons) == 0:
            return -1
        return self.exons[0][0]
    
    def getNextEnd(self):
        if len(self.exons) == 0:
            return -1
        return self.exons[0][1]
    
    def popExon(self):
        if len(self.exons) == 0:
            return -1
        return self.exons.pop(0)
    
    def putExon(self, exon):
        self.exons.append(exon)
        if self.start == -1:
            self.start = exon[0]
    
    def empty(self):
        return not self.exons
    
    def __str__(self):
        return self.name + " " + str(self.exons)
    
if __name__ == "__main__":   
    __main__()
