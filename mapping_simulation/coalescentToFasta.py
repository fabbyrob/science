doc = '''
Takes in a fasta of sequences, and 3 coalescent simulations to make.

It uses the first coalescent file to simulate divergence between duplicates of each sequence in the fasta.

Next it uses the other two coalescent simulation files to simulate alleles at each duplicate of each gene in the fasta.
One file specifies the diversity at the first duplicate of every locus, the other file specifies the diversity at the other duplicate.
'''

import sys
import getopt
import random

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 5:
        usage()
    
    processArgs(5)
    
    fasta = open(sys.argv[1],"r")  
    divergenceFile = open(sys.argv[2],"r")
    diversity1File = open(sys.argv[3],"r")
    diversity2File = open(sys.argv[4],"r")

    #read in the coalescent info
    #((positions),(haplotype 1), (haplotype 2), ... (haplotype n))
    divergence = readMS(divergenceFile)
    diversity1 = readMS(diversity1File)
    diversity2 = readMS(diversity2File)

    #read in gene sequences
    name = ""
    seqs = []
    for line in fasta:
        line = line.rstrip()
        if line[0] == ">":
            name = line
        else:
            seqs.append((name, line))
            
    #print(seqs)
    
    #make alleles for each gene sequence based on the coalescence
    if len(seqs) != len(divergence) or len(seqs) != len(diversity1) or len(seqs) != len(diversity2):
        use = min(len(seqs), len(divergence), len(diversity1), len(diversity2))
        sys.stderr.write("fasta %s\ndivergence %s\ndiversity 1 %s\ndiversity 2 %s\n" % (len(seqs), len(divergence), len(diversity1), len(diversity2)))
        sys.stderr.write("Mismatched numbers of sequences and coalescence simulations. Using the smaller number (%s).\n" % use)
        seqs = seqs[0:use]
        divergence = divergence[0:use]
        diversity1 = diversity1[0:use]
        diversity2 = diversity2[0:use]
        
    #print(len(divergence))
    #print(len(seqs))
    
    #print(divergence)
    duplicates = makeSequences(seqs, divergence)
    
    #interweave the two diversity files
    totalDiversity = [d for d in flat(zip(diversity1, diversity2))]
    #make the alleles for each duplicate
    makeSequences(duplicates, totalDiversity, pr = True)
    
def flat(ls):
    for i,j in ls:
        yield i
        yield j
    
def makeSequences(seqs, divergence, pr = False):
    newSeqs = []
    for s, c in zip(seqs, divergence):
        #sys.stderr.write("Assigning %s alleles to sequence %s.\n" % (len(c[1]), s[0]))
        name, seq = s
        pos, haplos = c
        
        pos = [int(x*len(seq)) for x in map(float, pos)]
        #print(pos)
        pos = correctPositions(pos, len(seq)-1)
        alts = pickAlt(seq, pos)
        #print(pos)
        #print(alts)
        ctr = 0
        for h in haplos:#generate each haplotype and output it
            currentHaplo = list(seq[:])
            
            for i, p in enumerate(pos):
                alleles = alts[i]
                currentHaplo[p] = alleles[int(h[i])]
                
            currentHaplo = "".join(currentHaplo)
            
            #if we're not printing, then we're making duplicates, give each a unique name
            if not pr:
                myName = name + "_" + str(ctr)
            else:
                myName = name
                print("%s\n%s" % (name, currentHaplo))
                
            newSeqs.append((myName, currentHaplo))
            
            ctr += 1
    return newSeqs
        
def readMS(MSfile):
    #read in the coalescent info
    #((positions),(haplotype 1), (haplotype 2), ... (haplotype n))
    found = False
    pos = []
    haplo = []
    divergence = []
    for line in MSfile:
        line = line.rstrip()
        if "positions:" in line:
            line = line.split()
            pos = line[1:]
            haplo = []
        elif "//" in line:
            #blank line
            if found:
                found = False
                divergence.append((pos, haplo))
        elif "segsites:" in line:
            found = True
            pos = []
        elif found and len(line):
            haplo.append(line)
            
    divergence.append((pos, haplo))
    
    return divergence
        
def correctPositions(ls, max):
    newls = []
    for i, x in enumerate(ls):
        if i == 0:
            val = x
        elif x <= newls[i-1]:
            val = newls[i-1] + 1
        else:
            val = x
            
        if val > max:
            val = max
            
        if newls and newls[-1] == max and val > max:
            #skipping, reached the end already
            continue
        newls.append(val)
    return newls
            
def pickAlt(seq, pos):
    alts = []
    vals = "ACTG"
    #sys.stderr.write("seq: %s\n\tlen: %s\n" % (seq, len(seq)))
    for i in pos:
        #sys.stderr.write("\t%s ," % i)
        alt = random.choice(vals)
        while alt == seq[i]:
            alt = random.choice(vals)
        alts.append((seq[i], alt))
    sys.stderr.write("\n")
    return alts

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        continue
#        if opt == "-o":
#            global _o 
#            _o = arg
#        else:
#            print ("Unrecognized option: "+opt+"\n")
#            usage()

    sys.stderr.write("fasta: %s  divergence: %s  diversity 1:%s  diversity 2:%s\n" % (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))
   
use = "python "+__file__.split("/")[-1]+" fastaFile divergenceFile diversity1File diversity2File"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print(doc)

    
if __name__ == "__main__":   
    __main__()