doc = """
Takes a summary file with individual genotypes and a gene annotation and outputs fastas for each gene.
"""

import sys
import getopt
import summary
import annotation

_o = "./"

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    summaryReader = summary.Reader(open(sys.argv[1],"rb"))
    
    annotationReader = annotation.Reader(open(sys.argv[2],"rb"))  
    annotIter = annotationReader.__iter__()
        
    myGene = getNextGene(annotIter)
    
    if myGene == None:
        sys.stderr.write("No genes listed in annotation. Exiting.\n")
        sys.exit(0)
        
    if len(summaryReader.summary.Samples) == 0:
        sys.stderr.write("No individual genotype info in this Summary file. Cannot make fastas.\n")
        sys.exit(0)
        
    seqs = {}
    for samp in summaryReader.summary.Samples:
        seqs[samp] = [[],[]]
        
    doneScafs = []
        
    genoDict = {v:k for k, v in summaryReader.summary.Genotypes.items()}
    pScaf = ""
    pSite = 0
    #read the infile file...
    for site in summaryReader:
        #sys.stderr.write("Processing site %s\n" % site.prettyStr())
        if site.CHROM != myGene.scaf:
            while myGene.scaf in doneScafs:
                sys.stderr.write("Missed gene %s on scaf %s. Make sure file is sorted correctly.\n" % (myGene.name, myGene.scaf))
                myGene = getNextGene(annotIter)
                
            for samp in summaryReader.summary.Samples:
                seqs[samp] = [[],[]]
                
            if site.CHROM != myGene.scaf:#haven't reached this scaf yet, skip along in the summary
                continue
            
        if pScaf == "" or pScaf != site.CHROM:
            doneScafs.append(pScaf)
            pScaf = site.CHROM
            pSite = 0
            
        myGene = addNs(seqs, pSite, site.POS, myGene, annotIter)
        
        if myGene == None:
            #sys.stderr.write("main 71: ran out of genes. Exiting.\n")
            break
            
        if site.POS < myGene.exons[0][0]:
            pass
        
        elif site.POS <= myGene.exons[0][1]:
            for samp, geno in site.Genotypes.items():
                if geno == genoDict['heterozygote']:
                    seqs[samp][0].append(site.REF)
                    seqs[samp][1].append(site.ALT)
                elif geno == genoDict['homozygote reference']:
                    seqs[samp][0].append(site.REF)
                    seqs[samp][1].append(site.REF)
                elif geno == genoDict['homozygote alternate']:
                    seqs[samp][0].append(site.ALT)
                    seqs[samp][1].append(site.ALT)
                else:
                    seqs[samp][0].append("N")
                    seqs[samp][1].append("N")
                    
        if site.POS >= myGene.exons[0][1]:
            myGene.exons.pop(0)
            
        if not myGene.exons:
            #print the seqs to a file
            outputFasta(myGene, seqs)
            #grab new gene
            myGene = getNextGene(annotIter)
            #check it
            if myGene == None:
                break
            #reset seqs
            for samp in summaryReader.summary.Samples:
                seqs[samp] = [[],[]]
                
        pSite = site.POS

def addNs(seq, prevSite, currSite, gene, iter):
    #sys.stderr.write("addNs(%s, %s, %s, %s, %s)\n" % (seq[seq.keys()[0]], prevSite, currSite, gene, iter))
    while currSite != prevSite+1:
        #sys.stderr.write("addNs: Missed site %s. Processing...\n" % (prevSite+1))
        if prevSite + 1 < gene.exons[0][0]:
            pass
        elif prevSite + 1 <= gene.exons[0][1]:
            for samp in seq.keys():
                seq[samp][0].append("N")
                seq[samp][1].append("N")
        
        if prevSite + 1 >= gene.exons[0][1]:
            gene.exons.pop(0)
            
        if not gene.exons:
            #print("addNs: ran out of exons. Grabbing a new gene.\n")
            #print the seqs to a file
            outputFasta(gene, seq)
            #grab new gene
            gene = getNextGene(iter)
            #check it
            if gene == None:
                return None
            #reset seqs
            for samp in seq.keys():
                seq[samp] = [[],[]]
                
        prevSite += 1
                
    return gene
                    
def outputFasta(gene, seqs):
    myFile = open(_o+gene.name+".fasta", "w")
    
    if gene.direction == "-":
        reverseCompliment(seqs)
    
    for samp, seq in seqs.items():
        myFile.write(">%s_1 %s\n" % (samp, gene.name))
        myFile.write("".join(seq[0])+"\n")
        
        myFile.write(">%s_2 %s\n" % (samp, gene.name))
        myFile.write("".join(seq[1])+"\n")
    
    myFile.close()                    

def reverseCompliment(seqs):
    for samp, seq in seqs.items():
        for strand in seq:
            strand.reverse()
            for i, base in enumerate(strand):
                if base == "A":
                    strand[i] = "T"
                elif base == "T":
                    strand[i] = "A"
                elif base == "G":
                    strand[i] = "C"
                elif base == "C":
                    strand[i] = "G"

def getNextGene(iter):
    try:
        next = iter.next()
    except StopIteration:
        next = None
    
    return next

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"o:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-o":
            global _o 
            if arg[-1] != "/":
                arg += "/"
            _o = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("summary: %s annotation: %s -o %s\n" % (sys.argv[1], sys.argv[2], _o))
   
use = "python "+__file__.split("/")[-1]+" summaryFile.txt geneAnnotation.txt [options]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("o - STR - %s - output directory. WARNING: this program will make a lot of files. Be careful where you put this." % _o)

    
if __name__ == "__main__":   
    __main__()
