import sys
import re
import getopt
import cPickle as pickle
import vcf
from math import sqrt

_g = ""
_f = 500
_w = 100
_s = 50
_i = 0
_d = 20
_D = 60
_q = 90
_L = 40

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    annot = pickle.load(open(sys.argv[1], "rb"))
    #print(annot)
    vcf_reader = vcf.Reader(open(sys.argv[2],"rb"))
    names = vcf_reader.samples
    
    ######TAJ D constants  
    N = len(names)*2  
    a1 = 0
    a2 = 0
    for i in range(1,N):
        a1 += 1/float(i)
        a2 += 1/(float(i)*float(i))
        
    b1 = float(N+1)/(3*(N-1))#see Tajima 1989
    b2 = float(2*(N*N+N+3))/(9*N*(N-1))

    c1 = b1 - 1/a1
    c2 = b2-float(N+2)/(a1*N)+a2/(a1*a1)

    e1 = c1/a1
    e2 = c2/(a1*a1+a2)
    
    consts = (len(names), a1, e1, e2)
    ########
    
    flankers = None
    #if _g != "" do a subset, read them into a list
    myGenes = []
    if _g:
        for line in open(_g, "r"):
            line = line.rstrip()
            line = line.split()
            myGenes.append((line[0], int(line[1])))
        flankers = {}
        
        #modify annot to only have those genes
        #flankers = Add the flanking genes to a new dictionary, gene_of_interest:(pre-gene, post-gene)
        #so that I can test if the region overlaps
        new_annot = []
        p_gene = None
        for i in range(len(annot)):
            start, gene = annot[i]
            if gene['name'] in myGenes:
                new_annot.append((start, gene))
                if p_gene:
                    prev = p_gene
                else:
                    prev = None
                    
                if i+1 < len(annot):
                    next = annot[i+1][1]
                else:
                    next = None
                    
                flankers[gene['name']] = (prev, next)
            p_gene = gene
        annot = new_annot
    
    p_gene = None
    gene, name, start, end, dir, pre_range, post_range = getNextRanges(annot, p_gene, flankers)
    
    pre_seq = []
    seq = []
    post_seq = []
    
    #for each gene, step through the file
    #build the total sequence for each individual
    #invert sequence for - strand genes
    #calculate the windows in each of the 3 zones
    #print the windows as we calculate them
    print("Gene,Midpoint,pi,theta,TajD")
    for record in vcf_reader:
        if gene == None:#if we're out of genes we're done
            break
        
        if record.POS >= pre_range[0] and record.POS <= pre_range[1]:#in pre flanking region
            processSite(pre_seq, record)
        elif record.POS >= start and record.POS <= end: #in gene
            processSite(seq, record)
        elif record.POS >= post_range[0]and record.POS <= post_range[1]: #in post flanking region
            processSite(post_seq, record)
        
        if record.POS >= post_range[1]:#done with gene
            if dir == "-":
                pre_seq, post_seq = post_seq, pre_seq
                pre_seq.reverse()
                post_seq.reverse()
                seq.reverse()
            
            #calculate the windows in each sequence
            sys.stderr.write("Process windows in pre for gene: "+name+"\n")
            processWindows(name, pre_seq, -1, consts)
            sys.stderr.write("Process windows in gene: "+name+"\n")
            processWindows(name, seq, 0, consts)
            sys.stderr.write("Process windows in post for gene: "+name+"\n")
            processWindows(name, post_seq, 1, consts)
            
            p_gene = gene
            gene, name, start, end, dir, pre_range, post_range = getNextRanges(annot, p_gene, flankers)
            
            #reset our sequences
            pre_seq = []
            seq = []
            post_seq = []
            

def getNextRanges(annot, p_gene, flankers = None):
    #flankers is the dict of genes flanking our genes of interest, if no -g then no flankers, use the adjacent genes in annotation
    if len(annot) == 0:
        return None, None, None, None, None, None, None
    
    start, gene = annot.pop(0)
    end = gene['end']
    name = gene['name']
    dir = gene['direction']
    
    if flankers:
        pre_gene = flankers[name][0]
        post_gene = flankers[name][1]
    else:
        pre_gene = p_gene
        if len(annot) > 0:
            post_gene = annot[0][1]
        else:
            post_gene = None
    
    #calc the pre_gene flanking region
    if pre_gene:
        pre_min = max(pre_gene['end']+1, start-_f)
    else:
        pre_min = start-_f
    pre_range = (pre_min, start-1)
    
    #calc the post_gene flanking region
    if post_gene:
        post_max = min(end+_f, post_gene['start']-1)
    else:
        post_max = end+_f
    post_range = (end+1, post_max)
        
    return gene, name, start, end, dir, pre_range, post_range

def processSite(seq, record):
    if 'Dels' in record.INFO.keys() and float(record.INFO['Dels']) > _i:
        sys.stderr.write("Site is an indel.\n\t"+str(record.POS)+"\n")
        seq.append(-1)
    elif record.QUAL < _q:
        sys.stderr.write("Site below quality threshold.\n\t"+str(record.POS)+"\n")
        seq.append(-1)
    else:
        seq.append(callSite(record))
 
#return the minor allele frequency at a given site, or -1 if the site fails filters for any individual
def callSite(record):
    freq = -1
    if(record.ALT[0] not in ['A','T','C','G','.']):
        #not a biallelic site
        #ambig
        sys.stderr.write("Non-standard base for alternate: "+str(record.ALT)+" at "+str(record.CHROM)+", "+str(record.POS)+". Calling ambiguous\n")
        return (freq)
    
    freq = 0
    freqFail = False

    for sample in record.samples:
        if freqFail:
            break
        
        if 'PL' not in sample.keys() or 'DP' not in sample.keys():#one individual is missing data, throw out the site
            sys.stderr.write("One sample: "+ str(sample['name']) +"lacks depth or likelihoods at\n\t"+str(record)+"\n")
            freq = -1
            freqFail = True
            break
         
        myDepth = sample['DP'][0]
            
        if myDepth < _d or myDepth > _D:
            #did not pass depth filters
            sys.stderr.write("One sample: "+ str(sample['name']) +" does not pass depth filters\n\t"+str(record)+"\n")
            freqFail = True
            break
        
        myPLs = map(int, sample['PL'])
        (minQ, midQ, maxQ) = getMinMaxIndex(myPLs)
        
        if (myPLs[minQ] != 0):
            #site is ambiguous, no max likelihood
            sys.stderr.write("One sample: "+ str(sample['name']) +" does not pass likelihood filters\n\t"+str(record)+"\n")
            freqFail = True
            break
        
        if (myPLs[midQ] < _L):
            #dont pass quality cut off
            sys.stderr.write("One sample: "+ str(sample['name']) +" does not pass likelihood filters\n\t"+str(record)+"\n")
            freqFail = True
            break
        else:
            freq += minQ#increment freq by the number of ALT alleles
            continue

    if freqFail:#if some filter failed return -1
        freq = -1
        return (freq)

    return (min(freq, 2*len(record.samples)-freq))#return minor allele frequency
 
def processWindows(name, seq, type, consts):
    #type = -1 for pre
    #type = 0 for gene
    #type = 1 for post
    #do shit to each window in this seq, print out the result
    N, a1, e1, e2 = consts
    
    for s in range(0, len(seq), _s):
        myWin = seq[s:min(s+_w, len(seq))]
        #sys.stderr.write(str(myWin)+"\n")
        if type == -1:
            mid = (s+int(len(myWin)/2))-len(seq)
        elif type == 1:
            mid = _f+s+int(len(myWin)/2)+1
        else:
            mid = _f*(s+int(len(myWin)/2))/float(len(seq))
        
        calcWindow(name, myWin, mid, N, a1, e1, e2)
        
    return
        
def calcWindow(name, afs, mid, numInds, a1, e1, e2):
    (pi, sitesUsedPi) = piWindow(afs, numInds*2)
    (tajd, theta, sitesUsedTheta) = tajDWindow(afs, pi, numInds*2, a1, e1, e2)
    if sitesUsedPi == 0 or sitesUsedTheta == 0:
        sys.stderr.write("No viable sites in window: "+str(mid)+"\n")
    elif pi == 0 and theta == 0:
        sys.stderr.write("No diversity in window: "+str(mid)+"\n")
    else:
        print(name+","+str(mid)+","+str(float(pi)/sitesUsedPi)+","+str(float(theta)/sitesUsedTheta)+","+str(tajd))#+","+str(filter(lambda a: a != -1, afs)))

def getMinMaxIndex(ls):
    min = 0
    minV = ls[0]
    max = 0
    maxV = ls[0]
    for i in range(0, len(ls)):
        if ls[i] < minV:
            min = i
            minV = ls[i]
        if ls[i] > maxV:
            max = i
            maxV = ls[i]
            
    return (min, max, list(set([1,2,0])-set([min, max]))[0])

def piWindow(afs, sampSize):
    total = 0
    numSites = 0
    for f in afs:
        if f == -1:#skip
            continue
        minorfreq = (float(f)/sampSize)
        #print("minorFreq ",minorfreq)
        majorfreq = 1-minorfreq
        #print("majorFreq ",majorfreq)
        #print("1-(majorFreq^2+minorFreq^2) ",1-(minorfreq*minorfreq+majorfreq*majorfreq))
        
        total += (1-(minorfreq*minorfreq+majorfreq*majorfreq))
        #print("total ", total)
        numSites += 1


    total = (float(sampSize)/(sampSize-1))*total
    #print("total*n/n-1 ", total)

    return (total, numSites)
                
def tajDWindow(afs, pi, sampSize, a1, e1 ,e2):
    #print(a1)
    #print(afs)
    S = 0#number of polymorphic sites
    nS = 0 #number non-polymorphic
    for s in afs:
        if s == -1:
            continue
        elif s == 0 or s == sampSize:
            nS += 1
        else:
            S += 1
         
    #print("S ",S," nS ",nS)
    if (S+nS) == 0:#no viable sites
        return (".",".",0)
            
    theta = float(S)/a1
    #print("theta ", theta)
    thetaPerSite = theta/(S+nS)

    d = pi-theta
    if(pi == 0 or theta == 0):
        D = '.'
    else:
        p1 = e1*S
        p2 = e2*S*(S-1)
        #sys.stderr.write(mystr)
        D = d/sqrt(p1+p2)
    
    return (D, theta, S+nS)
 
def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"f:w:s:i:d:D:q:L:g:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        continue
        if opt == "-f":
            global _f 
            _f = int(arg)
        elif opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-s":
            global _s 
            _s = int(arg)
        elif opt == "-i":
            global _i 
            _i = float(arg)
        elif opt == "-d":
            global _d 
            _d = int(arg)
        elif opt == "-D":
            global _D 
            _D = int(arg)
        elif opt == "-q":
            global _q 
            _q = int(arg)
        elif opt == "-L":
            global _L 
            _L = int(arg)
        elif opt == "-g":
            global _g 
            _g = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" pickled_annotation VCF [-g gene_list -f flanking_size -w window_size -s step_size -i indel_cutoff -d min_depth -D max_depth -q min_qual -L min_likelihood]"
def usage():
    print (use)
    sys.exit()
    
def details():
    descrip = "This program takes an annotation, and a VCF and analyzes diversity around all genes present in that annotation. \nIt prints out \
Watterson's theta, pi, and Tajima's D for sliding windows around and within genes. \nIf a list of genes is provided it does this analysis only for a subset of genes.\n\
The annotation should be pre-processed using pickleAnnotation.py.\n\
The program outputs data for each window and each gene formatted like:\n\n\
    gene_name,midpoint,pi,theta,TajD\n\n\
    Midpoint is between -1*flanking size and 0 for windows upstream of the gene, \n\
    0 to 100 for windows within the gene (representing a % of the length of the gene). So midpoint of 50 would be halfway across the gene.\n\
    Or between 101 and 100+flanking size for windows downstream of the gene.\n"
    print(descrip)
    print (use)
    print("\nOption details below are printed in the form:\n\t OPTION_NAME ARGUMENT_TYPE (DEFAULT OPTION) - DESCRIPTION\n")
    print("-g STR ("+str(_g)+") - Used to provide a list containing a subset of the genes to analyze.")
    print("-f INT ("+str(_f)+") - The number of bases on either side of a gene to consider. If two genes are closer together than this number, then the full distance is used.")
    print("-w INT ("+str(_w)+") - Size of the sliding window to analyze sequence in.")
    print("-s INT ("+str(_s)+") - Step size for the sliding window.")
    print("-i FLOAT ("+str(_i)+") - Maximum 'Dels' value a site can have and be included in analysis.")
    print("-d INT ("+str(_d)+") - Minimum depth to include a site in the analysis, this is the individual, not total site, depth.")
    print("-D INT ("+str(_D)+") - Maximum depth to include a site in the analysis, this is the individual, not total site, depth.")
    print("-q INT ("+str(_q)+") - Minimum quality to include a site in the analysis.")
    print("-L INT ("+str(_L)+") - Minimum likelihood of a particular individual's genotype call to include a site in the analysis.")

    
if __name__ == "__main__":   
    __main__()