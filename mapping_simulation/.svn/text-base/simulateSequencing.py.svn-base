import sys
sys.path.append("/opt/python27/lib")
sys.path.append("/usr/lib64/python2.6/site-packages/") 
from numpy.random import exponential, poisson, beta
import random
import getopt
import time

_len = 100 #read length
_exp = 1 #exponential shape parameter
_err = 10**-7 #per base error rate in sequencing
_N = 1000#number of reads to sample
_paired = False#sample paired end reads?
_pair_size = 300#distance between paired end reads

_plow_tail = 0.25#probability the tail will decrease
_low_tail_len = 0.2#how much of the tail gets a hit in quality score
_low_tail_cost = 0.2#how much lower the tail is relative to the rest of the read

_plow = 0.06#probability of getting the lowest quality score
_alpha = 2#shape parameter of beta dist for sampling qualities
_beta = 30
_max = 75 #maximum value beta distribution can take, affects how truncated the beta is

_exp_file = "expression.txt"
_v = False#verbose mode

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    infile = open(sys.argv[1],"r")
    
    sortedLoci = []
    lociAlleles = {}
    lociExpression = {}
    expression = open(_exp_file, "w")
    
    totalExpression = 0.0
    #read in the sequences for each gene
    verbose("Reading in sequences")
    name = ""
    seq = ""
    for line in infile:
        line = line.rstrip()
        if line[0] == ">":
            if name != "":
                if name in lociAlleles.keys():
                    lociAlleles[name].append(seq)
                else:
                    lociAlleles[name] = [seq]
                    #assign an expression level to each gene, sample from exponential
                    myExpression = exponential(_exp)
                    lociExpression[name] = myExpression
                    
                    totalExpression += myExpression
            #reset variables
            name = line[1:]
            seq = ""
        else:
            seq += line
            
    #grab the last sequence
    if name in lociAlleles.keys():
        lociAlleles[name].append(seq)
    else:
        lociAlleles[name] = [seq]
        #assign an expression level to each gene, sample from exponential
        myExpression = exponential(_exp)
        lociExpression[name] = myExpression
        
        totalExpression += myExpression
        
    #normalize the expression of each gene by the "total expression"
    verbose("Normalizing expression levels")
    for locus in lociExpression.keys():
        lociExpression[locus] = lociExpression[locus]/totalExpression
        expression.write("%s %s\n"% (locus, lociExpression[locus]))
        #remove short genes
        safe = True
        if _paired:
            if len(lociAlleles[locus][0]) < _pair_size:
                sys.stderr.write("Alleles for locus %s are shorter than fragment size (locus: %s, fragments: %s), removing it from pool.\n" % (locus, len(lociAlleles[locus][0]), _pair_size))
                safe = False
        if safe:
            sortedLoci.append((lociExpression[locus], locus))
    expression.close()
    sortedLoci.sort()

    #simulate extraction
    #sample randomly from all genes/alleles, with probability equal to their expression
        #sample a set number? _N
            #for each one, simulate some sequencing error # = poisson(len(seq)*_err)
            #randomly pick loci to mutate, randomly assign one of the other 3 bases
            #randomly pick a stretch of the sequence that is _len bases long
        
            #read quality? what do we?
        
            #print the sequence
    i = 0
    while i < _N:
        verbose("Sampling Read", i)
        myLocus, myAllele = pickLocus(sortedLoci, lociAlleles)
        
        myReads = pickRead(myAllele)
        if not _paired:
            myRead = myReads[0]
            myErroredRead, myErrorLocs = errorSequence(myRead)
            myQualities = assignQualities(myErroredRead)
            printRead(myErroredRead, myQualities, i)
        else:
            myFirstRead = myReads[0]
            myErroredRead1, myErrorLocs1 = errorSequence(myFirstRead)
            myQualities1 = assignQualities(myErroredRead1)
            printRead(myErroredRead1, myQualities1, i, 1)
            
            mySecondRead = myReads[1]
            myErroredRead2, myErrorLocs2 = errorSequence(mySecondRead)
            myQualities2 = assignQualities(myErroredRead2)
            printRead(myErroredRead2, myQualities2, i, 2)
        i += 1
        

def pickLocus(loci, alleles):
    verbose("pickLocus")
    p = random.random()
    tot = 0
    i = 0
    locus = ""
    while tot < p and i < len(loci)-1:
        tot += loci[i][0]
        locus = loci[i][1]
        i += 1
        
    #locus = loci[pick[i%4]][1]
    allele = random.choice(alleles[locus])
    return (locus, allele)

bases = ["A","C","G","T"]
def errorSequence(seq):
    verbose("errorSequence")
    numErrors = poisson(_err*len(seq))
    errors = []
    newSeq = list(seq)
    for i in range(numErrors):
        #pick place to error
        loc = random.randint(0,len(seq)-1)
        while loc in errors:#dont error the same site twice
            loc = random.randint(0,len(seq)-1)
            
        #pick error
        mut = random.choice(bases)
        while mut == newSeq[loc]:#don't pick the same base as was already there
            mut = random.choice(bases)
            
        newSeq[loc] = mut
            
    newSeq = "".join(newSeq)
    return (newSeq, errors)

def pickRead(seq):
    verbose("pickRead")
    if not _paired:
        start = random.randint(0, (len(seq)-_len)+1)
        return [seq[start:start+_len]]
    else:
        start = random.randint(0, (len(seq)-_pair_size)+1)
        fragment = seq[start:start+_pair_size]
        return [fragment[0:_len], fragment[-1*_len-1:-1]]

def assignQualities(seq):#does it need to know where errors were?
    verbose("assignQualities")
    quals = ""
    myScore = -1
    up = 0
    down = 0
    
    #sample starting quality
    if random.random() < _plow:
        qual = 35
    else:
        qual = _max-beta(_alpha, _beta)*_max
        while qual < 36 or qual >= 74:
            qual = _max-beta(_alpha, _beta)*_max
    myScore = int(qual)
    quals = [myScore]*len(seq)#set entire read to this quality
    
    #roll to see if tail decreases
    if random.random() < _plow_tail:
        #decrease tail
        low_start = int((1-_low_tail_len)*len(seq))
        lowScore = int(min(73, max(35, myScore*(1-_low_tail_cost))))
        quals[low_start:] = [lowScore]*(len(quals[low_start:]))
    
    quals = map(chr, quals)
    quals = "".join(quals)

    return quals

def printRead(read, quals, name, pair=False):
    verbose("printRead")
    if not pair:
        myStr = "@%s:1:1:1:1:#0\n%s\n+%s:1:1:1:1:#0\n%s" % (name, read, name, quals) #TODO check this
        print(myStr)
    else:
        myStr = "@%s:1:1:1:1:#0/%s\n%s\n+%s:1:1:1:1:#0/%s\n%s" % (name, pair, read, name, pair, quals) #TODO check this
        print(myStr)

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"l:e:r:N:s:L:A:B:m:t:z:c:x:pv")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _len 
            _len = int(arg)
        elif opt == "-e":
            global _exp 
            _exp = float(arg)
        elif opt == "-r":
            global _err 
            _err = float(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-s":
            global _pair_size 
            _pair_size = int(arg)
        elif opt == "-L":
            global _plow 
            _plow = float(arg)
        elif opt == "-A":
            global _alpha 
            _alpha = float(arg)
        elif opt == "-B":
            global _beta 
            _beta = float(arg)
        elif opt == "-m":
            global _max 
            _max = float(arg)
        elif opt == "-c":
            global _low_tail_cost 
            _low_tail_cost = float(arg)
        elif opt == "-t":
            global _plow_tail 
            _plow_tail = float(arg)
        elif opt == "-z":
            global _low_tail_len 
            _low_tail_len = float(arg)
        elif opt == "-x":
            global _exp_file 
            _exp_file = arg
        elif opt == "-p":
            global _paired 
            _paired = True
        elif opt == "-v":
            global _v 
            _v = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
            
    sys.stderr.write("infile %s expression file %s -l %s -e %s -r %s -N %s -p %s -s %s -L %s -A %s -B %s -m %s -t %s -z %s -c %s\n" % (sys.argv[1], _exp_file, _len, _exp, _err, _N, _paired, _pair_size, _plow, _alpha, _beta, _max, _plow_tail, _low_tail_len, _low_tail_cost))

def verbose(name, other = ""):
    if _v:
        sys.stderr.write("%s\t%s\n" % (name, other))

use = "python "+__file__.split("/")[-1]+" sequence_file.fasta [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("")
    print("Takes in a fasta of sequences from many loci. If sequences are from the same locus, then the gene name should be identical in the fasta label line.\
     Samples from those loci to simulate library prep, and sequencing.\
    Loci are assigned an expression level randomly from an exponential distribution. Each read has a random number of errors sampled from a poisson distribution.\n")
    print("Note on exponential: numpy uses the exponential parameterization such that you input beta where lambda = 1/beta.\n")
    print("option - type - default - description")
    print("N - INT - %s - the number of reads to sample for sequencing" % _N)
    print("l - INT - %s - the length of reads to be sampled" % _len)
    print("e - FLOAT - %s - shape parameter for the exponential distribution for gene expression" % _exp)
    print("r - FLOAT - %s - per base sequencing error rate" % _err)
    print("p - NONE - %s - use this option to output paired reads to stdout" % _paired)
    print("s - INT - %s - specifies size of fragments to sample paired end reads from (only used with -p)" % _pair_size)
    print("L - FLOAT - %s - probability of getting the lowest quality score" % _plow)
    print("A - FLOAT - %s - alpha parameter for beta distribution from which quality scores are sampled" % _alpha)
    print("B - FLOAT - %s - beta parameter for beta distribution from which quality scores are sampled" % _beta)
    print("m - FLOAT - %s - determines the maximum value of the beta distribution. Affects how truncated the distribution is. A higher value truncates more of the distribution." % _beta)
    print("t - FLOAT - %s - the probability that the tail of a read decreases in quality" % _plow_tail)
    print("z - FLOAT - %s - the fraction of the read that is decreased in quality if has a low quality tail" % _low_tail_len)
    print("c - FLOAT - %s - the decrease in quality score the tail takes relative to the rest of the read (i.e. if -c 0.5 then the tail has a 50 percent lower quality score than the rest of the read)" % _low_tail_cost)
    print("x - STR - %s - filename that expression levels are output to" % _exp_file)
    print("v - NONE - %s - flag that turns on verbose mode for debugging" % _v)


if __name__ == "__main__":
    t0 = time.clock()
    __main__()
    sys.stderr.write("elapsed time: %.2f seconds\n" % (time.clock()-t0))

