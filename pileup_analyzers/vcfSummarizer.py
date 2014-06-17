'''
This program takes in a VCF and an annotation then outputs a summary of the VCF
Each line includes information about one site, what kind of site it is based on the annotation and the allele frequencies at that site
and column indicating whether the site is diverged or not (0 no, 1 yes, -1 no data):

#site types: {'0fold': 3, 'stop': 8, 'intergene': 0, 'downstream': 6, 'upstream': 5, 'exon': 2, 'intron': 1, 'istop': 7, '4fold': 4, 'unknown': 9}
#CHROM    POS    REF    ALT    REF_NUMBER    ALT_NUMBER    TOTAL    SITE_TYPE DIVERGENCE
scaffold_1    3    A    T    7    6    13    2    0

To calculate allele counts at each site individuals who do not pass filters are excluded
or if the entire site does not pass the quality or indel filters all individuals are excluded

usage:
python vcfSummarizer.py mySamples.vcf myAnnotation.txt pickled_divergence [OPTIONS] > myOutput.txt 

run the program with no arguments for a full list of options and descriptions of each.

Requires the VCF parser, vcf.py.
'''

import sys
import getopt
import vcf
import annotation
import cPickle as pickle

_q = 40 #min quality
_d = 20 #min depth
_D = 60 #max depth
_L = 40 #GQ cutoff
_i = 0.0 #Dels cutoff
_v = False #verbose mode
_G = False #genotype mode
_a = False #gene mode
_t = "" #file for divergence

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    reader = vcf.Reader(open(sys.argv[1], 'rb'))
    names = reader.samples
        
    sites_file = open(sys.argv[2],"r")  
        
    if (sites_file == None):
        print("Bad sites file name: "+sys.argv[2])
        sys.exit()
        
    if _t:
        divergence = pickle.load(open(_t, "rb"))

    types = {'intergene':0 , 'intron':1, 'exon':2, '0fold':3, '4fold':4, '3utr':5, '5utr':6, 'istop':7, 'stop':8, 'unknown':9, 'cnc':10}
        
    #print("#site types: "+str(types))
    for k,v in types.items():
        print("#TYPE %s %s" % (k, v))
        
    if _G:
        print("#GENOTYPE;homozygote reference;R\n#GENOTYPE;homozygote alternate;A\n#GENOTYPE;heterozygote;H\n#GENOTYPE;unknown;N")
    header = "#CHROM\tPOS\tREF\tALT\tREF_NUMBER\tALT_NUMBER\tTOTAL\tSITE_TYPE\tDIVERGENCE"
    if _a:
        header += "\tGENE\tDIR"
    if _G:
        for samp in reader.samples:
            header += "\t%s" % samp
    print(header)
    
    #read the infile file...
    scaf, pos, type, gene, dir = getNextSite(sites_file)#first annotated site
    for record in reader:
        if scaf == None:
            if _v:
                sys.stderr.write("Ran out of sites.\n")
            break
        
        ref, alt, total, genos = processSite(record)
        
        if _t and record.POS < len(divergence.get(record.CHROM)):
            div = divergence[record.CHROM][record.POS]
        else:
            div = (None, None, None)
       
        #sys.stderr.write("%s\n" % div)
        if div[2] == None:
            div = -1
        elif div[2]:
            div = 0
        else:
            div = 1
        
        if record.CHROM < scaf:#we haven't reached the next annotated line yet, carry on
            if _v:
                sys.stderr.write("Site not annotated, looking for next chromosome."+record.CHROM +" "+str(record.POS)+"\n")
            printSummary(record, ref, alt, total, -1, div, genos, "-", "0")
            continue
        
        if record.CHROM > scaf:# we missed some sites or the sites file is not sorted correctly
            if _v:
                sys.stderr.write("Missed sites, or sites file not sorted. Looking for site on chrom: "+scaf+" at chrom: "+record.CHROM+" in VCF. Skipping through sites... ")
            ctr = 0
            while scaf != None and scaf < record.CHROM:
                scaf, pos, type, gene, dir = getNextSite(sites_file)
                ctr += 1
            if _v:
                sys.stderr.write("skipped "+str(ctr)+" sites.\n")
        
        if scaf == None:
            if _v:
                sys.stderr.write("Ran out of sites.\n")
            break
        
        #CHROM is fine, now just need to deal with pos
        
        if record.POS < pos:#haven't reached our next site yet
            if _v:
                sys.stderr.write("No type data for site: "+record.CHROM+" "+str(record.POS)+"\n")
            printSummary(record, ref, alt, total, -1, div, genos, "-", "0")
            continue
        
        if record.POS > pos:#we missed our last site, grab a new one
            if _v:
                sys.stderr.write("Missed sites, or sites file not sorted. Looking for site at: "+scaf+ " "+ str(pos)+" at: "+record.CHROM+ " " + str(record.POS)+" in VCF. Skipping through sites... ")
            ctr = 0
            while pos < record.POS:
                scaf, pos, type, gene, dir = getNextSite(sites_file)
                ctr += 1
            if _v:
                sys.stderr.write("skipped "+str(ctr)+" sites.\n")
                
        if scaf == None:
            if _v:
                sys.stderr.write("Ran out of sites.\n")
            break
        
        if record.POS == pos:#at a good site, print the data
            printSummary(record, ref, alt, total, type, div, genos, gene, dir)
            scaf, pos, type, gene, dir = getNextSite(sites_file)
            
    if _v:
        sys.stderr.write("Ran out of VCF.\n")
        
def printSummary(record, ref, alt, total, type, div, genos, gene, dir):
    myStr = record.CHROM+"\t"+str(record.POS)+"\t"+str(record.REF)+"\t"+str(record.ALT[0])+"\t"+str(ref)+"\t"+str(alt)+"\t"+str(total)+"\t"+str(type)+"\t"+str(div)
    if _a:
        myStr += "\t%s\t%s" % (gene, dir)
    if _G:
        for g in genos:
            myStr += "\t%s" % g 
    print(myStr)

def getNextSite(sitefile):
    line = sitefile.readline()
    if line == "":
        return (None, None, None, None, None)
    
    line = line.rstrip()
    
    while line[0] == "{":#dictionary line
        return getNextSite(sitefile)
    
    sline = line.split()
    if len(sline) < 6:
        if _v:
            sys.stderr.write("Line not formatted correctly\n\t"+str(line)+"\n")
        return getNextSite(sitefile)
    
    scaf = sline[0]
    pos = int(sline[1])
    gene = sline[3]
    dir = sline[4]
    type = sline[5]
    return (scaf, pos, type, gene, dir)

def getMinMax(ls):
    min = 0
    mid = 0
    max = 0
    
    for i in range(1, len(ls)):
        if ls[i] <= ls[min]:
            mid = min
            min = i
        elif ls[i] >= ls[max]:
            mid = max
            max = i
        else:
            mid = i
    return (min, mid, max)

def processSite(record):
    ref = 0
    alt = 0
    
    #site quality too low
    if record.QUAL < _q:
        if _v:
            sys.stderr.write("Site lacks sufficient quality at\n\t"+str(record)+"\n")
        return (ref, alt, ref+alt, ["N"]*len(record.samples))
    
    #site indel rate too high
    if 'Dels' in record.INFO.keys() and float(record.INFO['Dels']) > _i:
        if _v:
            sys.stderr.write("Site lacks Dels info, or Dels is too high at\n\t"+str(record)+"\n")
        return (ref, alt, ref+alt, ["N"]*len(record.samples))
    
    genos = []#R = reference homo, A = alt homo, H = het, U = unknown
    for sample in record.samples:
        
        if not ('DP' in sample.keys()):
            if _v:
                sys.stderr.write("One sample: "+ str(sample['name']) +"lacks depths at\n\t"+str(record)+"\n")
            genos.append("N")
            continue
        
        
        depth = sample['DP'][0]
        if depth < _d or depth > _D:#not in depth range
                if _v:
                    sys.stderr.write("'N' because depth outside range ("+str(_d)+", "+str(_D)+") \n\t"+str(sample['name'])+" "+str(record)+"\n")
                genos.append("N")
                continue
        
        #deprecated ML cutoff code
#        if 'PL' in sample.keys() and record.ALT[0] != ".":
#            myPLs = map(int, sample['PL'])
#            (minQ, midQ, maxQ) = getMinMax(myPLs)
#            
#            if myPLs[minQ] != 0:#no minimum quality N
#                if _v:
#                    sys.stderr.write("'N' because min PL != 0\n\t"+str(sample['name'])+" "+str(record)+"\n")
#                genos.append("N")
#                continue
#            elif myPLs[midQ] < _L:#min quality not high enough
#                if _v:
#                    sys.stderr.write("'N' because mid PL < "+str(_L)+"\n\t"+str(sample['name'])+" "+str(record)+"\n")
#                genos.append("N")
#                continue
#            else:
#                if minQ == 0:#homo REF
#                    ref += 2
#                    genos.append("R")
#                elif minQ == 1:#heterozygote
#                    ref += 1
#                    alt += 1
#                    genos.append("H")
#                else:#homo ALT
#                    alt += 2
#                    genos.append("A")
            
        if 'GQ' in sample.keys() and sample['GQ'] >= _L:
            if sample['GT'] == "0/0":
                ref += 2
                genos.append("R")
            elif sample['GT'] == "0/1":
                ref += 1
                alt += 1
                genos.append("H")
            elif sample['GT'] == '1/1':
                alt += 2
                genos.append("A")
            elif sample['GT'] == './.':
                genos.append("N")
        elif record.ALT[0] != ".":
            if _v:
                sys.stderr.write("One sample: "+ str(sample['name']) +"lacks likelihoods at\n\t"+str(record)+"\n")
            genos.append("N")
            continue
        else:#its a homozygote, no alternate qualities
            ref += 2
            genos.append("R")
        
    return(ref, alt, ref+alt, genos)

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"q:d:D:L:i:t:aGv")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-q":
            global _q
            _q = int(arg)
        elif opt == "-d":
            global _d
            _d = int(arg)
        elif opt == "-D":
            global _D
            _D = int(arg)
        elif opt == "-L":
            global _L 
            _L = int(arg)
        elif opt == "-i":
            global _i 
            _i = int(arg)
        elif opt == "-G":
            global _G 
            _G = True
        elif opt == "-v":
            global _v 
            _v = True
        elif opt == "-a":
            global _a 
            _a = True
        elif opt == "-t":
            global _t 
            _t = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" VCF ANNOTATION [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("This program takes in a VCF and an annotation and a pickled divergence dictionary then outputs a summary of the VCF\n\
Each line includes information about one site, what kind of site it is based on the annotation and the allele frequencies at that site:\n\
\n\
#site types: {'0fold': 3, 'stop': 8, 'intergene': 0, 'downstream': 6, 'upstream': 5, 'exon': 2, 'intron': 1, 'istop': 7, '4fold': 4, 'unknown': 9}\n\
#CHROM    POS    REF    ALT    REF_NUMBER    ALT_NUMBER    TOTAL    SITE_TYPE\n\
scaffold_1    3    A    T    7    6    13    2\n\
\n\
To calculate allele counts at each site individuals who do not pass filters are excluded\n\
or if the entire site does not pass the quality or indel filters all individuals are excluded")
    print()
    print("option - argument type - default - description")
    print("q - INT - "+str(_q)+" - the minimum quality for a site to be considered")
    print("d - INT - "+str(_d)+" - the minimum coverage for an individual to be considered")
    print("D - INT - "+str(_D)+" - the maximum coverage for an individual to be considered")
    print("L - INT - "+str(_L)+" - The minimum GQ for an individual's genotype to be considered valid.")
    print("i - float - "+str(_i)+" - the maximum Dels value for a site to be considered")
    print("G - None - "+str(_G)+" - this flag turns on output of individual genotypes for every sample")
    print("v - None - "+str(_v)+" - this flag indicated that verbose mode should be used, the reason for excluding each site and individual will be printed to stderr")
    print("a - None - "+str(_a)+" - if this argument is used then the gene names from the given annotation are added to the summary")
    print("t - Str - "+str(_t)+" - this option takes a pickled file name where divergence info is located, if used divergence info is added to the summary.")
if __name__ == "__main__":   
    __main__()
