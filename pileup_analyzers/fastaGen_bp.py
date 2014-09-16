import sys
import annotation
import vcf
import pileup
import getopt

_o = ""#output directory
_v = False#verbose mode?
_L = 40#minimum likelihood to call a site
_d = 20#minimum depth
_D = 100#maximum depth
_i = 0#indel threshold for Dels value
_r = False#reference mode?
_q = 40#min quality

def __main__():
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    processArgs(3)

    #initialize my parsers
    annotation_reader = annotation.Reader(open(sys.argv[1], 'rb'))
    if not _r:
        reader = vcf.Reader(open(sys.argv[2], 'rb'))
        names = reader.samples
    else:
        reader = pileup.Reader(open(sys.argv[2], 'rb'))
        names = ["outgroup"]
        
    annotation_iter = annotation_reader.__iter__()
    
    #pull my first gene
    next_gene = getNextGene(annotation_iter)
    
    if next_gene == None:#No genes = done
        sys.stderr.write("No valid genes in annotation. Exiting.\n")
        sys.exit(0)
    
    exons = next_gene.exons
    last_site = -1
    myFasta = Fasta(names)
    
    for record in reader:
        indel = 0
        #if we're out of exons get the next gene
        if len(exons) == 0:
            myFasta.writeFasta(next_gene.name, next_gene.start, next_gene.end, next_gene.direction)
            next_gene = getNextGene(annotation_iter)
            while next_gene != None and next_gene.start < record.POS:
                sys.stderr.write("Gene missed, read in gene at POS = "+str(record.POS)+"\n\t"+str(next_gene)+"\n")
                next_gene = getNextGene(annotation_iter)
                
            if next_gene == None:#no more genes
                break
            exons = next_gene.exons
            
            myFasta = Fasta(names)#reset our fasta data
        
        if not _r and 'Dels' in record.INFO.keys() and record.INFO['Dels'] > _i:#if this site has some probability of having an INDEL
            indel = 1
            #TODO check reference for deletions
        
        if record.POS < next_gene.start:#if we haven't reached the next gene yet
            continue
        
        if record.POS >= exons[0][0]:#if we're past the beginning of the next exon
            if (last_site != -1 and last_site != record.POS-1) or (last_site == -1 and exons[0][0] != record.POS):
                #if we missed some sites in the middle
                #OR
                #if we missed sites at the beginning of a gene
                num = myFasta.fillRemainder(max(last_site, exons[0][0]), record.POS)#fill in any missing sites with N
                #WARNING: if an exon starts right after another ended then this will N the first base in the new exon
                if _v and num > 0:
                    sys.stderr.write("'N'ed "+str(num)+" sites because we missed the beginning of an exon, or something in the middle\n\t"+str(record))
            if record.POS <= exons[0][1]:#if we're within the current exon
                processSite(record, myFasta, indel, next_gene.direction)
                last_site = record.POS
            else:#this exon is done
                num = myFasta.fillRemainder(max(last_site, exons[0][0]), exons[0][1])#fill in any missing sites with N
                if _v and num > 0:
                    sys.stderr.write("'N'ed "+str(num)+" sites because we have passed the end of the exon\n\t"+str(record))
                last_site = -1#we finished this exon
                exons.pop(0)

    if next_gene != None:#finish the last gene with N's
        for exon in exons:
            myFasta.fillRemainder(last_site, exon)
        myFasta.writeFasta(next_gene.name, next_gene.start, next_gene.end, next_gene.direction)


def getNextGene(iter):
    try:
        next = iter.next()
    except StopIteration:
        next = None
    
    return next

def processRefSite(record, fasta, dir):
    if record.REF == "*":#INDEL
        #remove the previous one
        #put an N
        fasta.removePrevious()
        base = ('N','N')
    elif record.QUAL < _q:
        base = ('N','N')
    elif record.ALT not in ['A','T','C','G']:#there should be no heterozygous sites
        base = ('N','N')
    else:
        base = (record.REF, record.ALT)
        
    if dir == '-':#if we are a - direction gene then take the compliment 
        sense = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
        if base[0] not in sense.keys() or base[1] not in sense.keys():
            sys.stderr.write("KEY ERROR: invalid base.\n\tBases: "+str(base)+"\n\trecord:"+str(record)+"\n\tsamples: "+str(sample))
            base = ('N','N')
        else:
            base = (sense[base[0]], sense[base[1]])
        
    fasta.callSite('outgroup', base[0], base[1])

def processSite(record, fasta, indel, dir):
    if _r:
        processRefSite(record, fasta, dir)
        return
    
    hets = 0
    for sample in record.samples:
        #site quality too low
        if record.QUAL < _q:
            base = ('N', 'N')
            fasta.callSite(sample['name'], base[0], base[1])
            continue
        
        if 'PL' in sample.keys():
            (minQ, midQ, maxQ) = getMinMax(sample['PL'])
        else:
            base = ('N', 'N')
            if _v:
                sys.stderr.write("'N' because no 'PL' in record\n\t"+str(record)+"\n")
            fasta.callSite(sample['name'], base[0], base[1])
            continue
        
        if indel:
            base = ('N', 'N')
            if _v:
                sys.stderr.write("'N' because in INDEL site\n\t"+str(record)+"\n")
        elif not ('DP' in sample.keys() or 'PL' in sample.keys):
            sys.stderr.write("One sample: "+ sample +"lacks depth or likelihoods at\n\t"+str(record)+"\n")
            base = ('N','N')
        else:
            if sample['PL'][minQ] != 0:#no minimum quality N
                base = ('N','N')
                if _v:
                    sys.stderr.write("'N' because min PL != 0\n\t"+str(record)+"\n")
            elif sample['PL'][midQ] < _L:#min quality not high enough
                base = ('N','N')
                if _v:
                    sys.stderr.write("'N' because mid PL < "+str(_L)+"\n\t"+str(record)+"\n")
            elif sample['DP'][0] < _d or sample['DP'][0] > _D:#not in depth range
                base = ('N','N')
                if _v:
                    sys.stderr.write("'N' because depth outside range ("+str(_d)+", "+str(_D)+") \n\t"+str(record)+"\n")
            else:
                if minQ == 0:#homo REF
                    base = (record.REF, record.REF)
                elif minQ == 1:#heterozygote
                    base = (record.REF, record.ALT[0])
                    hets += 1
                else:#homo ALT
                    base = (record.ALT[0], record.ALT[0])
              
        if dir == '-':#if we are a - direction gene then take the compliment 
            sense = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
            if base[0] not in sense.keys() or base[1] not in sense.keys():
                sys.stderr.write("KEY ERROR: invalid base.\n\tBases: "+str(base)+"\n\trecord:"+str(record)+"\n\tsamples: "+str(sample))
                base = ('N','N')
            else:
                base = (sense[base[0]], sense[base[1]])
                    
        fasta.callSite(sample['name'], base[0], base[1])
    
    if hets == len(record.samples):#if all samples are heterozygous
        fasta.removePrevious()#take off the calls we made
        for sample in record.samples:#make all of them N
            fasta.callSite(sample['name'], 'N', 'N')
                    
 
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
 
class Fasta():
    def __init__(self, samples):
        self.names = samples
        self.Seqs = {}
        
        for ind in samples:
            self.Seqs[ind] = ([],[])
            
    def callSite(self, ind, base, base2):
        if ind not in self.names:
            sys.stderr.write("Invalid sample name: "+str(ind))
            sys.exit(0)
            
        self.Seqs[ind][0].append(base)
        self.Seqs[ind][1].append(base2)
        
    def fillRemainder(self, last_set, current):
        for ind in self.names:
            i = 0
            if not isinstance(current, int):
                if last_set < current[0]:
                    last_set = current[0]
                current = current[1]+1
                    
            while i < max(0, current - last_set):
                self.callSite(ind, 'N', 'N') 
                i += 1

        return (current - last_set)
      
    def removePrevious(self, n = 1):
        for ind in self.names:
            self.Seqs[ind][0].pop(-1)
            self.Seqs[ind][1].pop(-1)
                
    def writeFasta(self, gene, start, end, direction, mode = "a"):
        if _r:#if in reference mode then overwrite files, otherwise append
            mode = "w"
            
        file_name = _o+"".join(gene.split(':'))+"_"+str(start)+"_"+str(end)+".fasta"
        #print(file_name)
        myFile = open(file_name, mode)
        
        for ind in self.names:
            if direction == "-":
                self.Seqs[ind][0].reverse()
                self.Seqs[ind][1].reverse()
            
            if not _r:
                myFile.write(">"+ind+"_1\n")
                myFile.write(self.prettyList(self.Seqs[ind][0])+"\n")
            
            myFile.write(">"+ind+"_2\n")
            myFile.write(self.prettyList(self.Seqs[ind][1])+"\n")
            
        myFile.close()
            
    def prettyList(self, list):
        str = ""
        for n in list:
            str+=n
        return str

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"o:L:d:D:i:q:vr")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-o":
            global _o 
            _o = arg
        elif opt == "-L":
            global _L 
            _L = int(arg)
        elif opt == "-d":
            global _d 
            _d = int(arg)
        elif opt == "-D":
            global _D 
            _D = int(arg)
        elif opt == "-i":
            global _i 
            _i = float(arg)
        elif opt == "-q":
            global _q 
            _q = int(arg)
        elif opt == "-v":
            global _v 
            _v = True
        elif opt == "-r":
            global _r 
            _r = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
            
use = "python "+__file__.split("/")[-1]+" <ANNOT> <VCF> [options]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print("Takes in an annotation, and a VCF file and generates one fasta for each gene. \n~~~Assumes all sites/genes are on the same scaffold.~~~ \n~~~Assumes that the VCF is sorted.~~~")
    print("")
    print("python fastaGen.py scaffold1_annotation.txt scaffold1.vcf -o scaffold1_fastas/ 2> fastaGen.log &")
    print("")
    print (use)
    print("options (defaults in brackets):")
    print("    ANNOT - an annotation file \n\
    VCF - a VCF file \n\
    o - STR -  output directory for fasta files, must end in \'/\' [\'\'] \n\
    L - INT -  minimum likelihood to call a genotype ["+str(_L)+"]\n\
    d - INT -  minimum depth to call a site ["+str(_d)+"]\n\
    D - INT -  maximum depth to call a site ["+str(_D)+"]\n\
    i - FLOAT -  maximum value of 'Dels' to consider a site ["+str(_i)+"]\n\
    r - flag that indicates we are outputting outgroup values into the fastas ["+str(_r)+"]\n\
    v - verbose mode, if on will tell you the reason for every 'N' call ["+str(_v)+"]\n\
    ")

    
if __name__ == "__main__":   
    __main__()