import sys, argparse
from collections import namedtuple

codes = {'2fold':2, '3fold':2,'intergenic':0 , 'intron':1, 'exon':2, '0fold':3, '4fold':4, '3utr':5, '5utr':6, 'istop':7, 'stop':8, 'unknown':9}

def main():
    print(codes)#output so that the codes are stored in the file
  
    #read in the annotation
    #throw out anything not CDS or UTR
    parser = GFFparser(open(args.annotation, 'r'))
    annotation = {"CDS":{}, "three_prime_UTR":{}, "five_prime_UTR":{}}#type: {scaffold: [GFF items]}
    for item in parser.generator():
        if item.type not in ["CDS", "three_prime_UTR", "five_prime_UTR"]:
            continue
        
        if item.scaf not in annotation[item.type].keys():
            annotation[item.type][item.scaf] = []
        
        annotation[item.type][item.scaf].append(item)
    
    #sort them on each scaffold so they appear in the right order
    #since GFFs aren't sorted for some reason
    for k in annotation.keys():
        for s in annotation[k].keys():
            annotation[k][s] = sorted(annotation[k][s])
      
    sys.stderr.write("Checking overlaps. If you have any you may want to completely remove these regions downstream.\n")
    #check for overlaps, print errors for problem regions, and eliminate one of them
    for k in annotation.keys():#for each type
        for s in annotation[k].keys():#for each scaf
            removeOverlaps(annotation[k][s])
                
    #read in the reference
    #pop off the next item in the GFF as we go
    #if we've missed all or part of it output an error
    scaf = None
    seq = ""
    for line in open(args.reference, 'r'):
        if line.startswith(">"):
            #deal with the old scaffold first
            if scaf and seq:
                processSeq(scaf, seq, annotation)
            #deal with the new scaffold
            #assumes the header is of the form:
                #>CHROM_NAME other stuff by white space
            scaf = line[1:].split()[0]
            seq = ""
        else:
            seq += line.rstrip()
    #get the last one
    if scaf and seq: processSeq(scaf, seq, annotation)

'''removes overlapping sequences...
'''
def removeOverlaps(ls):
    i = 0
    while i < len(ls)-1:
        item1 = ls[i]
        item2 = ls[i+1]
        if item2.regions[0][0] < item1.regions[-1][1]:
            sys.stderr.write("Overlap between %s (%s, %s) and %s (%s, %s). Throwing out %s.\n" % \
                             (item1.name, item1.regions[0][0], item1.regions[-1][1],\
                              item2.name, item2.regions[0][0], item2.regions[-1][1], item2.name))
            ls.pop(i+1)
            continue
        i += 1                 
        
def processSeq(scaf, seq, annotation):  
    #this assumes there is at least one region of each type on each scaffold
    items = {}
    regions = {}
    for k in annotation.keys():
        if scaf not in annotation[k].keys() or not annotation[k][scaf]:
            sys.stderr.write("Warning: no annotated sites of type %s on scaf %s.\n" % (k, scaf))
            items[k] = None
            regions[k] = []
            continue
        items[k] = annotation[k][scaf].pop(0)
        regions[k] = sorted(items[k].regions)
    
    gene = []#for holding the cds bases and positions, in order
    introns = []#same for introns
    
    for i, base in enumerate(seq):
        pos = i+1
        utr3, items["three_prime_UTR"], regions["three_prime_UTR"] = checkType(items["three_prime_UTR"], pos, annotation["three_prime_UTR"], regions["three_prime_UTR"])
        utr5, items["five_prime_UTR"], regions["five_prime_UTR"] = checkType(items["five_prime_UTR"], pos, annotation["five_prime_UTR"], regions["five_prime_UTR"])
        cds, newCDS, regions["CDS"] = checkType(items["CDS"], pos, annotation["CDS"], regions["CDS"], gene, introns)
        
        my_types = []
        if utr3:
            my_types.append(codes['3utr'])
        if utr5:
            my_types.append(codes['5utr'])
        if cds:
            gene.append([pos, base, my_types])
            
        #check if we finished a gene
        if items["CDS"] and (not newCDS or items["CDS"].name != newCDS.name):
            #output the old gene
            processGene(items["CDS"], gene, introns)
            gene = []
            introns = []
            items["CDS"] = newCDS
            
        #check if we were in an intron
        if not cds and gene:
            my_types.append(codes['intron'])
            introns.append([pos, base, my_types])
            
        #if no other types, then it was intergenic
        if not my_types and not cds:
            my_types.append(codes['intergenic'])
        
        if not cds and not gene:
            processBase(scaf, pos, base, my_types)
            
    if gene:
        sys.stderr.write("Warning: There was a gene right at the end of a scaffold, you should make sure none of it was truncated. (%s%s)\n" % (items["CDS"], ""))
        processGene(items["CDS"], gene, introns)
            
def parseArgs():
  parser = argparse.ArgumentParser(description="This takes a gff annotion and a reference genome in fasta format and lines them up annotating every site in the genome.")
  parser.add_argument("reference", type=str, help="the reference genome in fasta format")
  parser.add_argument("annotation", type=str, help="the gff annotation; this file MUST be sorted by position, or many regions will be skipped.")
  return parser.parse_args()

args = parseArgs()
sys.stderr.write("%s\n" % args)

'''takes an item, a position, a type dictionary, and a list of sorted regions
checks if the current position overlaps within the next item in the given type
updates the next item and regions as appropriate
'''
def checkType(item, pos, type, regions, gene=None, introns=None):
    if not type or not regions:
        return False, item, regions
    
    while item and regions[0][1] < pos:
        sys.stderr.write("We missed the beginning of an item. Are you sure that your annotation has no overlapping items?\
         \nSkipping GFFitem %s at pos = %s region = %s.\n" % (item, pos, regions[0]))
        if type[item.scaf]:
            item = type[item.scaf].pop(0)
            regions = sorted(item.regions)
        else:
            item = None
            regions = []
            return False, item, regions
    
    if pos < regions[0][0]:
        return False, item, regions
    
    if pos >= regions[0][0] and pos <= regions[0][1]:
        if pos == regions[0][1]:
            regions.pop(0)
            if not regions:
                if type[item.scaf]:
                    item = type[item.scaf].pop(0)
                    regions = sorted(item.regions)
                else:
                    item = None
                    regions = []
        return True, item, regions
    return False, item, regions

#outputs the given base
def processBase(scaf, pos, base, code, gene = "N", dir = "0"):
    print("%s %s %s %s %s %s" % (scaf, pos, base, gene, dir, ",".join(map(str, code))))

#processes a gene and outputs all bases in it
def processGene(item, cds, introns):
    if len(cds) % 3 != 0:
        sys.stderr.write("Gene length not a multiple of codons (%s len= %s) annotation error?\n" % (item, len(cds)))
        
    #deal with reverse genes, get the reverse compliment
    if item.dir == "-":
        cds.reverse()
        cds = sense(cds)
    
    bases = []
    for i in range(0, len(cds), 3):
        if i+2 >= len(cds):
            sys.stderr.write("Codon not a multiple of 3, setting it to unknown.\n")
            for j in range(i, len(cds)):
                bases.append([cds[i][0], cds[i][1], cds[i][2]+[codes["unknown"]]])
            continue
        
        #grab the codon
        codon = cds[i][1]+cds[i+1][1]+cds[i+2][1]
        #get the appropriate codes
        codon_codes = degeneracy(codon)
        #add them to the list 
        for j in range(3):
            bases.append([cds[i+j][0], cds[i+j][1], cds[i][2]+[codon_codes[j]]])           
            
    #un-compliment the bases if needed
    if item.dir == "-":
        bases = sense(bases)
        
    #add the introns in
    bases += introns
    #sort the bases by position
    bases.sort()
        
    #output the bases
    for pos, allele, code in bases:
        processBase(item.scaf, pos, allele, code, item.name, item.dir)

#returns the opposite strand of a sequence
def sense(bases):
  for i, b in enumerate(bases):
    if b[1] == "A":
      bases[i][1] = "T"
    elif b[1] == "T":
      bases[i][1] = "A"
    elif b[1] == "G":
      bases[i][1] = "C"
    elif b[1] == "C":
      bases[i][1] = "G"
    elif b[1] == "N":
      bases[i][1] = "N"
  return bases

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

  #if the codon isnt 3 bases callit unknown
  if len(codon) < 3:
    res = []
    for i in range (0, len(codon)):
      res.append(codes['unknown'])
    sys.stderr.write("Codon not long enough %s.\n" % codon)
    return res

  #if the codon has an N call it all unknown
  if ("N" in codon):
    return (codes['unknown'],codes['unknown'],codes['unknown'])
  if codon in _codonDict:#otherwise call the codon appropriately
    return _codonDict[codon]
  else:#oh no, the file seems to have non-standard bases
    sys.stderr.write("Unknown codon %s (does your reference have non-standard bases?)\n" % codon)
    #sys.exit()
    return(codes['unknown'],codes['unknown'],codes['unknown'])


'''
Pulls items out of a GFF, each line is pulled and returned as a tuple, except for CDS. CDS from the same gene are all 
put together in the same tuple (each one is added to the regions item in the tuple). 
Note that this parser will not return items in the exact order they appear in the GFF,
CDS will always be the last item returned for a given gene after the UTRs etc.
'''
class GFFparser:
  def __init__(self, file):
    self.file = file
  
  def generator(self):
    GFFitem = namedtuple("GFFitem", "regions type scaf name dir")
    currentItem = None
    for line in self.file:
      if line.startswith("#") or not line:
        continue
      
      sline = line.split()
      scaf = sline[0]
      type = sline[2]
      start = int(sline[3])
      end = int(sline[4])
      dir = sline[6]
      #this is a really gross way of pulling out only the gene number
      #this likely wont work for other genome out of the box
      name = sline[-1].split(";")[0].replace("ID=","").split(".")[0]
      
      #if we've reached a new item we should yield teh old one first
      if currentItem and name != currentItem.name:
        yield currentItem
        currentItem = None
        
      if not currentItem and type == "CDS":
        #if we have no current item, make this the start of the new one
        currentItem = GFFitem([(start, end)], type, scaf, name, dir)
      elif currentItem and currentItem.type == "CDS" and currentItem.type == type: 
        #if this is a gene, and we already have one going add to it
        currentItem.regions.append((start, end))
      else:
        yield(GFFitem([(start, end)], type, scaf, name, dir))
        
    #yield whatever is left over
    if currentItem:
      yield currentItem
if __name__ == "__main__":
  main()
  
  sys.exit(0)
  #tests
  #Testing the GFF parser
  gff = '''##gff-version 3
scaffold_1   JGI_gene    gene  750912 765975 .    -    .    ID=Carubv10008059m.g;Name=Carubv10008059m.g;
scaffold_1   JGI_gene    mRNA  750912 765975 .    -    .    ID=PAC:20892461;Name=Carubv10008059m;pacid=20892461;Parent=Carubv10008059m.g;
scaffold_1   JGI_gene    three_prime_UTR 750912 751458 .    -    .    ID=PAC:20892461.three_prime_UTR.1;Parent=PAC:20892461;pacid=20892461;
scaffold_1   JGI_gene    CDS   751459 752673 .    -    0    ID=PAC:20892461.CDS.19;Parent=PAC:20892461;pacid=20892461;
scaffold_1   JGI_gene    CDS   764128 764263 .    -    1    ID=PAC:20892461.CDS.4;Parent=PAC:20892461;pacid=20892461;
scaffold_1   JGI_gene    CDS   764478 764563 .    -    0    ID=PAC:20892461.CDS.3;Parent=PAC:20892461;pacid=20892461;
scaffold_1   JGI_gene    CDS   764711 764780 .    -    1    ID=PAC:20892461.CDS.2;Parent=PAC:20892461;pacid=20892461;
scaffold_1   JGI_gene    CDS   765203 765372 .    -    0    ID=PAC:20892461.CDS.1;Parent=PAC:20892461;pacid=20892461;
scaffold_1   JGI_gene    five_prime_UTR 765373 765560 .    -    .    ID=PAC:20892461.five_prime_UTR.2;Parent=PAC:20892461;pacid=20892461;
scaffold_1   JGI_gene    five_prime_UTR 765910 765975 .    -    .    ID=PAC:20892461.five_prime_UTR.1;Parent=PAC:20892461;pacid=20892461;
scaffold_1   JGI_gene    gene  8601147 8612398 .    -    .    ID=Carubv10008061m.g;Name=Carubv10008061m.g;
scaffold_1   JGI_gene    mRNA  8601147 8612398 .    -    .    ID=PAC:20890324;Name=Carubv10008061m;pacid=20890324;Parent=Carubv10008061m.g;
scaffold_1   JGI_gene    CDS   8612320 8612397 .    -    0    ID=PAC:20890324.CDS.1;Parent=PAC:20890324;pacid=20890324;
scaffold_1   JGI_gene    five_prime_UTR 8612398 8612398 .    -    .    ID=PAC:20890324.five_prime_UTR.1;Parent=PAC:20890324;pacid=20890324;
'''
  expectedItems = [([(750912, 765975)], 'gene', 'scaffold_1', 'Carubv10008059m', '-'),
           ([(750912, 765975)], 'mRNA', 'scaffold_1', 'PAC:20892461', '-'),
           ([(750912, 751458)], 'three_prime_UTR', 'scaffold_1', 'PAC:20892461', '-'),
           ([(765373, 765560)], 'five_prime_UTR', 'scaffold_1', 'PAC:20892461', '-'),
           ([(765910, 765975)], 'five_prime_UTR', 'scaffold_1', 'PAC:20892461', '-'),
           ([(751459, 752673),(764128, 764263),(764478, 764563),(764711, 764780),(765203, 765372)], 'CDS', 'scaffold_1', 'PAC:20892461', '-'),
           ([(8601147, 8612398)], 'gene', 'scaffold_1', 'Carubv10008061m', '-'),
           ([(8601147, 8612398)], 'mRNA', 'scaffold_1', 'PAC:20890324', '-'),
           ([(8612398, 8612398)], 'five_prime_UTR', 'scaffold_1', 'PAC:20890324', '-'),
           ([(8612320, 8612397)], 'CDS', 'scaffold_1', 'PAC:20890324', '-')]
           
  parser = GFFparser(gff.split("\n"))
  tries = 0
  miss = 0
  for i, item in enumerate(parser.generator()):
    tries += 1
    #print(item)
    if item != expectedItems[i]:
      miss += 1
      sys.stderr.write("FAIL. expected %s \n\t got: %s\n\n" % (expectedItems[i], item))
    else:
      sys.stderr.write("PASS. expected %s\n\t got: %s\n" % (expectedItems[i], item))
      
  sys.stderr.write("Pass/Tests: %s/%s" % ((tries-miss), tries))
    