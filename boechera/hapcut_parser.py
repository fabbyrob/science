details = '''
This is a parser for output from HapCut. The reader reads in blocks individually from the file and
stores them as Block objects.

Basic Usage:
>from hapcut_parser import Reader, Block
>
>myReader = Reader(open("myOutput.hapcut"))
>
>for block in myReader:
>    print(block)
 
Reader:
    Reader.hapcut - the hapcut file being read
    
Block:
    Block.offset
    Block.len
    Block.phased        These are all the metrics from the BLOCK line in the hapcut file
    Block.SPAN
    Block.MECscore
    Block.Fragments
    
    Block.chrom - the chromosome the block is on
    Block.start - the position of the first SNP in the block 
    Block.end - the position of the last SNP in the block
    
    Block.SNPs - a list of positions of all the SNPs that define the block
    Block.genotypes - a dictionary mapping Positions onto the Ref and Alt alleles for each SNP
        Block.genotypes[576].REF or Block.genotypes[765].ALT
        REF and ALT are parts of a named tuple, and are immutable
    Block.quals - a dictionary mapping Positions onto the genotype and quality columns from the file
        TODO: expand this
        
    Block.haplotype1 / Block.haplotype2 - these attributes are lists specifying each haplotype (i.e. ["A","T","T"..])
    
    Block.difference(hap1, hap2) - this function counts and returns the number of differences between two haplotypes. By default it returns the differences of
        the two haplotypes of self
    Block.getOverlap(start, end) - returns 2 lists of the haplotype nucleotides that fall within start and stop
'''

from collections import namedtuple
import sys

#sets up the named tuple for genotypes later
Alleles = namedtuple("Alleles", ["REF","ALT"])

#this class reds in teh files, and builds the Blocks
class Reader:
    def __init__(self, hapcut):
        self.hapcut = hapcut
        
    def __iter__(self):
        myBlock = None
        for line in self.hapcut:
            if line.startswith("*"):
                yield myBlock
                myBlock = None
                continue
            sline = line.split()

            if line.startswith("BLOCK"):
                if myBlock:
                    raise(UnfinishedBlock("A new Block started in the file before the end of the previous one was delimited (Old block: %s).\n\tDo your blocks have \'********\' on the lines between them?" % myBlock))
                myBlock = Block(sline[2], sline[4], sline[6], sline[8], sline[10], sline[12])
            else:
                if not myBlock:
                    raise(MissingBlock("Trying to add to a block that hasn't been started. Is the file formatted correctly? Line: %s" % (line)))
                myBlock.addSNP(sline[3], int(sline[4]), sline[5], sline[6], int(sline[1]), int(sline[2]), " ".join(sline[7:]))
            
        yield myBlock

#this class holds all the data from a single block
class Block:
    def __init__(self, offset, length, phased, SPAN, MECscore, fragments):
        self.offset = int(offset)
        self.length = int(length)
        self.phased = int(phased)
        self.SPAN = int(SPAN)
        self.MECscore = float(MECscore)
        self.fragments = int(fragments)
        
        self.chrom = None
        self.SNPs = []
        self.start = None
        self.end = None
        self.genotypes = {}
        self.quals = {}
        self.haplotype1 = []
        self.haplotype2 = []
        
    #adds a SNP to this haplotype block
    def addSNP(self, chrom, pos, REF, ALT, hap1, hap2, qual):
        if not self.chrom:
            self.chrom = chrom
        
        #make sure we're on the same chromosome
        if self.chrom != chrom:
            raise ChromosomeMissmatch("Chromosomes do not match within Block (%s and %s)." % (self.chrom, chrom))
            
        #make sure this position hasn't been done already
        if pos in self.SNPs:
            raise DuplicatePosition("The position %s is already present in this block." % pos)
        
        self.SNPs.append(pos)
        self.SNPs.sort() # This could be slow for big haplotypes...
        i = self.SNPs.index(pos)
        
        self.start = self.SNPs[0]
        self.end = self.SNPs[-1]
                
        self.genotypes[pos] = Alleles(REF, ALT)
        
        self.quals[pos] = qual
        
        if hap1 == 1:
            self.haplotype1.insert(i, ALT)
        elif hap1 == 0:
            self.haplotype1.insert(i, REF)
        else:
            self.haplotype1.insert(i, "N")
            
        if hap2 == 1:
            self.haplotype2.insert(i, ALT)
        elif hap2 == 0:
            self.haplotype2.insert(i, REF)
        else:
            self.haplotype2.insert(i, "N")
    
    def getOverlap(self, oStart, oEnd):
        overlap1 = []
        overlap2 = []
        for i, s in enumerate(self.SNPs):
            if s > oEnd:#we've reached the end
                break
            
            if s >= oStart:
                overlap1.append(self.haplotype1[i])
                overlap2.append(self.haplotype2[i])
        
        return overlap1, overlap2
        
    def difference(self, hap1, hap2):
        diffs = 0.0
        
        for h1, h2 in zip(hap1, hap2):
            if h1 != h2 and "N" not in [h1, h2]:
                diffs += 1
        
        return diffs
    
    def __str__(self):
        mystr = "BLOCK(chrom: %s start: %s offset %s len: %s SPAN: %s)" % (self.chrom, self.start, self.offset, self.length, self.SPAN)
        return mystr
    
    def __lt__(self, other):
        return self.start < other.start
    
    def __gt__(self, other):
        return self.start > other.start
    
####Some custom exceptions
class ChromosomeMissmatch(Exception):
    pass

class DuplicatePosition(Exception):
    pass

class UnfinishedBlock(Exception):
    pass

class MissingBlock(Exeption):
    pass

if __name__ == "__main__":
    print("Testing Parser")
    
    myReader = Reader(open("test.hapcut"))
    
    for block in myReader:
        print("%s\n\t%s\n\t%s" % (block, "".join(block.haplotype1), "".join(block.haplotype2)))
        