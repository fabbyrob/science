'''
Requires python 2.7 or earlier

A parser for samtools generated pileup files. Makes a reader that returns one line at a time. 

    >>> import pileup
    >>> pileup_reader = pileup.Reader(open('example.pileup','r'))
    >>> for record in pileup_reader:
    >>>    print record
    
Record has the following attributes:

    Record.CHROM
    Record.POS
    Record.REF
    Record.ALT
    Record.QUAL
    Record.SNPQUAL
    Record.MAPQUAL
    Record.DEPTH
    Record.READS
    REcord.READ_QUALS #this item will be a list, rather than with a string, since INDEL lines have extra stuff
'''
class Record(object):
    def __init__(self, chrom, pos, ref, alt, qual, snpqual, mapqual, depth, reads, read_quals):
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALT = alt
        self.QUAL = qual
        self.SNPQUAL = snpqual
        self.MAPQUAL = mapqual
        self.DEPTH = depth 
        self.READS = reads
        self.READ_QUALS = read_quals
        
    def __str__(self):
        return "Record(chrom="+self.CHROM+", pos="+str(self.POS)+", ref="+self.REF+", alt="+self.ALT+")"
    
class Reader(object):
    def __init__(self, filename):
        self.reader = filename
        
    def __iter__(self):
        return self
    
    def next(self):
        line = self.reader.next().split()
        
        chrom = line[0]
        pos = int(line[1])
        ref = line[2]
        alt = line[3]
        qual = int(line[4])
        snpqual = int(line[5])
        mapqual = int(line[6])
        depth = int(line[7])
        reads = line[8]
        read_quals = line[9:]
        
        record = Record(chrom, pos, ref, alt, qual, snpqual, mapqual, depth, reads, read_quals)
        
        return record
    
#TEST
#file = open('io/pile.txt', 'rb')
#myReader = Reader(file)
#print(myReader.next())
#for r in myReader:
#    print(r)
        
    