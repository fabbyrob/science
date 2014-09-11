details = '''
Takes in 2 hapcut files and a VCF, and the sample names the two hapcut files come from. Finds overlap in haplotypes, and calculates differences in the overlap haplotypes.
ONly calculates differences between the 1st haplotype in file 1 and the haplotypes from file 2.
'''

import sys
sys.path.append("/data/robert.williamson/bin")

import vcf
from hapcut_parser import Reader, Block
from collections import defaultdict

def __main__():
    hap1 = Reader(open(sys.argv[1]))
    hap2 = Reader(open(sys.argv[2]))
    
    file1Blocks = defaultdict(list)
    file2Blocks = defaultdict(list)
    
    #read in the blocks
    for block in hap1:
        file1Blocks[block.chrom].append(block)
    for block in hap2:
        file2Blocks[block.chrom].append(block) 
        
    #make sure they are sorted
    for k in file1Blocks.keys():
        file1Blocks[k].sort()
    for k in file2Blocks.keys():
        file2Blocks[k].sort()
        
    #add in all the homozygous data
    currf1 = 0
    currf2 = 0
    currChrom = ""
    vReader = vcf.Reader(open(sys.argv[3]))
    
    f1Samp = vReader.sampled.index(sys.argv[4])
    f2Samp = vReader.sampled.index(sys.argv[5])
    
    for site in vReader:
        if not currChrom:
            currChrom = site.CHROM
            
        if currChrom != site.CHROM:
            currChrom = site.CHROM
            currf1 = 0
            currf2 = 0
            
        if file1Blocks[currChrom][currf1].end < site.POS:
            currf1 = min(currf1+1, len(currf1)-1)
            
        if file1Blocks[currChrom][currf2].end < site.POS:
            currf2 = min(currf2+1, len(currf2)-1)
            
        insertSite(site, file1Blocks[currChrom][currf1], f1Samp)
        insertSite(site, file2Blocks[currChrom][currf2], f1Samp)        
        
    #find overlaps
    print("START\tEND\tDiff1\tDiff2")
    for chrom in file1Blocks.keys():
        if not file2Blocks.has_key(chrom) or not file2Blocks[chrom]:
            continue
        
        for b1 in file1Blocks[chrom]:
            if not file2Blocks[chrom]:
                break
            while file2Blocks[chrom][0].start < b1.start and file2Blocks[chrom][0].end < b1.start:
                file2Blocks[chrom].pop() 
                
            for b2 in file2Blocks[chrom]:
                if b2.start > b1.end:
                    break
                
                over = overlap(b1, b2)
                if over:
                    print("%s\t%s\t%s\t%s" % over)
  
#finds the overlap between 2 blocks and the # diffs between haplotypes
#returns the overlap and diffs or None if no overlap  
def overlap(B1, B2):
    start = max(B1.start, B2.start)
    end = min(B1.end, B2.end)
    b1Hap1, = B1.getOverlap(start, end)
    b2Hap1, b2Hap2 = B2.getOverlap(start, end)
    
    diff1 = B1.difference(b1Hap1, b2hap1)
    diff2 = B1.difference(b1Hap1, b2Hap2)
    
    return (start, end, diff1, diff2) if end - start > 0 else None

def insertSite(site, block, samp):
    if site.POS >= block.start and site.POS <= block.end:
        if site.POS in block.SNPs:
            #already in there, presumably a het ignore this site
            return
        geno = site.samples[samp]
        
        if "GT" in geno.keys():
            g1, g2 = map(int, geno["GT"].split("/"))
            if g1 == g2:#dont add heterozygotes... dont know the haplotype
                block.addSNP(site.CHROM, site.POS, site.REF, site.ALT, g1, g2, "")
                return
        #don't know the site genotype add "N"
        block.addSNP(site.CHROM, site.POS, site.REF, site.ALT, -1, -1, "")

use = "python haplo_overlap.py myHap1.hapcut myHap2.hapcut myPooled.vcf Sample1 Sample2"

if __name__ == "__main__":
    __main__()