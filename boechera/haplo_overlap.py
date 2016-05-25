details = '''
Takes in 2 hapcut files and a VCF, and the sample names the two hapcut files come from. Finds overlap in haplotypes, and calculates differences in the overlap haplotypes.
ONly calculates differences between the 1st haplotype in file 1 and the haplotypes from file 2.
'''

import sys
sys.path.append("/data/robert.williamson/bin")

import vcf
import time
from hapcut_parser import Reader, Block
from collections import defaultdict

def __main__():
    start_time = time.time()
    #print(sys.argv)
    hap1 = Reader(open(sys.argv[1]))
    hap2 = Reader(open(sys.argv[2]))
    
    file1Blocks = defaultdict(list)
    file2Blocks = defaultdict(list)
    
    #read in the blocks
    sys.stderr.write("Reading in haplotypes...\n")
    for block in hap1:
        file1Blocks[block.chrom].append(block)
    for block in hap2:
        file2Blocks[block.chrom].append(block) 
        
    sys.stderr.write("Sorting haplotypes 1 (Run time: %.2f)...\n" % ((time.time()-start_time)/60))
    #make sure they are sorted
    for k in file1Blocks.keys():
        file1Blocks[k].sort()
    sys.stderr.write("Sorting haplotypes 2 (Run time: %.2f)...\n" % ((time.time()-start_time)/60)) 
    for k in file2Blocks.keys():
        file2Blocks[k].sort()

    #read in the annotation, if it is there
    mySites = None
    if len(sys.argv) == 7:
        mySites = defaultdict(dict)
        annot = open(sys.argv[6])
        for line in annot:
            if line.startswith("{"):
                continue
            sline = line.split()
            mySites[sline[0]][sline[1]] = sline[-1]

    sys.stderr.write("Processing VCFs (Run time: %.2f)...\n" % ((time.time()-start_time)/60)) 
    ctr = 0
    #add in all the homozygous data
    currf1 = 0
    currf2 = 0
    currChrom = ""
    vReader = vcf.Reader(open(sys.argv[3]))
    
    f1Samp = vReader.samples.index(sys.argv[4])
    f2Samp = vReader.samples.index(sys.argv[5])
    
    for site in vReader:
        ctr += 1
        if ctr % 50000 == 0:
            sys.stderr.write("\tProcessed %s VCF lines (Run time: %.2f)...\n" % (ctr, (time.time()-start_time)/60))
        if not currChrom:
            currChrom = site.CHROM
            
        if currChrom != site.CHROM:
            currChrom = site.CHROM
            currf1 = 0
            currf2 = 0
            
        if file1Blocks[currChrom] and file1Blocks[currChrom][currf1].end < site.POS:
            currf1 = min(currf1+1, len(file1Blocks[currChrom])-1)
            
        if file2Blocks[currChrom] and file2Blocks[currChrom][currf2].end < site.POS:
            currf2 = min(currf2+1, len(file2Blocks[currChrom])-1)
            
        if site.ALT[0] != "." and file1Blocks[currChrom]:#skip non-SNP sites, they wont be different
            #sys.stderr.write("%s\n%s\n%s\n" % (site, file1Blocks[currChrom][currf1], file1Blocks[currChrom][currf1].haplotype1))    
            try:
                insertSite(site, file1Blocks[currChrom][currf1], f1Samp)
                insertSite(site, file2Blocks[currChrom][currf2], f1Samp)  
            except IndexError:
                sys.stderr.write("Index error. Chrom: %s Block 1 len: %s Block 2 len: %s. Index 1: %s Index 2: %s.\n" % (currChrom, len(file1Blocks[currChrom]), len(file2Blocks[currChrom]), currf1, currf2)) 
                sys.stderr.write("Skipping line...")
            #sys.stderr.write("%s\n" % file1Blocks[currChrom][currf1].haplotype1)     
        
    #find overlaps
    sys.stderr.write("Finding overlaps (Run time: %.2f)...\n" % ((time.time()-start_time)/60)) 
    print("CHROM\tSTART\tEND\tDiff1\tDiff2\tDiff3\tDiff4\tDiff5\tDiff6")
    #diffs are calculated on line 107
    for chrom in file1Blocks.keys():
        if not file2Blocks.has_key(chrom) or not file2Blocks[chrom]:
            continue
        
        for b1 in file1Blocks[chrom]:
            if not file2Blocks[chrom]:
                break
            while file2Blocks[chrom] and file2Blocks[chrom][0].start < b1.start and file2Blocks[chrom][0].end < b1.start:
                file2Blocks[chrom].pop() 
                
            for b2 in file2Blocks[chrom]:
                if b2.start > b1.end:
                    break
                
                over = overlap(b1, b2, mySites)
                if over and not mySites:
                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % over)
  
#finds the overlap between 2 blocks and the # diffs between haplotypes
#returns the overlap and diffs or None if no overlap  
def overlap(B1, B2, myAnnot):
    start = max(B1.start, B2.start)
    end = min(B1.end, B2.end)
    scaf = B1.chrom
    b1Hap1, b1Hap2 = B1.getOverlap(start, end)
    b2Hap1, b2Hap2 = B2.getOverlap(start, end)
    
    if myAnnot:
        perTypeDiff(b1Hap1, b1Hap2, scaf, start, end, sys.argv[4]+"1", sys.argv[4]+"2", myAnnot)
        perTypeDiff(b1Hap1, b2Hap1, scaf, start, end, sys.argv[4]+"1", sys.argv[5]+"1", myAnnot)
        perTypeDiff(b1Hap1, b2Hap2, scaf, start, end, sys.argv[4]+"1", sys.argv[5]+"2", myAnnot)
        perTypeDiff(b1Hap2, b2Hap1, scaf, start, end, sys.argv[4]+"2", sys.argv[5]+"1", myAnnot)
        perTypeDiff(b1Hap2, b2Hap2, scaf, start, end, sys.argv[4]+"2", sys.argv[5]+"2", myAnnot)
        perTypeDiff(b2Hap1, b2Hap2, scaf, start, end, sys.argv[5]+"1", sys.argv[5]+"2", myAnnot)
        return

    diff1 = B1.difference(b1Hap1, b2Hap1)
    diff2 = B1.difference(b1Hap1, b2Hap2)
    diff3 = B1.difference(b1Hap1, b1Hap2)
    diff4 = B1.difference(b1Hap2, b2Hap1)
    diff5 = B1.difference(b1Hap2, b2Hap2)
    diff6 = B1.difference(b2Hap1, b2Hap2)
    
    return (B1.chrom, start, end, diff1, diff2, diff3, diff4, diff5, diff6) if end - start > 0 else None

def perTypeDiff(b1, b2, scaf,  start, end, n1, n2, myAnnot):
    diffs = defaultdict(int)
    i = start
    for x, y in zip(b1, b2):
        d = 1 if x == y else 0 
        for t in myAnnot[scaf][i]:
            diffs[t] += d
        i += 1

    for t in diffs.keys():
        print("%s %s %s %s %s %s %s" % (scaf, start, end, n1, n2, t, diffs[t]))

def insertSite(site, block, samp):
    if site.POS >= block.start and site.POS <= block.end:
        if site.POS in block.SNPs:
            #already in there, presumably a het ignore this site
            return
        geno = site.samples[samp]
        
        #sys.stderr.write("%s\n%s\n%s\n" % (site, block, block.haplotype1))
        
        if "GT" in geno.keys() and "." not in geno["GT"]:
            g1, g2 = map(int, geno["GT"].split("/"))
            if g1 == g2:#dont add heterozygotes... dont know the haplotype
                block.addSNP(site.CHROM, site.POS, site.REF, site.ALT[0], g1, g2, "")
                #sys.stderr.write("%s\n" % (block.haplotype1))
                return
        #don't know the site genotype add "N"
        block.addSNP(site.CHROM, site.POS, site.REF, site.ALT, -1, -1, "")
        #sys.stderr.write("%s\n" % (block.haplotype1))

use = "python haplo_overlap.py myHap1.hapcut myHap2.hapcut myPooled.vcf Sample1 Sample2"

if __name__ == "__main__":
    __main__()
