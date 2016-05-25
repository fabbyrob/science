doc = '''
Takes a hapcut file that specifies haplotypes in an individual, and a summary file. Makes a fasta with one entry for every haplotype.
'''

import sys, argparse
sys.path.append("/data/robert.williamson/bin")
import hapcut_parser as hapcut
import summary

_args = None

def __main__():
    processArgs()
    hapcutReader = hapcut.Reader(open(_args.hapcut))
    summaryReader = summary.Reader(open(_args.summary))
    
    myGenerator = getNextBlock(hapcutReader)
    block = myGenerator.next()
    pSite = None
    #print(summaryReader.Genotypes)
    for site in summaryReader:
        if site.CHROM != block.chrom:
            sys.stderr.write("Not on the correct chromosome. Skipping ahead. (Summary: %s Haplotype: %s)\n" % (site.CHROM, block.chrom))
            continue

        if site.POS >= block.start and site.POS <= block.end and _args.name in site.Genotypes.keys():#within the block
            #add the site to the haplotypes
            #sys.stderr.write("WITHIN BLOCK %s %s.\n" % (site.CHROM, site.POS))
            #sys.stderr.write("%s %s\n" % (block.start, block.end))
            #sys.stderr.write("%s %s %s\n" % (site.CHROM, site.POS ,"".join(block.haplotype1)))
            #sys.stderr.write("%s\n" % summaryReader.Genotypes[site.Genotypes[_args.name]])
            if site.POS in block.SNPs:
                #site already phased NOTE the same qual cutoffs may NOT have happened here
                continue


            if len(site.REF) > 1 or len(site.ALT) > 1:
                block.addSNP(site.CHROM, site.POS, "N", "N", 1, 0, 0)
            if summaryReader.Genotypes[site.Genotypes[_args.name]] == "heterozygote":
                block.addSNP(site.CHROM, site.POS, "N", "N", 1, 0, 0)
            elif summaryReader.Genotypes[site.Genotypes[_args.name]] == "homozygote reference":
                block.addSNP(site.CHROM, site.POS, site.REF, site.ALT, 0, 0, 0)
            elif summaryReader.Genotypes[site.Genotypes[_args.name]] == "homozygote alternate":
                block.addSNP(site.CHROM, site.POS, site.REF, site.ALT, 1, 1, 0)
            else:#N
                continue

            #sys.stderr.write("%s %s %s\n" % (site.CHROM, site.POS, "".join(block.haplotype1)))

        if site.POS == block.end:#output this block
            outputBlock(block)
            block = myGenerator.next()
            if block == None:
                break

        if site.POS > block.end:#we missed something
            #add N's to fill it out
            fillN(block, pSite, site.POS)
            outputBlock(block)
            block = myGenerator.next()
            if block == None:
                break

        pSite = site.POS

def getNextBlock(reader):
    for block in reader:
        yield block

def fillN(block, pSite, pos):
    if len(block.haplotype1) == block.end - block.start + 1:
        return

    for i in range(block.start, block.end):
        if i not in block.SNPs:
            block.addSNP(block.chrom, i, 'N', 'N', 1, 1, 0)

def outputBlock(block):
    header = ">%s %s:%s-%s " % (_args.name, block.chrom, block.start, block.end)
    print(header + "haplotype1")
    print("".join(block.haplotype1)+"\n")
    print(header + "haplotype2")
    print("".join(block.haplotype2)+"\n")

    if len(block.haplotype1) != block.end-block.start + 1 or len(block.haplotype2) != block.end-block.start + 1:
        sys.stderr.write("Warning: haplotype at %s %s is not the correct length. Were some sites missed?\n" % (block.chrom, block.start))

def processArgs():
    parser = argparse.ArgumentParser(description=doc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("hapcut", help="The hapcut output file.")
    parser.add_argument("summary", help="The summary file with individual genotypes.")
    parser.add_argument("name", help="The sample name as specified in the summary file.")
    global _args
    _args = parser.parse_args()

if __name__ == '__main__':
    __main__()
