doc = '''
Takes a hapcut file that specifies haplotypes in an individual, and a summary file. Makes a fasta with one entry for every haplotype.
'''

import sys, argparse
import hapcut, summary

_args = None

def __main__():
    processArgs()
    hapcutReader = hapcut.Reader(open(_args.hapcut))
    summaryReader = summary.Reader(open(_args.summary))
    
    block = hapcutReader.next()
    pSite = None
    for site in summaryReader:
        if site.CHROM != block.chrom:
            sys.stderr.write("Not on the correct chromosome. Skipping ahead. (Summary: %s Haplotype: %s)\n" % (site.CHROM, block.chrom))
            continue

        if site.POS >= block.start and site.POS <= block.end:#within the block
            #add the site to the haplotypes
            if summaryReader.Genotypes[site.Genotypes[_args.name]] == "heterozygous":
                if site.POS not in block.SNPs:#an unphased het, make it N in the fasta
                    block.addSNP(site.CHROM, site.POS, "N", "N", 1, 0, 0)
            elif summaryReader.Genotypes[site.Genotypes[_args.name]] == "homozygous reference":
                block.addSNP(site.CHROM, site.POS, site.REF, site.ALT, 0, 0, 0)
            elif summaryReader.Genotypes[site.Genotypes[_args.name]] == "homozygous alternate":
                block.addSNP(site.CHROM, site.POS, site.REF, site.ALT, 1, 1, 0)
            else:#N
                continue

        if site.POS == block.end:#output this block
            outputBlock(block)
            block = hapcutReader.next()
            if block == None:
                break

        if site.POS > block.end:#we missed something
            #add N's to fill it out
            fillN(block, pSite, site.POS)
            outputBlock(block)
            block = hapcutReader.next()
            if block == None:
                break

        pSite = site.POS

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
    parser = argparse.ArgumentParser(decription=doc, formatter_class=argparse.ArgumentDefailtsHelpFormatter)
    parser.add_argument("hapcut", help="The hapcut output file.")
    parser.add_argument("summary", help="The summary file with individual genotypes.")
    parser.add_argument("name", help="The sample name as specified in the summary file.")
    global _args
    _args = parser.parse_args()

if __name__ == '__main__':
    __main__()
