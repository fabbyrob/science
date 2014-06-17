import vcf
import sys

vcf_reader = vcf.Reader(open(sys.argv[1], 'rb'))

print("CHROM\tPOS\tSAMP\tGENO")
for record in vcf_reader:
    for samp in record.samples:
        if samp['GT'] == "0/0":
            geno = record.REF[0]*2
        elif samp['GT'] == "0/1":
            geno = record.REF[0]+record.ALT[0]
        else:
            geno = record.ALT[0]*2
            
        print("%s\t%s\t%s\t%s" % (record.CHROM, record.POS, samp['name'], geno))