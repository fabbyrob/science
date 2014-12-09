import sys
import summary
import argparse

doc = '''takes a summary file with individual genotypes and a text file listing a subset of the samples. Splits the file into two summaries with different individuals.'''

_args = None

def __main__():
    processArgs()
    
    reader = summary.Reader(open(_args.summary))
    samples = []
    for line in open(_args.samples):
        line = line.rstrip()
        if line not in reader.Samples:
            sys.stderr.write("Sample '%s' not found in the Summary. Did you misspell something? Ignoring this sample name.\n" % line)
            continue
        samples.append(line)

    other_samples = []
    for s in reader.Samples:
        if s not in samples:
            other_samples.append(s)

    all_samples = reader.Samples[:]
    file_one = open(_args.file_one, "w")
    file_two = open(_args.file_two, "w")
    
    #generate the headers
    header = reader.Header
    sheader = header.split("\n")
    new_header = "\n".join(sheader[0:-1])
    file_one.write(new_header)
    file_two.write(new_header)
    fields = sheader[-1].split()
    fields = "\t".join(fields[0:len(fields)-len(reader.Samples)])
    file_one.write(fields+"\t"+"\t".join(samples)+"\n")
    file_two.write(fields+"\t"+"\t".join(other_samples)+"\n")
    
    #outputs the lines
    for record in reader:
        filtered = not record.TOTAL

        if len(record.Genos) != len(reader.Samples):#weird site
            sys.stderr.write("ERROR: missing genotypes from site. Skipping. (%s %s)\n" % (record.CHROM, record.POS))
            continue

        #get the allele counts
        allele_counts = [0,0,0]#REF, ALT, TOTAL
        
        #output the lines
        newgenos = []
        for s in samples:
            newgenos.append(record.Genotypes[s])
            if record.summary.Genotypes[newgenos[-1]] == "heterozygote":
                allele_counts[0] += 1
                allele_counts[1] += 1
            elif record.summary.Genotypes[newgenos[-1]] == "homozygote reference":
                allele_counts[0] += 2
            elif record.summary.Genotypes[newgenos[-1]] == "homozygote alternate":
                allele_counts[1] += 2

            if not filtered:
                allele_counts[2] = allele_counts[0]+allele_counts[1]
       
        record.Genos = newgenos 
        record.REF_NUM, record.ALT_NUM, record.TOTAL = allele_counts

        file_one.write(str(record)+"\n")

        #get the allele counts
        allele_counts = [0,0,0]#REF, ALT, TOTAL

        #output the lines
        newgenos = []
        for s in other_samples:
            newgenos.append(record.Genotypes[s])
            if record.summary.Genotypes[newgenos[-1]] == "heterozygote":
                allele_counts[0] += 1
                allele_counts[1] += 1
            elif record.summary.Genotypes[newgenos[-1]] == "homozygote reference":
                allele_counts[0] += 2
            elif record.summary.Genotypes[newgenos[-1]] == "homozygote alternate":
                allele_counts[1] += 2

            if not filtered:
                allele_counts[2] = allele_counts[0]+allele_counts[1]
       

        record.REF_NUM, record.ALT_NUM, record.TOTAL = allele_counts
        record.Genos = newgenos
        file_two.write(str(record)+"\n")

    file_one.close()
    file_two.close()

def processArgs():
    parser = argparse.ArgumentParser(description=doc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("summary", help="summary file that will be split up")
    parser.add_argument("samples", help="test file listing samples to split out, one sample name per line")
    parser.add_argument("-o","--file_one", help="the file name where the new summary, with samples from the input, is saved", type=str, default="samples.summary")
    parser.add_argument("-t","--file_two", help="the file where the the new summary with the remaining samples is saved", type=str, default="other_samples.summary")
    global _args
    _args = parser.parse_args()


if __name__ == '__main__':
    __main__()
