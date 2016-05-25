import sys, argparse
#you may need to add your path to summary.py here
sys.path.append("/data/robert.williamson/bin")
import summary, AFSstats, foldAFS

def main():
    reader = summary.Reader(open(args.summary,'rb'))
    if  reader.summary.Samples:
        N = reader.summary.Ploidy*len(reader.summary.Samples)
    elif not args.sample_size:
        sys.stderr.write("No samples listed in this summary. You must provide a sample size with -N.\n")
        sys.exit()
    else:
        N = args.sample_size

    genAFS = lambda : [0]*(N+1)

    if "GENE" not in reader.summary.Fields:
        sys.stderr.write("The gene names are not in this summary. All gene names will be blank if you use this file.\n")
        

    try:
        intronCode = reader.summary.typeToCode['intron']
    except KeyError:
        sys.stderr.write("Summary has no 'intron' codes, you need to add them or use a different summary.\n")
        sys.exit()

    #output header
    outputIntron(N=N)

    currentIntron = None
    for record in reader:
        if intronCode in record.Types:
            if not currentIntron:
                #new one
                currentIntron = Intron(record.GENE, record.CHROM, record.POS, genAFS())
            elif currentIntron.gene != record.GENE:
                #end of previous one, we're on a new one
                #this should really not happen
                sys.stderr.write("Reached a new intron adjacent to an old one at %s %s, you may want to check your annotation.\n" % (record.CHROM, record.POS))
                outputIntron(N, currentIntron)
                currentIntron = Intron(record.GENE, record.CHROM, record.POS, genAFS())
    #        print(record.POS, currentIntron.afs)
            #doesn't count AF at sites with too few alleles
            if record.TOTAL == N:
                try:
                    currentIntron.afs[record.ALT_NUM] += 1 
                except IndexError:
                    sys.stderr.write("At %s %s there is a problem with the ALT_NUM. It is out of acceptable ranges. Is there a problem with your summary? Or with the provided sample size?\n" % (record.CHROM, record.POS))
            currentIntron.end = record.POS
        elif currentIntron:
            #we reached the end on this one
            outputIntron(N, currentIntron)
            currentIntron = None
    #get the last one
    if currentIntron:
        outputIntron(N, currentIntron)

class Intron:
    def __init__(self, gene, scaf, start, afs):
        self.gene = gene
        self.scaf = scaf
        self.start = start
        self.end = start
        self.afs = afs        

def outputIntron(N, intron=None):
    if intron:
        s = "%s %s %s %s " % (intron.gene, intron.scaf, intron.start, intron.end)
        folded = foldAFS.foldAFS(intron.afs)
        theta, pi, d = AFSstats.stats(folded, N)
        s += "%s %s %s " % (pi, theta, d)
        s += " ".join(map(str, intron.afs))
    else:
        s = "Gene Scaffold Start End Pi Theta_W Tajimas_D "
        s += " ".join(["AFS_%s" % i for i in range(N+1)])
    print(s)

def parseArgs():
    parser = argparse.ArgumentParser(description="Takes a summary with genes listed, outputs an AFS for each Intron listed in the file, as well as summary statistics for each one (pi, theta, D).\nThis program skips sites where the number of alleles is not equal to the samplesize.")
    parser.add_argument("summary", help="the summary file to read data from")
    parser.add_argument("-N", "--sample-size", type=int, help="The sample size (in number of alleles), only needed if the summary does not list samples individually.") 
    return parser.parse_args()


if __name__=="__main__":
    args = parseArgs()
    main()

