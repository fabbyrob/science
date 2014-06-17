import sys

infile = open(sys.argv[1])

def getStartStop(lab, line):
    sline = line.split()
    
    chrom = sline[0]
    start = sline[3]
    end = sline[4]
    dir = sline[6]
    sname = sline[-1].split("=")[-1][:-1]
    name = "PAC:"+sname
    
    print("%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, dir, name, lab))

for line in infile:
    line = line.rstrip()
    
    if "three_prime_UTR" in line:
        getStartStop("3utr", line)
    elif "five_prime_UTR" in line:
        getStartStop("5utr", line)
        