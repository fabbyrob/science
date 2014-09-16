import sys

infile = open(sys.argv[1])
infile.readline()

type = sys.argv[2]

for line in infile:
    sline = line.split()
    
    scaf = sline[0]
    pos = sline[1]
    cat = sline[-1]
    
    if cat == type:
        print(scaf+"\t"+pos)