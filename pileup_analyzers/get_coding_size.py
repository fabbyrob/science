import sys
from collections import defaultdict

infile = open(sys.argv[1], "r")

if len(sys.argv) < 3:
    counts = defaultdict(int)
    for line in infile:
        line = line.split()
        if len(line) != 5:
            continue
        scaf = line[0]
        length = int(line[2])-int(line[1])+1
        counts[scaf] += length
        
    for k in counts.keys():
        #if k.split("_")[1] in ["1","2","3","4","5",'6',"7",'8']:
        print(k+"\t"+str(counts[k]))
else:
    lengths = {}
    for line in infile:
        line = line.split()
        if len(line) < 3:
            continue
        length = line[2].split("LN:")
        scaf = line[1].split("SN:")
        if len(length) == 1 or len(scaf) == 1:
            continue
        lengths[scaf[1]] = int(length[1])
    for k in lengths.keys():
        print(k+"\t"+str(lengths[k]))