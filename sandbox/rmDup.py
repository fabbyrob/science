import sys

infile = open(sys.argv[1])
prev = None

for line in infile:
    line=line.rstrip()
    if not prev:
        print(line)
        prev = line.split()
        continue
    else:
        sline = line.split()
        if sline[0] == prev[0] and sline[1] == prev[1]:
            continue
        print(line)
        prev = sline