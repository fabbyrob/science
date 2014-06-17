import sys

infile = open(sys.argv[1])

for line in infile:
    line = line.rstrip()
    sline = line.split()
    
    if len(sline) < 9:
        sys.stderr.write("weird line: "+line+"\n")
        continue
    
    scaf = sline[0]
    type = sline[2]
    start = sline[3]
    end = sline[4]
    dir = sline[6]
    rest = sline[8:]
    
    if type != "CDS":
        continue
    
    name = rest[1].replace("\"", "").replace(";","")
    
    print("%s\t%s\t%s\t%s\t%s" % (scaf, start, end, dir, name))