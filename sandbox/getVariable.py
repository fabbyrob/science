
for file in sys.argv:
    infile = open(file, "r")
    
    seqs = []

    for line in infile:
        if line[0] == ">":
            continue
        seqs.append(line)
        
    diffs = []
    for i in range(0, len(seqs[0])):
        same = True
        first = seqs[0][i]
        
        for j in range(1, len(seqs)):
            if seqs[j][i] != first:
                same = False
        
        diffs.append(i)
        
    
    outfile = open(file.split("/")[-1]+".variable", "w")
    
    outfile.write(str(diffs)+"\n")
    outfile.close()

