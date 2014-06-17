import sys

if ".sam" in sys.argv[1]:#sam file
    samfile = open(sys.argv[1],"r")
    
    for line in samfile:
        line = line.rstrip()
        if line.startswith("@SQ"):
            sline = line.split()
            name = sline[1].split("|")[0]
            sline[1] = name
            line = "\t".join(sline)
            print(line)
elif ".fasta" in sys.argv[1]:#reference
    reference = open(sys.argv[1],"r")
    
    for line in reference:
        line = line.rstrip()
        if line.startswith(">"):
            sline = line.split("|")
            line = sline[0]
        print(line)
else:
    sys.stderr.write("Incorrect file type.")
    sys.exit(0)
            