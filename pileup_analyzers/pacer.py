infile = open("io/sporo.pac.txt", "r")
outfile = open("io/sporo.full.pac.txt", "w")

for line in infile:
    line = line.rstrip()
    outfile.write("PAC:%s\n" % line)