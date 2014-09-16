import sys

infile = open(sys.argv[1],"r")

quals = {}
for i in range(35,74):
    quals[i] = 0

ready = False
for line in infile:
    line = line.rstrip()
    if len(line) == 0 :
        continue
    
    if line[0] == "+":
        ready = True
        continue
    
    if ready:
        for q in line:
            q = ord(q)
            if q in quals.keys():
                quals[q] += 1
            else:
                quals[q] = 1
        ready = False
        
keys = quals.keys()
keys.sort()
for k in keys:
    print("%s %s" % (k, quals[k]))