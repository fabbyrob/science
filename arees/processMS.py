import sys

infile = open(sys.argv[1])

N = 26

harmonic = 0

for i in range(1,N):
    harmonic += 1.0/i

sys.stderr.write("Harmonic mean: %.4f\n" % harmonic)

print("theta\tdivergence\tratio")

cmd = ""
seg = -1
out = None
ref = None
length = -1
for line in infile:
    sline = line.split()
    if len(sline) == 0 :
        continue
    
    if 'ms' in line:
        cmd = line
    elif sline[0] == "segsites:":
        ###process the previous set
        if seg != -1:
            divergence = sum([1 if a != b else 0 for a,b in zip(out, ref)])
            theta = seg/harmonic
            
            if divergence == 0:
                sys.stderr.write("Zero divergence. %s\n" % cmd)
                continue
            
            print("%s\t%s\t%s" % (theta/length, divergence/length, theta/(divergence)))
        
        
        ###set up variables for new set
        seg = int(sline[1])
        out = None
        ref = None
        length = -1
    elif sline[0] == "positions:":
        continue
    elif sline[0] == "length:":
        length = float(sline[1])
    elif out == None and seg != -1:
        out = line
    elif ref == None and seg != -1:
        ref = line
        
###get the last one
if seg != -1 and out != None and ref != None and length != -1:
    divergence = sum([1 if a != b else 0 for a,b in zip(out, ref)])
    theta = seg/harmonic

    print("%s\t%s\t%s" % (theta/length, divergence/length, theta/divergence))