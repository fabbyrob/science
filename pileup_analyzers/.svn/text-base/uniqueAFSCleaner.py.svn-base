import sys
from math import ceil

myFile = open(sys.argv[1],"r")

afs = []
for line in myFile:
    line = line.split()
    afs.append(int(line[2]))
    
folded = []

numFolded = int(ceil(len(afs)/2.0))
for i in range(0,numFolded):
    if i+1 == numFolded and numFolded % 2 == 1:
        folded.append(afs[i])
    else:
        folded.append(afs[i]+afs[-1*(i+1)])

myStr = ""
for n in folded:
    myStr += str(n)+"\t"
    
print(myStr)