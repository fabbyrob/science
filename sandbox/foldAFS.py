import sys

afs = map(float, sys.argv[1:])

N = len(afs)
newAFS = []

if N % 2 == 0:
    even = True
else:
    even = False

for i, freq in enumerate(afs):
    if not even and i == int(N/2):
        newAFS.append(freq)
        break
    
    if even and i == int(N/2):
        break
    
    newAFS.append(freq+afs[N-i-1])
    
    
print("\t".join(map(str, newAFS)))