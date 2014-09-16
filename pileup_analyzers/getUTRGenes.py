import sys

infile = open(sys.argv[1],"r")

def getStart(ls):
    start = None
    for l in ls:
        if start == None or l[0] < start:
            start = l[0]
    return start

def getEnd(ls):
    end = None
    for l in ls:
        if end == None or l[1] > end:
            end = l[1]
    return end

#infile = open("/Users/wiliarj/Desktop/temp/annotation.gff3", "r")

myGene = None
myScaf = None
myDir = None
five = []
three = []
cds = []

for line in infile:
    line = line.rstrip()
    sline = line.split()
    
    if len(sline) != 9:
        continue
    
    scaf = sline[0]
    type = sline[2]
    dir = sline[6]
    start = int(sline[3])
    end = int(sline[4])
    
    info = sline[8]
    sinfo = info.split(";")
    
    name = sinfo[0].split("ID=")[1]
    name = name.split(".")[0]
    
    if "Carubv" in name:#skip Carubv lines
        continue
    
#    if myGene != name:
#        #process old stuff
#        if myGene:     
#            bestStart = getStart(three+five+cds)
#            bestEnd = getEnd(three+five+cds)
#            
#            print("%s %s %s %s" % (myGene, myScaf, bestStart, bestEnd))
#            
#        #reset
#        myGene = name
#        myScaf = scaf
#        myDir = dir
#        five = []
#        three = []
#        cds = []
    
#    if type == "five_prime_UTR":
#        five.append((start, end))
#    elif type == "three_prime_UTR":
#        three.append((start, end))
#    elif type == "CDS":
#        cds.append((start, end))
    if type == "mRNA":
        print("%s %s %s %s" % (name, scaf, start, end))
            
#if myGene:
#    bestStart = getStart(starts)
#    bestEnd = getEnd(ends)
    
#    print("%s %s %s %s" % (myGene, myScaf, bestStart, bestEnd))