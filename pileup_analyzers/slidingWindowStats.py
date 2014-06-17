import sys
import re
import getopt
import vcf
from math import sqrt

_w = 10000
_s = 1000
_W = 5
_q = 40
_d = 20
_D = 80
_l = 40 
_S = True#flag to do snp window, false == do bp window; currently not in use
_N = False

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[2:],"l:w:s:q:d:D:W:N:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-s":
            global _s 
            _s = int(arg)
        elif opt == "-W":
            global _W 
            _W = int(arg)
        elif opt == "-q":
            global _q 
            _q = int(arg)
        elif opt == "-d":
            global _d 
            _d = int(arg)
        elif opt == "-D":
            global _D 
            _D = int(arg)
        elif opt == "-l":
            global _l 
            _l = int(arg)
        elif opt == "-S":
            global _S 
            _S = True
        elif opt == "-N":
            global _N
            _N = True
            numInds = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    
    vcf_reader = vcf.Reader(open(sys.argv[1], "rb"))
     
    #grab all of the individual IDs from the second line in the vcf
    if not _N:
        numInds = len(vcf_reader.samples)
        indIDs = vcf_reader.samples
    
    sys.stderr.write("Detected "+str(numInds)+" individuals\n")
    
    #initialize our sequence dictionaries
    seq = {}#per individual sequences
    for i in range(0, numInds):
        seq[i] = []

    #TAJ D constants  
    N = numInds*2  
    a1 = 0
    a2 = 0
    for i in range(1,N):
        a1 += 1/float(i)
        a2 += 1/(float(i)*float(i))
        
    b1 = float(N+1)/(3*(N-1))#see Tajima 1989
    b2 = float(2*(N*N+N+3))/(9*N*(N-1))

    c1 = b1 - 1/a1
    c2 = b2-float(N+2)/(a1*N)+a2/(a1*a1)

    e1 = c1/a1
    e2 = c2/(a1*a1+a2)
        
    seq = []    
    snps = []
    skip = 0
    indelCheck = 0
    indelCount = 0
    
    print("MidPoint,pi,theta,TajD")
    
    for record in vcf_reader:
        if len(snps) == _w:#have a window
            start = snps[0][1]
            end = snps[_w-1][1]
            subFreqs = seq[snps[0][0]:(snps[_w-1][0]+1)]

            mid = (start+end)/2.0
            
            calcWindow(subFreqs, mid, numInds, a1, e1, e2)
            
            snps = snps[_s+1:len(snps)]
            
            if _w < _s:
                skip = _w-_s
       
        indelCount += 1
        indelFlag = 0
        
        scaf = record.CHROM
        base = record.POS
        ref = record.REF
        alt = record.ALT[0]
        qual = record.QUAL
        
        if 'Dels' in record.INFO.keys():
            indelAmount = record.INFO['Dels']
        
            if indelAmount > 0:
                sys.stderr.write("INDEL found at\n\t %s %s\n" % (scaf, base))
                indelFlag = 1
                indelCount = -1*_W
        
        #make window ambig
        if indelFlag:
            makeWindowAmbig(seq, -1*indelCount, False)
        
        if indelCount <= 0 :
            #allSet(seq)
            seq.append(-1)
        else:
            af = callSite(scaf, base, ref, alt, qual, record)
            seq.append(af)
            if (af > 0 and af < numInds*2):
                if skip >= 0:
                    snps.append((len(seq)-1, base))
                skip += 1

        
    #get last window   
    if len(snps) > 0: 
        mid = (snps[0][1]+snps[-1][1])/2.0
        subFreq = seq[snps[0][1]:snps[-1][1]]
        calcWindow(subFreq, mid, numInds, a1, e1, e2)
              
def calcWindow(afs, mid, numInds, a1, e1, e2):
    (pi, sitesUsedPi) = piWindow(afs, numInds*2)
    (tajd, theta, sitesUsedTheta) = tajDWindow(afs, pi, numInds*2, a1, e1, e2)
    if sitesUsedPi == 0 or sitesUsedTheta == 0:
        sys.stderr.write("No viable sites in window: "+str(mid)+"\n")
    else:
        print(str(mid)+","+str(float(pi)/sitesUsedPi)+","+str(float(theta)/sitesUsedTheta)+","+str(tajd))#+","+str(filter(lambda a: a != -1, afs)))
    
def callSite(scaf, base, ref, alt, qual, record):
    freq = -1
    if(alt not in ['A','T','C','G','.']):
        #not a biallelic site
        #ambig
        sys.stderr.write("Non-standard base for alternate: "+alt+" at "+str(scaf)+", "+str(base)+". Calling ambiguous\n")
        return (freq)
    
    if qual < _q:
        sys.stderr.write("Site below quality threshold.\n\t"+str(base)+"\n")
        return (freq)
    
    freq = 0
    freqFail = False

    c = 0
    for sample in record.samples:
        if 'DP' in sample.keys():
            myDepth = sample['DP']
        else:
            freqFail = True
            continue
            
        if myDepth < _d or myDepth > _D:
            #did not pass depth filters
            c += 1
            freqFail = True
            continue
        
        if record.ALT[0] == '.' or 'PL' not in sample.keys():
            freqFail = True
            continue
        elif record.ALT[0] == '.':
            freq = 0
        else:
            data = map(int, sample['PL'])
            
            (min, max, mid) = getMinMaxIndex(data)
            
            if (data[min] != 0):
                #site is ambiguous, no max likelihood
                c += 1
                freqFail = True
                continue
            
            if (data[mid] < _l):
                #dont pass quality cut off
                freqFail = True
            else:
                freq += min
    #            if min == 1:#heterozygote
    #                freq += 1
    #            elif min == 0:#homozygous ref
    #                continue
    #            elif min == 2:#homozygous alt
    #                freq += 2
        c += 1

    if freqFail:
        freq = -1

    return (freq)

     
def makeWindowAmbig(seq, size, p = True):
#    for k in seq.keys():
#        if (len(seq[k]) > 0 and p):
#            seq[k].pop()#removes last one because of double INDEL lines
#        if len(seq[k]) < size:
#            seq[k] = ["N"]*len(seq[k])
#        else:
#            for i in range(-1*size, 0):
#                seq[k][i] = "N"
    if (len(seq) > 0 and p):
        seq.pop()#removes last one because of double INDEL lines
    if len(seq) < size:
        seq = [-1]*len(seq)
    else:
        for i in range(-1*size, 0):
            seq[i] = -1 

def makeInts(ls):
    for i in range(0, len(ls)):
        ls[i] = int(ls[i])

def getMinMaxIndex(ls):
    min = 0
    minV = ls[0]
    max = 0
    maxV = ls[0]
    for i in range(0, len(ls)):
        if ls[i] < minV:
            min = i
            minV = ls[i]
        if ls[i] > maxV:
            max = i
            maxV = ls[i]
            
    return (min, max, list(set([1,2,0])-set([min, max]))[0])

def countMinor(seqs):
    counts = []
    snps = []
    for i in range(0, len(seqs[0])):#for each position
        ref = ''
        refcount = 0
        alt = ''
        altcount = 0
        for myInd in seqs.keys():
            ind = seqs[myInd]
            if ind[i] == 'N':
                refcount = -1
                break
            if ref == '' and ind[i] in ['A','T','C','G']:
                ref = ind[i]
                refcount += 2
            elif ref == ind[i]:
                refcount += 2
            elif alt == ''and ind[i] in ['A','T','C','G']:
                alt = ind[i]
                altcount += 2
            elif alt == ind[i]:
                altcount += 2
            else:#heterozygote
                refcount += 1
                altcount += 1
        
        if refcount < altcount:
            counts.append(refcount)
        else:
            counts.append(altcount)
        if counts[-1] > 0:
            snps.append(i)
    return (counts, snps)

def piWindow(afs, sampSize):
    total = 0
    numSites = 0
    for f in afs:
        if f == -1:#skip
            continue
        minorfreq = (float(f)/sampSize)
        #print("minorFreq ",minorfreq)
        majorfreq = 1-minorfreq
        #print("majorFreq ",majorfreq)
        #print("1-(majorFreq^2+minorFreq^2) ",1-(minorfreq*minorfreq+majorfreq*majorfreq))
        
        total += (1-(minorfreq*minorfreq+majorfreq*majorfreq))
        #print("total ", total)
        numSites += 1


    total = (float(sampSize)/(sampSize-1))*total
    #print("total*n/n-1 ", total)

    return (total, numSites)
                
def tajDWindow(afs, pi, sampSize, a1, e1 ,e2):

    #print(a1)
    #print(afs)
    S = 0#number of polymorphic sites
    nS = 0 #number non-polymorphic
    for s in afs:
        if s == -1:
            continue
        elif s == 0 or s == sampSize:
            nS += 1
        else:
            S += 1
         
    #print("S ",S," nS ",nS)
    if (S+nS) == 0:#no viable sites
        return (".",".",0)
            
    theta = float(S)/a1
    #print("theta ", theta)
    thetaPerSite = theta/(S+nS)

    d = pi-theta
    if(pi == 0 or theta == 0):
        D = '.'
    else:
        p1 = e1*S
        p2 = e2*S*(S-1)
        #sys.stderr.write(mystr)
        D = d/sqrt(p1+p2)
    
    return (D, theta, S+nS)
            
def mean(numberList):
    if len(numberList) == 0:
        return float('nan')
 
    floatNums = [float(x) for x in numberList]
    return sum(floatNums) / len(numberList)
        
def processWindow(win):
    data = []
    for i in win:
        data.append(mean(i))
    return data
   
use = "python "+__file__.split("/")[-1]+" <VCF> [options]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("w - the width of the sliding window to analyze (10000)")
    print ("s - the step size the window will increment by (1000)")
    print ("q - the minimum quality to include a site in the analysis (40)")
    print ("d - the minimum depth to include a site in the analysis (20)")
    print ("D - the maximum depth to include a site in the analysis (60)")
    print ("l - the minimum likelihood to call a site as homozygous, or heterozygous rather than N (40)")
    print ("W - the number of sites on either side of INDELs to call ambiguous (5)")
    print ("N - sample size (If the VCF has a header this will be automatically detected)")
    
    
if __name__ == "__main__":   
    __main__()
#    numInds = 10
#    a1 = 0
#    a2 = 0
#    for i in range(1,numInds):
#        a1 += 1/float(i)
#        a2 += 1/(float(i)*float(i))
#        
#    b1 = float(numInds+1)/(3*(numInds-1))#see Tajima 1989
#    b2 = float(2*(numInds*numInds+numInds+3))/(9*numInds*(numInds-1))
#
#    c1 = b1 - 1/a1
#    c2 = b2-float(numInds+2)/(a1*numInds)+a2/(a1*a1)
#
#    e1 = c1/a1
#    e2 = c2/(a1*a1+a2)
#    
#    fakeData = [0,0,0,0,0]
#    myPi,numSites = piWindow(fakeData, 10)
#    print(myPi)
#    myD, myT, numSites2 = tajDWindow(fakeData, myPi, 10, a1, e1, e2)
#    print(myT)
#    print(myD)
#    print()
#    fakeData = [10,10,10,10,10]
#    myPi,numSites = piWindow(fakeData, 10)
#    print(myPi)
#    myD, myT, numSites2 = tajDWindow(fakeData, myPi, 10, a1, e1, e2)
#    print(myT)
#    print(myD)
#    print()
#    fakeData = [5,5,5,5,5,-1]
#    myPi,numSites = piWindow(fakeData, 10)
#    print(myPi)
#    myD, myT, numSites2 = tajDWindow(fakeData, myPi, 10, a1, e1, e2)
#    print(myT)
#    print(myD)
#    print()
#    print()
    
    