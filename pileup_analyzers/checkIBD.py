doc = """
Takes in a summary file, with individual genotype calls, and outputs the number of homozygous and heterozygous calls for each individual in windows across the genome.
The Fis for each window is also output.
The program will optionally output all regions in which the Fis is above a given threshold that defines a region of IBD.
"""

import sys
import getopt

_w = 10000
_o = ""#output file name for regions of high Fis
_n = 1#minimum number of continuous windows needed to define a region of IBD
_Fis = 0.9#cutoff for region of IBD
_N = 26#sample size

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    if _o:
        regionsOutput = open(_o,"w")  
        regionsOutput.write("CHROM\tSTART\tEND\tSAMPLE\n")
        if (regionsOutput == None):
            print("Bad regions file name: "+_o)
            sys.exit()
            
    global _N
    regions = {}
    window = {}
    start = -1
    prevPos = 0
    currScaf = ""
    afs = [0]*_N
    samples = ""
    
    print("CHROM\tMIDPOINT\tSAMPLE\tHOMOZYGOUS\tHETEROZYGOUS\tALL_SITES\tFis")
    #read the infile file...
    for line in infile:
        line = line.rstrip()
        sline = line.split()
        
        if len(line) > 0 and line[0] == "#":
            #header line
            if sline[0] == "#CHROM":
                samples = sline[9:]
                for ind in samples:
                    regions[ind] = []
                    window[ind] = []
                
                if len(samples)*2 != _N:
                    sys.stderr.write("Number of samples in file (%s) does not match input parameter -N (%s). Using the number from the file.\n" % (len(samples), _N))
                    _N = len(samples)
                    afs = [0]*_N
            continue
        
        scaf = sline[0]
        alt = int(sline[5])
        pos = int(sline[1])
        genos = sline[9:]
        
        if not currScaf:
            currScaf = sline[0]
            
        if start == -1:
            start = pos
        
        if pos >= start+_w or currScaf != scaf:
            regions = calcWindow(window, currScaf, afs, regions, start, prevPos, regionsOutput)
            
            #reset variables for next window
            if currScaf != scaf:
                outputRegions(regions, currScaf, regionsOutput)
                
            for ind in samples:
                if currScaf != scaf:
                    regions[ind] = []
                window[ind] = []
                
            afs = [0]*_N
            currScaf = scaf
            start = pos
            
        if len(genos) != len(samples):
            sys.stderr.write("Line with not enough genotypes (%s, %s), skipping:\n\t%s\n" % (len(genos), len(samples), line))
            continue
        
        afs.append(min(_N-alt, alt))
        #convert fixed sites into F's so later we can count homozygous SNPs separately
        if alt == 0:
            genos = ["F" if x == "R" else x for x in genos]
        elif alt == _N:
            genos = ["F" if x == "A" else x for x in genos]
            
        for ind, geno in zip(samples, genos):
            window[ind].append(geno)
            
        prevPos = pos

    #last window
    calcWindow(window, currScaf, afs, regions, start, prevPos, regionsOutput)
    outputRegions(regions, currScaf, regionsOutput)
    
    infile.close()
    regionsOutput.close()
  
#output: SCAF MIDPOINT IND HOMO REF TOTAL FIS       
def calcWindow(window, scaf, afs, regions, start, stop, output):
    betweenPi = calcPi(afs)
    indivPis = calcPiIndiv(window)
    Fis = calcFis(betweenPi, indivPis)
    regions = checkIBD(regions, Fis, scaf, (start, stop), output)
    midpoint = (start+stop)/2.0
    
    for ind in window.keys():
        homo = window[ind].count("R")+window[ind].count("A")
        homoFixed = window[ind].count("F")
        het = window[ind].count("H")
        total = homo+het+homoFixed
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (scaf, midpoint, ind, homo, het, total, Fis[ind]))
        
    return regions
     
#output: SCAF START END IND
def checkIBD(regions, Fis, scaf, limits, output):
    for ind in regions.keys():
        if Fis[ind] != "NA" and Fis[ind] > _Fis:#count it!
            regions[ind].append(limits)
        else:
            if len(regions[ind]) >= _n:#have enough regions of IBD, output it
                start = regions[ind][0][0]
                end = regions[ind][-1][1]
                if _o:
                    output.write("%s\t%s\t%s\t%s\n" % (scaf, start, end, ind))
            #clear the region
            regions[ind] = []
    return regions
 
def outputRegions(regions, scaf, output):
    for ind in regions.keys():
        if regions[ind]:
            start = regions[ind][0][0]
            end = regions[ind][-1][1]
            if _o:
                output.write("%s\t%s\t%s\t%s\n" % (scaf, start, end, ind))
        
#N = unknown, H = ref homo, A = alt homo, H = het
def calcFis(betweenPi, indivPis):
    myFis = {}
    for ind in indivPis.keys():
        if indivPis[ind] != "NA":
            if betweenPi > 0:
                myFis[ind] = float(betweenPi-indivPis[ind])/betweenPi
            else:
                myFis[ind] = 0
        else:
            myFis[ind] = "NA"
        
    return myFis

def calcPiIndiv(window):
    myPis = {}
    for ind in window.keys():
        if (len(window[ind])-window[ind].count("N")) > 0:
            myPis[ind] = float(window[ind].count("H"))/(len(window[ind])-window[ind].count("N"))
        else:
            myPis[ind] ="NA"
            
    return myPis

def calcPi(window):
    total = 0
    numSites = 0.0
    for frequency in window:
        if frequency == -1:#no data
            continue
        
        minorfreq = float(frequency)/_N
        majorfreq = 1-minorfreq
        
        total += 2*minorfreq*majorfreq
        numSites += 1.0
        
    total *= _N/(_N-1.0)
    if numSites == 0:
        return "NA"
    return total/numSites

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:o:n:N:f:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-o":
            global _o 
            _o = arg
        elif opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-n":
            global _n 
            _n = int(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-f":
            global _Fis 
            _Fis = float(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -o %s -w %s -n %s -N %s -f %s\n" % (sys.argv[1], _o, _w, _n, _N, _Fis))
   
use = "python "+__file__.split("/")[-1]+" VCFSUMAMRY [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("w - INT - %s - window size in which to calculate statistics for each individual." % _w)
    print("N - INT - %s - sample size, number of chromosomes in your full sample." % _N)
    print("o - STR - %s - file that regions of IBD should be output to. If this option is not provided no regions will be output." % _o)
    print("n - INT - %s - number of continuous windows of high Fis needed to define a region of IBD." % _n)
    print("f - FLOAT - %s - any windows with Fis higher than this cutoff are defined as IBD." % _Fis)

    
if __name__ == "__main__":   
    __main__()