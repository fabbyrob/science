doc ="""
Reads in a summary file, then computes diversity statistics in windows across the genome. 
Calculates:
    *Pi (per site)
    *Waterson's theta (per site)
    *Tajima's D
    *% of sites diverged
in each window. 
Can use windows based on either based on physical distance (use the -w option) or based on number of snps (use the -s option).

OUTPUT:
chromosome    midpoint    pi    theta    tajd    divergence
scaffold_1    10000    0.01    0.009    -0.03    0.15
scaffold_1    30000    0.02    0.011    -0.05    0.27
"""


import sys
import getopt

_w = 10000
_s = 0
_t = 1000
_N = 26
_c = None
_d = False
_h = False

params = {}

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    calcParams()#calculates all the Tajima's D constants
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    if _d:
        print("chromosome\tmidpoint\tpi\ttheta\ttajd\tdivergence\tCDSdensity")
    else:
        print("chromosome\tmidpoint\tpi\ttheta\ttajd\tdivergence")

    #read the infile file...
    window = []
    window_divergence = []
    window_coding = []
    snps = []
    curr_chrom = ""
    window_start = -1
    window_end = -1
    for line in infile:
        line = line.rstrip()
        sline = line.split()
        
        if line[0] == "#":#comment line, skip
            continue
        
        if len(sline) < 9:
            sys.stderr.write("Incorrectly formatted line:\n\t%s\n" % line)
            continue
        
        chrom = sline[0]
        pos = int(sline[1])
        ref = int(sline[4])#number of reference alleles at site
        tot = int(sline[6])#total number of alleles at site, if this is not _N throw out the site
        if _h:
            ref /= 2
            tot /=2
        type = sline[7].split(",")
        div = int(sline[8])
        
        if not curr_chrom:
            curr_chrom = chrom
        
        if window_start == -1:
            window_start = pos
            window_end = window_start + _w
        
        if chrom != curr_chrom:
            #moved to new chromosome, finish last window and reset the windows
            #print("curr chrom")
            calcWindow(curr_chrom, window_start, window_end, window, window_divergence, window_coding)
            window = []
            window_divergence = []
            window_coding = []
            snps = []
            window_start = pos
            window_end = pos+_w
            curr_chrom = chrom
        
        if _s and len(snps) == _s:
            #print("snp window")
            calcWindow(curr_chrom, window_start, window_end, window, window_divergence, window_coding)
            if _s == _t:
                window = []
                window_divergence = []
                window_coding = []
                snps = []
                window_start = window_end + 1
                window_end = window_start
            else:
                window = window[snps[_t]:]
                window_divergence = window_divergence[snps[_t]:]
                window_coding = window_coding[snps[_t]:]
                window_start += snps[_t]
                window_end = window_start
                last = snps[_t-1]
                sps = snps[_t:]
                snps = [x-last-1 for x in sps]
        elif not _s and pos >= window_end:
            #print("size window")
            calcWindow(curr_chrom, window_start, window_end-1, window, window_divergence, window_coding)
            if _w == _t:
                window = []
                window_divergence = []
                window_coding = []
            else:
                window = window[_t:]
                window_divergence = window_divergence[_t:]
                window_coding = window_coding[_t:]
            window_start += _t
            window_end = window_start + _w 
        
        ####### save info for current site ##########
        if _c != None and _c in type:
            continue
        
        minorfreq = min(ref, _N-ref)
        
        coding = 0
        if "4" in type or "3" in type or "2" in type:
            coding = 1
        
        if tot != _N:
            window.append(-1)
            window_divergence.append(-1)
        else:
            window.append(minorfreq)
            window_divergence.append(div)     
            if minorfreq:
                snps.append(len(window)-1)   
            if _s:
                window_end = pos  
        window_coding.append(coding) 

    #last window, midpoint may be slightly off if it was -w
    calcWindow(curr_chrom, window_start, window_end, window, window_divergence, window_coding)

def calcParams():
    a1 = 0
    a2 = 0
    for i in range(1,_N):
        a1 += 1.0/i
        a2 += 1.0/(i**2)
        
    params["a1"] = a1
    params["a2"] = a2
    
    b1 = float(_N+1)/(3*(_N-1))#see Tajima 1989
    b2 = float(2*(_N**2+_N+3))/(9*_N*(_N-1))

    params["b1"] = b1
    params["b2"] = b2

    c1 = b1 - 1/a1
    c2 = b2-float(_N+2)/(a1*_N)+a2/(a1**2)

    params["c1"] = c1
    params["c2"] = c2

    e1 = c1/a1
    e2 = c2/(a1**2+a2)
    
    params["e1"] = e1
    params["e2"] = e2
    
    sys.stderr.write("Tajima's D parameters, named as in Tajima 1989 Genetics:\n\t%s\n" % params)

def calcWindow(chrom, start, end, window, window_divergence, window_coding):
    midpoint = (start+end)/2.0
    #print("start: %s end: %s" % (start, end))
    #print(window)
    pi, piSites = calcPi(midpoint, window)
    theta, S, nS = calcTheta(midpoint, window)
    tajd = calcTajD(midpoint, pi, theta, S)
    div = calcDivergence(midpoint, window_divergence)
    try:
        out = "%s\t%s\t%s\t%s\t%s\t%s" % (chrom, midpoint, float(pi)/piSites, float(theta)/nS, tajd, div)
        if _d:
            out += "\t"+str(sum(window_coding)/len(window))
        print(out)
    except ZeroDivisionError:
        out = "%s\t%s\t%s\t%s\t%s\t%s" % (chrom, midpoint, "NA", "NA", tajd, div)
        if _d:
            out += "\t"+str(sum(window_coding)/len(window))
        print(out) 
    except ValueError:
        out = "%s\t%s\t%s\t%s\t%s\t%s" % (chrom, midpoint, "NA", "NA", tajd, div)
        if _d:
            out += "\t"+str(sum(window_coding)/len(window)) 
        print(out)    

def calcPi(midpoint, window):
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
        return "NA", 0
    return (total, numSites)

def calcTheta(midpoint, window):
    S = 0#num polymorphic
    nS = 0#num fixed
    for frequency in window:
        if frequency == -1:#no data
            continue
        elif frequency == 0 or frequency == _N:
            nS += 1
        else:
            S += 1
            
    theta = S/params["a1"]
    if S+nS == 0:
        return ("NA", 0, 0)
    return (theta, S, nS)

def calcTajD(midpoint, pi, theta, S):
    if S == 0:#no segregating sites
        sys.stderr.write("No segregating sites for window at %s. Returning NA for Taj D.\n" % midpoint)
        return "NA"
    if pi == "NA" or theta == "NA":#no pi or theta
        sys.stderr.write("Bad pi (%s) or theta (%s) values. Returning NA for Taj D.\n" % (pi, theta))
        return "NA"
    
    
    d = pi - theta
    d /= (params["e1"]*S+params["e2"]*S*(S-1))**0.5
    return d

def calcDivergence(midpoint, window):
    if len(window) > 0:
        div = 0.0
        S = 0.0
        for site in window:
            if site == -1:#missing data
                continue
            if site:
                div += 1
            S += 1
        if S == 0:
            sys.stderr.write("No divergence info for window %s. Returning NA.\n" % midpoint)
            return "NA"
        return div/S
    else:
        sys.stderr.write("No divergence info for window %s. Returning NA.\n" % midpoint)
        return "NA"

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:s:t:N:c:dh")
    except getopt.GetoptError:
        usage()
    
    global _w 
    global _s
    global _t
    global _N
    global _c
    global _d
    global _h
    
    for opt, arg in opts:
        if opt == "-w":
            _w = int(arg)
        elif opt == "-s":
            _s = int(arg)
            _w = 0
        elif opt == "-t":
            _t = int(arg)
        elif opt == "-N":
            _N = int(arg)
        elif opt == "-c":
            _c = arg
        elif opt == "-d":
            _d = True
        elif opt == "-h":
            _h = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -w %s -s %s -t %s -N %s -c %s -d %s\n" % (sys.argv[1], _w, _s, _t, _N, _c, _d))
   
use = "python "+__file__.split("/")[-1]+" SummaryFile [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print (doc)
    print("______________________________________")
    print("OPTION - TYPE - DESCRIPTION (default)")
    print("w - INT - window size to calculate stats in, number of bases (%s). This is the default window mode." % _w)
    print("s - INT - number of snps to define a window (%s). Overrides any input for -w." % _s)
    print("t - INT - step size for window. For fixed windows make this the same as -w or -s. (%s)" % _t)
    print("N - INT - sample size, number of alleles. (%s)" % _N)
    print("c - STR - if given a site type with this option, then diversity statistics are only calculated for these sites. (%s)" % _c)
    print("d - NONE - if this flag is used then the coding density in each window is also output (%s)" % _d)

    
if __name__ == "__main__":   
    __main__()
