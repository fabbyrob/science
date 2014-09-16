import sys
import re
import getopt

_log = __file__.split("/")[-1]+"."+sys.argv[1].split("/")[-1]+".log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
#    try: 
#        opts, args = getopt.getopt(sys.argv[2:],"l:")
#    except getopt.GetoptError:
#        usage()
#    
#    for opt, arg in opts:
#        if opt == "-l":
#            global _log 
#            _log = str(arg)
#        else:
#            print ("Unrecognized option: "+opt+"\n")
#            usage()
    logfile = open(_log,"w")
        
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
    
    sitefile = open(sys.argv[1], "r")
    
    if (sitefile == None):
        print("Bad sitefile name: "+sys.argv[1])
        sys.exit()
    
    site_pat = re.compile("\D+(\d+)\s+(\S+)")
    
    scaf = -1
    start = end = 0
    totalSites = 0
    totalDif = 0
    seqs = []
    last = (0,0)
    first = 0
    bases = {'A':0,'T':0,'C':0,'G':0}
    
    for file in sys.argv[2:]:
        infile = open(file,"r")
        
        if (infile == None):
            print("Bad infile name: "+file)
            sys.exit()
            
        
            
        for line in infile:
            line = line.rstrip()
            if len(line) == 0 or line[0] == "#":
                continue
            
            if line[0] == "a" and scaf != -1:
                #if all 4 samples are listed
                (t, d, last, baseComp) = process((scaf, start, end), seqs, last, sitefile, site_pat, logfile)
                if last == None:#if we have no more sites to look for
                    break
                totalSites += t
                totalDif += d
                seqs = []
                first = 0
                for b in bases.keys():
                    bases[b] += baseComp[b]
                continue
            
            if line[0] == "s":
                items = line.split()
                #print (items[-1])
                if first == 2:#makes sure we get data only from first line (rubella)
                    #scaf = int(items[1].split("_")[1])
                    scaf = int (items[1][-1])#for use with AT references
                    start = int(items[2])+1
                    end = start+int(items[3])-1
                    first = False
                seqs.append((items[1],items[-1]))
                first += 1
                continue
            
        (t, d, last, baseComp) = process((scaf, start, end), seqs, last, sitefile, site_pat, logfile)
        for b in bases.keys():
            bases[b] += baseComp[b]
        totalSites += t
        totalDif += d
        seqs = []
      
    print ("Total\tDivergence")  
    print (str(totalSites)+"\t"+str(totalDif))
    if totalSites > 0:
        print(float(totalDif)/float(totalSites))
    print ("Base composition")

    t = 0
    for k in bases.keys():
        t += bases[k]
    
    if t > 0:
        for k in bases.keys():
            print(k+": "+str(bases[k])+"\t"+str(float(bases[k])/t))
 
def process(range, seqs, last, sites, site_pat, logfile):
    tot = dif = 0
    bases = {'A':0,'T':0,'C':0,'G':0}
    if len(seqs) == 3:
        #(tot, dif, last) = roughCompare(seqs[1][1], seqs[2][1])
        (tot, dif, last, bases) = compare(range, seqs[2][1], seqs[1][1], last, sites, site_pat, logfile)
    return (tot, dif, last, bases)
   
def roughCompare(seq1, seq2):
    tot = dif = 0
    i = 0
    
    for i in range(0, len(seq1)):
        if seq1[i] == "-" or seq2[i] == "-":
            continue
        
        if seq1[i] != seq2[i]:
            dif += 1
            
        tot += 1
        
    return (tot, dif, (0,0))
     
def compare(range, seq1, seq2, last, sites, site_pat, logfile):
    #print (seq1, seq2)
    tot = dif = 0
    i = 0
    scaf = range[0]
    start = range[1]
    end = range[2]
    bases = {'A':0,'T':0,'C':0,'G':0}
    
    indices = correctIndex(start, seq1)
    
    #get new site
    if last == (0,0):
        logfile.write("site not found: "+str(last)+"\n")
        last = getNewSite(sites, site_pat, logfile)
    
    while last != None and last[0] < scaf:#move us up to the right scaffold
        logfile.write("site not found: "+str(last)+"\n")
        last = getNewSite(sites, site_pat, logfile)
        
    if last != None and last[0] > scaf:#nothing on this scaffold
        return (0,0, last, bases)    
    
    while last != None and last[1] < start:#move us along the chromosome to the correct location
        logfile.write("site not found: "+str(last)+"\n")
        last = getNewSite(sites, site_pat, logfile)
    
    while last != None and last[0] == scaf and last[1] >= start and last[1] <= end:
        i = indices[last[1]]
        
        if seq2[i] == "-" or seq1[i] == "-":
            logfile.write("no neslia data for: "+str(last)+"\n")
        elif seq1[i].upper() != seq2[i].upper():
            tot += 1
            dif += 1
            logfile.write("DIVERGED: "+str(last)+"\n")
        else:
            bases[seq1[i].upper()] += 1
            tot += 1
            logfile.write("NOT DIVERGED: "+str(last)+"\n")
            
        last = getNewSite(sites, site_pat, logfile)
            
    return (tot, dif, last, bases)
  
def correctIndex(start, seq):
    x = start
    dict = {}
    for i in range(0, len(seq)):
        if seq[i] == "-":
            continue
        dict[x] = i
        x += 1
    
    return dict
  
#grabs a new site too look for from the sites input file
def getNewSite(sites, site_pat, logfile):
    next_site = None
    while(next_site == None):
        new_site = sites.readline()
        
        if (new_site == ""):
            return None
        
        m = site_pat.match(new_site)
        
        if (m == None):
            logfile.write("Non-standard line (SITE FILE):\n\t"+new_site+"\n")
            continue
        next_site = (int(m.group(1)), int(m.group(2)))
        
    return next_site
   
use = "python "+__file__.split("/")[-1]+" SITEFILE ALIGNMENT [ALIGNMENT2 ... ALIGNMENTn]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()
