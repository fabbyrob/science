import sys
import re
import getopt

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
    
    site_pat = re.compile("\D+(\d+)\s+(\S+)\s+(\S+)\s+(\S+).+")
    
    scaf = -1
    start = end = 0
    totalSites = 0
    totalDif = 0
    seqs = []
    last = (0,0,0,0,0)
    first = 0
    #bases = {'A':0,'T':0,'C':0,'G':0}
    acc = (0,0)
    
    print("Total_Site_Num\tSites_Diverged\tPercent_Diverged\tIntron_Length")
    
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
                (last, acc) = process((scaf, start, end), seqs, last, acc, sitefile, site_pat, logfile)
                if last == None:#if we have no more sites to look for
                    break
                seqs = []
                first = 0
                continue
            
            if line[0] == "s":
                items = line.split()
                #print (items[-1])
                if first == 0:#makes sure we get data only from first line (rubella)
                    scaf = int(items[1].split("_")[1])
                    #scaf = int (items[1][-1])#for use with AT references
                    start = int(items[2])+1
                    end = start+int(items[3])-1
                    first = False
                seqs.append((items[1],items[-1]))
                first += 1
                continue
            
        (last, acc) = process((scaf, start, end), seqs, last, acc, sitefile, site_pat, logfile)

 
def printData(last, acc):
    tot = acc[0]
    dif = acc[1]
    string = str(tot)+"\t"+str(dif)+"\t"
    if tot > 0:
        string += str(float(dif)/tot)
    else:
        string += "."
    string += "\t"+str(last[4])
    print(string)#print out the divergence and length of this intron
    
def process(range, seqs, last, acc, sites, site_pat, logfile):
    #bases = {'A':0,'T':0,'C':0,'G':0}
    if len(seqs) == 2:
        #(tot, dif, last) = roughCompare(seqs[1][1], seqs[2][1])
        (last, acc) = compare(range, seqs[0][1], seqs[1][1], last, acc, sites, site_pat, logfile)
    return (last, acc)
   
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
     
def compare(myRange, seq1, seq2, last, acc, sites, site_pat, logfile):
    #print (seq1, seq2)
    tot = dif = 0
    i = 0
    scaf = myRange[0]
    start = myRange[1]
    end = myRange[2]
    #bases = {'A':0,'T':0,'C':0,'G':0}
    
    indices = correctIndex(start, seq1)
    
    #get new site
    if last == (0,0,0,0,0):
        logfile.write("site not found: "+str(last)+"\n")
        if acc != (0,0):
            printData(last, acc)
        last = getNewRange(sites, site_pat, last, logfile)
        acc = (0,0)
    
    while last != None and last[1] < scaf:#move us up to the right scaffold
        logfile.write("site not found: "+str(last)+"\n")
        if acc != (0,0):
            printData(last, acc)
        last = getNewRange(sites, site_pat, last, logfile)
        acc = (0,0)
        
    if last != None and last[1] > scaf:#nothing on this scaffold
        return (last, acc)    
    
    while last != None and last[2] < start:#move us along the chromosome to the correct location
        logfile.write("site not found: "+str(last)+"\n")
        if acc != (0,0):
            printData(last, acc)
        last = getNewRange(sites, site_pat, last, logfile)
        acc = (0,0)
    
    #I'm sorry future me, you probably just want to completely re-write this script. It's horrible.
    while last != None and last[1] == scaf and ((last[2] >= start and last[2] <= end) or (last[3] >= start and last[3] <= end)):
        tot = dif = 0
        if last[3] <= end and last[2] >= start:#then we have the whole thing
            whole = True
            ran = range(last[2], last[3]+1)
        else:#we need to add to the accumulator 
            if last[2] >= start: #begining okay
                ran = range(last[2], end+1)
            else:#end okay
                ran = range(start, last[3]+1)
            whole = False
        
        for j in ran:
            i = indices[j]
            
            if seq2[i] == "-" or seq1[i] == "-":
                logfile.write("no neslia data for: "+str(last)+"\n")
            elif seq1[i].upper() != seq2[i].upper():
                tot += 1
                dif += 1
                logfile.write("DIVERGED: "+str(last)+"\n")
            else:
                #bases[seq1[i].upper()] += 1
                tot += 1
                logfile.write("NOT DIVERGED: "+str(last)+"\n")
            
        if whole:
            printData(last, (acc[0]+tot, acc[1]+dif))
        else:
            last = (last[0], last[1], end+1, last[3], last[4])
            return (last, (tot,dif))#return the accumulator  
            
        last = getNewRange(sites, site_pat, last, logfile)
        acc = (0,0)
            
    return (last, acc)
  
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
def getNewRange(sites, site_pat, last, logfile):
    new_site = None
    while(new_site != ""):#will stop if it reaches the end of the file
        new_site = sites.readline()
        
        if (new_site == ""):
            return None
        
        m = site_pat.match(new_site)
        
        if (m == None):
            logfile.write("Non-standard line (SITE FILE):\n\t"+new_site+"\n")
            continue
        
        scaf = int(m.group(1))
        start = int(m.group(2))
        end = int(m.group(3))
        mod = m.group(4)
        
        return((mod, scaf, start, end, start-end-1))
        
    return None
   
use = "python "+__file__.split("/")[-1]+" SITEFILE ALIGNMENT [ALIGNMENT2 ... ALIGNMENTn]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()
