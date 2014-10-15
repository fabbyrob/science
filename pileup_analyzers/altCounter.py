import sys
import getopt
import vcf

_w = 100
_i = 0
_d = 20
_D = 60
_q = 90
_L = 40
_N = 26
_Q = True

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    vcf_reader = vcf.Reader(open(sys.argv[1],"rb"))
    names = vcf_reader.samples
    _N = len(names)
    window = resetWindow({}, names)
    dwindow = resetWindow({}, names)
    start = 0
    end = start + _w
    currChrom = ""
    prevPos = -1
    
    
    print("ind\tpos\tmidpoint\tnum_homo\tnum_het\tnum_homo_ref\ttotal_depth")
    for record in vcf_reader:
        if record.CHROM != currChrom:
            if currChrom != "":
                processWindow(currChrom, start, prevPos, window, dwindow)
                start = record.POS
                end = start + _w
                window = resetWindow(window, names)
                dwindow = resetWindow(window, names)
            currChrom = record.CHROM
        elif record.POS > end:
            processWindow(currChrom, start, end, window, dwindow)
            start = end + 1
            end = start + _w
            window = resetWindow(window, names)
            dwindow = resetWindow(window, names)
                
        window, dwindow = processSite(record, window, dwindow)
        prevPos = record.POS
            
    #get the last window
    processWindow(currChrom, start, prevPos, window, dwindow)
        

def processWindow(chrom, start, end, window, dwindow):
    midpoint = (start+end)/2
    for n in window.keys():
        homos = window[n].count(2)
        hets = window[n].count(1)
        rhomos = window[n].count(0)
        tdepth = sum(dwindow[n])
        print(n+"\t"+chrom+"\t"+str(midpoint)+"\t"+str(homos)+"\t"+str(hets)+"\t"+str(rhomos)+"\t"+str(tdepth))
 
def resetWindow(window, names):
    window = {}
    for n in names:
        window[n] = []
    return window
        
def processSite(record, window, dwindow):
    if 'Dels' in record.INFO.keys() and float(record.INFO['Dels']) > _i:
        sys.stderr.write("Site is an indel.\n\t"+str(record.POS)+"\n")
        return window, dwindow
    elif record.QUAL < _q:
        sys.stderr.write("Site below quality threshold.\n\t"+str(record.POS)+"\n")
        return window, dwindow
    else:
        return(callSite(record, window, dwindow))
 
#return the minor allele frequency at a given site, or -1 if the site fails filters for any individual
def callSite(record, window, dwindow):
    freq = -1
    if(record.ALT[0] not in ['A','T','C','G','.']):
        #not a biallelic site
        #ambig
        sys.stderr.write("Non-standard base for alternate: "+str(record.ALT)+" at "+str(record.CHROM)+", "+str(record.POS)+". Calling ambiguous\n")
        return (window, dwindow)
    
    freq = 0
    freqFail = False

    for sample in record.samples:
        if freqFail:
            break
        
        if (_Q and 'PL' not in sample.keys()) or 'DP' not in sample.keys():#one individual is missing data, throw out the site
            sys.stderr.write("One sample: "+ str(sample['name']) +"lacks depth or likelihoods at\n\t"+str(record)+"\n")
            freq = -1
            freqFail = True
            break
         
        myDepth = sample['DP'][0]
            
        if myDepth < _d or myDepth > _D:
            #did not pass depth filters
            sys.stderr.write("One sample: "+ str(sample['name']) +" does not pass depth filters\n\t"+str(record)+"\n")
            freqFail = True
            break
        
        if (_Q or sample['GT'] != "0/0") and (record.ALT != "."):
            myPLs = sample['PL']
            if type(myPLs) == list:
                myPLs = map(int, myPLs)
            else:
                myPLs = map(int, myPLs.split(","))
            (minQ, midQ, maxQ) = getMinMaxIndex(myPLs)
            
            if (myPLs[minQ] != 0):
                #site is ambiguous, no max likelihood
                sys.stderr.write("One sample: "+ str(sample['name']) +" does not pass likelihood filters\n\t"+str(record)+"\n")
                freqFail = True
                break
            
            if (myPLs[midQ] < _L):
                #dont pass quality cut off
                sys.stderr.write("One sample: "+ str(sample['name']) +" does not pass likelihood filters\n\t"+str(record)+"\n")
                freqFail = True
                break
            else:
                window[sample['name']].append(minQ)
                dwindow[sample['name']].append(myDepth)
                #freq += minQ#increment freq by the number of ALT alleles
                continue
        else:
            if sample['GT'] == "0/0":
                window[sample['name']].append(0)
                dwindow[sample['name']].append(myDepth)

#    if freqFail:#if some filter failed return -1
#        freq = -1
#        return (freq)

    return (window, dwindow)#return minor allele frequency
 
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


def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:i:d:D:q:L:Q")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-i":
            global _i 
            _i = float(arg)
        elif opt == "-d":
            global _d 
            _d = int(arg)
        elif opt == "-D":
            global _D 
            _D = int(arg)
        elif opt == "-q":
            global _q 
            _q = int(arg)
        elif opt == "-L":
            global _L 
            _L = int(arg)
        elif opt == "-Q":
            global _Q 
            _Q = False
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" VCF [-w window_size -i indel_cutoff -d min_depth -D max_depth -q min_qual -L min_likelihood]"
def usage():
    print (use)
    sys.exit()
    
def details():
    descrip = "This program takes a VCF and outputs the number of heterozygotes and homo-non-ref for each individual in windows across the genome."
    print(descrip)
    print (use)
    print("\nOption details below are printed in the form:\n\t OPTION_NAME ARGUMENT_TYPE (DEFAULT OPTION) - DESCRIPTION\n")
    print("-w INT ("+str(_w)+") - Size of the window to analyze sequence in.")
    print("-i FLOAT ("+str(_i)+") - Maximum 'Dels' value a site can have and be included in analysis.")
    print("-d INT ("+str(_d)+") - Minimum depth to include a site in the analysis, this is the individual, not total site, depth.")
    print("-D INT ("+str(_D)+") - Maximum depth to include a site in the analysis, this is the individual, not total site, depth.")
    print("-q INT ("+str(_q)+") - Minimum quality to include a site in the analysis.")
    print("-L INT ("+str(_L)+") - Minimum likelihood of a particular individual's genotype call to include a site in the analysis.")
    print("-Q NONE ("+str(_Q)+") - Indicates that we should not throw out sites with no genotype qualities (for use with VCF v4.1+)")

    
if __name__ == "__main__":   
    __main__()
