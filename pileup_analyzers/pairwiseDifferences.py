'''
2012 Robert williamson
calculates the differences between pairs of samples from a VCF. Run the program with no options for full decription of input.

python pairwiseDifferences.py

'''

import sys
import getopt
import vcf

_i = 0 #deletetion cutoff
_d = 20 #min depth
_D = 60 #max depth
_l = 40 #likelihood cutoff
_q = 90 #min quality

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    reader = vcf.Reader(open(sys.argv[1], 'rb'))
    
    samps = reader.samples
    combs = len(samps)*(len(samps)-1)/2
    diffs = [0]*combs#our list of differences
    
    for record in reader:
        #determine genotype of each sample, after filtering
        #0 = homo ref, 0.5 = het, 1 = homo alt, -1 = no info
        genos = []
        if "Dels" not in record.INFO.keys() or float(record.INFO['Dels']) > _i:
            #site spans a deletion, skip it 
            continue
        
        if record.QUAL < _q:
            #site quality too low, skip it
            continue
        
        for samp in record.samples:
            #check the sample depth
            if 'DP' not in samp.keys() or samp['DP'][0] < _d or samp['DP'][0] > _D:
                #depth too low, or not present, no data available
                genos.append(-1)
                continue
            #check the sample likelihood
            if 'PL' not in samp.keys():
                #no likelihood data
                genos.append(-1)
                continue
            
            myPLs = map(int, samp['PL'])
            (minQ, midQ, maxQ) = getMinMax(myPLs)
            
            if myPLs[minQ] != 0:
                #lowest likelihood not zero, skip
                genos.append(-1)
                continue
            
            if myPLs[midQ] < _l:
                #next closest is too likely, skip
                genos.append(-1)
                continue
            
            #everything is good, store the genotype
            genos.append(minQ*0.5)
            
           
        if -1 in genos:#some individual did not pass filters, skip site
            continue
        
        #add the differences for this site to the totals 
        comp = 0
        for i in range(len(genos)):
            for j in range(i+1, len(genos)):
                diff = abs(genos[i]-genos[j])
                diffs[comp] += diff
                comp += 1

    #print out the total differences
    print("SAMP1\tSAMP2\tDIV")
    comp = 0
    for i in range(len(samps)):
        for j in range(i+1, len(samps)):
            print("%s\t%s\t%.1f" % (samps[i], samps[j], diffs[comp]))
            comp += 1

def getMinMax(ls):
    min = 0
    mid = 0
    max = 0
    
    for i in range(1, len(ls)):
        if ls[i] <= ls[min]:
            mid = min
            min = i
        elif ls[i] >= ls[max]:
            mid = max
            max = i
        else:
            mid = i
    return (min, mid, max)

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"d:D:i:l:q:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-d":
            global _d 
            _d = int(arg)
        elif opt == "-D":
            global _D 
            _D = int(arg)
        elif opt == "-i":
            global _i 
            _i = float(arg)
        elif opt == "-l":
            global _l 
            _l = int(arg)
        elif opt == "-q":
            global _q 
            _q = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" VCF [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("Takes a VCF and calculates the difference between each pair of samples. \
(Difference is incremented by 1 for homozygous differences, and 0.5 when one sample is heterozygous).")
    print("Sites where any individual does not pass filters (see below), are ignored.")
    print("Output consists of three columns: the two samples being compared, and the total difference:\n")
    print("SAMP1\tSAMP2\tDIV\nsample_korfu\tsample_greence\t12001.5\nsample_korfu\tsample_albania\t15001.0\nsample_greece\tsample_albania\t10231.0")
    print("")
    print("option - TYPE - default - description")
    print("d - INT - %d - minimum individual read depth" % _d)
    print("D - INT - %d - maximum individual read depth" % _D)
    print("l - INT - %d - minimum genotype likelihood of second most likely genotype (most likely must be 0)" % _l)
    print("q - INT - %d - minimum site quality score" % _q)
    print("i - FLOAT - %.2f - maximum fractions of reads with a spanning deletion at a site (\'Dels\' value in VCF\'s INFO column)" % _i)

    
if __name__ == "__main__":   
    __main__()