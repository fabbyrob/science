doc = """
Takes in a global alignment. Calculates lineage specific divergence at a given set of sites.
"""

import sys
import getopt
import cPickle as pickle

_w = 5
_p = "divergence.pickle"
_m = False

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    alignment = open(sys.argv[1],"r")

    goodScafs =  ["scaffold_1","scaffold_2","scaffold_3","scaffold_4","scaffold_5","scaffold_6","scaffold_7","scaffold_8"]

    align = {}
    
    divergence = {}

    #read the infile file...
    for line in alignment:
        line = line.rstrip()
        sline = line.split()
        
        if len(line) == 0:#blank line; process old alignment
            if align.has_key('rubella') and align.has_key('arabidopsis') and align.has_key('neslia'):
                processDivergence(align, divergence)

            align = {}
            continue
        
        if line[0] == "a":
            continue
        elif line[0] == "s":
            if line.startswith("s Scaffold_"):#CR
                if sline[1].lower() not in goodScafs:#skip extra scafs
                    continue
                
                align['scaf'] = sline[1]
                align['start'] = int(sline[2])+1
                align['end'] = align['start']+int(sline[3])-1
                align['rubella'] = sline[6].upper()
                
                if not divergence.has_key(align['scaf']):  
                    divergence[align['scaf']] = [None]+[-1]*int(sline[5])
                
            elif line.startswith("s Scaffold"):#AT
                align['arabidopsis'] = sline[6]
            elif line.startswith("s Chr"):#NP
                align['neslia'] = sline[6]
    
    if align.has_key('rubella') and align.has_key('arabidopsis') and align.has_key('neslia'):
        processDivergence(align, divergence)

    pickle.dump(divergence, open(_p, "wb"))
            
            
def processDivergence(align, divergence):
    offset = 0
    
    if not _m:
        align['rubella'] = align['rubella'].upper()
        align['neslia'] = align['neslia'].upper()
        align['arabidopsis'] = align['arabidopsis'].upper()
    
    for i in range(0, len(align['rubella'])):
        if align['rubella'][i] not in ['A','T','C','G']:
            continue
        
        base = align['start']+offset
        offset += 1
            
        minSafe = max(0, i-_w)   
        maxSafe = min(len(align['rubella'])-1, i+_w)+1     
            
        try:
            if "-" in align['rubella'][minSafe:maxSafe] or "-" in align['neslia'][minSafe:maxSafe] or "-" in align['arabidopsis'][minSafe:maxSafe]:
                #indel too close to site in one of the sequences
                divergence[align['scaf']][base] = -1
                continue
                
            if align['rubella'][i] == align['neslia'][i] and align['rubella'][i] == align['arabidopsis'][i]:
                divergence[align['scaf']][base] = 0
            elif align['rubella'][i] != align['neslia'][i] and align['neslia'][i] == align['arabidopsis'][i]:
                divergence[align['scaf']][base] = 1
            else:
                divergence[align['scaf']][base] = -1
        except IndexError:
            sys.stderr.write("Tried to assign divergence for site %s but there are only %s sites on chromosome %s. Alignment start: %s. Skipping.\n" % (base, len(divergence[align['scaf']]), align['scaf'], align['start']))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:p:m")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-p":
            global _p 
            _p = arg
        elif opt == "-m":
            global _m 
            _m = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("alignment: %s -w %s -p %s -m %s\n" % (sys.argv[1], _w, _p, _m))
   
use = "python "+__file__.split("/")[-1]+" alignmentFile.txt [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("w - INT - %s - number of bases on either side of an indel to ignore." % _w)
    print("p - STR - %s - output file to save the divergence info to in a binary format." % _p)
    print("m - NONE - %s - using this flag will make the script assign a divergence value of unknown for any sites that have been masked (i.e. the bases in thealignment are lower case)." % _m)

    
if __name__ == "__main__":   
    __main__()