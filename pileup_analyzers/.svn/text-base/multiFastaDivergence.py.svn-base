doc = """
Given a file containing a list of many fastas containing alignments in many files and produces a pickled dictionary\
of lists indicating divergence states at each site.

Header lines are expected to be formatted:
>CG_scaffold3517__184728-185030__strand+
>CO_scaffold-61__4438-4740__strand+

>SPECIES_CHROM__START-END__STRAND
"""

import sys
import getopt
import cPickle

_v = False

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 5:
        usage()
    
    spp1 = sys.argv[1]
    spp2 = sys.argv[2]
    refSpp = sys.argv[3]
    output = sys.argv[4]
    
#    if sys.argv[5] == "n":
#        _v = False
#    else:
#        _v = True
    
    divergence = {}
    
    for f in open(sys.argv[5], "r"):
        infile = open(f.rstrip(), "r")
        
        alignment = {}#spp:{scaf, start, stop, seq}
        
        spp = ""
        seq = ""
        for line in infile:
            if line.startswith(">"):
                if spp:
                    alignment[spp]["seq"] = seq
                    seq = ""
                    
                sline = line[1:].split('_')
                spp = sline[0]
                scaf = sline[1]
                #sys.stderr.write("%s\n" % sline[2])
                start, end = map(int, sline[3].split("-"))
                alignment[spp] = {"scaf":scaf, "start":start, "end":end}
            else:
                seq += line.upper()
                
        localDiv = []
    
        if _v: 
            tot = 0
            diff = 0
            
        for i in range(len(alignment[spp1]["seq"])):
            if alignment[spp1]["seq"][i] == "N" or alignment[spp2]["seq"][i] == "N":
                localDiv.append(-1)
            elif alignment[spp1]["seq"][i] == "-" or alignment[spp2]["seq"][i] == "-":
                localDiv.append(-1)
            elif alignment[spp1]["seq"][i] == alignment[spp2]["seq"][i]:
                localDiv.append(0)
                if _v: 
                    tot += 1
            else:
                localDiv.append(1)
                if _v: 
                    tot += 1
                    diff += 1
        
        if _v:
            sys.stderr.write("For alignment %s %s-%s, %s sites had data %s were divgerged.\n" % (alignment[refSpp]["scaf"], alignment[refSpp]["start"], alignment[refSpp]["end"], tot, diff))
        
        #extend the divergence list for the scaffold if needed
        if alignment[refSpp]["scaf"] not in divergence.keys():
            divergence[alignment[refSpp]["scaf"]] = [-1]*(alignment[refSpp]["end"]+1)
        elif len(divergence[alignment[refSpp]["scaf"]])-1 < alignment[refSpp]["end"]:
            divergence[alignment[refSpp]["scaf"]].extend([-1]*(alignment[refSpp]["end"]-len(divergence[alignment[refSpp]["scaf"]])+1))
        
        #add the divergence info to the dictionary
        for pos, div in zip(range(alignment[refSpp]["start"], alignment[refSpp]["end"]+1), localDiv):
            #sys.stderr.write("%s %s %s\n"%(len(divergence[alignment[refSpp]["scaf"]]), alignment[refSpp]["end"], pos))
            divergence[alignment[refSpp]["scaf"]][pos] = div
        
    #pickle the dictionary 
    cPickle.dump(divergence, open(output, 'wb'))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"")
    except getopt.GetoptError:
        usage()
    
#    for opt, arg in opts:
#        if opt == "-v":
#            global _v 
#            _v = True
#        else:
#            print ("Unrecognized option: "+opt+"\n")
#            usage()

    sys.stderr.write("species to compare: %s %s reference species: %s output file: %s  file list: %s\n" % (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]))
   
use = "python "+__file__.split("/")[-1]+" Species1 Species2 refSpecies outputFile fileList"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")

    
if __name__ == "__main__":   
    __main__()