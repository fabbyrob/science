import sys
import getopt

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
        
    processArgs(3)
    
    vcf_fastas = open(sys.argv[1],"r")
    
    if (vcf_fastas == None):
        print("Bad vcf name: "+sys.argv[1])
        sys.exit()
        
    sanger_fastas = open(sys.argv[2],"r")  
        
    if (sanger_fastas == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()
        
    vfastas = {}
    
    #read in all the fastas from the VCF
    #will exist in name:Sequence pairs, where Sequence is a list of two items, one for each chromosome
    for line in vcf_fastas:
        line = line.rstrip()
        if line == "":
            continue
        
        if line[0] == ">":#new sequence!
            name = line[1:].split("_")[0]
        else:
            if name in vfastas.keys():
                vfastas[name].append(line)
            else:
                vfastas[name] = [line]
    
    sfastas = {}
    #read in all the sanger data of the format name:Sequence, where Sequence is a string of the sequence
    for line in sanger_fastas:
        line = line.rstrip()
        if line == "":
            continue
        
        if line[0] == ">":#new sequence!
            name = line[1:].split(" ")
            if name [0] == "":
                name = name[1]
            else:
                name = name[0]
        else:
            if name in sfastas.keys():
                sys.stderr.write("Duplicate sample in sanger data, ignoring second sequence for: "+ind+"\n")
            else:
                sfastas[name] = line
    
    #print(sfastas)
    
    #nucleotide codes
    nucs = {"-":["N","N"],"N":["N","N"],"A":["A","A"],"T":["T","T"],"C":["C","C"],"G":["G","G"],"R":["A","G"], "Y":["C","T"], "K":["G", "T"], "M":["A", "C"], "S":["C", "G"], "W":["A", "T"]}
    matches = {}
    size = 0
    tSize = 0
    totals = [0,0,0,0,0,0]#missmatches, matches, false positive (VCF shows vairant sanger doesnt), false negative (VCF doesn't show variant sanger does), Ns (sites where VCF or sanger is N), # snps
    for ind in vfastas.keys():
        if ind not in sfastas.keys():
            sys.stderr.write("Missing ind from sanger data: "+ind+"\n")
            continue
        

        if size == 0:
            size = len(sfastas[ind])
 
        if len(vfastas[ind][0]) != len(sfastas[ind]):
            sys.stderr.write("VCF and sanger sequences are not the same length for: "+ind+"\n")
            continue
        
        if len(sfastas[ind]) != size:
            sys.stderr.write("Sequences not the correct length for: "+ind+"\n")
            continue
        
        tSize += size
        v = vfastas[ind]
        s = sfastas[ind]
        
        if len(v) != 2:
            sys.stderr.write("Incorrect number of sequences from VCF ("+str(len(v))+")\n")
            continue
        
        match = 0
        missmatch = 0
        fp = 0
        fn = 0
        ns = 0
        snps = 0
        for i in range(len(v[0])):
            v1 = v[0][i]
            v2 = v[1][i]
            si = s[i]
            
            #missing data don't count it
            if v1 == "N" or v2 == "N" or si == "N" or si == "-":
                ns += 1
                continue
            
            s1 = nucs[si][0]
            s2 = nucs[si][1]
            
            if (v1, v2) == (s1, s2) or (v1, v2) == (s2, s1):
                #we're good!
                match += 1
                if v1 != v2:
                    snps += 1
                continue
            else:
                sys.stderr.write("Sites do not match "+ind+" site# "+str(i)+"\n\tVCF: "+str((v1, v2))+"\n\tsanger: "+str((s1, s2))+"\n")
                missmatch += 1
                
                if v1 == v2 and s1 != s2:#false neg
                    fn += 1
                elif v1 != v2 and s1 == s2:#false pos
                    fp += 1
                continue
            
        matches[ind] = (missmatch, match, fp, fn, ns, snps)
        totals[0] += missmatch
        totals[1] += match
        totals[2] += fp
        totals[3] += fn
        totals[4] += ns
        totals[5] += snps
        
    if sum(totals) != 0:
        print("Ind\tMatch\tMiss\t%Match\tFalse_pos\t%FP\tFalse_neg\t%FN\tWrong_Homo\t%WH\tMissing_data/Total_data\tSNPs")
        printData("Total", totals, tSize)
        for ind in matches.keys():
            vals = matches[ind]
            printData(ind, vals, size)
    else:
        print("No individuals had correct data.")

    

def printData(name, data, size):
    miss = data[0]
    match = data[1]
    fp = data[2]
    fn = data[3]
    ns = data[4]
    snps = data[5]
    wrongHomo = miss-fp-fn
    if miss+match != 0:
        pMatch = 100*float(match)/(miss+match)
        pfp = 100*float(fp)/(miss+match)
        pfn = 100*float(fn)/(miss+match)
        pwrongHomo = 100*float(wrongHomo)/(miss+match)
    else:
        pMiss=pfp=pfn=pwrongHomo=0.0
    print(str(name+"\t"+str(match)+"\t"+str(miss)+"\t%.2f\t"+str(fp)+"\t%.2f\t"+str(fn)+"\t%.2f\t"+str(wrongHomo)+"\t%.2f\t"+str(ns)+"/"+str(size)) % (pMatch, pfp, pfn, pwrongHomo)+"\t"+str(snps))
def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        continue
#        if opt == "-o":
#            global _o 
#            _o = arg
#        else:
#            print ("Unrecognized option: "+opt+"\n")
#            usage()
   
use = "python "+__file__.split("/")[-1]+" VCF_fasta Sanger_fasta"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()