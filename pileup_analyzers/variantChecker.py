import sys
import re
import getopt

_l = 40
_r = False
_s = False

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"l:sr")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _l 
            _l = int(arg)
        elif opt == "-s":
            global _s 
            _s = True
        elif opt == "-r":
            global _r 
            _r = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    vcf = open(sys.argv[1],"r")
    
    if (vcf == None):
        print("Bad vcf name: "+sys.argv[1])
        sys.exit()
        
    fasta = open(sys.argv[2],"r")  
        
    if (fasta == None):
        print("Bad fasta file name: "+sys.argv[1])
        sys.exit()
    
    #fetch all the samples from the fasta
    fastaSamps = {}
    samp = ""
    seq = ""
    for line in fasta:
        line = line.rstrip() 
        if line != "" and line[0] == ">":
             if (samp != ""):
                fastaSamps[samp] = seq
             #print (line)
             samp = line.split()[0][1:]#TODO
             #print (samp)
             seq = ""
        else:
            seq += line
            
    fastaSamps [samp] = seq
        
    
        
    #reg ex
    vcf_pat = re.compile("\w+_(\S+)\s(\d+)\s\S+\s(\S+)\s(\S+)\s\S+\s\S+\s(\S+)\s\S+\s(.+)")
    
    #figure out which column is which sample
    line = vcf.readline()
    while (line.split()[0] != '#CHROM'):
        line = vcf.readline()
    header = line.split()[9:]
    
    vcfSamps = {}
    sampIndexs = {}
    for i in range(0, len(header)):
        tmp = header[i].split(".sort")[0]
        sys.stderr.write("Name found: "+tmp+"\n")
        #fixed = ".".join(tmp)
        if tmp in fastaSamps.keys():
            sampIndexs[i] = tmp
            vcfSamps[tmp] = ""
    
    #print (sampIndexs.keys())
    
    pAdded = ""
    #read the sequences for each individual from the vcf 
    for line in vcf:
        line = line.rstrip()
        m = vcf_pat.match(line)
        
        if (m == None):# if we find a weird line...
            sys.stderr.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            scaf = m.group(1)
            base = m.group(2)
            ref = m.group(3)
            alt = m.group(4)
            stuff = m.group(5)
            genos = m.group(6).split()
            
            #sys.stderr.write(scaf+" "+base+"\n")
            
            if pAdded == base:#repeat base
                sys.stderr.write("Previous line is the same as this one! "+ pAdded + " " + base +"\n")
                continue
            
            deletion = False
            if stuff.split("=")[0]=="Dels" or stuff[0]=="." or stuff.split(";")[0]=="INDEL":#at deletion
                sys.stderr.write("Deletion detected: "+str(scaf)+" "+str(base)+"\n")
                deletion = True
            
            for i in sampIndexs.keys():
                if not deletion:
                    vcfGeno = genotypeCall(ref, alt, genos[i])
                else:
                    vcfGeno = "N"
                vcfSamps[sampIndexs[i]] += vcfGeno
                
            pAdded = base
            #sys.stderr.write(base+"\n")
            continue

    if _r:#reverse vcf sequences
        for i in vcfSamps.keys():
            vcfSamps[i] = vcfSamps[i].reverse()

    totalMatch = totalFalsePositive = totalMis = totalMisHom = totalN = totalFalseNegative = 0


    results = {}
    for samp in fastaSamps.keys():
        if samp not in vcfSamps.keys():
            sys.stderr.write("Samp missing from VCF: "+samp+"\n")
            continue
        
        data = eval(fastaSamps[samp], vcfSamps[samp])
        
        if data != None:
            (match, fp, fn, miscalled, misHom, nCount) = data
            results[samp] = data
            totalMatch += match
            totalFalsePositive += fp
            totalFalseNegative += fn
            totalMis += miscalled
            totalMisHom += misHom
            totalN += nCount
        else:
            sys.stderr.write("Eval returned None for: "+str(samp)+"\n")
        
    total = totalMatch+totalFalsePositive+totalFalseNegative+totalMis+totalN+totalMisHom
    print ("S\tM\tFP\tFN\tMisHet\tMisHom\tN\tT")
    print ("Tot\t"+str(totalMatch)+"\t"+str(totalFalsePositive)+"\t"+str(totalFalseNegative)+"\t"+str(totalMis)+"\t"+str(totalMisHom)+"\t"+str(totalN)+"\t"+str(total))
    print("-------------------------------------------------------------")
    for samp in results.keys():
        data = results[samp]
        print(samp+"\t"+str(data[0])+"\t"+str(data[1])+"\t"+str(data[2])+"\t"+str(data[3])+"\t"+str(data[4])+"\t"+str(data[5])+"\t"+str(sum(data)))
  
def eval(FastaL, VcfL):
    Match = 0
    FalseNeg = 0
    FalsePos = 0
    HetCount = 0
    MissCalledHet = 0
    MissCalledHom = 0
    NCount = 0
    if len(FastaL) != len(VcfL):
        sys.stderr.write("Fasta Len: "+str(len(FastaL))+" VCF Len: "+str(len(VcfL))+"\n")
        return None #confirms that seq are the same length
    else:
        #loop through each base in the sequence
        for i in range(len(FastaL)):
            #find anything with any sort of het call going on
            if (VcfL[i] in ["N","-"]) or (FastaL[i] in ["N","-"]):
                NCount = NCount + 1
            elif VcfL[i] == FastaL[i]:
                Match = Match + 1
                if VcfL[i] not in ["A","C","G", "T"]:
                    HetCount += 1
            else:
                #count false positives (VCF has het that the Fasta doesn't)
                if (FastaL[i] in ["A","C","G","T"]) and (VcfL[i] not in ["A","C","G","T"]):
                    FalsePos = FalsePos + 1
                #count false neg (Fasta has het that VCF doesn't)
                elif (VcfL[i] in ["A","C","G","T"]) and (FastaL[i] not in ["A","C","G","T"]):
                    FalseNeg = FalseNeg + 1
                #both have different hets
                elif (VcfL[i] not in ["A","C","G","T"]) and (FastaL[i] not in ["A","C","G","T"]):
                    MissCalledHet += 1
                #both have different homos
                else:
                    MissCalledHom += 1
            
    #count number of missed Het calls (both files have a Het, but have different ideas of what base it is -- should we say that this is 'correct' anyway?
    #MissCalledHet = HetCount - Match - FalsePos - FalseNeg
    result = (Match, FalsePos, FalseNeg, MissCalledHet, MissCalledHom, NCount)
    return(result)

def genotypeCall(ref, alt, vcf):
    #sys.stderr.write(str(vcf)+"\n")
    #samtools or gatk? split each
    if _s:
        if vcf == "0":
            return(ref)
        
        data = vcf.split(",")
    
    else:
        data = vcf.split(":")
        data = data[-1].split(",")
    
    #make dictionary of fasta codes for hets
    key = dict([('GA','R'),('AG','R'),('AC','M'),('CA','M'),('GC','S'),('CG','S'),('TC','Y'),('CT','Y'),('GT','K'),('TG','K'),('AT','W'),('TA','W')])
    
    #sys.stderr.write(str(data)+"\n")
    if (len(ref) > 1 or len(alt) > 1):
        return "N"
    
    if len(data) > 3:
        return "N"
    
    #make calls
    if (float(data[1]) == 0) and (float(data[0]) > _l) and (float(data[2]) > _l) and alt != ".":#if het but alt is . then we are N
        _f = ref + alt
        return(key[_f]) #het
    elif (float(data[0]) == 0) and (float(data[1]) > _l) and (float(data[2]) > _l):
        return(ref) #hom reference
    elif (float(data[2]) == 0) and (float(data[1]) > _l) and (float(data[0]) > _l):
        return(alt) #hom alt
    else:
        return("N") #unclear
   
use = "python "+__file__.split("/")[-1]+" VCF Fasta [options]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("-l - (40) - Likelihood cutoff for calling a site from vcf genotype calls")
    print("-s - (False) - flag indicating that the input is a vcf generated by samtools")
    print("-r - (False) - flag indicating that gene we are checking is on the negative strand (i.e. reverse the sequence we get from the vcf)")

    
if __name__ == "__main__":   
    __main__()