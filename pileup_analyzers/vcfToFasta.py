#!/usr/bin/python
import sys
import re
import getopt

_l = 60
_q = 40
_d = 20
_D = 60
_w = 5
_log = "vcfToFasta.log"
_o = ""
_input = ""
_R = False
_s = False
_G = False

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"q:d:D:L:w:l:o:RsG")
    except getopt.GetError:
        usage()
    
    for opt, arg in opts:
        if opt == "-q":
            global _q
            _q = int(arg)
        elif opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-d":
            global _d 
            _d = int(arg)
        elif opt == "-D":
            global _D 
            _D = int(arg)
        elif opt == "-L":
            global _l 
            _l = int(arg)
        elif opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-o":
            global _o 
            _o = str(arg)
        elif opt == "-R":
            global _R 
            _R = True
        elif opt == "-s":
            global _s 
            _s = True
        elif opt == "-G":
            global _G 
            _G = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    vcf = open(sys.argv[2],"r")
    
    if (vcf == None):
        print("Bad vcf name: "+sys.argv[2])
        sys.exit()
        
    global _input 
    _input = sys.argv[1]
    annot = open(sys.argv[1],"r")  
        
    if (annot == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    logfile.write("Running with: -q "+str(_q)+" -d "+str(_d)+" -D "+str(_D)+" -L "+str(_l)+" -w "+str(_w)+"\n")
    logfile.write("Outputting files to: "+_o+"\n")
    if(_R):
        logfile.write("Running in outgroup mode (-R)\n")
        
    #reg ex
    #1 = scaf, 2 = pos
    if (not _R):
        if not _G:
            indel = re.compile("\w+_(\d+)\s(\S+).+INDEL.+")
        else:
            indel = re.compile("\w+_(\d+)\s(\S+).+\tDels=.+")
    else:
        indel = re.compile("\w+_\d+\s\S+\s\*.+")
    #1 = scaf, 2 = pos, 3 = ref, 4 = alt, 5 = qual, 6 = depth, 7 = inds
    if (not _R):
        if not _G:
            vcf_pat = re.compile("\w+_(\d+)\s(\d+)\s\S+\s(\S+)\s(\S+)\s(\S+)\s\S+\s.*DP=(\d+);\S+\s\S+\s(.+)")
        else:
            vcf_pat = re.compile("\w+_(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(.+)")
    else:
        #1 = scaf, 2 = pos, 3 = ref, 4 = alt, 5 = snpqual, 6 = depth
        vcf_pat = re.compile("\w+_(\d+)\s(\d+)\s(\S+)\s(\S+)\s(\S+)\s\S+\s\S+\s(\S+).+")
        
    #1 = scaffold, 2 = start, 3 = end, 4 = model, 5 = +/-
    if (not _s):
        annotpat = re.compile("\w+_(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(.*)")
    else:
        annotpat = re.compile("\w+_(\d+)\s+(\S+)")
    
    #grab all of the individual IDs from the second line in the vcf
    if not _R:
        line = vcf.readline().split()
        while line[0] != "#CHROM":
            line = vcf.readline().split()
        
        numInds = len(line)-9
        indIDs = line[9:len(line)]
    else:
        numInds = 1
        indIDs = ["outgroup"]
    
    logfile.write("Detected "+str(numInds)+" individuals\n")
    
    #initialize our sequence dictionaries
    seq1 = {}
    seq2 = {}
    for i in range(0, numInds):
        seq1[i] = []
        seq2[i] = []
    
    #initialize other variables
    indelCount = 0
    
    lowest = -1
    maxest = -1
    modScaf = -1
    modMin = -1
    modMax = -1
    mod = ""
    modDir = ""
    modP = ""
    modBlast = ""
    dirP = ""
    scafP = ""
    addedP = -1
    genos = "" 
    
    #grab our first model to search for
    if (not _s):
        (mod, modScaf, modMin, modMax, modDir, modBlast) = getNextModel(annot, annotpat)
    else:
        (modScaf, modMin) = getNextSite(annot, annotpat)
        modMax = modMin
    lowest = modMin
    maxest = modMax
    
    #read the vcf file...
    for line in vcf: 
        indelCount += 1
        indelFlag = 0
        line = line.rstrip()
        m = vcf_pat.match(line)
        m2 = indel.match(line)
        
        if(m2 != None):
            logfile.write("INDEL found at\n\t"+line+"\n")
            indelFlag = 1
            indelCount = -1*_w
            
        if (m == None):
            logfile.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            #grab the relevant data
            scaf = int(m.group(1))
            base = int(m.group(2))
            ref = m.group(3)
            alt = m.group(4)
            qual = m.group(5)
            if qual == ".":
                qual = 0
            else:
                qual = float(qual)
            if not _G:
                depth = int(m.group(6))
                if (not _R):
                    genos = m.group(7)
            else:
                pat = m.group(6).split(":")
                if 'DP' in pat:
                    depthIndex = pat.index('DP')
                    genos = m.group(7)
                    #depth = getTotalDepth(genos, depthIndex)
                    depth = depthIndex
                else:
                    depth = -1
            
            #grab next model if we reached the end or skipped a site
            if(base > modMax or scaf > modScaf):
                if(base > modMax and addedP < modMax):
                    if (not _s):
                        if addedP == -1:
                            addedP = modMin-1
                        #if the VCF skipped some lines mid intron, then add N's to the sequences
                        for i in range(addedP, modMax):
                            logfile.write("VCF skipped to base: "+str(scaf)+", "+str(base)+". Calling ambiguous\n")
                            allSet(seq1, seq2)
                    else:
                        logfile.write("VCF skipped base: "+str(modScaf)+", "+str(modMin)+". Calling ambiguous\n")
                        allSet(seq1, seq2)
                        
                if not _s:
                    modP = mod
                    dirP = modDir
                    scafP = modScaf
                    vals = getNextModel(annot, annotpat)
                    addedP = -1
                    if vals == None:
                        break
                    (mod, scaf, start, end, dir, blast) = vals
                    
                    if mod != modP:#done with gene output data
                        printFasta(seq1, seq2, indIDs, modP, scafP, lowest, maxest, modBlast, dirP)
                        reset(seq1, seq2)
                        lowest = start
                        modScaf = scaf
                        modMin = start
                        modMax = end
                        maxest = end
                        modDir = dir
                        modBlast = blast
                    else:#on next exon
                        maxest = end
                        modMin = start
                        modMax = end
                        modBlast += blast
                else:#grab the nest site in the sitelist
                    vals = getNextSite(annot, annotpat)
                    if vals == None:
                        break
                    (modScaf, modMin) = vals
                    modMax = modMin
            
            #if we are within a model range
            if base >= modMin and (_s or base <= modMax) and scaf == modScaf:
                if(not _s and base > modMin and base > addedP + 1):
                    if addedP == -1:
                        addedP = modMin
                        end = base
                    else:
                        end = base-1
                    #if the VCF skipped some lines mid intron, then add N's to the sequences
                    for i in range(addedP, end):
                        logfile.write("VCF skipped base: "+str(scaf)+", "+str(base)+". Calling ambiguous\n")
                        allSet(seq1, seq2)
                        
                addedP = base
                
                if (not _s) and indelFlag:#if we are at an indel
                    makeWindowAmbig(seq1, -1*indelCount)
                    makeWindowAmbig(seq2, -1*indelCount)
                    allSet(seq1, seq2)
                    indelSite = base
                    continue
                
                if (not _s) and indelCount <= 0:
                    #if we are in a site list, and still WITHIN the ambig window
                    #OR if we are doing an annotation
                   # if (_s and abs(base - indelSite) < _w) or not _s:
                        #we are in an INDEL window print ambig
                    allSet(seq1, seq2)
                    continue
                
                #decide if the site is het/homo/ambig for each individual
                callSite(scaf, base, ref, alt, qual, depth, genos, seq1, seq2, line, logfile)
                    
                #get the next site from the site file
                if _s:
                    vals = getNextSite(annot, annotpat)
                    if vals == None:
                        break
                    (modScaf, modMin) = vals
                    modMax = modMin    
        
    #Need to look at _w more sites to make sure the end of our final model is not in an indel window
    i = 0
    indelFlag = 0
    for line in vcf:
        line = line.rstrip()
        m = vcf_pat.match(line)
        m2 = indel.match(line)
        
        if(m2 != None):
            logfile.write("INDEL found at\n\t"+line+"\n")
            indelFlag = 1
            indelCount = -1*_w
            
        if (m == None):
            logfile.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            #grab the relevant data
            base = int(m.group(2))
                
            if base > modMax + _w:
                break
            
            if indelFlag:
                indelCount += base-modMax-1
                makeWindowAmbig(seq1, -1*indelCount, False)
                break
                        
    #print last chuck of data
    printFasta(seq1, seq2, indIDs, modP, scafP, lowest, maxest, modBlast, dirP)
        
    vcf.close()
    annot.close()
    logfile.close()
 
def callSite(scaf, base, ref, alt, qual, depth, genos, seq1, seq2, line, logfile):
    if(alt not in ['A','T','C','G','.']):
        #not a biallelic site
        #ambig
        logfile.write("Non-standard base for alternate: "+alt+" at "+str(scaf)+", "+str(base)+". Calling ambiguous\n")
        allSet(seq1, seq2)
        return 
    
    if qual < _q:
        logfile.write("Site below quality threshold.\n\t"+line+"\n")
        allSet(seq1, seq2)
        return 
                    
    if not _G and depth < _d or depth > _D:
        logfile.write("Site outside depth range.\n\t"+line+"\n")
        allSet(seq1, seq2)
        return 
    
    if _R and qual >= _q:
        #nonreference
        #logfile.write("Possible heterozygote.\n\t"+line+"\n")
        allSet(seq1, seq2, alt)
        return 
    
    if (not _R):
        inds = genos.split()
    
        #everybody homozygous reference
        if not _G and (inds[0].split(":") == [inds[0]]):
            allSet(seq1, seq2, ref)
        else:
            #other
            c = 0
            for i in inds:
                if _G and i == "./.":#missing data for this individual TODO: WHOLE SITE to N
                    bothN(seq1, seq2, c)
                    c += 1
                    continue
                    
                data = i.split(":")
                if _G:
                    if depth == -1:
                        myDepth = 0
                    else:
                        myDepth = int(data[depth])
                    if myDepth < _d or myDepth > _D:
                        #did not pass depth filters
                        bothN(seq1, seq2, c)
                        c += 1
                        continue
                
                if not _G:
                    data = data[0].split(",")
                else:
                    data = data[-1].split(",")
                makeInts(data)
                
                (min, max, mid) = getMinMaxIndex(data)
                
                if (data[min] != 0):
                    #site is ambiguous, no max likelihood
                    bothN(seq1, seq2, c)
                    c += 1
                    continue
                
                if (data[mid] < _l):
                    #dont pass quality cut off print N's
                    bothN(seq1, seq2, c)
                else:
                    #print(str(i)+" "+str(scaf)+" "+str(base))
                    if min == 1:#heterozygote
                        if(not _R):
                            seq1[c].append(ref)
                            seq2[c].append(alt)
                        else:#outgroup should be homozygous for everything
                            bothN(seq1, seq2, c)
                        #print ("het" + str(seq1[c]))
                    elif min == 0:#homozygous ref
                        seq1[c].append(ref)
                        seq2[c].append(ref)
                        #print ("hom r"+ str(seq1[c]))
                    elif min == 2:#homozygous alt
                        seq1[c].append(alt)
                        seq2[c].append(alt)
                        #print ("hom a"+ str(seq1[c]))
                c+=1  
        return 
  
def getNextSite(file, pat):
    dat = None
    while(dat == None):
        aLine = file.readline()
        
        if (aLine == ""):
            return
        
        m = pat.match(aLine)
        
        if(m == None):
            continue
        scaf = int(m.group(1))
        site = int(m.group(2))
        dat = (scaf, site)
    return dat

def getNextModel(file, pat):
    mod = None
    while (mod == None):
        aLine = file.readline()
        
        if (aLine == ""):
            return
        
        m = pat.match(aLine)
        
        if(m == None):
            continue
        modP = mod
        
        scaf = int(m.group(1))
        start = int(m.group(2))
        end = int(m.group(3))
        mod = m.group(4)
        dir = m.group(5)
        blast = m.group(6)
            
    return (mod, scaf, start, end, dir, blast)
  
def reset(seq1, seq2):
    for k in seq1.keys():
        seq1[k] = []
        seq2[k] = []
    
def printFasta(seq1, seq2, ids, model, scaf, min, max, blast, dir):
    if not _s:
        filename = model+"_"+str(min)+"_"+str(max)
    else:
        filename = _input.split("/")[-1]
    file = open(_o+filename+".fasta", "a")
    if file == None:
        logfile.write("Could not open fasta file: \n\t"+_o+model+".fasta"+"\noutput here instead:\n")

    blast = list(set(blast.split()))
    cblast = len(blast)
    blast = " ".join(blast)
    
    for k in seq1.keys():
        myStr = ""
        myStr += ">"+ids[k]+"_1 "+model+": "+"super_"+str(scaf)+" "+str(min)+" "+str(max)+" "+str(cblast)
        if (not _R):
            myStr += " "+blast+"\n"
        else:
            myStr += "\n"
            
        code = "".join(seq1[k])
        if(dir == "-"):
            code = sense(code)
        myStr += code
        
        if (file):
            file.write(myStr+"\n\n")
        else:
            logfile.write(myStr+"\n\n")
            
        if (not _R):
            myStr = ""
            myStr += ">"+ids[k]+"_2 "+model+": "+"super_"+str(scaf)+" "+str(min)+" "+str(max)+" "+str(cblast)+" "+blast+"\n"
            
            code = "".join(seq2[k])
            if(dir == "-"):
                code = sense(code)
            myStr += code
            
            if (file):
                file.write(myStr+"\n\n")
            else:
                logfile.write(myStr+"\n\n")

def sense(codon):
    ncodon = ""
    for b in codon:
        if b == "A":
            ncodon += "T"
        elif b == "T":
            ncodon += "A"
        elif b == "G":
            ncodon += "C"
        elif b == "C":
            ncodon += "G"
        elif b == "N":
            ncodon += "N"

    ncodon = ncodon[::-1]#reverse the string

    return ncodon
    
def bothN(seq1, seq2, ind):
    seq1[ind].append('N') 
    seq2[ind].append('N')
    #print(str(seq1[ind]))
    
def allSet(seq1, seq2, v = "N"):
    for k in seq1.keys():
        seq1[k].append(v)
        seq2[k].append(v)
        #print(str(seq1[k]))
     
def makeWindowAmbig(seq, size, p = True):
    for k in seq.keys():
        if (len(seq[k]) > 0 and p):
            seq[k].pop()#removes last one because of double INDEL lines
        if len(seq[k]) < size:
            seq[k] = ["N"]*len(seq[k])
        else:
            for i in range(-1*size, 0):
                seq[k][i] = "N"
     
def makeInts(ls):
    for i in range(0, len(ls)):
        ls[i] = int(ls[i])
   
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
   
def getTotalDepth(inds, dIndex):
    tot = 0
    inds = inds.split("\t")
    for i in inds:
        if i == "./.":
            continue
        tot += int(i.split(":")[dIndex])
    return tot
   
use = "usage: vcfToFasta.py <INPUT> <VCF> [-q <MIN QUALITY> -d <MIN DEPTH> -D <MAX DEPTH> -L <likelihood cutoff> -w <WINDOW SIZE> -o <DIRECTORY> -s -R -G -l <LOG FILE NAME>]\n"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("<INPUT> - an annotation file or a list of sites")
    print ("<VCF> - a vcf file")
    print ("q - the minimum quality to include a site in the analysis (40)")
    print ("d - the minimum depth to include a site in the analysis (20)")
    print ("D - the maximum depth to include a site in the analysis (60)")
    print ("L - the minimum likelihood to call a site as homozygous, or heterozygous rather than N (60)")
    print ("w - the number of sites on either side of INDELs to call ambiguous (5)\n\t -s Mode does NOT ignore windows around INDELs, it assumes you have removed INDEL windows from your site file already.")
    print ("o - the output directory for all the fasta files (\"\")")
    print ("s - flag that indicates we are outputting the fasta information for the a list of sites not an annotation")
    print ("R - flag that indicates we are outputting the fasta information for the outgroup")
    print ("G - flag that indicates we are reading data from a GATK generated vcf")
    print ("l - the name of the log file (vcfToFasta.log)")
    
if __name__ == "__main__":   
    __main__()
