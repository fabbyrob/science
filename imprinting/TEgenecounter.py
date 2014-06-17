#!/usr/bin/python
import sys
import re
import getopt

_log = __file__.split("/")[-1]+".log"#log file name

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"l:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    tefile = open(sys.argv[2],"r")  
        
    if (tefile == None):
        print("Bad TE file name: "+sys.argv[1])
        sys.exit()
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    TEs = ['DNA',
            'DNA/En-Spm',
            'DNA/Harbinger',
            'DNA/hAT',
            'DNA/hAT-Ac',
            'DNA/hAT-Tag1',
            'DNA/hAT-Tip100',
            'DNA/MuDR',
            'DNA/TcMar-Pogo',
            'DNA/TcMar-Stowaway',
            'LINE/L1',    
            'LTR',
            'LTR/Copia',
            'LTR/ERVK',
            'LTR/Gypsy',
            'RC/Helitron',
            'SINE']
        
    dict = {}
    for t in TEs:
        dict[t] = 0    
    
    mystr = "gene"
    for k in dict.keys():
        mystr += ","+k
    print(mystr)
    
    #reg ex
    annot_pat = re.compile("\w+_(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(.*)")
    te_pat = re.compile("\s+\d+\s+\S+\s+\S+\s+\S+\s+\w+_(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+).*")
    
    modP = None
    scafP = None
    startP = None
    endP = None
    TEP = None
    #read the vcf file...
    for line in infile:
        line = line.rstrip()
        m = annot_pat.match(line)
        
        if (m == None):# if we find a weird line...
            logfile.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            scaf = int(m.group(1))+1
            start = int(m.group(2))
            end = int(m.group(3))
            mod = m.group(4)
            
            if mod != modP:
                if modP != None:
                    #print (scafP, modP, startP, endP)
                    dict = {}
                    for t in TEs:
                        dict[t] = 0
                    #find all the TEs
                    while TEP == None or (TEP[0] == scafP and (TEP[1] < endP)):
                        if TEP != None and (TEP[1] > startP or TEP[2] > startP):
                            if TEP[3] in dict.keys():
                                dict[TEP[3]] += 1
                            else:
                                logfile.write("Unknown TE class: "+str(TEP))
                        TEP = getNextTE(tefile, te_pat)
                        if TEP == None:
                            break
                    
                    myStr = str(modP)
                    
                    print(str(modP)+","+ str(list(dict.values()))[1:-1]) 
                
                #reset main variables
                modP = mod
                scafP = scaf
                startP = start
                endP = end
            else:
                #not at the end yet...
                endP = end
            

def getNextTE(file, pat):
    teClass = None
    while (teClass == None):
        aLine = file.readline()
        
        if (aLine == ""):
            return
        
        m = pat.match(aLine)
        
        if(m == None):
            continue

        
        scaf = int(m.group(1))
        start = int(m.group(2))
        end = int(m.group(3))
        teClass = m.group(4)
        
        #if it is an unidentified or non-TE class toss it
        if (teClass in ['Low_complexity','Simple_repeat', 'rRNA', 'snRNA', 'Satellite', 'Other/Composite', 'Unknown', 'buffer']):
            teClass = None
            
    #print (aLine)
    return (scaf, start, end, teClass)

use = "python "+__file__.split("/")[-1]+"<Annotation> <TE locations>"

def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("<Annotation> - the annotation file for gene locations")
    print ("<TE locations> - file listing TE locations and types")

    
if __name__ == "__main__":   
    __main__()
