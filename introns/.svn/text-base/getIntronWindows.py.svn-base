doc = """
Reads in a VCF summary with gene names and outputs AFS and CNC data from windows along every intron.
Windows are numbered from 1 on the 5' to the 3' end.
The window immediately in the middle of the intron may be thrown out if it is smaller than the window size.

python mySummary.txt -w 5 -N 6

#window size 5
GENE DIR IntronLen IntronNum WindowNum CNCNum afs_0 afs_1 afs_2 afs_3
ATG1 - 100 1 1 2 4 0 0 1
ATG1 - 100 1 2 1 3 1 1 0
ATG1 - 100 1 3 0 5 0 0 0
ATG1 - 100 1 4 2 2 1 1 1

the above example file contains info from 1 intron, the first intron in gene ATG1.
There are 4 5 bp windows in this intron.

"""

import sys
sys.path.append("/data/robert.williamson/bin")

import getopt
import summary
from math import floor, ceil

_w = 5
_N = 26
_intron = "intron"
_cnc = "allCNC"
_v = False #verbose mode

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    reader = summary.Reader(open(sys.argv[1], 'rb'))
    
    if not "GENE" in reader.summary.Fields:
        sys.stderr.write("Your summary does not have gene names. Please add gene info to this summary and try again.\n")
        sys.exit()
    
    myHeader = "GENE DIR IntronNum IntronLen Window CNCnum "
    for i in range(int(floor(_N/2.0)+1)):
        myHeader += "afs_%s "%i
    
    print(myHeader)
    
    currentGene = ""
    dir = ""
    introns = []
    currentIntronAF = []
    currentIntronCNC = []
    if _v: sys.stderr.write("Reading in Summary...\n")
    for record in reader:
        if record.GENE != "N":#we are in a gene
            if record.GENE != currentGene and currentGene != "N":#finish off the previous gene
                if currentIntronAF: introns.append(zip(currentIntronAF, currentIntronCNC))
                processGene(currentGene, dir, introns)
                currentGene = record.GENE
                dir = record.DIR
                introns = []
                currentIntronAF = []
                currentIntronCNC = []
            elif not currentGene:#start this new gene
                currentGene = record.GENE
                dir = record.DIR
                introns = []
                currentIntronAF = []
                currentIntronCNC = []
            
            if reader.TypeToCode[_intron] in record.Types:#in an intron
                if record.TOTAL == _N:
                    currentIntronAF.append(min(record.REF_NUM, record.ALT_NUM))
                else:#no data
                    currentIntronAF.append(-1)
                    
                if reader.TypeToCode[_cnc] in record.Types:#also a CNC site
                    currentIntronCNC.append(1)
                else:#not a CNC
                    currentIntronCNC.append(0)
            
            else:#not in an intron
                if currentIntronAF:#we have passed the previous intron
                    introns.append(zip(currentIntronAF, currentIntronCNC))
                    currentIntronAF = []
                    currentIntronCNC = []
                
        else:#we arent in a gene
            if currentGene and currentGene != "N":
                if currentIntronAF: introns.append(zip(currentIntronAF, currentIntronCNC))
                processGene(currentGene, dir, introns)
                currentGene = ""
                dir = ""
                introns = []
                currentIntronAF = []
                currentIntronCNC = []
                
    if introns:#get the last set of introns, if any
        if currentIntronAF: introns.append(zip(currentIntronAF, currentIntronCNC))
        processGene(currentGene, dir, introns)

def processGene(gene, dir, introns):
    if _v: sys.stderr.write("Processing gene data - %s...\n" % gene)
    #number introns
    num = 1
    if dir == "-":
        num = len(introns)
    
    for intron in introns:
        intronLen = len(intron)
        toss = intronLen % _w #the left over sites in the middle of the window that are gonna be excluded
        tossSite = int(ceil(floor(intronLen/_w)/2))#where the sites are tossed from
        
        windowNum = 0
        #get all the windows up to the toss point
        for i in range(0, tossSite, _w):
            window = intron[i:i+_w]
            processWindow(gene, dir, num, intronLen, windowNum, window)
            windowNum += 1
            
        #get all the windows after the tossed sites
        for i in range(tossSite+toss, intronLen, _w):
            window = intron[i:i+_w]
            processWindow(gene, dir, num, intronLen, windowNum, window)
            windowNum += 1
            
        if dir == "+":
            num += 1
        else:
            num -= 1
            
def processWindow(gene, dir, intronNum, intronLen, windowNum, window):
    if _v: sys.stderr.write("Processing introns...\n")
    cncNum = 0
    myAFS = [0]*int(floor(_N/2.0)+1)
    #build AFS
    #count CNCs
    for af, cnc in window:
        cncNum += cnc
        if af != -1: myAFS[af] += 1
    
    #output line
    myStr = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene, dir, intronNum, intronLen, windowNum, cncNum, "\t".join(map(str, myAFS)))
    print(myStr)

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:N:i:c:v")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _o = int(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-i":
            global _intron 
            _intron = arg
        elif opt == "-c":
            global _cnc 
            _cnc = arg
        elif opt == "-v":
            global _v 
            _v = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("infile: %s -N %s -w %s -i %s -c %s -v %s\n" % (sys.argv[1], _N, _w, _intron, _cnc, _v))
   
use = "python "+__file__.split("/")[-1]+" Summary [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("w - INT - %s - the window size to split individual introns into." % _w)
    print("N - INT - %s - the sample size in number of alleles. Sites without this many total alleles will not be included in the output." % _N)
    print("i - STR - %s - the name used to identify introns in the summary file." % _w)
    print("c - STR - %s - the name used to identify CNCs in the summary file." % _w)
    print("v - None - %s - verbose mode." % _w)
    
if __name__ == "__main__":   
    __main__()