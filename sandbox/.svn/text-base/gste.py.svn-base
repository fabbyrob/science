import sys
import re
import getopt

_log = __file__.split("/")[-1]+".log"#log file name
_w = 50
_b = 1000

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 7:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[7:],"b:w:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w
            _w = int(arg)
        elif opt == "-b":
            global _b
            _b = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    ortho = open(sys.argv[1],"r")
    
    if (ortho == None):
        print("Bad vcf name: "+sys.argv[1])
        sys.exit()
        
    teposR = open(sys.argv[2],"r")  
        
    if (teposR == None):
        print("Bad sites file name: "+sys.argv[2])
        sys.exit()

    teposL = open(sys.argv[3],"r")  
        
    if (teposL == None):
        print("Bad sites file name: "+sys.argv[3])
        sys.exit()
        
    teposT = open(sys.argv[4],"r")  
        
    if (teposT == None):
        print("Bad sites file name: "+sys.argv[4])
        sys.exit()

    orthLength = open(sys.argv[5],"w")  
        
    if (orthLength == None):
        print("Bad sites file name: "+sys.argv[5])
        sys.exit()
        
    teCounts = open(sys.argv[6],"w")  
        
    if (teCounts == None):
        print("Bad sites file name: "+sys.argv[6])
        sys.exit()

     
    #get the length and midpoint of each section   
    rWindows = []
    lWindows = []
    tWindows = []

    numRead = 0
    starts = []
    prev = []
    end = 0
    orthLength.write("midpoint\trubella\tlyrata\tthaliana\n")
    ortho.readline()
    for line in ortho:
        line = line.rstrip()
        vals = line.split(",")
        
        fail = False
        for i in range(4,9):
            if vals[i] == "":
                #sys.stderr.write("Missing values in line:\n\t"+str(line)+"\n")
                fail = True
        
        if fail:
            continue
        
        rStart = int(vals[1])
        rEnd = int(vals[2])
        
        lStart = int(vals[4])
        lEnd = int(vals[5])
        
        tStart = int(vals[7])
        tEnd = int(vals[8])
        foundBreak = False
        
        
        if not starts:
            starts = [rStart, lStart, tStart]
        
        
        foundBreak = False
        if (prev) and (lStart-prev[1] > _b or tStart-prev[2] > _b):
#            print("break! "+str(starts)+str(rStart))
            foundBreak = True
        
        
        if numRead == _w or foundBreak:
            if foundBreak:
#                print("BREAK!"+str(starts)+str(rStart))
                midpoint = (starts[0]+prev[0])/2
                rLength = prev[0] - starts[0]
                lLength = prev[1] - starts[1]
                tLength = prev[2] - starts[2]
                
                rWindows.append((starts[0], prev[0], midpoint))
                lWindows.append((starts[1], prev[1], midpoint))
                tWindows.append((starts[2], prev[2], midpoint))
                
                starts = [rStart, lStart, tStart]
            else:
                midpoint = (starts[0]+rEnd)/2
                rLength = rEnd - starts[0]
                lLength = lEnd - starts[1]
                tLength = tEnd - starts[2]
            
                rWindows.append((starts[0], rEnd, midpoint))
                lWindows.append((starts[1], lEnd, midpoint))
                tWindows.append((starts[2], tEnd, midpoint))
                
                starts = []
            
#            if rLength < 0 or lLength < 0 or tLength < 0:
#                print (str(starts)+"\t"+str(prev)+"\t"+str(rEnd)+" "+str(rStart))
            orthLength.write(str(midpoint)+"\t"+str(rLength)+"\t"+str(lLength)+"\t"+str(tLength)+"\n")
            
            numRead = 0
            
        prev = [rEnd, lEnd, tEnd]
        numRead += 1
        
    if starts:
        midpoint = (starts[0]+prev[0])/2
        rLength = prev[0] - starts[0]
        lLength = prev[1] - starts[1]
        tLength = prev[2] - starts[2]
        
        rWindows.append((starts[0], prev[0], midpoint))
        lWindows.append((starts[1], prev[1], midpoint))
        tWindows.append((starts[2], prev[2], midpoint))
            
        orthLength.write(str(midpoint)+"\t"+str(rLength)+"\t"+str(lLength)+"\t"+str(tLength)+"\n")
     
     #get the number of TEs in each window   
    #print("DONE WITH LENGTHS")
    orthLength.close()
    families = ['DNA/Helitron','LINE?','LTR/ERVK','DNA/TcMar-Stowaway','DNA/TcMar-Mariner','DNA/TcMar-Pogo','DNA/hAT-Tip100','DNA/hAT-Tag1','DNA/hAT-Ac','DNA/Harbinger','DNA/MuDR','DNA','DNA/En-Spm','DNA/hAT','LINE','LINE/L1','LTR','LTR/Copia','LTR/Gypsy','LTR/DIRS', 'LTR/ERV1','RC/Helitron','SINE','SINE?', 'DNA?']

    rubTEs, rubFams = countTEs(rWindows, teposR, families)
    lyrTEs, lyrFams = countTEs(lWindows, teposL, families)
    thTEs, thFams = countTEs(tWindows, teposT, families)
    
    #print("TES COUNTED")
    i = 0 
    myStr = "midpoint\trubella"
    for f in families:
        myStr += "\trub."+f
    myStr += "\tlyrata"
    for f in families:
        myStr += "\tlyr."+f
    myStr += "\tthaliana"
    for f in families:
        myStr += "\tth."+f
        
    myStr += "\n"
    teCounts.write(myStr)
    for i in range(0, len(rubTEs)):
        myStr = str(rubTEs[i][0])+"\t"+str(rubTEs[i][1])
        for f in families:
            myStr += "\t"+str(rubFams[i][f])
            
        myStr += "\t"+str(lyrTEs[i][1])
        for f in families:
            myStr += "\t"+str(lyrFams[i][f])
            
        myStr += "\t"+str(thTEs[i][1])
        for f in families:
            myStr += "\t"+str(thFams[i][f])
             
        myStr += '\n'
        teCounts.write(myStr)#+"\t"+str(lyrTEs[i][1])+"\t"+str(thTEs[i][1])+"\n")
    
    teCounts.close()

def countTEs(windows, file, families):
    lineStart = -1
    lineEnd = -1
    counts = []
    fcounts = []
    pline = ""
    for win in windows:
        print(win)
        famCounts = {}
        for f in families:
            famCounts[f] = 0
            
        #print(win)
        start = win[0]
        end = win[1]
        count = 0
        
        if start <= lineEnd and start >= lineStart:
            #print("\t"+str(pline))
            count += 1
            famCounts[family] += 1
        
        while lineEnd <= end:
            line = file.readline()
            if line == "":
                break
            pline = line
            if len(line) == 0:
                continue
            line = line.split(",")
            lineStart = int(line[1])
            lineEnd = int(line[2])
            family = line[-1].split("#")[-1].rstrip()
            
            if lineStart >= start and lineStart <= end:
                #print("\t"+str(line))
                famCounts[family] += 1
                count += 1
                
        #print(famCounts) 
        counts.append((win[2], count))
        fcounts.append(famCounts)
        
    return (counts, fcounts)
            
        
        
    return
   
use = "python "+__file__.split("/")[-1]+" OrthologFile TEpositionFile -w <WINDOW SIZE> -b <BREAK SIZE>"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()