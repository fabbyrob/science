import sys
import re
import getopt

_log = __file__.split("/")[-1]+".log"#log file name

_w = 10000
_s = 1000

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[2:],"l:w:s:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        elif opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-s":
            global _s 
            _s = int(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    logfile = open(_log,"w")
    
    if (logfile == None):
        print("Bad logfile name: "+_log)
        sys.exit()
        
    #reg ex
    vcf_pat = re.compile("(\S+)\s(\d+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)")
    
    myWindowStarts = []
    myWindowVals = [[],[],[],[]]
    currentWindowStart = 0
    currentWindowEnd = 0
    
    #read the infile file...
    for line in infile:
        line = line.rstrip()
        m = vcf_pat.match(line)
        
        if (m == None):# if we find a weird line...
            sys.stderr.write("Non-standard line (VCF):\n\t"+line+"\n")
            continue
        else:
            scaf = m.group(1)
            base = int(m.group(2))
            col1 = float(m.group(3))
            col2 = float(m.group(4))
            col3 = float(m.group(5))
            col4 = float(m.group(6))
            
            if currentWindowStart == 0:#first line
                currentWindowStart = base
                myWindowStarts.append(base)
                myWindowVals[0].append(col1)
                myWindowVals[1].append(col2)
                myWindowVals[2].append(col3)
                myWindowVals[3].append(col4)
                currentWindowEnd = base + _w
                continue
            
            
            
            
            while base > currentWindowEnd:
                #calc window
                result = processWindow(myWindowVals)
                #print results: scaffold, start, v1,v2,v3,v4
                print (scaf+"\t"+str(currentWindowStart)+"\t"+str(result[0])+"\t"+str(result[1])+"\t"+str(result[2])+"\t"+str(result[3]))
                #update start and end
                currentWindowStart += _s
                currentWindowEnd += _s
                #remove elements as necessary
                pops = 0
                for site in myWindowStarts:
                    if site < currentWindowStart:
                        pops += 1
                    else:
                        break
                    
                #remove any elements that are earlier than the new window
                while pops > 0:
                    myWindowStarts.pop(0)
                    myWindowVals[0].pop(0)
                    myWindowVals[1].pop(0)
                    myWindowVals[2].pop(0)
                    myWindowVals[3].pop(0)
                    pops -= 1

            #add last base to the window
            myWindowStarts.append(base)
            myWindowVals[0].append(col1)
            myWindowVals[1].append(col2)
            myWindowVals[2].append(col3)
            myWindowVals[3].append(col4)
            
            #finish the last windows
    while len(myWindowStarts) > 0:
        #calc window
        result = processWindow(myWindowVals)
        #print results: scaffold, start, v1,v2,v3,v4
        print (scaf+"\t"+str(currentWindowStart)+"\t"+str(result[0])+"\t"+str(result[1])+"\t"+str(result[2])+"\t"+str(result[3]))
        #update start and end
        currentWindowStart += _s
        currentWindowEnd += _s
        #remove elements as necessary
        pops = 0
        for site in myWindowStarts:
            if site < currentWindowStart:
                pops += 1
            else:
                break
            
        #remove any elements that are earlier than the new window
        while pops > 0:
            myWindowStarts.pop(0)
            myWindowVals[0].pop(0)
            myWindowVals[1].pop(0)
            myWindowVals[2].pop(0)
            myWindowVals[3].pop(0)
            pops -= 1
            
def mean(numberList):
    if len(numberList) == 0:
        return float('nan')
 
    floatNums = [float(x) for x in numberList]
    return sum(floatNums) / len(numberList)
        
def processWindow(win):
    data = []
    for i in win:
        data.append(mean(i))
    return data
   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print ("w - the width of the sliding window to analyze (10000)")
    print ("s - the step size the window will increment by (1000)")

    
if __name__ == "__main__":   
    __main__()