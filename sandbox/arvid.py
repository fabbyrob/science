import sys
import re
import getopt

def __main__():
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad vcf name: "+sys.argv[1])
        sys.exit()

    #reg ex
    pat = re.compile(" (\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")

    counts = {}
    pTE = ""
    
    for i in range(0,10):
        counts[i] = []

    for line in infile:
        line= line.rstrip()
        m = pat.match(line)
        
        if (m == None):
            print ("NO!")
            continue
        else:
        
            TE = m.group(1)
            
            if (pTE != "" and pTE != TE):
                for i in range(0,10):
                    myMax = max(counts[i])
                    print(pTE+"\t"+str(i)+"\t"+str(myMax))
                    counts[i] = []
            
            for i in range(0,10):
                counts[i].append(m.group(i+2))
            pTE = TE
          
    for i in range(0,10):
        myMax = max(counts[i])
        print(pTE+"\t"+str(i)+"\t"+str(myMax))
        counts[i] = []      
    
if __name__ == "__main__":   
    __main__()      