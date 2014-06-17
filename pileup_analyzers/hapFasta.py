doc = """

"""

import sys
import getopt
import vcf

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(3)
    
    reader = vcf.Reader(open(sys.argv[1],'rb'))

    cols = {"A":"r","C":"b","G":"g","T":"y","bold":"_","black":"black"}
    weird = False
    weirds = []

    fastas = {}
    for s in reader.samples:
        fastas[s] = [[],[]]

    i = 0
    #read the infile file...
    for record in reader:
        genos = [record.REF, record.ALT]
        i += 1
        for samp in record.samples:
            if "AD" in samp.keys():
                fastas[samp['name']][0].append(genos[int(samp['GT'][0])][0])
                fastas[samp['name']][1].append(genos[int(samp['GT'][2])][0])
                
                if samp["AD"] == ['.'] or (('0' in samp["GT"]) and (samp["AD"][0] == 0)) or (('1' in samp["GT"]) and (samp["AD"][1] == 0)):
                    #fastas[samp['name']][0].append(cols["bold"])
                    #fastas[samp['name']][1].append(cols["bold"])
                    print("%s %s %s" %(samp['name'], i, record.POS))
#                 elif record.ALT != ".":
#                     fastas[samp['name']][0].append(cols[genos[int(samp['GT'][0])][0]])
#                     fastas[samp['name']][1].append(cols[genos[int(samp['GT'][2])][0]])
                    
            else:
                fastas[samp['name']][0].append("N")
                fastas[samp['name']][1].append("N")
                
    for s in fastas.keys():
        print(">%s"%s)
        print("".join(fastas[s][0]))
        print(">%s"%s)
        print("".join(fastas[s][1]))
                

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

    sys.stderr.write("infile: %s -arg %s\n" % (sys.argv[1], 1))
   
use = "python "+__file__.split("/")[-1]+""
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