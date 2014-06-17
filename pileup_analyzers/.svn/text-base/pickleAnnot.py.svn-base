import sys
import re
import getopt
import pickle

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-l":
            global _log 
            _log = str(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    annot = open(sys.argv[1],"r")
    
    if (annot == None):
        print("Bad vcf name: "+sys.argv[1])
        sys.exit()
        
    annotation = []
        
    name = ""
    start = -1
    end = -1
    direction = ""
    scaffold = ""
    exons = []
    pGene = None
    #read the annotation.
    for line in annot:
        line = line.rstrip()
        line = line.split()
        
        if len(line) != 5:
            sys.stderr.write("Incorrectly formatted line: "+str(line)+"\n")
            continue
        
        (line_scaf, line_start, line_end, line_name, line_dir) = line
        
        if name != line_name:
            #process previous
            if name != "":
                if pGene:
                    annotation.append(pGene)
                pGene = (start, {'start':start, 'end':end, 'exons':exons[:], 'scaffold':scaffold, 'name':name, 'direction':direction})
            #store new
            name = line_name
            start = int(line_start)
            end = int(line_end)
            exons = [(start, end)]
            direction = line_dir
            scaffold = line_scaf  
        else:
            #update end
            end = int(line_end)
            #add new exon
            exons.append((int(line_start), int(line_end)))
        
    #get the last two, if applicable  
    if pGene:  
        annotation.append(pGene)
    annotation.append((start, {'start':start, 'end':end, 'exons':exons[:], 'scaffold':scaffold, 'name':name, 'direction':direction}))
        
    annotation.sort()#sort the annotation by base start position
       
    #remove genes that overlap (assumes only 2 genes can overlap in a row) 
    i = 0
    while i < len(annotation)-2:
        g1 = annotation[i][1]
        g2 = annotation[i+1][1]
        if g1['end'] >= g2['start']:
            sys.stderr.write("Overlapping genes, excluding from annotation:\n\t"+g1['name']+"\t"+g2['name']+"\n")
            annotation.pop(i)
            annotation.pop(i)
        i += 1
            
    pickle.dump(annotation, open(sys.argv[2], "wb"))
   
use = "python "+__file__.split("/")[-1]+" AnnotationFile OutFileName"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()