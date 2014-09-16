import sys
import re
import getopt
import annotation
import pickle

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    f = "first_introns.sites"
    n = "nonfirst_introns.sites"
    
    try: 
        opts, args = getopt.getopt(sys.argv[3:],"f:n:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-f":
            f = arg
        elif opt == "-n":
            n = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
    
    first = open(f, "w")
    nonfirst = open(n, "w")
    annot = annotation.Reader(open(sys.argv[1], "rb"))
      
    divergence = pickle.load(open(sys.argv[2], "rb"))
        
    print("Gene\tIntronNum\tIntronLen\tSitesTot\tDiverged")
    for gene in annot:
        #print gene_name intron_number intron_length site_count divergence_count
        #sys.stderr.write(gene)
        introns = gene.makeIntrons()
        for j, intron in enumerate(introns):
            #sys.stderr.write(intron)
            siteCount = 0
            divCount = 0
            for i in range(intron[0], intron[1]+1):
                div = divergence[i][2]#get the flag indicating divergence
                if div == None:
                    sys.stderr.write("NO DATA: "+str(i)+"\n")
                    continue
                siteCount += 1
                if not div:
                    divCount += 1
                if j == 0:#a first intron
                    first.write(gene.scaf+"\t"+str(i)+"\n")
                else:
                    nonfirst.write(gene.scaf+"\t"+str(i)+"\n")
                    
            print(gene.name+"\t"+str(j)+"\t"+str(intron[1]-intron[0]+1)+"\t"+str(siteCount)+"\t"+str(divCount))
    

   
use = "python "+__file__.split("/")[-1]+" ANNOTATION PICKLED_DIVERGENCE [-f <FIRSTINTRON_SITES> -n <NONFIRST_SITES>]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()