import sys
import re
import getopt
import vcf

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    reader = vcf.Reader(open(sys.argv[1],'rb'))

    print (reader.metadata)
    
    #read the vcf file...
    for record in reader:

        for sample in record.samples:
            print(sample['DP'][0])

   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)

    
if __name__ == "__main__":   
    __main__()