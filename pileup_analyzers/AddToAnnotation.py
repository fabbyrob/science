'''
'Takes in an annotation, a list of sites, and a site type, and modifies the annotation so all sites in the site list now have the new indicated code.'

If the Annotation file includes a line specifying the site types at the top of the file, those codes are used rather than the defaults.

Used to add UTRs and CNCs to annotated references. 

python AddToAnnotation.py myAnnotation.txt mySites.txt -t 3utr

----- myAnnotation.txt -----
scaffold_1    3756819    C    PAC:20888905    0    1
scaffold_1    3756820    A    PAC:20888905    0    1
scaffold_1    3756821    C    PAC:20888905    0    1
scaffold_1    3756822    A    PAC:20888905    0    1
scaffold_1    3756823    A    PAC:20888905    0    1
scaffold_1    3756824    G    PAC:20888905    0    1
scaffold_1    3756825    G    PAC:20888905    0    1
scaffold_1    3756826    T    PAC:20888905    0    1
scaffold_1    3756827    G    PAC:20888905    0    1

----- mySites.txt -----
scaffold_1 3756820
scaffold_1 3756821
scaffold_1 3756822

----- output -----
scaffold_1    3756819    C    PAC:20888905    0    1
scaffold_1    3756820    A    PAC:20888905    0    1,5
scaffold_1    3756821    C    PAC:20888905    0    1,5
scaffold_1    3756822    A    PAC:20888905    0    1,5
scaffold_1    3756823    A    PAC:20888905    0    1
scaffold_1    3756824    G    PAC:20888905    0    1
scaffold_1    3756825    G    PAC:20888905    0    1
scaffold_1    3756826    T    PAC:20888905    0    1
scaffold_1    3756827    G    PAC:20888905    0    1

'''

import sys
import getopt

_t = '3utr'
_r = False
codes = {'intergene':0 , 'intron':1, 'exon':2, '0fold':3, '4fold':4, '3utr':5, '5utr':6, 'istop':7, 'stop':8, 'unknown':9, 'cnc':10}

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    addCode(_t)
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
        
    sites = open(sys.argv[2],"r")  
        
    if (sites == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()

    next_site = getNextSite(sites)
    codesPrinted = False
    #read the infile file...
    for line in infile:
        line = line.rstrip()
        sline = line.split()
        
        if line[0] == "#":
            if codesPrinted:
                sys.stderr.write("Got to a header line beginning with '#' after printing the codes, ignoring the line. If your codes look wrong make sure the first line in your file specifies the code dictionary.\n")
                continue
            global codes
            codes = strToDict(line[1:])
            addCode(_t)
            sys.stderr.write("Using codes from file: %s\n" % codes)
            print("#%s" % codes)
            codesPrinted = True
            continue
        
        if not codesPrinted:
            print("#%s" % codes)
            codesPrinted = True
        
        if next_site == (None, None):#print the remaining sites
            print(line)
            continue
        
        if len(sline) != 6:
            continue
        
        scaf, pos, ref, gene, dir, type = sline
        pos = int(pos)
        
        while next_site[0] != None and scaf > next_site[0]:#we missed this site
            sys.stderr.write("Missed a site (CHROM): "+str(next_site)+"\n")
            next_site = getNextSite(sites)
            
        if next_site == (None, None):
            sys.stderr.write("Ran out of sites\n")
            print(line)
            continue
        
        if scaf < next_site[0]:#we havent found the right scaf yet, go to the next site. 
            print(line)
            continue
        
        #on the right scaf, now deal with pos info
        
        while next_site[1] != None and pos > next_site[1]:#we missed a site WARNING, this while loop could bring us to the next chromosome....
            sys.stderr.write("Missed a site (POS): "+str(next_site)+"\n")
            next_site = getNextSite(sites)
            
        if next_site[1] > pos:#haven't reached the site yet
            print(line)
            continue
        
        if next_site == (None, None):
            sys.stderr.write("Ran out of sites\n")
            print(line)
            continue
        
        if next_site[1] == pos:#at the site
            if type == str(codes['intergene']) and _r:
                new_code = str(codes[_t])
            else:
                new_code = str(type)+","+str(codes[_t])
                
            print(scaf+"\t"+str(pos)+"\t"+ref+"\t"+gene+"\t"+dir+"\t"+new_code)#print out the site with the new code
            next_site = getNextSite(sites)

def strToDict(aStr):
    newStr = aStr[1:-1]
    ls = newStr.split(", ")
    myDict = {}
    for entry in ls:
        key, val = entry.split(": ") 
        key = key[1:-1]
        val = int(val)
        myDict[key] = val
    return myDict         

def addCode(type):
    global codes
    if type in codes.keys():
        return codes
    
    codes[type] = max(codes.values())+1
    
    return codes
      
def getNextSite(infile):
    line = infile.readline()
    if line == "":#EOF
        return None, None
    
    sline = line.rstrip().split()
    return sline[0], int(sline[1])
        

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"t:r")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-t":
            global _t 
            _t = arg
        elif opt == "-r":
            global _r 
            _r = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+" AnnotationFile SitesFile [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print('Takes in an annotation, a list of sites, and a site type, and modifies the annotation so all sites in the site list now have the new indicated code.')
    print()
    print("option - argument type - default - description")
    print("t - STR - "+str(_t)+" - the default types of sites include: "+str(codes.keys())+"\nHowever, arbitrary site types can be added, a new code will be generated for any unrecognized site type.")
    print("r - NONE - "+str(_r)+" - If this flag is used then any site with a code of 0 (noncoding) is replaced with the new code. Otherwise, the code is appended")

    
if __name__ == "__main__":   
    __main__()