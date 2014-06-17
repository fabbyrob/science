doc = """
Much like AddToAnnotation.py adds site a given site type to the summary at specified sites.
"""

import sys
import getopt
import summary

_t = '3utr'
_r = False
_d = False

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    summaryReader = summary.Reader(open(sys.argv[1],"rb"))
    
    if _t not in summaryReader.Types:
        intCodes = map(int, summaryReader.Codes)
        summaryReader.addCode(max(intCodes)+1, _t)
        
    print(summaryReader.Header)
    
    sites = open(sys.argv[2],"r")  
        
    if (sites == None):
        print("Bad sites file name: "+sys.argv[1])
        sys.exit()

    next_site = getNextSite(sites)
    #read the infile file...
    for site in summaryReader:
        if next_site == (None, None):#print the remaining sites
            if _d and summaryReader.TypeToCode[_t] in site.Types:
                site.Types.remove(summaryReader.TypeToCode[_t])
                site.Types.append(0)
                
            print(site)
            continue
    
        while next_site[0] != None and site.CHROM > next_site[0]:#we missed this site
            sys.stderr.write("Missed a site (CHROM): "+str(next_site)+"\n")
            next_site = getNextSite(sites)
            
        if next_site == (None, None):
            sys.stderr.write("Ran out of sites\n")
            if _d and summaryReader.TypeToCode[_t] in site.Types:
                site.Types.remove(summaryReader.TypeToCode[_t])
                site.Types.append(0)
            print(site)
            continue
        
        if site.CHROM < next_site[0]:#we havent found the right scaf yet, go to the next site. 
            if _d and summaryReader.TypeToCode[_t] in site.Types:
                site.Types.remove(summaryReader.TypeToCode[_t])
                site.Types.append(0)
            print(site)
            continue
        
        #on the right scaf, now deal with pos info
        
        while next_site[1] != None and site.POS > next_site[1]:#we missed a site WARNING, this while loop could bring us to the next chromosome....
            sys.stderr.write("Missed a site (POS): "+str(next_site)+"\n")
            next_site = getNextSite(sites)
            
        if next_site[1] > site.POS:#haven't reached the site yet
            if _d and summaryReader.TypeToCode[_t] in site.Types:
                site.Types.remove(summaryReader.TypeToCode[_t])
                site.Types.append(0)
            print(site)
            continue
        
        if next_site == (None, None):
            sys.stderr.write("Ran out of sites\n")
            if _d and summaryReader.TypeToCode[_t] in site.Types:
                site.Types.remove(summaryReader.TypeToCode[_t])
                site.Types.append(0)
            print(site)
            continue
        
        if next_site[1] == site.POS:#at the site
            if _r:
                for r in _r:
                    if summaryReader.TypeToCode[r] in site.Types:
                        site.Types.remove(summaryReader.TypeToCode[r])
                        
            site.Types.append(summaryReader.TypeToCode[_t])
                
            print(site)
            next_site = getNextSite(sites)
      
def getNextSite(infile):
    line = infile.readline()
    if line == "":#EOF
        return None, None
    
    sline = line.rstrip().split()
    return sline[0], int(sline[1])

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"t:r:d")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-t":
            global _t 
            _t = arg
        elif opt == "-r":
            global _r 
            if _r:
                _r.append(arg)
            else:
                _r = [arg]
        elif opt == "-d":
            global _d 
            _d = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
            
    sys.stderr.write("python "+__file__.split("/")[-1]+" %s %s -t %s -r %s -d %s\n" % (sys.argv[1], sys.argv[2], _t, _r, _d))
   
use = "python "+__file__.split("/")[-1]+" AnnotationFile SitesFile [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print('Takes in an annotation, a list of sites, and a site type, and modifies the annotation so all sites in the site list now have the new indicated code.')
    print()
    print("option - argument type - default - description")
    print("t - STR - "+str(_t)+" - the type of site being added. Can be a type already annotated in tehile header or an arbitrary site types can be added, a new code will be generated for any unrecognized site type.")
    print("r - INT - "+str(_r)+" - If this option is used then any site with the given code is replaced with the new code. Otherwise, the code is simply appended to the list. This option can be used multiple times, if you want to replaceseveral different categories.")
    print("d - NONE - "+str(_r)+" - If this flag is used then any site with the given code given in -t that is not inluded on the input list has that code removed from its type list.")

   
if __name__ == "__main__":   
    __main__()