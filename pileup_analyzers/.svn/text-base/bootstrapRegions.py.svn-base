'''
This program takes in a vcf summary file with site types, and bootstraps region s of a given size to build 
AFS for each type of site. 

python bootstrapRegions.py vcfSummary.txt [OPTIONS]

run program with no arguments for a detailed description of each option. 
'''

import sys
import getopt
from math import ceil
from random import choice
import summary

_w = 10000#size of regions to bootstrap
_b = 200#number of bootstraps to perform
_N = 26#sample size
_wanted = ['4', '3', '1', '0', '10']#site types to output AFS for

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()

    processArgs(2)
    
    summaryReader = summary.Reader(open(sys.argv[1],"rb"))
    
    #dictionary of site codes
    #types = {'0fold': 3, 'stop': 8, 'intergene': 0, '5utr': 6, '3utr': 5, 'exon': 2, 'intron': 1, 'istop': 7, '4fold': 4, 'unknown': 9, 'cnc':10}
    #reverse_types = {v:k for k, v in types.items()}
        
    #Read in the file and build the regions of each size. Calculate the local AFS for each region.
    if _N % 2 == 0:
        folded_size = int(_N/2+1)
    else:
        folded_size = int(ceil(float(_N)/2))
        
    regions = {}
    i = -1
    currscaf = None
    afs = resetAFS({}, summaryReader.summary.typeToCode, folded_size)
    div = resetAFS({}, summaryReader.summary.typeToCode, 2)#just makes a dictionary for each type with a list of length 2, store num sites and num div sites
    for site in summaryReader:#"#CHROM\tPOS\tREF\tALT\tREF_NUMBER\tALT_NUMBER\tTOTAL\tSITE_TYPE")
        if currscaf == None:
            currscaf = site.CHROM
        
        if i == -1:#start the window at the first site we have
            i = site.POS
        
        if site.POS >= i+_w or currscaf != site.CHROM:
            regions[currscaf+"_"+str(i)] = (afs, div)
            afs = resetAFS({}, summaryReader.summary.typeToCode, folded_size)
            div = resetAFS({}, summaryReader.summary.typeToCode, 2)
            i += _w
            if currscaf != site.CHROM:
                i = site.POS
                currscaf = site.CHROM
        
        if site.TOTAL != _N or site.TOTAL != site.REF_NUM+site.ALT_NUM:#not enough alleles at this locus, or total is wrong (filtered data)
            sys.stderr.write("Skipping site, either not enough alleles or the total number of alleles does not match the sum of alleles.\n\t%s\n"%site)
            continue
        
        af = min(site.REF_NUM, site.ALT_NUM)
        
        for t in site.Types:
            if t not in afs.keys():
                sys.stderr.write("Weird site type encountered: \n\t"+line+"\n")
                continue
                
            afs[t][af] += 1
            
            #store divergence info for the window
            if site.DIVERGENCE != -1:
                div[t][0] += 1
            
            if site.DIVERGENCE == 1:
                div[t][1] += 1
    
    regions[currscaf+"_"+str(i)] = (afs, div)#store the last region
    
    #make sure all wanted types are actually present
    global _wanted
    safe = []
    for w in _wanted:
        if w in summaryReader.summary.Types.keys():
            safe.append(w)
        else:
            sys.stderr.write("Wanted type not present in codes: %s. Excluding it from analysis.\n" % w)
    _wanted = safe
    
    counts = resetAFS({}, summaryReader.summary.typeToCode, 1)
    for r in regions.values():
        for t in summaryReader.summary.Types.keys():
            counts[t][0] += sum(r[0][t])
            
    sys.stderr.write("Numbers of each site type in whole analysis:\n")
    for t, val in summaryReader.summary.typeToCode.items():
        sys.stderr.write("%s %s\n" % (t, counts[val][0]) )
        
    
    region_names = list(regions.keys())
    num_regions = len(region_names)
    
    #print out the non-bootstrapped afs
    output_bootstrap(summaryReader.summary.typeToCode, summaryReader.summary.Types, folded_size, region_names, regions, "real")

    #perform the bootstraps
    for i in range(_b):
        boot_regions = []
        while len(boot_regions) < num_regions:
           boot_regions.append(choice(region_names))
        
        #output the AFS for this bootstrap
        output_bootstrap(summaryReader.summary.typeToCode, summaryReader.summary.Types, folded_size, boot_regions, regions, str(i))

#takes in all the regions, and a list of region names (allows duplicates) to calculate an AFS from
#prints that AFS out in the format for dfealpha.pl
def output_bootstrap(types, reverse_types, folded_size, region_names, regions, name="real"):
    #calculate the afs and div for each type for this set of regions
    total_afs = resetAFS({}, types, folded_size)
    total_div = resetAFS({}, types, 2)
    for region in region_names:
        new_region = regions[region]
        for type in _wanted:
            total_afs[type] = [sum(values) for values in zip(*[new_region[0][type], total_afs[type]])]
            total_div[type] = [sum(values) for values in zip(*[new_region[1][type], total_div[type]])]
    
    #calculate the info for the neutral sites
    #assumes the first site type in _wanted is the neutral category
    neut = _wanted[0]
    neut_div = str(total_div[neut][0])+"\t"+str(total_div[neut][1])
    afs = total_afs[neut]
    myStr = ""
    for x in afs:
        myStr += str(x)+"\t"
    neut_afs = myStr+"0\t"*int(_N/2)
    
    #print out info for each type
    for type in _wanted[1:]:
        #calculate the selected info
        afs = total_afs[type]
        div = total_div[type]
        sel_div = str(div[0])+"\t"+str(div[1])
        sel_afs = ""
        for x in afs:
            sel_afs += str(x)+"\t"
        sel_afs += "0\t"*int(_N/2)
            
        #print the data
        print("%s_%s" % (reverse_types[type], name))
        print(sel_div)
        print(neut_div)
        print(str(_N))
        print(sel_afs)
        print(neut_afs)  
        print("")  
        
def resetAFS(dict, types, size):
    for key in types.keys():
        dict[types[key]] = [0]*size
    return dict        

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"w:b:N:s:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-w":
            global _w 
            _w = int(arg)
        elif opt == "-b":
            global _b 
            _b = int(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-s":
            global _wanted 
            _wanted = arg.split(",")
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
            
    sys.stderr.write("infile: %s -w %s -b %s -N %s -s %s\n" % (sys.argv[1], _w, _b, _N, _wanted))
   
use = "python "+__file__.split("/")[-1]+" vcfSummary [OPTIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("This program takes in a vcf summary file with site types, and bootstraps region s of a given size to build AFS for each type of site.")
    print()
    print("option - argument type - default - description")
    print("w - INT - "+str(_w)+" - the size of the regions to be bootstrapped in bases")
    print("b - INT - "+str(_b)+" - the number of bootstraps to perform")
    print("N - INT - "+str(_N)+" - the sample size a site must have to be included in analysis")
    print("s - STR (comma separated INTs) - "+str(",".join(map(str, _wanted)))+" - the site categories that should be output after bootstrapping - ASSUMES THE FIRST IS YOUR NEUTRAL CATEGORY")

    
if __name__ == "__main__":   
    __main__()