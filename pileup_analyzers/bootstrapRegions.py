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

def __main__():
    #check aruguments
    processArgs()
    
    summaryReader = summary.Reader(open(_args.summary,"rb"))
    
    #dictionary of site codes
    #types = {'0fold': 3, 'stop': 8, 'intergene': 0, '5utr': 6, '3utr': 5, 'exon': 2, 'intron': 1, 'istop': 7, '4fold': 4, 'unknown': 9, 'cnc':10}
    #reverse_types = {v:k for k, v in types.items()}
        
    #Read in the file and build the regions of each size. Calculate the local AFS for each region.
    if _args.sample_size % 2 == 0:
        folded_size = int(_args.sample_size/2+1)
    else:
        folded_size = int(ceil(float(_args.sample_size)/2))
        
    regions = {}
    i = -1
    currscaf = None
    afs = resetAFS({}, summaryReader.summary.typeToCode, folded_size)
    div = resetAFS({}, summaryReader.summary.typeToCode, 2)#just makes a dictionary for each type with a list of length 2, store num sites and num div sites
    added = 0
    for site in summaryReader:#"#CHROM\tPOS\tREF\tALT\tREF_NUMBER\tALT_NUMBER\tTOTAL\tSITE_TYPE")
        if currscaf == None:
            currscaf = site.CHROM
        
        #if i == -1:#start the window at the first site we have
        #    i = site.POS
        
        if added >= _args.window or currscaf != site.CHROM:
        #if site.POS >= i+_args.window or currscaf != site.CHROM:
            regions[currscaf+"_"+str(site.POS)] = (afs, div)
            afs = resetAFS({}, summaryReader.summary.typeToCode, folded_size)
            div = resetAFS({}, summaryReader.summary.typeToCode, 2)
         #   i += _args.window
            added = 0
            if currscaf != site.CHROM:
                currscaf = site.CHROM
        
        if site.TOTAL != _args.sample_size or site.TOTAL != site.REF_NUM+site.ALT_NUM:#not enough alleles at this locus, or total is wrong (filtered data)
            sys.stderr.write("Skipping site, either not enough alleles or the total number of alleles does not match the sum of alleles.\n\t%s\n"%site)
            continue
        
        af = min(site.REF_NUM, site.ALT_NUM)
       
        #if in SNP mode only add to the window size if we are not at a fixed site
        if af != 0 and af != 1:
            added += 1
        elif not _args.snps:
            added += 1
 
        for t in site.Types:
            if t not in afs.keys():
                sys.stderr.write("Weird site type encountered: \n\t"+str(site)+"\n")
                continue
                
            afs[t][af] += 1
            
            #store divergence info for the window
            if site.DIVERGENCE != -1:
                div[t][0] += 1
            
            if site.DIVERGENCE == 1:
                div[t][1] += 1
    
    #this is commented out because the last region is usually not long enough to be useful
    #regions[currscaf+"_"+str(i)] = (afs, div)#store the last region
    
    #make sure all wanted types are actually present
    safe = []
    for w in _args.site_types:
        if w in summaryReader.summary.Types.keys():
            safe.append(w)
        else:
            sys.stderr.write("Wanted type not present in codes: %s. Excluding it from analysis.\n" % w)
    _args.site_types = safe
    
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
    for i in range(_args.bootstraps):
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
        for type in _args.site_types:
            total_afs[type] = [sum(values) for values in zip(*[new_region[0][type], total_afs[type]])]
            total_div[type] = [sum(values) for values in zip(*[new_region[1][type], total_div[type]])]
    
    #calculate the info for the neutral sites
    #assumes the first site type in _args.site_types is the neutral category
    neut = _args.site_types[0]
    neut_div = str(total_div[neut][0])+"\t"+str(total_div[neut][1])
    afs = total_afs[neut]
    myStr = ""
    for x in afs:
        myStr += str(x)+"\t"
    neut_afs = myStr+"0\t"*int(_args.sample_size/2)
    
    #print out info for each type
    for type in _args.site_types[1:]:
        #calculate the selected info
        afs = total_afs[type]
        div = total_div[type]
        sel_div = str(div[0])+"\t"+str(div[1])
        sel_afs = ""
        for x in afs:
            sel_afs += str(x)+"\t"
        sel_afs += "0\t"*int(_args.sample_size/2)
            
        #print the data
        print("%s_%s" % (reverse_types[type], name))
        print(sel_div)
        print(neut_div)
        print(str(_args.sample_size))
        print(sel_afs)
        print(neut_afs)  
        print("")  
        
def resetAFS(dict, types, size):
    for key in types.keys():
        dict[types[key]] = [0]*size
    return dict        

import argparse
_args = None
def processArgs():
    parser = argparse.ArgumentParser(description="takes a summary file and bootstraps the AFS by either bp or snp windows. Outputs DFE formmatted files.")
    parser.add_argument("summary", help="the summary file to read data from")
    parser.add_argument("sample_size", type = int, help="the sample size, in number of alleles. Sites without this size are excluded from the AFS.")
    parser.add_argument("-w", "--window", default = 10000, type = int, help="the window size")
    parser.add_argument("-b", "--bootstraps", default= 200, type = int, help="the number of bootstraps to perform.")
    parser.add_argument("-s", "--site-types", type=str, help="a comma separated list of the site types you want to get AFS for. The first one should be your NEUTRAL sites.")
    parser.add_argument("-p", "--snps", action="store_true", help="If this option is turned on then windows are defined by number of SNPs not bp.")
    global _args
    _args = parser.parse_args()
    sys.stderr.write(str(_args)+"\n")
    
if __name__ == "__main__":   
    __main__()
