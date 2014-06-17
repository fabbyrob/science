# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 12:24:25 2013

@author: wiliarj
"""

doc = """
Takes two lists of ranges and pulls out any from the first list that overlap a range in the second list at a minimum of sites.

Assumes the files are sorted!
"""

cases = '''-------------=============------------- # = is the range in list 2
case 1 -----------------------********- # * is the range in list 1
2 --******-----------------------------
3 -------***********-------------------
4 --------------------**********-------
5 --------------********--------------- '''

import sys
import getopt

_o = 0.5
_v = False

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    file1 = open(sys.argv[1],"r")
        
    file2 = open(sys.argv[2],"r")  
    
    safeRanges = {}# {scaffold:[(start, end, name)]}
    
    if _v: 
        sys.stderr.write("The different possibe overlaps are:\n%s\n\nNumbers in parentheses below indicate which type of overlap was being tested, and whether or not the test failed.\n\n" % cases)    
    
    for line in file2:
        line = line.rstrip()
        sline = line.split()
        
        scaf = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        name = sline[3]

        if safeRanges.has_key(scaf):
            safeRanges[scaf].append((start, end, name))
        else:
            safeRanges[scaf] = [(start, end, name)]

    #read the infile file...
    for line in file1:
        line = line.rstrip()
        sline = line.split()
        
        scaf = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        name = sline[3]
        
        if not safeRanges.has_key(scaf):
            continue
        
        #-------------=============------------- # = is the range in list 2
        #case 1 -----------------------********- # * is the range in list 1
        #2 --******-----------------------------
        #3 -------***********-------------------
        #4 --------------------**********-------
        #5 --------------********---------------        
             
        for safeStart, safeEnd, safeName in safeRanges[scaf]:
            if end < safeStart:#case 1
                if _v:              
                    sys.stderr.write("No overlap (1): %s\t%s\n" % (name, safeName))
                break
            
            if start > safeEnd:#case 2
                if _v:
                    sys.stderr.write("No overlap (2): %s\t%s\n" % (name, safeName))
                continue
            
            if safeStart > start and end > safeStart:#case 3
                overlap = end - safeStart
                percentOverlap = float(overlap)/(end-start)
                if percentOverlap >= _o:
                    print("%s\t%s\t%s\t%s" % (scaf, start, end, safeName))
                    if _v:
                        sys.stderr.write("Accepted Region (3): %s\t%s\t%s\n" % (name, safeName, percentOverlap))
                    break
                else:
                    if _v:
                        sys.stderr.write("Insufficient overlap (3): %s\t%s\t%s\n" % (name, safeName, percentOverlap))
            
            if end > safeEnd and safeEnd > start:#case 4
                overlap = safeEnd - start
                percentOverlap = float(overlap)/(end-start)
                if percentOverlap >= _o:
                    print("%s\t%s\t%s\t%s" % (scaf, start, end, safeName))
                    if _v:
                        sys.stderr.write("Accepted Region (4): %s\t%s\t%s\n" % (name, safeName, percentOverlap))
                    break
                else:
                    if _v:
                        sys.stderr.write("Insufficient overlap (4): %s\t%s\t%s\n" % (name, safeName, percentOverlap))
            
            if safeStart >= start and end <= safeEnd:#case 5
                percentOverlap = 1
                if percentOverlap >= _o:
                    print("%s\t%s\t%s\t%s" % (scaf, start, end, safeName))
                    if _v:
                        sys.stderr.write("Accepted Region (5): %s\t%s\t%s\n" % (name, safeName, percentOverlap))
                    break
            else:
                if _v:
                    sys.stderr.write("Insufficient overlap (5): %s\t%s\t%s\n" % (name, safeName, percentOverlap))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"o:v")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-o":
            global _o 
            _o = float(arg)
        elif opt == "-v":
            global _v 
            _v = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("RangeList1: %s RangeList2 %s -o %s -v %s\n" % (sys.argv[1], sys.argv[2], _o, _v))
   
use = "python "+__file__.split("/")[-1]+" RangeList1.txt RangeList2.txt [OPIONS]"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("o - FLOAT - %s - the amount of overlap in ranges to include a range." % _o)
    print("v - NONE - %s - thi flag turns on verbose mode, the result of every overap test is printed to stderr." % _o)

    
if __name__ == "__main__":   
    __main__()