import sys
import getopt
import random

_o = ""
_matches = {}
_I = 95
_r = False
_t = 1

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 2:
        usage()
    
    processArgs(2)
    
    infile = open(sys.argv[1],"r")
    
    if (infile == None):
        print("Bad infile name: "+sys.argv[1])
        sys.exit()
    
    
    cols = ["alpha","omega","nes_1","nes1_10", "nes10_100", "nes_100"]
    
    samples, reals = processDFE(infile)
      
    if not _o:  
        sys.stderr.write("Calculating confidence intervals.\n")    
        print("type\talpha_min\talpha\talpha_max\tomega_min\tomega\tomega_max\tnes_1_min\tnes_1\tnes_1_max\tnes1_10_min\tnes1_10\tnes1_10_max\tnes10_100_min\tnes10_100\tnes10_100_max\tnes100_min\tnes_100\tnes100_max")
        for type in samples.keys():
            data = samples[type]
            myStr = type
            for i, stat in enumerate(cols):
                vals = [ls[i] for ls in data]
                ci = calculateCI(vals)
                if type in reals.keys():
                    myStr += "\t%s\t%s\t%s" % (ci[0], reals[type][i],ci[1])
                else:
                    sys.stderr.write("No point estimate for type %s. Make sure you name the point estimate as %s_real.\n" % (type, type))
                    myStr += "\t%s\t%s\t%s" % (ci[0], "NA",ci[1])
                    
            print(myStr)
        
    if _o:
        sys.stderr.write("Comparing two files to get p-values.\n")
        comparefile = open(_o)
        samples2, reals2 = processDFE(comparefile)
        
        print("type\talpha\tomega\tnes_1\tnes1_10\tnes10_100\tnes100")
        for type in samples.keys():
            s_data = samples[type]
            random.shuffle(s_data)
            
            otherType = type
            if type in _matches.keys():
                otherType = _matches[type]
                
            if otherType not in samples2.keys():
                sys.stderr.write("Compare file has no outputs of type %s.\n" % type)
                continue
            
            s2_data = samples2[otherType]
            random.shuffle(s2_data)
            
            if len(s_data) != len(s2_data):
                sys.stderr.write("Number of samples from two files not the same for type %s:\n\t%s - %s\n\t%s - %s\n" % (type, sys.argv[1], len(s_data), _o, len(s2_data)))
                sys.stderr.write("Using the smaller number.\n")
                num = min([len(s_data), len(s2_data)])
                s_data = s_data[0:num]
                s2_data = s2_data[0:num]
                
            paired_data = zip(s_data, s2_data)
            myStr = type
            for i, stat in enumerate(cols):
                greaters = [1 if ls[i] < ls2[i] else 0  for ls, ls2 in paired_data]
                p = _t*float(sum(greaters))/len(greaters)
                if _t == 2:
                    p = min(p, 2-p)
                myStr += "\t%s" % p
            print(myStr)

def processDFE(infile):
    samples = {}
    reals = {}

    #read the infile file...
    for line in infile:
        line = line.rstrip()
        sline = line.split()
        if len(sline) < 22:
            sys.stderr.write("weird line:\n\t %s" % line)
            continue
        type = sline[0].split("_")[0]
        #alpha, omega, nes < 1, nes 1-10, nes 10-100, nes 100+
        try:
            data = [float(sline[-6]),float(sline[-5]),float(sline[-4]),float(sline[-3]),float(sline[-2]),float(sline[-1])]
        except:
            sys.stderr.write("some none number (probably NA) values in line:\n\t %s\n" % line)
            continue
        
        if "real" in sline[0]:
            reals[type] = data 
        else:
            if type in samples.keys():
                samples[type].append(data)
            else:
                samples[type] = [data]
            
    for t in reals.keys():
        if t not in reals.keys():
            reals[t] = ["NA"]*len(cols)
            
    return samples, reals
    

def calculateCI(vals):
    throw = ((100.0-_I)/200)*len(vals)

    sys.stderr.write("throwing %s values from top and bottom.\n" % int(throw))
    vals.sort()
    
    return((vals[int(throw)], vals[-1*int(throw+1)]))

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"o:I:m:t:r")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-o":
            global _o 
            _o = arg
        elif opt == "-I":
            global _I 
            _I = float(arg)
        elif opt == "-t":
            global _t 
            if int(arg) not in [1,2]:
                sys.stderr.write("Invalid input to -t, needs to be either 1 for a single tailed test, or 2 for a two tailed test.\n")
                sys.exit()
            _t = int(arg)
        elif opt == "-m":
            global _matches
            info = arg.split(":")
            if len(info) != 2:
                sys.stderr.write("Malformed -m input %s, skipping.\n" % arg)
                continue
            _matches[info[0]] = info[1]
        elif opt == "-r":
            global _r 
            _r = True
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()
   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    print("")
    print("This program takes in an output file from dfealpha.pl, it will generate confidence intervals for each category of sites in the file.\n\
    If another file is provided with the -o option, then samples will be randomly paired, and the proportion that have a higher value in the -o file for each statistic\n\
    (4 Nes categories, alpha, and omega) will be reported for each category.\n")
    print("The categories should have a distinct name, individual samples can be labeled with an underscore followed by identifying information. \
    \nExample:\n\
intron_0        -1767601.2058   660     NA      540.77  NA      587.0000        NA      100.0000        -0.002999       0.782347656     NA      NA      NA      NA      NA      0.903239        1.256290        0.000004        0.999996        0.000000        0.000000\n\
intergene_0     -2008617.2470   660     NA      535.02  NA      559.0000        NA      100.0000        -0.003818       0.771542157     NA      NA      NA      NA      NA      0.952723        1.493557        0.000000        1.000000        0.000000        0.000000\n\
intron_1   -376738.4379    660     NA      510.75  NA      460.0000        NA      0.9403  -0.094031       0.786260336     NA      NA      NA      NA      NA      -nan    -nan    0.025125        0.176258        0.670813        0.127804\n")
    print("option - type - default - description")
    print("o - STR - NONE - the alternate file to test for significance")
    print("m - STR - NONE - for use with -o, if you want to compare types among files that are not the same name, use this to specify them in the form: -m file1Name:file2Name.\n\tSo if in my mail file the sites were called '5utr' and in the -o file they are called '5utrcns' I would use: -m 5utr:5utrcns.\nYou can use this option multiple times, once for each pair of site types.")
    print("I - FLOAT - %s - the size of the confidence interval to be calculated" % _I)

    
if __name__ == "__main__":   
    __main__()