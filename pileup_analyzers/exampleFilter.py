import vcf
import sys
import argparse
    
#process the input arguments I want to use for my filter
#you can see how argpars works here: https://docs.python.org/2/howto/argparse.html
parser = argparse.ArgumentParser(description="This is my filter for VCFs it filters based on ...") #fill in the ...
parser.add_argument("-v", "--verbose", action="store_true", help = "this flag turns on verbose mode, so more info is printed to stderr")
############ add your own arguments here ############
#before you use these on the command line you need to add in an extra -- to separate the vcfSummarizer.py options from your filter's options
parser.add_argument("-m", "--mapqual", type = int, default = 15, help = "The minimum map quality to call a site")#this is an example, feel free to delete it

#####################################################
#this only grabs the arguments after the -- for the filter
if "--" in sys.argv:
    i = sys.argv.index("--")+1
else:
    i = 0
args = parser.parse_args(sys.argv[i:])
sys.stderr.write("Arguments for filter: %s\n" % args)#this is just so you always have a record of what you ran

#set up the filter function
#this function always takes exactly 1 argument a vcf.Record object, containing all the info from one VCF line
#this function MUST return exactly 4 arguments
#ref_num is the number of reference alleles at this site, after filtering
#    e.g. 5
#alt_num is the number of alternate alleles at this site, after filtering
#    e.g. 3
#total is the total number of alleles at this site, after filtering
#    e.g. 8
#genos is a list of all the genotype calls after filtering, in the order that samples appear in the VCF
#    e.g. ["H","N","A","R","R","N"]
#    in this list H is a heterozygote
#    R is a reference homozygote
#    A is a alternate homozygote
#    N is an unknown genotype (this can mean both no data or it did not pass filters)
def filter(record):
    #this version of the filter is really dumb, 
    #it just checks the MQ and outputs the genotypes in the GT column
    ref_num = 0
    alt_num = 0
    total = 0
    genos = []
    
    #you can raise custom (or default) exceptions at any time to stop this site from being output
    #NO data from sites with exceptions will be output not even Ns
    if record.POS == 192:
        raise BadSiteException("I don't like site %s. Skipping it." % (record.POS))
    
    if record.INFO["MQ"] < args.mapqual:
        if args.verbose:
            sys.stderr.write("Site at %s %s has low MQ.\n" % (record.CHROM, record.POS))
        return (0, 0, 0, ["N"]*len(record.samples))
    
    for sample in record.samples:
        if sample["GT"] == "./.":
            genos.append("N")
        elif sample["GT"] == "0/0":
            genos.append("R")
            ref_num += 2
            total += 2
        elif sample["GT"] == "0/1":
            genos.append("H")
            ref_num += 1
            alt_num += 1
            total += 2
        elif sample["GT"] == "1/1":
            genos.append("A")
            alt_num += 2
            total += 2
            
    return (ref_num, alt_num, total, genos)

#A custom exception for whatever you're doing
class BadSiteException(Exception):
    pass

#this following if statement will never be true if you are using vcfSummariser.py and giving it this script
#it will only be true if you directly run *this* script on the command line
#     e.g. >python emptyFilter.py
#so it is a good place to test your filter
#give it some example Records and the expected output and it will test for you
if __name__=="__main__":
    from collections import OrderedDict
    tests_passed = 0
    tests_done = 0
    
    #You should modify (and add more) tests like these 3 examples
    #I HIGHLY encourage you to write these BEFORE you write the filter
    #If you are testing on real data, and find a line that breaks your code
    #Add it as a test here as you are fixing your code, so that you will always know that it wont
    #re-break in the future
    
    ##### TEST 1  ##### 
    #tests that we count reference homozygotes correctly
    #info is a dictionary form of the INFO line from the VCF
    info = {'AC':[0], 'AF':[0.0], 'MQ0':7, 'MQ':19.58}
    #info is an OrderedDict of info from the genotypes columns - it is an OrderedDict so that the samples stay in the right order
    #i'm spreading this out a lot so it is easy to see what all the elements are
    genos = OrderedDict()
    #add these one at a time or the order will be wrong
    genos['sample_1'] = {'name':'sample_1', 'GT':'0/0', 'DP': 20, 'GQ':63, 'PL':(0.0, 21.0, 55.0)} #the keys in this dictionary are the 'name' and all the items in the FORMAT column of the VCF
    genos['sample_2'] = {'name':'sample_2', 'GT':'0/0', 'DP': 20, 'GQ':63, 'PL':(0.0, 21.0, 55.0)}
    genos['sample_3'] = {'name':'sample_3', 'GT':'0/0', 'DP': 20, 'GQ':63, 'PL':(0.0, 21.0, 55.0)}            
    
    record_1 = vcf._Record(CHROM="scaffold_1", POS=100, ID='.', REF='A', ALT='.', QUAL='.', FILTER='.', INFO=info, FORMAT='GT:DP:GQ:PL', genotypes=genos)
    result = filter(record_1)
    
    tests_done += 1
    expected = (6, 0, 6, ["R", "R", "R"])
    if result == expected:
        sys.stderr.write("Test 1 passed.\n")
        tests_passed += 1
    else:
        sys.stderr.write("Test 1 FAILED.\n\tExpected: %s\b\t Got: %s\n\n" % (expected, result))
        
    ##### TEST 2  ##### 
    #tests that we count heterozygotes correctly
    info = {'AC':[0], 'AF':[0.5], 'MQ0':7, 'MQ':19.58}
    genos = OrderedDict()
    genos['sample_1'] = {'name':'sample_1', 'GT':'0/1', 'DP': 20, 'GQ':63, 'PL':(20.0, 0.0, 55.0)} #the keys in this dictionary are the 'name' and all the items in the FORMAT column of the VCF
    genos['sample_2'] = {'name':'sample_2', 'GT':'0/1', 'DP': 20, 'GQ':63, 'PL':(20.0, 0.0, 55.0)}
    genos['sample_3'] = {'name':'sample_3', 'GT':'0/1', 'DP': 20, 'GQ':63, 'PL':(20.0, 0.0, 55.0)}            
    
    record_1 = vcf._Record(CHROM="scaffold_1", POS=100, ID='.', REF='A', ALT='T', QUAL='.', FILTER='.', INFO=info, FORMAT='GT:DP:GQ:PL', genotypes=genos)
    result = filter(record_1)
    
    tests_done += 1
    expected = (3, 3, 6, ["H", "H", "H"])
    if result == expected:
        sys.stderr.write("Test 2 passed.\n")
        tests_passed += 1
    else:
        sys.stderr.write("Test 2 FAILED.\n\tExpected: %s\b\t Got: %s\n\n" % (expected, result))
        
    ##### TEST 3  ##### 
    #tests that we count alternate homozygotes and no data correctly
    info = {'AC':[0], 'AF':[0.5], 'MQ0':7, 'MQ':19.58}
    genos = OrderedDict()
    genos['sample_1'] = {'name':'sample_1', 'GT':'0/0', 'DP': 20, 'GQ':63, 'PL':(20.0, 10.1, 55.0)} #the keys in this dictionary are the 'name' and all the items in the FORMAT column of the VCF
    genos['sample_2'] = {'name':'sample_2', 'GT':'./.', 'DP': 20, 'GQ':63, 'PL':(20.0, 50.0, 0.0)}
    genos['sample_3'] = {'name':'sample_3', 'GT':'1/1', 'DP': 20, 'GQ':63, 'PL':(20.0, 0.0, 55.0)}            
    
    record_1 = vcf._Record(CHROM="scaffold_1", POS=100, ID='.', REF='A', ALT='T', QUAL='.', FILTER='.', INFO=info, FORMAT='GT:DP:GQ:PL', genotypes=genos)
    result = filter(record_1)
    
    tests_done += 1
    expected = (2, 2, 4, ["R", "N", "A"])
    if result == expected:
        sys.stderr.write("Test 3 passed.\n")
        tests_passed += 1
    else:
        sys.stderr.write("Test 3 FAILED.\n\tExpected: %s\b\t Got: %s\n\n" % (expected, result))
        
    ##### TEST 3  ##### 
    #tests that we filter on MQ correctly
    info = {'AC':[0], 'AF':[0.5], 'MQ0':7, 'MQ':5.0}
    genos = OrderedDict()
    genos['sample_1'] = {'name':'sample_1', 'GT':'0/0', 'DP': 20, 'GQ':63, 'PL':(20.0, 10.1, 55.0)} #the keys in this dictionary are the 'name' and all the items in the FORMAT column of the VCF
    genos['sample_2'] = {'name':'sample_2', 'GT':'./.', 'DP': 20, 'GQ':63, 'PL':(20.0, 50.0, 0.0)}
    genos['sample_3'] = {'name':'sample_3', 'GT':'1/1', 'DP': 20, 'GQ':63, 'PL':(20.0, 0.0, 55.0)}            
    
    record_1 = vcf._Record(CHROM="scaffold_1", POS=100, ID='.', REF='A', ALT='T', QUAL='.', FILTER='.', INFO=info, FORMAT='GT:DP:GQ:PL', genotypes=genos)
    result = filter(record_1)
    
    tests_done += 1
    expected = (0, 0, 0, ["N", "N", "N"])
    if result == expected:
        sys.stderr.write("Test 4 passed.\n")
        tests_passed += 1
    else:
        sys.stderr.write("Test 4 FAILED.\n\tExpected: %s\b\t Got: %s\n\n" % (expected, result))
        
        
    ##### Summary  ##### 
    sys.stderr.write("Summary:\n\t%s/%s tests passed.\n" % (tests_passed, tests_done))
    if tests_passed != tests_done:
        sys.stderr.write("\tWARNING NOT ALL YOUR TESTS PASSED!\n")
    
    
    
    
    
    