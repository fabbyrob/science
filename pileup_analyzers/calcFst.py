import argparse

def parseArgs():
    parser = argparse.ArgumentParser(description="takes two summary files and calculates Fst between the populations in windows across the genome. uses Hudson's estimator of Fst, see Bhatia et al 2013 Genome Res. As suggested by Bhatia the Fst of individual loci in windows is combined using an 'ratio of averages' where both the top half and bottom half of the Fst are averaged across loci then divided to get the estimate for the window.")
    parser.add_argument("summary1", type=str, help="the first summary file")
    parser.add_argument("summary2", type=str, help="the second summary file")
    parser.add_argument("n1", type=int, help="number of alleles in population 1")
    parser.add_argument("n2", type=int, help="number of alleles in population 2")

    parser.add_argument("-w","--window_size", default = 25, type=int, help="the window size in snps to report Fst in")
    parser.add_argument("-t","--test", action="store_true", help="this flag turns on testing of this module.")
    parser.add_argument("-v", "--verbose", action="store_true", help="this option turns on verbose error reporting")
    return parser.parse_args()

args = parseArgs()

def main():
    return

def calcFst(window1, window2):
    '''
    Given a list of allele frequencies for each population this function calculates the two components of Fst for each SNP.
    The result should be used as an argument to avgFst to get the whole-window average.
    >>> args.n1, args.n2 = 10, 10

    >>> calcFst([], [])
    []

    >>> calcFst([-1, 0, -1], [0, 0, 0])
    []

    >>> calcFst([2,3,0], [7,0,1])
    [(0.20888888888888882, 0.62), (0.06666666666666667, 0.3), (0.0, 0.1)]
 
    >>> args.n2 = 15
    >>> calcFst([2,3,0], [7,0,1])
    [(0.03555555555555555, 0.48000000000000004), (0.06666666666666667, 0.3), (0.0, 0.06666666666666667)]
    '''
    tops = []
    bottoms = []
    for p1, p2 in zip(window1, window2):
        if p1+p2 == 0 or p1+p2 == args.n1+args.n2 or -1 in (p1,p1):
            #not a SNP or unknown data
            continue
        p1 = float(p1)/args.n1
        p2 = float(p2)/args.n2

        t = (p1-p2)**2 - (p1*(1-p1))/(args.n1-1) - (p2*(1-p2))/(args.n2-1)
        b = p1*(1-p2) + p2*(1-p1)
        
        tops.append(t)
        bottoms.append(b)

    return zip(tops, bottoms)

def avgFst(window_fst):
    '''
    Return the ratio of averages of Fst across all loci in this window.

    >>> window = [(0.1,0.2),(0.5,0.4),(0.1,0.1),(0.01,0.2),(0.19,0.2)]
    >>> avgFst(window[0:3])
    0.9999999999999999

    >>> avgFst(window)
    0.818181818181818

    >>> avgFst([(0,1), (0.0,0.4)])
    0.0
    '''
    top = 0.0
    bottom = 0.0
    for t,b in window_fst:
        top += t
        bottom += b
    #note I dont divide by the length because I'd just do it to both top and bottom
    #likely introducing some rounding error and leading to the same result    

    return top/bottom

if __name__=="__main__":
    if args.test:
        import doctest
        doctest.testmod(verbose=args.verbose)
    else:
        main

