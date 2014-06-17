import sys
import re
import getopt
import time
from population import Population
#from collections import Counter

_N = 50 #pop size
_T = 10 #number of TEs per ind
_r = 0 #selfing rate
_S = 'T'#selection mode
_l = [0.1]#starting silencing rate for each class
_G = 10#number of loci per chromosome
_u = 0.1#transposition rate
_x = 0.001#selection coeff for TEs
_z = 0.001#selction ratefor silencing
_g = 'ALL'#fitness calc mode one of 'ALL' 'HOMO' or 'HET'
_t = 1 #selection exponent
_X = 2#number of crossovers per ind
_s = 0#mutation rate in silencing
_K = 100# number of generations to run
_h = False#help flag
_b = 1#max number of classes
_c = 0.01#probability of class mutation
_e = 0.01#mutation magnitude for silencing
_v = 0.01#excision rate

def __main__():
    try: 
        opts, args = getopt.getopt(sys.argv[1:],"hT:l:N:G:u:x:g:t:X:s:r:K:S:b:c:e:z:v:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-h":
            details()
            sys.exit(0)
        elif opt == "-T":
            global _T
            _T = int(arg)
        elif opt == "-l":
            global _l
            _l = map(float, arg.split(","))
        elif opt == "-N":
            global _N
            _N = int(arg)
        elif opt == "-G":
            global _G
            _G = int(arg)
        elif opt == "-u":
            global _u
            _u = float(arg)
        elif opt == "-x":
            global _x
            _x = float(arg)
        elif opt == "-t":
            global _t
            _t = float(arg)
        elif opt == "-X":
            global _X
            _X = int(arg)
        elif opt == "-s":
            global _s
            _s = float(arg)
        elif opt == "-r":
            global _r
            _r = float(arg)
        elif opt == "-K":
            global _K
            _K = int(arg)
        elif opt == "-c":
            global _c
            _c = float(arg)
        elif opt == "-b":
            global _b
            _b = int(arg)
        elif opt == "-e":
            global _e
            _e = float(arg)
        elif opt == "-z":
            global _z
            _z = float(arg)
        elif opt == "-v":
            global _v
            _v = float(arg)
        elif opt == "-g":
            global _g
            if arg not in ['ALL', 'HOMO', 'HET']:
                print("Invlaid selection calculation type: -g "+arg)
                sys.exit(0)
            _g = arg
        elif opt == "-S":
            global _S
            if arg not in ['R', 'T', 'O', 'P']:
                print("Invlaid selection type: -S "+arg)
                sys.exit(0)
            _S = arg
        else:
            sys.stderr.write("Unrecognized option: "+opt+"\n")
            usage()
    
    start = time.time()
    printParams()
        
    #initialize population    
    pop = Population(size = _N, selfing = _r, silence = _l, numCrosses = _X, selection = (_x, _z), t = _t, mode = _g, transRate = _u, maxClasses = _b, classMutation = _c, silenceMutation = _s, epsilon = _e, excision = _v)
    pop.firstGen(loci = _G, numTEs = _T)
     
    #phase 1 initial population stabilization   
    for gen in range(1, _K):
        pop.generation(_S, gen)
        
    sys.stderr.write("Run time (min): "+str((time.time()-start)/60.0)+"\n")

  
def prettyPrintList(ls):
    myStr = ""
    for i in ls:
        myStr += str(i)+","
    return myStr
  
def printParams():
  mystr = "-N "+str(_N)+" -T "+str(_T)+" -r "+str(_r)+" -S "+str(_S)+" -l "+str(_l)+" -G "+str(_G)+" -u "+str(_u)+" -x %.9f"%_x+" -g "+str(_g)+" -t "+str(_t)+" -X "+str(_X)+" -s "+str(_s)+" -K "+str(_K)+" -b "+str(_b)+" -c "+str(_c)+" -e %.9f"%_e+" -z %.9f"%_z+" -v %.9f"%_v+"\n"
  print(mystr)
  
use = "TEsim.py [options]" 
def usage():
    print (use)
    sys.exit(0)
    
def details():
    print (use)
    print ("T - initial number of TEs per genome (10)\
    \nl - initial silencing rate, list separated by commas [i.e. 0.1,0.2,0.3], all classes not specified will start at 0 ([0.1])\
    \nN - pop size (50)\
    \nG - number of loci per chromosome (10)\
    \nu - transposition rate (0.1)\
    \nx - selection coefficient for TEs (0.001)\
    \nz - selection coefficient for silencing (0.001)\
    \ng - how to calculate g for fitness (\'HOMO\', \'HET\', \'ALL\') ('ALL')\
    \nt - selection thingy in fitness function (1)\
    \nX - mean number of crossovers per individual (2)\
    \ns - mutation rate in silencing (0)\
    \nr - selfing rate (0)\
    \nK - number of generations to run (100)\
    \nc - mutation rate between TE classes (0.01)\
    \nb - maximum number of TE classes to allow (1)\
    \ne - mutation magnitude of a mutation in silencing (0.01)\
    \nv - excision rate (0.01)\
    \nS - selection type (\'T\' - tournament, \'R\' - random, \'O\' - roulette, \'P\' - proportional)(\'T\')")
    
if __name__ == "__main__":   
    __main__()
