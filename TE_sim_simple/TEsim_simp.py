from pop_simp import Population
import sys
import getopt
import time
#from collections import Counter

_N = 50 #pop size
_T = [10] #number of TEs per ind per class
_r = 0 #selfing rate
_S = 'R'#selection mode
_l = [0]#starting silencing rate for each class
_G = 200#number of loci per chromosome
_u = 0.1#transposition rate
_x = 0.001#selection coeff for TEs
_z = 0#selction ratefor silencing
_g = 'ALL'#fitness calc mode one of 'ALL' 'HOMO' or 'HET'
_t = 1 #selection exponent
_X = 200#number of crossovers per ind
_s = 0#mutation rate in silencing
_K = 100# number of generations to run
_h = False#help flag
_b = 0#number of classes
_c = 0#probability of class mutation
_e = 0#mutation magnitude for silencing
_v = 0#excision rate
_V = False#verbose mode for debugging
_I = None#number of generations to wait to introduce the silencing allele, if this is used silencing is never mutated
_R = False#use relative fitness to check for viability

def __main__():
    try: 
        opts, args = getopt.getopt(sys.argv[1:],"hT:N:G:u:x:g:t:X:r:K:v:e:l:s:z:c:S:I:R:V")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-h":
            details()
            sys.exit(0)
        elif opt == "-T":
            global _T
            _T = map(float,arg.split(","))
        elif opt == "-N":
            global _N
            _N = int(arg)
        elif opt == "-G":
            global _G
            _G = int(arg)
        elif opt == "-u":
            global _u
            _u = float(arg)
        elif opt == "-v":
            global _v
            _v = float(arg)
        elif opt == "-x":
            global _x
            _x = float(arg)
        elif opt == "-t":
            global _t
            _t = float(arg)
        elif opt == "-X":
            global _X
            _X = int(arg)
        elif opt == "-r":
            global _r
            _r = float(arg)
        elif opt == "-K":
            global _K
            _K = int(arg)
        elif opt == "-e":
            global _e
            _e = float(arg)
        elif opt == "-l":
            global _l
            _l = map(float,arg.split(","))
        elif opt == "-s":
            global _s
            _s = float(arg)
        elif opt == "-z":
            global _z
            _z = float(arg)
        elif opt == "-c":
            global _c
            _c = float(arg)
        elif opt == "-V":
            global _V
            _V = True
        elif opt == "-I":
            global _I
            _I = int(arg)
        elif opt == "-R":
            global _R
            _R = True
        elif opt == "-S":
            global _S
            _S = arg
            if _S not in ["R","O"]:
                sys.stderr.write("invalid selection type %s.\n" % _S)
                sys.exit(0)
        elif opt == "-g":
            global _g
            if arg not in ['ALL', 'HOMO', 'HET']:
                print("Invlaid selection calculation type: -g "+arg)
                sys.exit(0)
            _g = arg
        else:
            sys.stderr.write("Unrecognized option: "+opt+"\n")
            usage()
    
    if len(_T) != len(_l):
        sys.stderr.write("Silencing and starting TE values not the same length. \n\t-T %s -l %s\nExiting.\n" % (_T, _l))
        sys.exit(0)
    
    _b = len(_T)
    start = time.time()
    printParams()
        
    #initialize population    
    pop = Population(u = _u, v = _v, crosses = _X, loci = _G, startTEs = _T, size = _N, selection = _x, t = _t, r = _r, mode = _g, classes = _b, silenceU = _s, silenceMag = _e, silStart= _l, selectionSilence = _z, classU = _c, verbose = _V, selType = _S, silIntro = _I, relative = _R)
     
    #phase 1 initial population stabilization   
    for gen in range(1, _K):
        result = pop.generation(gen)
        if not result:
            break
        
    if _I:
        print("FREQ\t%s\t%s" % (gen, pop.freqs[0]))
        
    sys.stderr.write("Run time (min): "+str((time.time()-start)/60.0)+"\n")

  
def prettyPrintList(ls):
    myStr = ""
    for i in ls:
        myStr += str(i)+","
    return myStr
  
def printParams():
  mystr = "-N "+str(_N)+" -T "+str(_T)+" -r "+str(_r)+" -S "+str(_S)+" -l "+str(_l)+" -G "+str(_G)+" -u "+str(_u)+" -x %.9f"%_x+" -g "+str(_g)+" -t "+str(_t)+" -X "+str(_X)+" -s "+str(_s)+" -K "+str(_K)+" -b "+str(_b)+" -c "+str(_c)+" -e %.9f"%_e+" -z %.9f"%_z+" -v %.9f"%_v+" -I %s"%_I + "\n"
  print(mystr)
  
use = "TEsim.py [options]" 
def usage():
    print (use)
    sys.exit(0)
    
def details():
    print (use)
 
if __name__ == "__main__":   
    __main__()
