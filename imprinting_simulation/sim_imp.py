from pop_imp import Population
import getopt
import sys

_G = 200
_N = 500
_r = 0.5

def __main__():
    processArgs(1)
    
    pop = Population(_N, _r)
    
    g = 0
    while g < _G:
        pop.generation(g)
        g += 1
    return

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"G:N:r:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-G":
            global _G 
            _G = int(arg)
        elif opt == "-N":
            global _N 
            _N = int(arg)
        elif opt == "-r":
            global _r 
            _r = float(arg)
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

if __name__ == "__main__":   
    __main__()