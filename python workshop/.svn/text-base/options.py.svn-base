import sys
import getopt

print (sys.argv)

_d = 40
_D = 60
_f = False

try:
    opts, args = getopt.getopt(sys.argv[1:],"d:D:f")
except getopt.GetoptError:
    print ("Error")
    
print (opts)

print(_d, _D, _f)

for opt, arg in opts:
    if opt == "-d":
        _d = int(arg)
    elif opt == "-D":
        _D = int(arg)
    elif opt == "-f":
        _f = True
    else:
        print ("Unknown option entered")
        
print(_d, _D, _f)