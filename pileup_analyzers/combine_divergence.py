import sys
import getopt
import cPickle as pickle

div = {}

i = 1
for p in sys.argv[2:]:
    name = "scaffold_"+str(i)
    div[name] = pickle.load(open(p, "rb"))
    i +=1
    
#repickle as dict
pickle.dump(div, open(sys.argv[1], "wb"))