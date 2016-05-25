# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 13:14:09 2013

@author: wiliarj
"""

import sys

def __main__():
    if len(sys.argv) < 2:
        print("provide a sample size and a folded AFS:\n\tpython AFSstats.py 5 100 50 23")
        print("Outputs Taj D, pi, and watersons theta")
        sys.exit(0)
        
    N = int(sys.argv[1])
    if sys.argv[2] == "f":#reading many AFS from file
        print("Sample\ttheta\tpi\tTajD")
        for line in open(sys.argv[3]):
            sline = line.split()
            name = sline[0]
            afs = map(float, sline[1:])
            t, p, D = stats(afs,N)
            print("%s\t%s\t%s\t%s" % (name, t, p, D))
    
    else:
        afs = map(float, sys.argv[2:])
        t, p, D = stats(afs,N)
        print("theta: %s\npi: %s\nTajD: %s" % (t, p, D))

def calcParams(_N):
    params = {}
    a1 = 0
    a2 = 0
    for i in range(1,_N):
        a1 += 1.0/i
        a2 += 1.0/(i**2)
        
    params["a1"] = a1
    params["a2"] = a2
    
    b1 = float(_N+1)/(3*(_N-1))#see Tajima 1989
    b2 = float(2*(_N**2+_N+3))/(9*_N*(_N-1))

    params["b1"] = b1
    params["b2"] = b2

    c1 = b1 - 1/a1
    c2 = b2-float(_N+2)/(a1*_N)+a2/(a1**2)

    params["c1"] = c1
    params["c2"] = c2

    e1 = c1/a1
    e2 = c2/(a1**2+a2)
    
    params["e1"] = e1
    params["e2"] = e2
    return params

def stats(afs, N):
    params = calcParams(N)
    
    Ns = sum(afs)
    S = sum(afs[1:])
    theta = S/params["a1"]
    pi = 0 
    
    if S == 0:
        return 0, 0, 0   
 
    for i, numSites in enumerate(afs[1:]):
        i = float(i+1)
        j = 2*(i/N)*((N-i)/N)*numSites
        pi += j
    pi = pi * N/(N-1.0)
    D = (pi-theta)/(((params["e1"]*S)+(params["e2"]*S*(S-1)))**0.5)
    
    return(theta/Ns, pi/Ns, D)

if __name__ == "__main__":
    __main__()
