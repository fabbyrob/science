params = {}
'''
Used internally to calculate TajD stats if needed.
'''
def calcParams(_N):
    myparams = {}
    a1 = 0
    a2 = 0
    for i in range(1,_N):
        a1 += 1.0/i
        a2 += 1.0/(i**2)
        
    myparams["a1"] = a1
    myparams["a2"] = a2
    
    b1 = float(_N+1)/(3*(_N-1))#see Tajima 1989
    b2 = float(2*(_N**2+_N+3))/(9*_N*(_N-1))

    myparams["b1"] = b1
    myparams["b2"] = b2

    c1 = b1 - 1/a1
    c2 = b2-float(_N+2)/(a1*_N)+a2/(a1**2)

    myparams["c1"] = c1
    myparams["c2"] = c2

    e1 = c1/a1
    e2 = c2/(a1**2+a2)
    
    myparams["e1"] = e1
    myparams["e2"] = e2
    
    global params
    params[_N] = myparams

'''
Returns the per-site pairwise diversity from a given AFS and sample size.
'''
def Pi(afs, N):
    """
    >>> round(Pi([5,0,0,0],6),3)
    0.0
    
    >>> round(Pi([4,1,0,0],6),3)
    0.067
    
    >>> round(Pi([4,0,1,0],6),3)
    0.107
    
    >>> round(Pi([4,0,0,1],6),3)
    0.12
    
    >>> round(Pi([3,1,0,1],6),3)
    0.187
    
    >>> round(Pi([3,0,1,1],6),3)
    0.227
    
    >>> round(Pi([3,2,0,0],6),3)
    0.133
    
    >>> round(Pi([2,1,1,1],6),3)
    0.293
    """
    pi = 0
    total = float(sum(afs))
    N = float(N)
    for i, freq in enumerate(afs):
        if i == 0:
            continue
        pi += 2.0*(i/N)*((N-i)/N)*freq
    
    pi = pi * N/(N-1)
    
    return pi/total

'''
Returns the per-site watterson's theta for a given afs and sample size.
'''
def Theta(afs, N):
    """
    >>> round(Theta([5,0,0,0],6),3)
    0.0
    
    >>> round(Theta([4,1,0,0],6),3)
    0.088
    
    >>> round(Theta([4,0,1,0],6),3)
    0.088
    
    >>> round(Theta([4,0,0,1],6),3)
    0.088
    
    >>> round(Theta([3,1,0,1],6),3)
    0.175
    
    >>> round(Theta([3,0,1,1],6),3)
    0.175
    
    >>> round(Theta([3,2,0,0],6),3)
    0.175
    
    >>> round(Theta([2,1,1,1],6),3)
    0.263
    """
    S = sum(afs[1:])#num polymorphic sites
    
    if not params.has_key(N):
        calcParams(N)
        
    myparams = params[N]
    theta = float(S)/myparams["a1"]
    
    return (theta/sum(afs))

'''
Either calculates Tajima's D from an AFS and sample size 
or calculates it from a given pi, theta, number of polymorphic sites, and sample size.

Note if you provide pi and theta they should NOT be per-site estimates.

Returns None if there are no segregating sites.
'''
def tajD(afs=None, N=None, pi=None, theta=None, S=None):
    """
    >>> tajD(afs=[5,0,0,0],N=6)
    
    
    >>> tajD(N=6, pi=0.0, theta=0.0, S=0)
    
    
    >>> round(tajD([4,1,0,0],6),3)
    -1.12
    
    >>> round(tajD(N=6, pi=0.067*6, theta=0.088*6, S=1),3)
    -1.124
    
    >>> round(tajD([4,0,1,0],6),3)
    1.021
    
    >>> round(tajD(N=6, pi=0.107*6, theta=0.088*6, S=1),3)
    1.017
    
    >>> round(tajD([4,0,0,1],6),3)
    1.734
    
    >>> round(tajD(N=6, pi=0.12*6, theta=0.088*6, S=1),3)
    1.712
    
    >>> round(tajD([3,1,0,1],6),3)
    0.373
    
    >>> round(tajD(N=6, pi=0.187*6, theta=0.175*6, S=2),3)
    0.39
    
    >>> round(tajD([3,0,1,1],6),3)
    1.671
    
    >>> round(tajD(N=6, pi=0.227*6, theta=0.175*6, S=2),3)
    1.688
    
    >>> round(tajD([3,2,0,0],6),3)
    -1.358
    
    >>> round(tajD(N=6, pi=0.133*6, theta=0.175*6, S=2),3)
    -1.363
    
    >>> round(tajD([2,1,1,1],6),3)
    0.72
    
    >>> round(tajD(N=6, pi=0.293*6, theta=0.263*6, S=3),3)
    0.707
    """
    
    if not N:
        raise ValueError("You must always input a sample size.")
    
    if afs:
        S = float(sum(afs[1:]))
        pi = Pi(afs, N)*N
        theta = Theta(afs, N)*N
    else:
        if pi == None or theta == None or S == None:
            raise ValueError("If you provide one of pi, theta, or S you must provide them all. They should NOT be per-site estimates.") 
    
    if not params.has_key(N):
        calcParams(N)
    
    myparams = params[N]
    
    if S == 0:#no segregating sites
        return None
    
    d = pi - theta
    d /= (myparams["e1"]*S+myparams["e2"]*S*(S-1))**0.5
    return d

if __name__ == "__main__":
    import doctest
    doctest.testmod()