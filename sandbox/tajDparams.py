import sys

if len(sys.argv) < 2:
    sys.stderr.write("Please provide a sample size (chromosomes).\n")
    sys.stderr.write("Optionally provide pi, theta, and S to calculate Taj D for a sammple.\n")
    sys.stderr.write("example: python tajDParams.py 26 0.015 0.011 150\n")
    sys.exit()

def calcParams():
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
    
    print("Tajima's D parameters, named as in Tajima 1989 Genetics:\n\t%s\n" % params)
    
_N = int(sys.argv[1])
params = {}
calcParams()

print("Taj D = (pi - theta)/(e1*S+(e2*S*(S-1)))**0.5")
print("Where S is the number of polymorphic sites.")

if len(sys.argv) >= 5:
    pi = float(sys.argv[2])
    theta = float(sys.argv[3])
    S = float(sys.argv[4])
    
    tajD = (pi - theta)/((params["e1"]*S+(params["e2"]*S*(S-1)))**0.5)
    print("Taj D = %s" % tajD)