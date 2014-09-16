#temperature converter
def temperatureConverter(n, unit="C"):
    if unit == "C":
        temp = 9.0/5.0*n+32
        return temp
    elif unit == "F":
        temp = (n-32)*5.0/9.0
        return temp
    else:
        print("Invalid unit")
        return None
 
#factorial   
def factorial(n):
    fact = 1
    for i in range(1, n+1):
        fact = fact*i
    return fact

#fibonacci
def fibonacci(n):
    n_0 = 0
    n_1 = 1
    if n == 0:
        return n_0
    if n == 1:
        return n_1
    
    for i in range(2,n+1):
        temp = n_0+n_1
        n_0 = n_1
        n_1 = temp
        
    return n_1

#processing strings 1
filename = "annot.txt"
myFile = open(filename, "r")
for line in myFile:
    splitLine = line.split()
    print(splitLine[3]+": "+splitLine[1])
myFile.close() 


#processing strings 2
import math
filename = "annot.txt"
myFile = open(filename, "r")
for line in myFile:
    splitLine = line.split()
    codons = math.ceil((int(splitLine[2])-int(splitLine[1])+1)/3.0)
    print(splitLine[3]+": "+str(codons))
myFile.close() 