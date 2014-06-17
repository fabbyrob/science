import re

myPattern = "\w+_(\d+)\s+(\d+)"

regEx = re.compile(myPattern)

file = open("sites.txt", "r")

for line in file:
    line = line.rstrip()
    
    matchCheck = regEx.match(line)
    
    num1 = int(matchCheck.group(1))
    
    num2 = int(matchCheck.group(2))
    
    print(num1, num2)