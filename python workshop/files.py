#opening a file

myFile = open("text.txt","r")

#closing a file
myFile.close()

#reading a file
myFile = open("text.txt","r")

line = myFile.readline()
print("my line: "+line)
line = myFile.readline()
print("my line: "+line)
line = myFile.readline()
print("my line: "+line)

#getting rid of trailing newlines
line = myFile.readline()
line = line.rstrip()
print("my line: "+line)

line = myFile.readline()
line = line.rstrip()
print("my line: "+line)

line = myFile.readline()
line = line.rstrip()
print("my line: "+line)


myFile.close()
myFile = open("text.txt","r")

for line in myFile:
    line = line.rstrip()
    print (line)
myFile.close()
#define a function that returns the number of lines in that file that use the word "meat"
def meatCounter():
    myFile = open("text.txt","r")
    ctr = 0
    for line in myFile:
        if "meat" in line.lower():
            ctr+= 1
    return ctr        
    
    
print(meatCounter())

#define a function that prints the second word on each line of the given file
def printSecond(file):
    myFile = open(file, "r")
    for line in myFile:
        line = line.split()
        if len(line) > 1:
            print(line[1])
            
printSecond("text.txt")

#writing a file

myOut = open("out.txt", "w")
myOut.write("Hello!")
myOut.write("Hello!\n")
myOut.close()
#appending to a file
myOut = open("out.txt", "a")
myOut.write("Hello!")
myOut.write("Hello!\n")
myOut.close()

#using your factorial function, print the N! for N = 0..30 to a file
#read that file back in and print out the square of each of those numbers