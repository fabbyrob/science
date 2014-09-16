#none
None

#boolean
True
False


#strings
myStr = "apple"
myStr2 = "pear"

myStr3 = myStr + myStr2

print("myStr + myStr2 = " +myStr3+"\n")

#numbers
anInt = 4
anInt2 = 6

res = anInt + anInt2
print("anInt + anInt2 = " + str(res)+"\n")

res = anInt/anInt2
print("anInt/anInt2 = "+ str(res)+"\n")

#lists
myList = []
myLength = len(myList)

print("list len:")
print(myLength)

myList = [1,2,3]

myList[0]
myList[1]
myList[-1]
myList[-2]

print("list len:")
print(len(myList))

myList.append(4)
print(myList)

myList.pop()
print(myList)

myList.reverse()
print(myList)

myList.insert(1, 5)
print(myList)

mySubList = myList[0:2]
print(mySubList)

mySubList = myList[2:]
print(mySubList)

myList = range(10)
print(list(myList))

myList = range(4,10)
print(list(myList))

myList = range(0,10,3)
print(list(myList))

print(4 in myList)
print(3 in myList)

myStr = "Wright Lab does python"
print(myStr[0])

print("Lab" in myStr)

myList = myStr.split()
print(myList)



#dictionaries
myDict = {}

myDict = {'TTT':'Phe', 'TCT':'Ser'}
print(myDict)

myDict['CTT'] = 'Leu'
print(myDict)

keys = myDict.keys()
print(list(keys))

item = myDict['TTT']
print(item)