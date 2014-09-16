from LinkedNode import LinkedNode #change this to load in your version of LinkedNode

total = 16
passed = 0

print("Testing object creation")
try:
    node1 = LinkedNode(10)
    if node1.value == 10:
        print("Success!")
        passed += 1
    else:
        print("Test fail. Node1.value should be 10, it is: "+str(node1.value))
except Exception, e:
    print("Test broke.")
    print(e)
    
print("")
print("Testing node linkage")
try:
    node1 = LinkedNode(10)
    node2 = LinkedNode(20, node1)
    if node1.next is node2:
        print("Success!")
        passed += 1
    else:
        print("Test fail. node1.next should be node2 it is: "+str(node1.next))
except Exception, e:
    print("Test broke.")
    print(e)
    
print("")
print("Testing node addNode")
try:
    node1 = LinkedNode(10)
    node2 = LinkedNode(20, node1)
    node3 = LinkedNode(17)
    
    node2.addNode(node3)
    if node2.next is node3:
        print("Success!")
        passed += 1
    else:
        print("Test fail. node2.next should be node3, it is: "+str(node2.next))
        
    node1.addNode(node3)
    if node1.next is node3:
        print("Success!")
        passed += 1
    else:
        print("Test fail. node1.next should be node3, it is: "+str(node1.next))
except Exception, e:
    print("Test broke.")
    print(e)
    
print("")
print("Testing printing")
try:
    node1 = LinkedNode(10)
    newNode = LinkedNode(11, node1)
    newNode = LinkedNode(12, newNode)
    newNode = LinkedNode(14, newNode)
    newNode = LinkedNode(13, newNode)
    newNode = LinkedNode(15, newNode)
    result = "10, 11, 12, 14, 13, 15"
    if str(node1) == result:
        print("Success!")
        passed += 1
    else:
       print("Test fail. String should be \""+result+"\" but it is: \""+str(node1)+"\"")
except Exception, e:
    print("Test broke.")
    print(e)
    
print("")
print("Testing index")
try:
    node1 = LinkedNode(10)
    newNode = LinkedNode(11, node1)
    newNode = LinkedNode(12, newNode)
    newNode = LinkedNode(14, newNode)
    newNode = LinkedNode(13, newNode)
    newNode = LinkedNode(15, newNode)
    
    if node1.index(0) == 10:
        print("Success!")
        passed += 1
    else:
       print("Test fail. node1.index(0) should be 10, it is: "+str(node1.index(0)))
       
    if node1.index(2) == 12:
        print("Success!")
        passed += 1
    else:
       print("Test fail. node1.index(2) should be 12, it is: "+str(node1.index(0)))\
       
    if node1.index(10) == None:
        print("Success!")
        passed += 1
    else:
       print("Test fail. node1.index(10) should be None, it is: "+str(node1.index(0)))
except Exception, e:
    print("Test broke.")
    print(e)
    
print("")
print("Testing remove")
print("Note: this test counts on __str__ working above, if that test fails results from this one are not accurate.")
try:
    node1 = LinkedNode(10)
    newNode = LinkedNode(11, node1)
    newNode = LinkedNode(12, newNode)
    newNode = LinkedNode(14, newNode)
    newNode = LinkedNode(13, newNode)
    newNode = LinkedNode(15, newNode)
    
    result1 = "10, 12, 14, 13, 15"
    result2 = "10, 12, 14, 15"
    result3 = "10, 12, 14, 15"
    
    node1.remove(1)
    if str(node1) == result1:
        print("Success!")
        passed += 1
    else:
       print("Test fail. List should be \""+result1+"\", it is: \""+str(node1)+"\"")
       
    node1.remove(3)
    if str(node1) == result2:
        print("Success!")
        passed += 1
    else:
       print("Test fail. List should be \""+result2+"\", it is: \""+str(node1)+"\"")
       
    node1.remove(10)
    if str(node1) == result3:
        print("Success!")
        passed += 1
    else:
       print("Test fail. List should be \""+result3+"\", it is: \""+str(node1)+"\"")
except Exception, e:
    print("Test broke.")
    print(e)
    
print("")
print("Testing append")
print("Note: this test counts on __str__ working above, if that test fails results from this one are not accurate.")
try:
    node1 = LinkedNode(10)
    newNode = LinkedNode(11, node1)
    newNode = LinkedNode(12, newNode)
    newNode = LinkedNode(14, newNode)
    newNode = LinkedNode(13, newNode)
    newNode = LinkedNode(15, newNode)
    
    result1 = "10, 11, 12, 14, 13, 15, 72"
    result2 = "10, 11, 12, 14, 13, 15, 72, 165"
    
    node1.append(72)
    if str(node1) == result1:
        print("Success!")
        passed += 1
    else:
       print("Test fail. List should be \""+result1+"\", it is: \""+str(node1)+"\"")
       
    node1.append(165)
    if str(node1) == result2:
        print("Success!")
        passed += 1
    else:
       print("Test fail. List should be \""+result2+"\", it is: \""+str(node1)+"\"")

except Exception, e:
    print("Test broke.")
    print(e)
    
print("")
print("Testing insert")
print("Note: this test counts on __str__ working above, if that test fails results from this one are not accurate.")
try:
    node1 = LinkedNode(10)
    newNode = LinkedNode(11, node1)
    newNode = LinkedNode(12, newNode)
    newNode = LinkedNode(14, newNode)
    newNode = LinkedNode(13, newNode)
    newNode = LinkedNode(15, newNode)
    
    result1 = "10, 11, 72, 12, 14, 13, 15"
    result2 = "10, 11, 72, 12, 165, 14, 13, 15"
    result3 = "10, 11, 72, 12, 165, 14, 13, 15"
    
    node1.insert(72, 2)
    if str(node1) == result1:
        print("Success!")
        passed += 1
    else:
       print("Test fail. List should be \""+result1+"\", it is: \""+str(node1)+"\"")
       
    node1.insert(165, 4)
    if str(node1) == result2:
        print("Success!")
        passed += 1
    else:
       print("Test fail. List should be \""+result2+"\", it is: \""+str(node1)+"\"")
       
    node1.insert(720, 14)
    if str(node1) == result3:
        print("Success!")
        passed += 1
    else:
       print("Test fail. List should be \""+result3+"\", it is: \""+str(node1)+"\"")

except Exception, e:
    print("Test broke.")
    print(e)
    
print ("")
print("Passed "+str(passed)+" of "+str(total)+" tests.")
if passed == total:
    print("Great job!")