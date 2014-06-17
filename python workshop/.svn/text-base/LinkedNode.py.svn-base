class LinkedNode:
    def __init__(self, value, parent=None):
        self.value = value
        self.next = None
        if parent != None:
            parent.next = self
    
    def addNode(self, node):
        self.next = node
    
    def index(self, n):
        toCheck = self
        while n > 0 and toCheck != None:
            n -= 1
            toCheck = toCheck.next
        
        if toCheck:
            return toCheck.value
        else:
            return None
    
    #a more elegant implementation of index
#    def index(self, n):
#        if n == 0:
#            return self.value
#        elif self.next == None:
#            return None   
#        else:
#            return self.next.index(n-1)
                
    def remove(self, n):
        if n <= 0:
            return None
        previous = self
        toCheck = self.next
        while n > 1 and toCheck != None:
            previous = toCheck
            toCheck = toCheck.next
            n -= 1
        
        if toCheck == None:
            return None
        
        previous.next = toCheck.next

    #A more elegant implementation of remove
#    def remove(self, n):
#        if self.next == None:
#            return
#        
#        if n == 1:
#            self.next = self.next.next
#        else:
#            self.next.remove(n-1)
          
    def append(self, value):
        toCheck = self
        while toCheck.next != None:
            toCheck = toCheck.next
            
        LinkedNode(value, toCheck)
          
    #A more elegant implementation of append
#    def append(self, value):
#        if self.next == None:
#            LinkedNode(value, self)
#        else:
#            self.next.append(value)
             
    def insert(self, value, index):
        toCheck = self
        newNode = LinkedNode(value)
        while toCheck != None and index > 1:
            toCheck = toCheck.next
            index -= 1
        
        if toCheck != None:
            newNode.next = toCheck.next
            toCheck.next = newNode
             
    #a more elegant solution to insert
#    def insert(self, value, index):
#        if self.next == None:
#            return
#        
#        if index == 1:
#            remainder = self.next
#            newNode = LinkedNode(value, self)
#            newNode.next = remainder
#        else:
#            self.next.insert(value, index-1)
                             
    def __str__(self):
        myStr = str(self.value)
        toAdd = self.next
        while toAdd != None:
            myStr += ", "+str(toAdd.value)
            toAdd = toAdd.next
        return myStr
    
    #This is a much more elegant implementation of __str__ uncomment it out and give it a try
    #the output is the same as the other __str__ function
#    def __str__(self):
#        if self.next == None:
#            return str(self.value)
#        
#        return str(self.value)+", "+str(self.next)
            
