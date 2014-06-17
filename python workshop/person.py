import random

class Person:
    def __init__(self, myName, myBirthday, myHeight):
        self.name = myName
        self.birthday = myBirthday
        self.height = myHeight
        
        self.tiredness = random.randint(0,10)
        self.hunger = random.randint(0,10)#0 < hunger < 10
        
    def getHungry(self, amount = 1):
        self.hunger = self.hunger + amount
        if self.hunger > 10:
            self.hunger = 10
            
    def eat(self, amount = 1):
        self.hunger = self.hunger - amount
        if self.hunger < 0:
            self.hunger = 0
        
    def sleep(self, amount = 1):
        self.tiredness = self.tiredness - amount
        if self.tiredness < 0:
            self.tiredness = 0
            
    def work(self, amount = 1):
        self.tiredness = self.tiredness + amount
        if self.tiredness > 10:
            self.tiredness = 0
        
    def age(self, today):
        myAge = today[2] - self.birthday[2]
        
        if today[1] < self.birthday[1]:
            myAge = myAge - 1
        elif today[1] == self.birthday[1] and today[0] < self.birthday[0]:
            myAge = myAge -1
        elif today[1] == self.birthday[1] and today[0] == self.birthda[0]:
            print("Happy birthday!")
        
        return myAge
    
    def doSomething(self):
        do = random.randint(1,4)
        #1 == eat
        #2 == getHungry
        #3 == work
        #4 == sleep
        
        if do == 1:
            self.eat()
            print("Om nom nom!")
        elif do == 2:
            self.getHungry()
            print("*rumble*")
        elif do == 3:
            if self.tiredness == 10:
                print("Got no work done! :-(")
            else:
                self.work()
                print("Yay science!")
        elif do == 4:
            self.sleep()
            print("Zzzzz.")
        else:
            print("Wtf?")
    
    def __str__(self):
        myStr = "Person(name="+str(self.name)+", hunger="+str(self.hunger)+", tiredness="+str(self.tiredness)+")"
        return(myStr)
        