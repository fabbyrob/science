total = 0
x = -1

for i in range(0, 5):
    total = total + i
    x = x*x*i

    #note % is the modulus function, it gives you the remainder. i.e. 5 % 2 = 1, 11 % 3 = 2
    if x % 2 == 0:
        print(x)
        x = x + 1
    
print(total)
print(x)