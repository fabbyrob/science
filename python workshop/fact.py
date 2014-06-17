def mysteryFunction(value):
    final = 1
    for count in range(1, value+1):
        final = final * count
    return final

def mysteryFunction2(list):
    if list:
        return sum(list)/float(len(list)) 
    else:
        return 0