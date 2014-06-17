import multiprocessing
import time
import queue
 
def fact(n):
    tot = 1
    for i in range(1,n+1):
        tot *= i
   
    return tot
 
nums = [20000, 20000, 20000]
 
def noThreads(ls):
    for n in ls:
        print("%s! = %s" % (n, fact(n)))
 
def threads(ls):
    ts = []
    pool = multiprocessing.Pool()
    res = pool.map(fact, ls)
    for r in zip(ls, res):
        print("psh %s! = %s" % (r[0], r[1]))
 
if __name__ == '__main__':
    start = time.time()
    noThreads(nums)
    print("%.2e" % ((time.time()-start)))
     
    start = time.time()
    threads(nums)
    print("%.2e" % ((time.time()-start)))