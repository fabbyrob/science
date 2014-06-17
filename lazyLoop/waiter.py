import time
import sys

wait = 6

if len(sys.argv) > 1:
    wait = float(sys.argv[1])
    
num = 100
if len(sys.argv) > 2:
    num = int(sys.argv[2])

for i in range(0, num):
        print(i)
        sys.stdout.flush()
        time.sleep(wait)