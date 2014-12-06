import sys
import summary

reader = summary.Reader(open(sys.argv[1]))

_w = int(sys.argv[2])#this is a dumb way to do this!

print(reader.Header)

prev = []
skip = 0
for record in reader:
    if prev and record.POS == prev[-1].POS:#reached an indel
        #do stuff to fix it
        #clear the buffer, because we're skipping stuff
        for r in prev:
            r.TOTAL = 0
        #is it an deletion
        if len(record.REF) > 1:
            skip = len(record.REF)+_w
        else:#it is an insertion
            skip = _w
        #DONOT add this record to the list, because it is a duplicate
        #the previous non-indel version is already in there with a total of 0
    elif len(prev) == _w+1:#+1 because we want to keep one mroe than window size, to throw out the duplicate site too
        print(prev.pop(0))

    if not skip:
        prev.append(record)
    else:
        record.TOTAL = 0
        prev.append(record)
        skip -= 1
    
for record in prev:
    print(record)

