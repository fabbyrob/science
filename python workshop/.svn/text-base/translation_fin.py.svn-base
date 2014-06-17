def __main__():
    #input sequence to un-translate
    input = 'malwmrllpl lallalwgpd paaafvnqhl cgshlvealy lvcgergffy tpktrreaed lqvgqvelgg gpgagslqpl alegslqkrg iveqcctsic slyqlenycn'
    
    #TODO: put your dictionary here
    myDict = readCodon()

    #TODO: un-translate the sequence
    input = input.upper()
    mySeq = ""
    
    for acid in input:
        if acid in myDict.keys():
            mySeq += myDict[acid]

    #TODO: calculate GC content of this gene
    gcCont = getGC(mySeq)
    
    #print the DNA sequence and the GC content
    print("insulin sequence: "+mySeq)
    print("GC content: " + str(gcCont))
    
    #input sequence to translate
    input2 = 'ATGAGCTTCACGCTGACGAACAAGAACGTGATATTCGTGGCGGGGCTGGGGGGGATAGGGCTGGACACGAGCAAGGAGCTGCTGAAGAGGGACCTGAAGAACCTGGTGATACTGGACAGGATAGAGAACCCGGCGGCGATAGCGGAGCTGAAGGCGATAAACCCGAAGGTGACGGTGACGTTCTACCCGTACGACGTGACGGTGCCGATAGCGGAGACGACGAAGCTGCTGAAGACGATATTCGCGCAGCTGAAGACGGTGGACGTGCTGATAAACGGGGCGGGGATACTGGACGACCACCAGATAGAGAGGACGATAGCGGTGAACTACACGGGGCTGGTGAACACGACGACGGCGATACTGGACTTCTGGGACAAGAGGAAGGGGGGGCCGGGGGGGATAATATGCAACATAGGGAGCGTGACGGGGTTCAACGCGATATACCAGGTGCCGGTGTACAGCGGGACGAAGGCGGCGGTGGTGAACTTCACGAGCAGCCTGGCGAAGCTGGCGCCGATAACGGGGGTGACGGCGTACACGGTGAACCCGGGGATAACGAGGACGACGCTGGTGCACAAGTTCAACAGCTGGCTGGACGTGGAGCCGCAGGTGGCGGAGAAGCTGCTGGCGCACCCGACGCAGCCGAGCCTGGCGTGCGCGGAGAACTTCGTGAAGGCGATAGAGCTGAACCAGAACGGGGCGATATGGAAGCTGGACCTGGGGACGCTGGAGGCGATACAGTGGACGAAGCACTGGGACAGCGGGATA'
    
    #dictionary of items from codon->AA
    #this line of code automagically switched your dictionary from above for you
    invDict = {v:k for k, v in myDict.items()}

    #TODO: translate input2 into AA sequence
    myAAs = ''
    for i in range(0,len(input2),3):
        codon  = input2[i:i+3]
        if codon in invDict.keys():
            myAAs += invDict[codon]
    
    #print generated AA seq
    print("ADH protein: "+myAAs)
    
def readCodon():
    dict = {}
    file = open('codon.txt')
    for i in file:
            tuple = i.split()
            dict[tuple[1]] = tuple[0]
    return dict
    
def getGC(seq):
    #define a function that takes 'seq', a dna sequence, and returns the GC content
    gcCount = 0
    for i in seq:
        if i == 'G' or i == 'C':
            gcCount += 1
    
    return gcCount/len(seq)
    
if __name__ == "__main__":
    __main__()