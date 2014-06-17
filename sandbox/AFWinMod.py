#! /usr/bin/env python

# A list containing the valid alternative bases
validAltBase = ['.','A','G','C','T','A,G','G,A','A,C','C,A','A,T','T,A','G,C','C,G','G,T','T,G','C,T','T,C']

# A method to extract the depth of coverage (DP) in float form
def getDP(line):

        # Return false if the line is a newline, empty, or the splitLine is blank
        if line == '\n' or '':
                return False
        else:
                splitLine = line.split()
        if splitLine == [ ]:
                return False
        else:
                # Get the part of the line that starts at DP= and ends at AF1, which is
                # the next parameter.
                # The semicolon is at the beginning of AF1 because we dont' want it
                DP = line[line.find("DP="): line.find(";AF1")]

                # Get just teh DP number, so take everything after the first three c
                # haracters, which are "D" and "P" and "="
                DP = DP[3:]
                
                # Make DP a floating point, so the program knows it's a number
                DP = float(DP)
                
                return DP

def getCI(line):

                CI = line[line.find("CI95="): line.find(";DP4=")]
                CI = CI[5:]

                CLimit1 = CI[:CI.find(',')]
                CLimit2 = CI[CI.find(',')+1:]

                print ("CLimit1")
                print (CLimit1)
                print ("CLimit2")
                print (CLimit2)

                return CLimit1, CLimit2

def getCLimit1(line):

                CI = line[line.find(";CI95=")+1: line.find(";DP4=")]

                CI = CI[5:]
                CLimit1 = CI[:CI.find(',')]

                return float(CLimit1)


def getCLimit2(line):

                CI = line[line.find("CI95="): line.find(";DP4=")]
                CI = CI[5:]
                CLimit2 = CI[CI.find(',')+1:]

                return float(CLimit2)

def getAF1(line):
        # Return false if the line is a newline, empty, or the splitLine is blank
        if line == '\n' or '' or '':
                return False
        else:
                splitLine = line.split()
        if splitLine == [ ]:
                return false
        else:
                # Get the part of the line that starts at AF1 and ends at CI, which is
                # the next parameter.
                # The semicolon is at the beginning of CI because we dont' want it
                AF1 = line[line.find("AF1"): line.find(";CI")]

                # Get just the AF1 number, so chop the first four characters, which
                # are "A", "F", "1" and "=".
                AF1 = AF1[4:]

                # Make DP a floating point, so the program knows it's a number
                AF1 = float(AF1)
                return AF1

# A method to extract the quality score in float form
def getQUAL(line):

        # Return false if the line is a newline, empty, or the splitLine is blank
        if line == '\n' or '' or '':
                return False
        else:
                splitLine = line.split()
        if splitLine == [ ]:
                return false
        else:
                # The sixth delimitation in the line contains the quality score
                QUAL = splitLine[5]

                # Make it a floating point, so the program knows its a number
                QUAL = float(QUAL)
                return QUAL

def getScaffold(line):

        if line == '\n' or '' or '':
                return False
        else:
                splitLine = line.split()
        if splitLine == [ ]:
                return False
        else:
                scaffold = str(splitLine[0])
                return scaffold

def getPosition(line):
        if line == '\n' or '':
                return False
        else:
                splitLine = line.split()
        if splitLine == [ ]:
                return False
        else:
                position = str(splitLine[1])
                return position
                

# A method to extract the alternative base
def getALT(line):

        # Return false if the line is a newline, empty, or the splitLine is blank
        if line == '\n' or '' or '':
                return False
        else:
                line = line.split()
                if line == [ ]:
                        return False
                else:
                        # The alternative base is the fifth delimitation in the line
                        altBase = line[4]

                        # Make the alternative base into a string
                        altBase = str(altBase)
                        return altBase

        
# A function to check if the line in the file meets our criteria
def lineCriteria(currentLine):

        # Return false if the line is a return character, or blank
        if currentLine == '\n' or '' or '':
                return False
        elif currentLine.split() == [ ]:
                return False
        else:
                DP = getDP(currentLine)
                QUAL = getQUAL(currentLine)
                altBase = getALT(currentLine)

        # If the quality score and DP pass these criteria, return true
        #note that instead of 50, we want 160, and instead of 20, we want 40, this is just for test
        if QUAL < 40:
                return False
        if DP < 40 or DP > 600:
                return False
        else:
              if  altBase in validAltBase:
                      return True
              else:
                      return False

# A function to calculate the allele frequency, from 0 to 16

def alleleFreq(line):

        AF1 = getAF1(line)

        #If the confidence interval of AF1 concompasses 0.5,
        # change the AF value to 0.5. This is to make sure
        # we catch all differences between homeologs, and
        # do not count them as real SNPs.
        #C1 = getCLimit1(line)
        #C2 = getCLimit2(line)
        #if C1<=0.5 and C2 >= 0.5:
        #        AF1 = 0.5
        
        AF1 = AF1*26
        AF1 = round(AF1)
        AF1 = int(AF1)

        return AF1

# A function to count the number of alleles at each frequency
# by incrementing the count at allele frequency bins, where
# the frequency corresponds to the index in a list of bins.
def binAF(AF, bins):

        alleleFrequency = int(AF)
        bins[alleleFrequency] +=1
        return bins

        
