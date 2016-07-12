import numpy as np
import sys
from pprint import pprint
import random
import csv
import time

########################################
##         HAPLOTYPE ASSEMBLY         ##
########################################

X = "_" # the blank character in the read matrix, anything except [0,1]

# simulates the readMatrix based on number of SNPs and the read Length
def genHaps(numSNPs, readLen):
    length = numSNPs * 10
    readLength = readLen
    numReads = 15;
    seq1 = [X] * length

    randPos = np.random.randint(0, length - 1, numSNPs) # snps occur randomly in sequence

    hapDict = {}
    for i in randPos:
        # assign a random SNP, 0 or 1, for each random position
        hapDict[i] = random.randint(0,1)

    sortedKeys = sorted(hapDict.keys()) 

    # create the readMatrix
    readMatrix = []
    for i in range(0, numReads): # number of reads variable
        read = [X] * numSNPs
        startPos = np.random.randint(0, length) # random starting position
        flip = random.randint(0,1) # randonly choose which chromosome the read is from

        for i in range(0, readLength):
            try:
                index = sortedKeys.index(startPos + i)
                val = hapDict[startPos + i]
                if flip:           # select from other chromosome
                    val = 1 - val
                read[index] = val
            except Exception: # read goes past haplotype end (just truncate current read)
                continue
        readMatrix.append(read)

    return readMatrix


# helper function for merging two rows of readMatrix
# hap1, hap2 should be same length
def hapMerge(hap1, hap2):

    if hap1 == []:
        hap1 = [X] * len(hap2)

    if hap2 == []:
        hap2 = [X] * len(hap1)

    mergedHap = []

    for i in range(0, len(hap1)):
        #mismatch
        if (hap1[i] == 0 and hap2[i] == 1) or (hap1[i] == 1 and hap2[i] == 0):
            raise Exception("mismatch when merging haps")

        elif hap1[i] != X:
            mergedHap.append(hap1[i])
        elif hap2[i] != X:
            mergedHap.append(hap2[i])
        else:  
            mergedHap.append(X);

    return mergedHap

# helper, returns complement of a hap
def hapComp(hap):
    comp = []
    for i in hap:
        if i == X:
            comp.append(X)
        elif i == 0:
            comp.append(1)
        elif i == 1:
            comp.append(0)
    return comp


# greedy helper
# sort by first occuring SNP pos
def sortByPos(readMatrix):
    sortedRead = [];

    for j in range(0, len(readMatrix[0])):
        for i in readMatrix:
            if (i[j] == 0 or i[j] == 1):
                sortedRead.append(i)
                readMatrix.remove(i)

    for i in readMatrix: # add any blank reads
        sortedRead.append(i)

    return sortedRead

# the greedy algorithm for finding the haplotype given a readMatrix
def greedy(readMatrix): 
# the read matrix uses variable X for no info
# 2d array: columns = SNP positions
# rows  = reads at thsoe SNP positions

    hap1 = []
    hap2 = []
    numHap1 = 0
    numHap2 = 0

    # sort the readMatrix for the first occuring SNP?
    inOrder = sortByPos(readMatrix)

    # assign first read to hap1 to start
    # find longest hap
    lengths = []
    for i in readMatrix:
        snps = [1 for x in i if x != X]
        lengths.append(len(snps))

    startIndex = lengths.index(max(lengths)) # start with read with most SNPs

    hap1.extend(inOrder[startIndex])
    numHap1 += 1

    for j in range(0, len(inOrder) - (numHap1 + numHap2)): # need to run through array multiple times
        for i in range(0, len(inOrder)):
            readComp = zip(hap1, inOrder[i])
            for tup in readComp:
                # match with hap1
                if (tup[0] == 0 and tup[1] == 0) or (tup[0] == 1 and tup[1] == 1):
                    hap1 = hapMerge(hap1, inOrder[i])
                    numHap1 += 1
                    break

                # doesn't match with hap1 (so it matches with hap2)
                elif (tup[0] == 0 and tup[1] == 1) or (tup[0] == 1 and tup[1] == 0):
                    hap1 = hapMerge(hap1, hapComp(inOrder[i]))
                    numHap2 += 1
                    break

                # possibly something else here to deal with error rates
                else:
                    continue

    hap1 = hapMerge(hap1, hapComp(hap2))
    hap2 = hapMerge(hap2, hapComp(hap1))
    
    return [hap1, hap2]
    
def baseline(readMatrix): # slow baseline method, should be very slow ~20+ SNPs
    allBitStrings = genBitStrings(len(readMatrix[0])).split(" ")
    
    for i in allBitStrings: # every possible string
        bitStr = [int(j) for j in list(i)] # convert to format [0, 0, 0] ... [1, 1, 1] for hapMerge
        matchesAll = True
        for read in readMatrix:
            try:
                hapMerge(bitStr, read)
            except Exception:
                try:
                    hapMerge(bitStr, hapComp(read))
                except Exception:      # doesn't match read or comp(read)
                    matchesAll = False
                    break
                    
        if matchesAll:           
            return i # only reaches here if it matches all reads / comps(reads)
                

def genBitStrings(n): # all 2^N binary string combinations
    def genHelper(n, bitStr):
        if n == 0:
            return bitStr
            
        left = genHelper(n-1, bitStr + "0")
        right = genHelper(n-1, bitStr + "1")
        
        return str(left + " " + right) # returned as string with spaces, n = 2 -> "00 01 10 11"
    
    return genHelper(n, "")
    

# prints alligned elements of readmatrix
def allignedPrint(readMatrix):
    printMatrix = []
    for i in readMatrix:
        temp = []
        for j in range(0, len(i)):
            temp.append( str(i[j]).ljust(len(X)) )
        printMatrix.append(temp)
    with open('C:\\Users\\Neil\\Dropbox\\UCLA\\124\\Final Project\\out.txt', 'a') as out:
        pprint(printMatrix, out, indent=4, width=5000)
    print
    
def testPrint(readMatrix):
    for i in readMatrix:
        temp = []
        for j in range(0, len(i)):
            temp.append( str(i[j]).ljust(len(X)) )
        print temp


def writeCSV(data, nameIn): # write csv data, helper for the data gen functions
    fname = 'C:\\Users\\Neil\\Dropbox\\UCLA\\124\\Final Project\\' + nameIn
    with open(fname, 'w') as dataOut:
        wr = csv.writer(dataOut, quoting=csv.QUOTE_ALL)
        wr.writerows(data)
        
def snpSizeTest(): # simulates data, varying the # of snps
    graphData = []
    for snpSize in range(10, 200, 10):
        errSum = 0
        for numTrial in range(0, 50):  # take average of 50 trials
            readMatrix = genHaps(snpSize, 150)
            x = greedy(readMatrix)
            errSum += len([1 for i in x[0] if i == X])
            
        accuracy = 1 - (errSum/50.)/snpSize
        graphData.append([snpSize, accuracy])
    
    writeCSV(graphData, "sizeAccuracy.csv")
    return graphData
    
def readLenTest(): # simulates data, varying the readLen
    graphData = []
    for readLen in range(150, 351, 10):
        errSum = 0
        for numTrial in range(0, 50):  # take average of 50 trials
            readMatrix = genHaps(150, readLen)
            x = greedy(readMatrix)
            errSum += len([1 for i in x[0] if i == X])
            
        accuracy = 1 - (errSum/50.)/150
        graphData.append([readLen, accuracy])
    
    writeCSV(graphData, "readLenAccuracy.csv")
    return graphData
    
def greedyTimingTest(): # measures greedy time to completion, varies # of snps
    graphData = []
    for snpSize in range (5, 36, 2):
        readMatrix = genHaps(snpSize, 150)
        runTime = calcTime(greedy, readMatrix)
        graphData.append([snpSize, runTime])
        
    writeCSV(graphData, "greedyTiming.csv")    
    return graphData
    
def baselineTimingTest(): # measures baseline time to completion, varies # of snps
    graphData = []
    for snpSize in range (5, 24, 1):
        readMatrix = genHaps(snpSize, 150)
        runTime = calcTime(baseline, readMatrix)
        graphData.append([snpSize, runTime])
        
    writeCSV(graphData, "baselineTiming.csv")    
    return graphData
    
def calcTime(function, *args): # helper to calc function run time
    startTime = time.time()
    result = function(*args)
    return time.time() - startTime
    
    
def main():
    
    # readLenTest()
    # snpSizeTest()
    # print greedyTimingTest()
    # baselineTimingTest()
    numSNPS = 10
    readLen = 15

    # example problem, run above tests to write csv data to file
    print "Read Matrix"
    readMatrix = genHaps(numSNPS, readLen)
    testPrint(readMatrix)
    
    print "Hap solutions"
    solutions = greedy(readMatrix)
    testPrint(solutions)

    print "end main"

if __name__ == "__main__":
    main()
