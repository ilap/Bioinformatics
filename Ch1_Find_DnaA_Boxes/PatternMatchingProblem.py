__author__ = 'ilap'
# BioInf - Module 1.8
import os
import time

import time

def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0)
        return ret
    return wrap


kmer=["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT","GA", "GC", "GG", "GT","TA", "TC", "TG", "TT" ]
symbols = ["A", "C", "G", "T"]

def symbolToNumber (char):
    res = 0
    if char == "C":
        res = 1
    if char == "G":
        res = 2
    if char == "T":
        res = 3
    return res

'''
Input: Integers index and k.
Output: The string NumberToPattern(index, k).
'''
def numberToPattern (index, k):
    #print "K i: ", index, k
    if k == 1:
        return symbols[index]

    div = long (index / 4)
    mod = long (index % 4)

    # print " DIVMOD", div, mod
    result = numberToPattern(div, k - 1) + symbols[mod]
    return result

'''
Input: AGT
Output: 11
'''
def patternToNumber (pattern):
    if pattern == "":
        return 0

    length = len (pattern) - 1

    symbol = pattern[length] # last symbol of the pattern
    pattern = pattern [0: length] # remove last symbol from the pattern
    #print "symbol:", symbol, "pattern:", pattern
    return 4 * patternToNumber (pattern) + symbolToNumber (symbol)

'''
Hamming Distance Problem: Compute the Hamming distance between two strings.
     Input: Two strings of equal length.
     Output: The Hamming distance between these strings.
'''
def hammingDistance (p, q):
    mismatch = 0
    plen = len (p)
    qlen = len (q)
    if qlen != plen:
        return -1

    for i in range (0, plen):
       # print "PD: ", p[i], q[i]
        if p[i] != q[i]:
            mismatch += 1

    return mismatch

def patternMatching (pattern, dna, d):
    mismatches = []

    k = len (pattern)
    dlen = len (dna)

    for i in range (0, dlen-k+1):
        temp_pattern = dna[i:i+k]
        if hammingDistance(pattern, temp_pattern) <= d:
            mismatches.append (i)

    return mismatches

def computeCount (dna, pattern, d):
    count = 0

    k = len (pattern)
    dlen = len (dna)

    for i in range (0, dlen-k+1):
        temp_pattern = dna[i:i+k]

        if hammingDistance(pattern, temp_pattern) <= d:
            count += 1

    return count


'''
input:  "AAAACCCGGT"
output: "ACCGGGTTTT"
'''
def reverseComplement (dna):

    arr=[['A', 'C', 'G', 'T'],['T', 'G', 'C', 'A']]
    dlen = len(dna)

    res = ""
    #print "DNA LEN: ", dlen, dna
    for i in range (dlen-1, -1, -1):
        #print "III: ", i
        chr = dna[i:i+1]
        idx = arr[0].index (chr)
        comp_chr = arr[1][idx]
        res += comp_chr
    return res

'''
Sample Input:
TTTAGAGCCTTCAGAGG
GAGG
2
Sample Output:
4
'''

def approximatePatternCount (dna, pattern, d):

    count=0

    klen = len (pattern)
    dlen = len (dna)

    for i in range (0, dlen-klen+1):
        temp_pattern = dna[i:i+klen]

        if hammingDistance(temp_pattern, pattern) <= d:
            count +=1
    return count

nucleotides = ["A", "C", "G", "T"]

def immediateNeighbors (pattern):
    neighborhood = []
    plen = len (pattern)

    for i in range (1, plen+1):
        for j in range (0, len (nucleotides)):
            if nucleotides[j] != pattern[i-1]:
                neighbor =  pattern[0:i-1] + nucleotides[j] + pattern[i:]
                neighborhood.append (neighbor)

    return neighborhood

def neighbors_1 (pattern, d):
    neighborhood = []
    plen = len (pattern)

    for i in range (1, plen+1):
        for j in range (0, len (nucleotides)):
            if nucleotides[j] != pattern[i-1]:
                neighbor =  pattern[0:i-1] + nucleotides[j] + pattern[i:]
                neighborhood.append (neighbor)

    return neighborhood


def neighbors (pattern, d):
    if d == 0: return [pattern]
    elif len (pattern) == 1: return nucleotides

    neighborhood = []
    # means ACTG -> ACT
    suffix = pattern[1:] #0:len (pattern)-1]
    # firstSymbol means CTG -> CTG
    first_symbol = pattern [0]

    #print pattern, suffix, first_symbol
    #return
    suffix_neighbors = neighbors (suffix, d)

    for text in suffix_neighbors:
        ham = hammingDistance(suffix, text)
        #print "HAM:", ham, suffix, text
        if  ham < d:
            for j in nucleotides:
                neighborhood.append (j + text)
        else:
            neighborhood.append (first_symbol + text)

    return neighborhood


'''
Frequent Words with Mismatches Problem: Find the most frequent k-mers with mismatches in a string.
     Input: A string Text as well as integers k and d. (You may assume k <= 12 and d <=3.)
     Output: All most frequent k-mers with up to d mismatches in Text.
'''

#@timing
def frequentWordsWithMismatches (dna, k, d):

    length = len (dna) - k + 1

    count = [None]*length

    kmers = []
    for i in range (0, length):
        pattern = dna[i:i+k]
        count[i] = approximatePatternCount (dna, pattern, d)

    max_count = max (count)

    #print max_count

    freq_patterns = []
    for i in range (0, length):
        tstr = dna[i:i+k]
        if count[i] == max_count and not tstr in freq_patterns:
            freq_patterns.append (tstr)

    freq_patterns.sort ()

    result = []
    #for i in freq_patterns:
        #print "III:", i

        #for j in patternMatching(i, dna, d):
            #print "VVV:" , j
            #if str(dna[j:j+k]) == "GCACACAGAC":
            #    print "megvan"



    return freq_patterns


def frequentWordsBySorting (dna, k):
    frequent_patterns = []
    length = len (dna) - k + 1

    index = [0]*length
    count = [0]*length

    print len  (index), length, k
    print range (0, length)
    print index
    for i in range (0, length):
        pattern = dna[i:i+k]
        index[i] = patternToNumber(pattern)
        print pattern, index[i]
        count[i] = 1

    index = sorted(index)
    print "PRINT IDX", index
    for i in range (1, length):
        if index[i] == index[i-1]:
            count[i] = count[i-1] + 1

    max_count = max (count)
    for i in range (0, length):
        #print "III", i
        if count[i] == max_count:
            pattern = numberToPattern(index[i], k)
            frequent_patterns.append (pattern)

    return frequent_patterns

@timing
def frequentWordsWithMismatchesBySorting (dna, k, d):

    frequent_patterns = []
    neighborhoods = []
    neighborhood_array = []

    dlen = len (dna) - k


    for i in range (0, dlen +1):
        nh = neighbors(dna[i:i+k], d)
        neighborhoods += nh

    nlen = len (neighborhoods)
    index = [0]*nlen
    count = [0]*nlen

    for i in range(0, nlen):
        index[i] = patternToNumber (neighborhoods[i])
        count [i] = 1

    index = sorted (index)

    for i in range(0, nlen-1):
        #print "IND", nlen, i+1, i
        if index [i+1] == index[i]:
            count[i+1]=count[i]+1

    max_count = max (count)
    for i in range(0, nlen ):
        if count[i] == max_count:
            pattern = numberToPattern(index[i], k)
            frequent_patterns.append(pattern)

    return sorted (set (frequent_patterns))
    #   #neighborhood_array =


#@timing
def frequentWordsWithMismatchesBySortingAndComplement (dna, k, d):

    frequent_patterns = []
    neighborhoods = []
    complementaries = []

    dlen = len (dna) - k

    print "DNA length:", dlen
    for i in range (0, dlen +1):

        pattern = dna[i:i+k]
        #print "PATTERN:", pattern, i
        nh = neighbors(pattern, d)
        comp = reverseComplement(pattern)
        cm = neighbors (comp, d)

        neighborhoods += nh
        complementaries += cm


    clen = len (complementaries)

    neighborhoods += complementaries
    nlen = len (neighborhoods)

    #print "CLEN:", clen, nlen
    #print "NG", neighborhoods
    #print "CM", complementaries

    index = [0]*nlen
    count = [0]*nlen

    for i in range(0, nlen):
        index[i] = patternToNumber (neighborhoods[i])
        count [i] = 1

    index = sorted (index)

    for i in range(0, nlen-1):
        #print "IND", nlen, i+1, i
        if index [i+1] == index[i]:
            count[i+1]=count[i]+1

    max_count = max (count)
    print "MAX COUNT: ", max_count
    for i in range(0, nlen ):
        if count[i] == max_count:
            pattern = numberToPattern(index[i], k)
            frequent_patterns.append(pattern)

    return sorted (set (frequent_patterns))
    #   #neighborhood_array =


#@timing
def test():
    return neighbors("CACAGTAGGC", 2)


## MAIN
################################
gene_file = 'Salmonella_enterica_genome.txt'
dna = open(gene_file,'rU').read()
###########################################################################################


### 1.4

print  len ( frequentWordsWithMismatchesBySorting ("ACGT", 4, 3 ))
# print immediateNeighbors ("ACGTACGT")

#print '\r\n'.join (neighbors("GGCCCAGAG", 3))
#print '\r\n'.join (test())

print computeCount ("TACGCATTACAAAGCACA", "AA", 1)
print computeCount ("AACAAGCTGATAAACATTTAAAGAG", "AAAAA" ,1)
# print frequentWordsBySorting ("AACAAGCTGATAAACATTTAAAGAG", 3 )

#newdna = dna[3818620:3819000]
#print "Genome l: ", len (dna), len (newdna)

#print frequentWordsWithMismatchesBySortingAndComplement (newdna, 9, 2)
print hammingDistance ("CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA", "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG")

#print patternCount ("TTTAGAGCCTTCAGAGG", "GAGG", 2 )
#print patternMatching ("GCGCACACAC", "CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC", 2)
#print ' '.join (str(x) for x in mismatches)