__author__ = 'ilap'

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

def computeCount (dna, pattern):
    count = 0

    k = len (pattern)
    dlen = len (dna)

    for i in range (0, dlen-k+1):
        temp_pattern = dna[i:i+k]

        if hammingDistance(pattern, temp_pattern) <= 2:
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
    #print "DNA LEN: ", dlen
    for i in range (dlen-1, -1, -1):
        #print "III: ", i
        chr = dna[i:i+1]
        idx = arr[0].index (chr)
        comp_chr = arr[1][idx]
        res += comp_chr
    return res

'''
DNA: ACGTTGCATGTCGCATGATGCATGAGAGCT
FRQ: 4
OUT: CATG GCAT

Text: ACTGACTCCCACCCC
Count:2111211311133--
'''
def computingFrequencies1 (dna, k):

    length = len (dna)
    print "lenght: ", length
    array_length = (4**k)
    frequencies = [0]*array_length

    for i in range (0,  length - k + 1):
        # print "III: ", length, array_length, i
        pattern = dna[i:i+k]
        j = patternToNumber (pattern)
        frequencies [j] += 1

    return frequencies

def computingFrequencies (dna, k, c, t):

    length = len (dna)
    print "lenght: ", length
    array_length = (4**k)
    frequencies = [0]*array_length

    for i in range (0,  length - k + 1):
        # print "III: ", length, array_length, i
        pattern = dna[i:i+k]
        j = patternToNumber (pattern)
        frequencies [j] += 1
        if frequencies[j] >= t:
            c.append (j)

    return frequencies
'''
DNA: ACGTTGCATGTCGCATGATGCATGAGAGCT
FRQ: 4
OUT: CATG GCAT

Text: ACTGACTCCCACCCC
Count:2111211311133--
'''

def frequentWords (dna, k, t = 0):

    freq_patterns = []
    frequencies = computingFrequencies(dna,k)

    if t == 0:
        max_count = max (frequencies)
    else:
        max_count = t
    print "4KK: ", 4**k
    for i in range (0, 4**k):
        pattern = numberToPattern(i, k)
        if frequencies[i] == max_count and not pattern in freq_patterns:
            freq_patterns.append (pattern)

    freq_patterns.sort ()
    return freq_patterns


'''
Clump Finding Problem: Find patterns forming clumps in a string.
     Input: A string Genome, and integers k, L, and t.
     Output: All distinct k-mers forming (L, t)-clumps in Genome.

Sample Input:
     CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
     5 50 4

Sample Output:
     CGACA GAAGA
Steps are:
    1. compute frequencies for the kmers in the whole DNA sequence.
    2. then
'''
@timing
def findClumps (dna, k, L, t):

    length = len (dna)
    result = []
    c = []
    #frequencies = computingFrequencies(dna[0:L], k, c, t)
    frequencies = computingFrequencies1(dna[0:L], k)

    print "DNA legth: ", length, "Freq length: ", len (frequencies), L

    print "1st for: "
    for i in range (0, 4**k):
        if frequencies[i] >= t:
            c.append (i)

    print "2nd for: ", length - L + 1
    for i in range (1, length - L + 1):

        first_pat = dna[i-1:i-1+k]
        j = patternToNumber(first_pat)
        frequencies[j] -= 1


        last_pat = dna[i+L-k:i+L]
        j = patternToNumber(last_pat)
        frequencies[j] += 1


        if frequencies [j] >= t:
            c.append (j)

    result = []
    for i in range (0, len (c)):
        pat = numberToPattern(c[i], k)
        result.append (pat)

    return sorted(set(result))

@timing
def findClumpsOld (dna, k, L, t):

    length = len (dna)
    result = []
    c = [0]*(4**k-1)

    frequencies = computingFrequencies1(dna[0:L], k)

    print "DNA legth: ", length, "Freq length: ", len (frequencies), L

    print "1st for: "
    for i in range (0, 4**k):
        if frequencies[i] >= t:
            c[i] = 1

    print "2nd for: ", length - L + 1
    for i in range (1, length - L + 1):
        if i % 10000 == 0:
            print "Si", i
        first_pat = dna[i-1:i-1+k]
        j = patternToNumber(first_pat)
        frequencies[j] -= 1


        last_pat = dna[i+L-k:i+L]
        j = patternToNumber(last_pat)
        frequencies[j] += 1


        if frequencies [j] >= t:
            c.append (j)

    result = []
    for i in range (0, 4**k):
        if c[i] == 1:
            pat = numberToPattern(i, k)
            result += pat

    return sorted(set(result))

#
# MAIN
#
# Global variables
###########################################################################################
#gene_file = 'E-coli.txt'
gene_file = 'Salmonella_enterica_genome.txt'
dna = open(gene_file,'rU').read()

dna = open(gene_file,'rU').read()
# pattern = "AAAAACGA"
dna="GCACAAGGCCGACAATAGGACGTAGCCTTGAAGACGACGTAGCGTGGTCGCATAAGTACAGTAGATAGTACCTCCCCCGCGCATCCTATTATTAAGTTAATT"

k=4
L=30
t=3
###########################################################################################


### 1.4

f=findClumps (dna, k, L, t)
print "FF:", f
print ''.join (f)

#f=findClumpsOld (dna, k, L, t)
#print ' '.join (f)
# print patternToNumber ("ATGCAA")
# print numberToPattern (5437, 8)
# print computingFrequencies("ACCGGACACGGCTCGTACCCGGCGTAAGTGTCTCGAAGCGCCTATACTCTCGGGCAGATTATCAGGCAGTCGTCATGCGTGGGGCGGACAACAGTGGCCGGTCTTAGGTTGAAAAAGCTCAGGCCTTCCATTTAATAGAATCGGGAAAGCCCGTTAGCGTGCACAGCTTTCCGAATGCAAGAATGTCGACTTTTGCACGTTTAAGGATGACCTTGAGACGAGGTGGTAGCGTACGCAGCCGCCCGTGTTGGGACGGCACGGCAGGCCACTGCCCACTTCTGCACCGTATCTATAAGGCCCGACAACATCCCTCTTGACTAAAGGAACTACGGACCTATCGATAGGCTTGTCGAGACCGATCAACTATGGCCTAAGATAAAAGTGCTCGGTCGTTGCCACGACCTTTACCGCCGTATCTGACGATGGCAGTCCATAGCAGCATAAATTCCGGAAGACCTGTGGCGCAGGCTCCTTAACGGGGTCCTTCATTTATTCCGGGGTACTTCTCGTTACTTCATAGCGTTACCCGAATCTTAAGTGCGGAGAAGTACTATAGCTTCCTGCTTATAACGCTTTACTTTAACCCGGAATCGTCAAACATT", 5)

