__author__ = 'ilap'

'''
Hamming Distance Problem: Compute the Hamming distance between two strings.
     Input: Two strings of equal length.
     Output: The Hamming distance between these strings.
     Desc.: The Hamming distance is the mismatches between two string
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

'''
3.4 Equivalent Motif Finding Problem: Given a collection of strings, find a collection of k-mers (one from each string) that minimizes the distance between all possible patterns and all possible collections of k-mers.
     Input: A collection of strings Dna and an integer k.
     Output: A k-mer Pattern and a collection of k-mers, one from each string in Dna, minimizing
     d(Pattern, Motifs) among all possible choices of Pattern and Motifs
'''
import sys
def d1 (pattern, dna):
    plen = len (pattern)

    min = sys.maxint
    for i in range (0, len (dna) - plen + 1):
        kmer =  dna[i:i+plen]

        val = hammingDistance(pattern, kmer)
        if min > val:
            min = val


    return min

def d (pattern, dnas):

    result = 0

    for dna in dnas:
        val =    d1 (pattern, dna)
        #print val
        result += val
    return result

def motif (pattern, DNA):
    plen = len (pattern)

    min = sys.maxint
    for i in range (0, len (DNA) - plen + 1):
        kmer =  DNA[i:i+plen]

        val = hammingDistance(pattern, kmer)
        if min > val:
            min = val
            result = kmer

    return result

def generateAllKmers (pattern, k):

    if k == 0: return pattern

    result = ""

    for N in ["A", "C", "G", "T"]:
        tval = generateAllKmers(pattern + N, k-1)
        if k == 1:
            result += tval + " "
        else:
            result += tval

    return result

def medianString (dnas, k):
    distance = sys.maxint

    median = ""
    medians = []
    for kmer in [str (x) for x in generateAllKmers ("", k).split ()]:
        #print "KMER", kmer
        val = d (kmer, dnas)
        if distance >= val:
            distance = val
            median = kmer
            medians.append (kmer)

    print medians

    return median
######################################################
# #### Main
######################################################

print d1 ("GATTCTCA", "GCAAAGACGCTGACCAA")

print motif ("GATTCTCA", "GCAAAGACGCTGACCAA")
print motif ("AAG", "GCAATCCTCAGC")

print d ("AAA", ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"])

T = "TACCGCTAGATCATTAGCACCCCACGTCTGCACTTACAGAAG CGGACTTAACCTGTGGGCTTATCGCTGGTTGTTAGCACAATT ACTTTCTCTAAGATCTGTTTACGCACTAATCTTAGCGGACGG GAGCGGTAGCTTCGGGTCTACCCTGAAAAGCTTAGCGTTTGC CGCAGCGTACGACTCAGTTCAGGATTGCCTTTTAGCTCTTAA TGCTTGAGTTACACTAGTTGGCCGTCTTGGCGTAGATTTAGC TATCATGGGAACGTTAGCTGTTAACTTGACTTGCCTAACCTA GTATCAAACAAACCCGCACTTAGCGAACGAGTTGTGCTCCCC GTTAGCGTCTTTAATGGCGCCGTCCCTTGGTCTTCATTAATT TCAGTTATTAGCATAATGAGGAGTGCGCCCTCCTCCTGTTCT"

T ='''CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC
GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC
GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG'''

DNA = T.split ()

print "DNA", DNA
print "MS:", medianString(DNA, 7)

