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
 Implement Neighbors to find the d-neighborhood of a string.
     Input: A string Pattern and an integer d.
     Output: The collection of strings Neighbors(Pattern, d).
'''
NUCLEOTIDES = ["A", "C", "G", "T"]
def neighbors (pattern, d):

    if d == 0: return [pattern]
    elif len (pattern) == 1: return NUCLEOTIDES

    neighborhood = []

    suffix = pattern[1:]
    first_symbol = pattern [0]

    suffix_neighbors = neighbors (suffix, d)

    for text in suffix_neighbors:
        ham = hammingDistance(suffix, text)
        if  ham < d:
            for j in NUCLEOTIDES:
                neighborhood.append (j + text)
        else:
            neighborhood.append (first_symbol + text)

    return neighborhood

'''
Implanted Motif Problem: Find all (k, d)-motifs in a collection of strings.
     Input: A collection of strings Dna, and integers k and d.
     Output: All (k, d)-motifs in Dna.
'''
def motifEnumeration (dna_array, k, d):
    patterns =  []

    kmers = []
    for dna in dna_array:
        #print "DNA:", dna
        dlen = len (dna)

        for i in range (0, dlen - k + 1):
            kmers.append (dna[i:i+k])

    #print "KMER:", kmers

    for kmer in kmers:
        for pattern in neighbors(kmer, d):
            is_in = True
            for dna in dna_array:
                dlen = len (dna)
                tmp_in = False
                for i in range (0, dlen - k + 1):
                    tkmer = dna[i:i+k]
                    if hammingDistance(pattern, tkmer) <= d:
                        tmp_in = True
                is_in &= tmp_in
            if is_in == True:
                 patterns.append(pattern)

    return list(set (patterns))

#@@@@@ 3.3 Scoring Motifs
import math
import numpy
def entropy (matrix):

    res = 0.0

    for col_idx in range (0, len (matrix[0])):
        col = matrix[:, col_idx]
        ent = 0.0
        for val in col:
            if val != 0.0:
                ent += -(val * math.log (val,2))
        res += ent
        #print "VAL", ent
    return res

from collections import Counter

def countMatrix (motifs):
    col_len = len (motifs)
    row_len = len (motifs[0])

    count_matrix = [[0]*row_len]*4
    matrix = []
    for row in motifs:
        matrix.append([str (x) for x in row])

    nmatrix = np.array(matrix)
    cmatrix = np.array(count_matrix)

    for col_idx in range(0, row_len):
        col = nmatrix[:, col_idx]
        #print "COL", col
        counter = Counter (col.tolist())
        cmatrix[0, col_idx] = counter['A']
        cmatrix[1, col_idx] = counter['C']
        cmatrix[2, col_idx] = counter['G']
        cmatrix[3, col_idx] = counter['T']


    return cmatrix

def profileMatrix (matrix):
    return matrix/10.

def scoreMatrix (motifs):
    result = 0

    col_len = len (motifs)
    row_len = len (motifs[0])

    matrix = []
    for row in motifs:
        matrix.append([str (x) for x in row])

    nmatrix = np.array(matrix)

    for col_idx in range (0, row_len):
        col = nmatrix[:, col_idx]
        counter = Counter (col.tolist())
        (N, maxv) = counter.most_common()[0]
        #print "MAX", maxv, col_len

        result += col_len - maxv


    return result

def consensusMatrix (motifs):
    result = ""

    col_len = len (motifs)
    row_len = len (motifs[0])

    matrix = []
    for row in motifs:
        matrix.append([str (x) for x in row])

    nmatrix = np.array(matrix)

    for col_idx in range (0, row_len):
        col = nmatrix[:, col_idx]
        counter = Counter (col.tolist())
        (N, maxv) = counter.most_common()[0]
        #print "MAX", N

        result += N

    return result

######################################################
# #### Main
######################################################
gene_file = 'Implanted_Motif_2.txt'
fd = open(gene_file,'rU')

dnas = []
for line in fd:
    dnas.append(line)
######################################################
#@@@@@ 3.4 Fro Motif Finding to Finding a Median String


#@@@@@ 3.3 Scoring Motifs


import numpy as np

ENT = np.array ([
    [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]], np.float32)

print entropy(ENT)

T= "TCGGGGGTTTTT CCGGTGACTTAC ACGGGGATTTTC TTGGGGACTTTT AAGGGGACTTCC TTGGGGACTTCC TCGGGGATTCAT TCGGGGATTCCT TAGGGGAACTAC TCGGGTATAACC"

motifs = T.split(' ')


#print scoreMatrix(motifs)
#print "CONS", consensusMatrix (motifs)
#matrix = countMatrix(motifs)

#print matrix
#print profileMatrix(matrix)

#@@@@@ 3.2 Finding Motif - brute force - PASSED
'''
T= "GACCATCAGCACGCACTGAATGTGA TCTACCTTGTACCGATCGGGGGTAA TGTAATTGGTAGTAGCGTGAAATGC TGTGACCGCCGATTTAATATCGGGT AGACCTCAGGCGTAACCGCATAGAC ACGCAACTAATCCACAGTCAGTTAG"
DNA = T.split(' ')

k = 5
d = 2
RES = motifEnumeration (DNA, k , d)
print ' '.join (RES)
#print len (RES)
'''

print profileMatrix(countMatrix(["AAGTGA"]))