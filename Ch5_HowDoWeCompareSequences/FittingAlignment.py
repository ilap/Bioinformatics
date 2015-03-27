__author__ = 'ilap'

from SequenceCompareLib import *

import numpy as np
def fittingAlignmentGraph (v, w, ro=0, mu=0):

    inc=1
    v_len = len (v)+1
    w_len = len (w)+1



    # use simple scores matrix....
    s = np.zeros (shape=(v_len, w_len))
    #for i in range (v_len):
    #    s[i,0] = -i
    for i in range (w_len):
        s[0,i] = -i


    (row, col, matrix) = generateScoreMatrixFromFile("PAM250.txt")
    # set it for simple score,,, overwrite...
    # match +1 ...

    matrix[:,:] = -1
    for i in range (len(matrix)):
        matrix[i,i] = 1
    #print matrix

    max_row = -sys.maxint
    max_col = -sys.maxint
    max_v = -sys.maxint

    for i in range (1, v_len):
        for j in range (1, w_len):
            v_char = v[i-1]
            w_char = w[j-1]
            v_idx = row.index(v_char)
            w_idx = col.index(w_char)

            down = s[i-1,j] - ro
            left = s[i, j-1] - ro
            downleft = s[i-1, j-1] + matrix[v_idx, w_idx]


            max_t = max (down, left, downleft)
            s[i,j] = max_t

            if j == 1 and max_t > max_row:
                max_row = max_t
                max_start = i

            if j == w_len-1 and max_t > max_col:
                max_col = max_t
                max_end = i


            print i, j, v[i-1], w[j-1], down, left, downleft, matrix[v_idx, w_idx]


    return (s, row, col, matrix, max_start, max_end)


global _v_alignment
_v_alignment = str ()

def backTrackLocalScoreMatrix (backtrack, v, w, s_r, s_c, i, j, ro, row, col, matrix):


    global _v_alignment
    v_len = len (v)
    w_len = len (w)

    #print "IJSR", i,j, s_r, s_c


    if i < s_r:
            print _v_alignment
            return ""

    if j < s_c:
            print _v_alignment
            return ""

    v_idx= row.index(v[i-1])
    w_idx= col.index(w[j-1])

    diag = backtrack[i-1,j-1]+ matrix[v_idx, w_idx]
    down = backtrack[i-1,j] - ro
    right = backtrack[i,j-1] - ro

    max_v = max (0, down, right, diag)

    #print "(", i,",",j,") ",max_v, " DRG", down, right, diag,

    if v_idx == w_idx:
        #print "diag"
        _v_alignment = v[i-1] + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, s_r, s_c, i-1, j-1, ro, row, col, matrix) + w[j-1]
    elif max_v == down:
        #print "down"
        _v_alignment = v[i-1] + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, s_r, s_c, i-1, j, ro, row, col, matrix) + "-"
    elif max_v == right:
        #print "right"
        _v_alignment = "-" + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, s_r, s_c, i, j-1, ro, row, col, matrix) + w[j-1]
    else:
        #print "diag"
        _v_alignment = v[i-1] + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, s_r, s_c, i-1, j-1, ro, row, col, matrix) + w[j-1]


v="GTAGGCTTAAGGTTA"
w="TAGATA"

v="CCCATAATTTCGAAGGGGCTCACGTCTGACTATCGGTTCCTCTCAGGTTTCCAGGTAAATGGGGTGCATGCATACACTCGCAGGAATATCGGTTACTCCGTATGGAATCGCTACACCTAAACCCCAGAATGTCCATTTGATAGACGGATGGATCGACTAGTTGGCAGTCGAAGGTAAAAGTCCCCGAGAGCTTCACCATCAAAGTTAATCAAGGGTTTTTGCGGCCTTGGCCAAAGAGTGACCTCCAGTTCTGCCGAAAGGATCTGGCAGTCCCCTTCATTGATCAAAAAAGTAATTAAGTACCAATGGCCACGACCAACGGCCATCTGCACTGCACAGGGGCGAATCCTATCTGTGAGCCGCAGGATCAAAATACCATGAGTGACGTCAAACCGCGGAGCTATTCAACGAAAAATACAGTCGTGGTTGAAAAGTTTGTCAACGTTGGTGCCGGTTAGTTGGTGTGTTGACTTGCGTCGCAATTAGATAGCATAAATCGAATACCACTAACAGTGCTGACAACCCTAAGGTTAAACGCGCGCCTTTTCCAGTCAGCTTCTGACGCTCGCCTCGGCCATTAATATGGTCTCAAACTAAGCGCTTCCGCCCGTCCGATACCGGCAAGAATATCGGCGGCTGTTCTGCGCTTCGCATTCTGTTATCGTCGCATCGCTCCAAACTGGTCGAGGCCCCTCAAGGGAAATTTACGGCTTAGGATAAATACAGATCGGGTCGGGCCCAAAACGTTTAACTGCGTCTAGACAAGTACGGTACAAGGGGCCATTTGCGTCATCTCGTAAAAGAGGACGTTCAACTCTGGCAACACAACGGTATTCAATATTGTGATTCACGGGGTCAGTTGAGGGTAACTCCGATTGCACATGGTAGGATCGGACCGTTGTGTAAGGTTTGGTATGGAAGGCATTTTTCTCCATTCCACGCTCCTGCGGCTAATTGTGGGGCTAGAAT"


w="GGCGCGGAATGAAGCTGTAGTCGTCTATATAATGGCACAGTCGCAGCTATCGTCGTGAATGTTGAGTCTGCGAAGCTCAATCTG"
ro = 1
(b, row, col, matrix, max_start, max_end) =  fittingAlignmentGraph(v, w, ro)
print b, row, col
print int (max_start), int(max_end), b[max_end,len(w)]

print ''.join (backTrackLocalScoreMatrix(b, v, w, max_start, 1, max_end, len(w), ro, row, col, matrix))

