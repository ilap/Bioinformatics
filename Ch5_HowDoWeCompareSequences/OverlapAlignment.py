__author__ = 'ilap'

from SequenceCompareLib import *

import numpy as np
def overlapAlignmentGraph (v, w, ro=0, mu=0):

    inc=1
    v_len = len (v)+1
    w_len = len (w)+1
    print "LEM", v_len, w_len

    s = np.zeros (shape=(v_len, w_len))

    (row, col, matrix) = generateScoreMatrixFromFile("PAM250.txt")
    # set it for simple score,,, overwrite...
    # match +1 ...
    print col

    matrix[:,:] = -2
    for i in range (len(matrix)):
        matrix[i,i] = 1

    s_i = -sys.maxint
    s_j = -sys.maxint
    e_i = -sys.maxint
    e_j = -sys.maxint
    max_start = -sys.maxint
    max_end = -sys.maxint


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

            max_v = max (down, left, downleft)
            s[i,j] = max_v


            if j == 1  and max_v > max_start:
                max_start = max_v
                s_i = i
                s_j = j

            if i == len (v) and max_v > max_end:
                max_end = max_v
                e_i = i
                e_j = j

            #print i, j, v[i-1], w[j-1], down, left, downleft, matrix[v_idx, w_idx]
    print "RSCS", s_i, e_i, s_j, e_j
    return (s, row, col, matrix, s_i, s_j, e_i, e_j)


global _v_alignment
_v_alignment = str ()
global idx
idx=0
def backTrackLocalScoreMatrix (backtrack, v, w, s_r, s_c, i, j, ro, row, col, matrix):


    global _v_alignment
    global idx
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

    max_v = max (down, right, diag)
    idx +=1
    print idx,"(", i,",",j,") ",max_v, " DRG", down, right, diag,

    if max_v == down:
        print "down"
        _v_alignment = v[i-1] + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, s_r, s_c, i-1, j, ro, row, col, matrix) + "-"
    elif max_v == right:
        print "right"
        _v_alignment = "-" + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, s_r, s_c, i, j-1, ro, row, col, matrix) + w[j-1]
    else:
        print "diag"
        _v_alignment = v[i-1] + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, s_r, s_c, i-1, j-1, ro, row, col, matrix) + w[j-1]

v="TCTGCATATAACGGGCCAGGCCCCCTATTGGCGAAGCTCGAAATAACTGTCTCCGGTTGTGGAGACCACGGCTAATCCCATGGCCCAACGGGAGACGTCAATTCAACTCTTGTATTATCATCCCATATCCGAACATAGAAGGGCAAGAAACCATCCGTTTGCACTATGTTTACATCATACTACACCGTCTAGCCAGGAATTGGCCGGTTGTCATGAATCGCGCGGCATCCTGGAGGCAGTAGGGCTACCCTTCTACGCTCTAATTTCTAAGCCTATATACTAGTCTGGAGGACCAGCACGTGACTCCCACGAATATCTGAGGTACACGTGTCCGCGGAGCTGGGCTACCAATATCCGAGACGATAGAGTACCCTTTCCAGCCCCGTATGGAACCGAAATTGTTCTGCGATGATGGAGTGGACCTCAATCTATCACTTTTCATATTTGCACGTCGGGCGTGACGTGTCATCGAACATTGCCTAAGCAATACGGTAAAGTTCGACATCCAGCGCCAGTACACGGGGACTCCTCATGGAACAAGCCAGTAGAACGTGAGACTCCGCACATAATCATGACCAAGTTCCTAGGTACCACTATGTTATCCCTCATTGGACCGGGGGATAGACGAAGTACATGTTACCTCTGGAGAACTTAGAAACTCATGGCCGTACCCCGCGTTGACTCCACCAGATTGCCGACCGATACGACCTGCCCTGCCTCACCGGTTCTCGAAGGGTTGTCTGAGTTTATCTTTGGCTGCGAGAGGTAGTAAGAGATACATATTCGATCATAAGTGACCCATAGGCAATTCAATTATGCACGCCACCTCTCCTCATCTGTTCTCGTAGTACGTGGTACTACTTCACGGCATGCATGGACGTAATACGTTTCATGTCAAATGACATCGTCGTGCATGCTCCCAGAGTAGCCAGAGTCTCTATCACCGCTACCGCTCCACTTTTCGCTTGCTACG"
w="TGCAGAGAGGCATAGGAGAAATATACACCATAAGCGACCTCTAGGGGGATTTATTATGCACGCACCCTCCGGTTTCATCCTTCTATAGTATATGGTATGATTATACTGTGTAGGAAGTAATAGCGTTAGCATGCTAGTAACATTGGCCGTTAATTCACCACGGTACGCCGAGATAGGTTTTATCACAGCTGACCGATCCACCTATTACCCTTGCAACGCTGCCTGATACTCACGTTGCTTCCACCTGAAAACAACACGGTACCACGGCCCCTGTGTGACTAAGTAGTTCTCACTCTTGCTTCATAATGTGACATACGCATTCGCATTCTTGGAAATTACACAACTTTTCGTGGGGCCCTCTAAAACATAATTGTCCGACTGAGAATCTCTTGAGACTGCCTTCAAATAACCCGAGGCCTTTGTGATTAGGATTCGTACAATATTTTCAGTTCCGCGTCTAAACAATATAGATACTGATATATATGGCCGGGTCTCCCACACTACCCTAGATAACACTGTGATCCAATAGAAACTCCGCCATTGTGATAAAGCACGGACCATAGTGCCTGCCGCTATTTCTAACTCCGCCTTAGAGGTACACTTCAAACCGACGTGTCGGCCCAAAAGGTGTTTCCCGCCTCCTGTTGGAGGGGCGGGTAGTGTGAAACACAGAGGCCGTACAAGGCCCCTGTGATTTTTTCGTATTAAGGACTCCCTCCATTGCAGAATGTTCACCTAGGACATCTGAAAGCATCGATCTAATCCGGAAGTGCTCACCACAGTGGACGGCACGTGCCCGTCAGTCATGCCGCAAACGTTGTGTCGAGTAGCTGGGGAATGCCAGTAA"

vk="PAWHEAE"
wk="HEAGAWGHEE"
ro = 2
(b, row, col, matrix, r_s, r_e, c_s, c_e) =  overlapAlignmentGraph(v, w, ro)

print "S", b
print "STARTEND", int (r_s), int(r_e), int (c_s), int (c_e)
print "MOD", b[r_s, r_e], b[c_s, c_e]

print ''.join (backTrackLocalScoreMatrix(b, v, w, r_s, r_e, c_s, c_e, ro, row, col, matrix))
