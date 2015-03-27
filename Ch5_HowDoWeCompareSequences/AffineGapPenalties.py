__author__ = 'ilap'

from SequenceCompareLib import *
import numpy as np

np.set_printoptions(edgeitems=20)
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=300)
def affineGapPenalties (v, w, ro, ep):

    v_len = len (v)+1
    w_len = len (w)+1



    middle_diag = np.zeros (dtype=int, shape=(v_len, w_len))
    upper_right = np.zeros (dtype=int,shape=(v_len, w_len))
    lower_down = np.zeros (dtype=int, shape=(v_len, w_len))

    profile = ['']*(v_len*w_len)
    bm = np.array (profile)
    bm = bm.reshape(v_len,w_len)

    bu = np.array (profile)
    bu = bu.reshape(v_len,w_len)
    bu[0,:] = "x"
    bu[0,0] = "d"


    bl = np.array (profile)
    bl = bl.reshape(v_len,w_len)
    bl[:,0] = "y"
    bl[0,0] = "d"

    for i in range (1, v_len):
        lower_down[i,0] = -(ro + i-1)

    for i in range (1, w_len):
        upper_right[0,i] = -(ro + i-1)



    middle_diag[0,:] = -1000
    middle_diag[:,0] = -1000
    middle_diag[0,0] = 0

    upper_right[:,0] = -1000
    lower_down[0,:] = -1000


    '''print middle_diag
    print upper_right
    print lower_down'''

    (row, col, matrix) = generateScoreMatrixFromFile("BLOSUM62.txt")

    '''matrix[:,:] = 0
    for i in range (len(matrix)):
        matrix[i,i] = 1'''

    #print matrix
    for i in range (1, v_len):
        for j in range (1, w_len):
            v_char = v[i-1]
            w_char = w[j-1]
            v_idx = row.index(v_char)
            w_idx = col.index(w_char)


            upper_right[i, j] = max (upper_right[i,j-1] - ep, middle_diag[i,j-1]-ro)
            if upper_right[i, j] == middle_diag[i,j-1]-ro: bu[i,j-1] = "d"
            else: bu[i,j-1] = "x"

            lower_down[i, j] = max (lower_down[i-1,j] - ep, middle_diag[i-1,j]-ro)
            if lower_down[i, j] == middle_diag[i-1,j]-ro: bl[i-1,j] = "d"
            else: bl[i-1,j] = "y"

            middle_diag[i, j] = max (lower_down[i-1,j-1],middle_diag[i-1,j-1], upper_right[i-1,j-1])+matrix[v_idx,w_idx]
            if middle_diag[i, j] == (middle_diag[i-1,j-1] + matrix[v_idx,w_idx]): bm[i-1,j-1] = "d"
            if middle_diag[i, j] == (lower_down[i-1,j-1] + matrix[v_idx,w_idx]): bm[i-1,j-1] = "y"
            if middle_diag[i, j] == (upper_right[i-1,j-1] + matrix[v_idx,w_idx]): bm[i-1,j-1] = "x"




    score = int (max (middle_diag[i,j],lower_down[i,j],upper_right[i,j]))

    v_str = ""
    w_str = ""

    print
    print middle_diag
    print bm
    print
    print upper_right
    print bu
    print
    print lower_down
    print bl

    print
    print

    id_str = 0
    tstr = ""
    if middle_diag[i,j] == score: t = "d"; B=bm
    if upper_right[i,j] == score: t = "x"; B=bu
    if lower_down[i,j] == score: t = "y"; B=bl


    while i and j:
            v_char = v[i-1]
            w_char = w[j-1]
            v_idx = row.index(v_char)
            w_idx = col.index(w_char)


            m = middle_diag[i,j]
            l = lower_down[i,j]
            u = upper_right[i,j]

            max_v = max (m,u,l)
            #print "MAX", max_v, m, u, l

            v_chr = v[i-1]
            w_chr = w[j-1]
            '''
            print "IJ", id_str%10, i-1, j-1, v_chr, w_chr,
            if max_v == u:
                print "right", m, u, l

                v_chr = "-"
                j -= 1
            elif max_v == l:
                print "down", m, u, l
                w_chr = "-"
                i -= 1
            else:
                print "diag", m, u, l
                i -= 1
                j -= 1

            '''
            v_chr = v[i-1]
            w_chr = w[j-1]
            print "IJ", id_str%10, i-1, j-1, v_chr, w_chr,

            if t == "d":
                i -= 1
                j -= 1
                B=bm
                t = B[i,j]

            elif t == "y":
                w_chr = "-"
                i -= 1
                B=bl
                t = B[i,j]
            elif t == "x":
                v_chr = "-"
                j -= 1
                B=bu
                t = B[i,j]
            print t
            v_str = v_chr + v_str
            w_str = w_chr + w_str
            tstr = tstr + str (id_str%10)
            id_str += 1


    print "SCORE", score
    print
    print v_str
    print w_str

    return

### MAINC

v="PRTEINS"
w="PRTWPSEIN"
v="YHFDVPDCWAHRYWVENPQAIAQMEQICFNWFPSMMMKQPHVFKVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE"
w="YHEDVAHEDAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPIISATCARMRVRTVWE"

wl="RTFFVGGNFKLNTASIPENVEVVICPPATYLDY"
vl="RKFFVGGNWKMNGDKKSLNG"

v="SWNASTYQSMMVWSFTCQVAMDEMTQRHPRREHCVLVPDDYFILVETIIYLVYDKERHEAIQWWAKVTLDKCFRAEDYHKQ"
w="SWNASTYQSMMVWSFSYEAMCQVAMDEMTQRHPRRPHDYFILVETIMLLVYDKERHRAIQWWAKVTLDKCFRAEDYHKQ"

ro = 11
ep = 1

affineGapPenalties(v, w, ro, ep)


#print backtrackAffineGap (backtrack, v, w, len (v), len (w))