__author__ = 'ilap'

import math
import numpy as np



def score (v, w, u):

    return logScore(v, w, u)
def logScore (v, w, u):

    if v == w == u:
        return 1
    else:
        return 0
        if v == "-":
            if w == u: return 1
            else: return 0
        elif w == "-":
            if v == u: return 1
            else: return 0
        elif u == "-":
            if v == w: return 1
            else: return 0
        else:
            return 0

    symbols = ["-", "A", "C", "G", "T"]

    ent = 0.
    s_arr = [v, w, u]
    for s in symbols:
        c = s_arr.count(s)


        if c:
            val = c/3.
            ent -= (c * math.log (val,2))
    return ent

def multipleAlignment (v, w, u):

    v_len = len (v)+1
    w_len = len (w)+1
    u_len = len (u)+1



    s = np.zeros (shape=(v_len, w_len, u_len))
    s[:,:,:] = 0
    s[0,0,0] = 0

    for i in range (1, v_len):
        for j in range (1, w_len):
            for k in range (1, u_len):
                v_char = v[i-1]
                w_char = w[j-1]
                u_char = u[k-1]

                x = s[i-1,j,k] + score (v_char, "-", "-")
                y = s[i,j-1,k] + score ("-", w_char, "-")
                z = s[i,j,k-1] + score ("-", "-", u_char)

                xy = s[i-1,j-1,k] + score (v_char, w_char, "-")
                xz = s[i-1,j,k-1]+ score (v_char, "-", u_char)
                yz = s[i,j-1,k-1]+ score ("-", w_char, u_char)

                xyz = s[i-1,j-1,k-1]+ score (v_char, w_char, u_char)

                s[i,j,k] = max (x,y,z, xy,xz,yz, xyz)

                print "MAX", x,y,z, xy,xz,yz, xyz, "IJK", i, j, k, "MAAAAAXIZZZZZ:", s[i,j,k]
    return s


global v_align
global w_align

v_align = str ()
w_align = str ()

def backTrackMultipleScoreMatrix (s, v, w, u, i, j, k):

    global v_align
    global w_align

    v_len = len (v)
    w_len = len (w)
    u_len = len (u)

    min_v = min (v_len, w_len, u_len)

    if i == 0 or j == 0 or k == 0:
            print v_align
            print w_align
            return ""
    if min_v == v:
        if i == 0:
            print v_align
            print w_align
            return ""
    elif min_v == w:
        if j == 0:
            print v_align
            print w_align
            return ""
    elif k == 0:
        print v_align
        print w_align
        return ""

    v_char = v[i-1]
    w_char = w[j-1]
    u_char = u[k-1]

    x = s[i-1,j,k] + score (v_char, "-", "-")
    y = s[i,j-1,k] + score ("-", w_char, "-")
    z = s[i,j,k-1] + score ("-", "-", u_char)

    xy = s[i-1,j-1,k] + score (v_char, w_char, "-")
    xz = s[i-1,j,k-1]+ score (v_char, "-", u_char)
    yz = s[i,j-1,k-1]+ score ("-", w_char, u_char)

    xyz = s[i-1,j-1,k-1]+ score (v_char, w_char, u_char)

    max_v = max (x,y,z, xy,xz,yz, xyz)

    print "MAX", max_v, x,y,z, xy,xz,yz, xyz, "IJK", i, j, k, "MAAAAAXIZZZZZ:",

    if max_v == xyz:
        print "XYZ"
        v_align = v[i-1] + v_align
        w_align = w[j-1] + w_align
        return backTrackMultipleScoreMatrix(s, v, w, u, i-1, j-1, k-1) + u[k-1]
    elif max_v == yz:
        print "YZ"
        v_align = "-" + v_align
        w_align = w[j-1] + w_align
        return backTrackMultipleScoreMatrix(s, v, w, u, i, j-1, k-1) + u[k-1]
    elif max_v == xy:
        print "XY"
        v_align = v[i-1] + v_align
        w_align = w[j-1] + w_align
        return backTrackMultipleScoreMatrix(s, v, w, u, i-1, j-1, k) + "-"
    elif max_v == xz:
        print "XZ"
        v_align = v[i-1] + v_align
        w_align = "-" + w_align
        return backTrackMultipleScoreMatrix(s, v, w, u, i-1, j, k-1) + u[k-1]
    elif max_v == z:
        print "Z"
        v_align = "-" + v_align
        w_align = "-" + w_align
        return backTrackMultipleScoreMatrix(s, v, w, u, i, j, k-1) + u[k-1]
    elif max_v == y:
        print "Y"
        v_align = "-" + v_align
        w_align = w[j-1] + w_align
        return backTrackMultipleScoreMatrix(s, v, w, u, i, j-1, k) + "-"
    elif max_v == x:
        print "X"
        v_align = v[i-1] + v_align
        w_align = "-" + w_align
        return backTrackMultipleScoreMatrix(s, v, w, u, i-1, j, k) + "-"





def createScoreMatrix ():

    sm = s = np.zeros (shape=(21, 21, 21))


#### Main

vs=["ATGTTATA",
    "AGCGATCA",
    "ATCGTCTC"]

vs=["ATATCCG",
    "TCCGA",
    "ATGTACTG"]

vkk=["TA","CA","CG"]
v = vs[0]
w = vs[1]
u = vs[2]
v_len = len (v)
w_len = len (w)
u_len = len (u)

s =  multipleAlignment(v, w, u)

print int (s[v_len,w_len,u_len])
print s

print backTrackMultipleScoreMatrix(s, v, w, u, v_len, w_len, u_len)
