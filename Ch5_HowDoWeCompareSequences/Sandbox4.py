__author__ = 'ilap'
from SequenceCompareLib import *
import math

def logScore (v, w, u):

    if v == w == u:
        return 1

    symbols = ["-", "A", "C", "G", "T"]

    ent = 0.
    s_arr = [v, w, u]
    for s in symbols:
        c = s_arr.count(s)

        if c:
            val = c/3.
            print c, val
            ent += -(val * math.log (val,2))
    return ent


print logScore("A","A","G")

a = -(2 * math.log (2/3.,2))
print -a
print -(1 * math.log (1/3.,2))


longestCommonSubsequene ("TGTACG", "GCTAGT")