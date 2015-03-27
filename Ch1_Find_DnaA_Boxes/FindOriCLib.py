__author__ = 'ilap'

import os
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
