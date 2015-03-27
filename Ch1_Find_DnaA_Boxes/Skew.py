__author__ = 'ilap'

import sys

def skew (dna):
    dlen = len (dna)
    print "DNA length:", dlen
    skew_arr = [0]* (dlen+1)

    for i in range (1, len (dna)+1):
        if dna[i-1] == "G":
            skew_arr[i] = skew_arr[i-1] + 1
        elif dna[i-1] == "C":
            skew_arr[i] = skew_arr[i-1] - 1
        else:
            skew_arr[i] = skew_arr[i-1]

    return skew_arr

'''
Minimum Skew Problem: Find a position in a genome minimizing the skew.
Input: A DNA string Genome.
Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
'''
def minimumSkewProblem (data):
    dlen = len (data)

    mins = []
    skew_min = sys.maxint
    for i in range (1, dlen):
        if skew_min > data[i]:
            del mins[:]
            skew_min = data[i]
            mins.append (i)
            #print "MIN", i
        elif skew_min == data[i]:
            mins.append (i)

    return mins


gene_file = 'Salmonella_enterica_genome.txt'
dna = open(gene_file,'rU').read()
# pattern = "AAAAACGA"
#dna = "ACGTTGCATGTCGCATGATGCATGAGAGCT"

dna ="GATACACTTCCCAGTAGGTACTG"

data = skew (dna)

from pylab import *

#print "DATA: len", len (data)

x = arange(0, len (data), 1)
#print x
y = data
plot(x, y)

#xlabel('time (s)')
#ylabel('voltage (mV)')
#title('About as simple as it gets, folks')
#grid(True)
#savefig("test.png")


# print ' '.join(str(x) for x in data)
print ' '.join (str (i) for i in minimumSkewProblem (data))

show()