__author__ = 'ilap'

from FragileRegionsLib import *


PSTR="5 1"
P=getArrayFromTuppleString(PSTR)

N=chromosomeToCycle(P)
print "("+' '.join (str (x) for x in N)+")"


print cycleToChromosome(P)

'''

PSTR="(+1 -2 -3)(+4 +5 -6)"
PERMS=getPermsFromTuppleStrings(PSTR)

print coloredEdges(PERMS)

GRAPH=[[(2, 4), (3, 6), (5, 1)],[(7, 9), (10, 12), (11, 8)]]
#print graphToGenome(GRAPH)
'''
