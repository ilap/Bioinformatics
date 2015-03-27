__author__ = 'ilap'

from AssembleGenomeLib import *

#######

READS=''.join ('''
CAT|GGA
AAT|CCA
ATG|GAT
GGA|GTT
GGG|TGT
TAA|GCC
ATG|CAT
TGG|ATG
TGC|ATG
GCC|TGG
CCA|GGG
''').split()


klen=3
d=1

g =  deBrujinReadPairs(READS, klen, d)
#printGraph(g)
(end_node, start_node, g) =  getUnbalancedNode(g)

cy = eucladianPath(g, end_node, start_node)

print "RESULT", stringSpelledByGappedPatterns(cy, d)


