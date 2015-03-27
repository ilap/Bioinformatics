__author__ = 'ilap'

from SequenceCompareLib import *


v="MEANLY"
w="PLEASANTLY"

#generateScoreMatrixFromFile("BLOSUM62.txt")
ro = 5
(b, row, col, matrix) =  globalAlignmentGraph(v, w, ro)
print int (b[len(v),len(w)])

print b
print ''.join (backTrackScoreMatrix(b, v, w, len(v), len(w), ro, row, col, matrix))
