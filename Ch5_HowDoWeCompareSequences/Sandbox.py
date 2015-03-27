__author__ = 'ilap'
from SequenceCompareLib import *


line= '''
2 -> 3:4, 2:1
'''

line1= '''
2 -> 3, 2
'''
tup = re.findall(r'(\d+):(\d+)', line)

print "TUP", tup
print [ (int (x),int(y)) for (x,y) in tup]
print "F3", re.findall(r'\d+', line1)


print '->'.join (str (x) for x in [0,1,2,3])


str = '''
0 -> 35,40,9
10 -> 20
11 -> 16,22
12 -> 31
13 -> 20,35,38
14 -> 27
15 -> 35,42
16 -> 29
17 -> 28,46
18 -> 31
2 -> 27
21 -> 44,46
22 -> 39
24 -> 39,41
25 -> 36
27 -> 34
28 -> 42
29 -> 46
3 -> 32,43
32 -> 45,46
35 -> 46
36 -> 38,42
38 -> 39,47
4 -> 23,38
40 -> 45
45 -> 48
5 -> 26,29,39
6 -> 38,45
7 -> 19,41
8 -> 14,18
9 -> 10
'''

str = '''
1 -> 2, 3, 4, 5, 6
2 -> 3, 6
3 -> 4
5 -> 4, 6
'''
graph =  convertToTupleDAG (str)

print "GRAPH", graph

print topologicalOrder(graph)