__author__ = 'ilap'

from SequenceCompareLib import *

# Main
############################################################################
# Ch 5.8 Step 7
############################################################################
'''
CODE CHALLENGE: Solve the Longest Path in a DAG Problem.
     Input: An integer representing the source node of a graph, followed by an integer representing the
     sink node of the graph, followed by a list of edges in the graph. The edge notation 0->1:7 indicates
     that an edge connects node 0 to node 1 with weight 7.
     Output: The length of a longest path in the graph, followed by a longest path.

Sample Input:
     0
     4
     0->1:7
     0->2:4
     2->3:2
     1->4:1
     3->4:3
'''
source = 0
sink = 4
graph = {
    0: [(1,7), (2,4)],
    2: [(3,2)],
    1: [(4,1)],
    3: [(4,3)]
}

graph2 = {
    0: [1, 2],
    2: [3],
    1: [4],
    3: [4]
}

tstr = '''
2->21:1
11->20:33
11->21:9
19->27:23
2->24:32
22->35:32
18->28:33
13->25:12
10->25:28
7->19:12
8->28:35
21->27:6
16->30:1
5->26:1
25->26:31
5->28:3
0->4:23
27->32:17
29->31:0
12->29:36
2->26:36
2->3:2
19->34:20
18->24:37
9->22:14
8->19:20
4->14:7
6->26:1
14->29:3
7->28:32
21->28:3
4->25:18
2->34:15
7->29:37
20->35:38
15->32:15
23->34:3
6->35:6
6->31:19
3->20:28
13->16:31
11->35:32
2->12:4
11->32:34
5->10:3
9->18:20
14->21:18
2->19:1
10->12:1
26->36:22
16->27:15
8->29:3
4->6:22
16->24:8
16->22:32
20->26:33
21->35:16
28->31:38
28->30:5
21->24:24
17->35:0
0->12:38
0->32:26
17->32:25
16->18:23
23->31:25
1->25:28
0->22:22
1->21:1
6->29:24
3->12:32
8->17:34
8->14:17
7->30:32
15->16:4
20->25:36
'''

tstr2 = '''
0->1:7
0->2:4
2->3:2
1->4:1
3->4:3'''

tstr = '''
1->2:5
1->3:6
1->4:5
2->3:2
2->6:4
3->5:4
3->6:3
3->7:5
4->5:6
4->6:8
5->7:2
6->7:1
'''
s=1
e=

graph =  convertToTupleDAG (tstr)
print "GRAPH", graph

print backtrackDAG(s, e, graph)
#print "FP", find_all_paths(graph,0, 4)

#print backtrackDAG (source, sink, graph)