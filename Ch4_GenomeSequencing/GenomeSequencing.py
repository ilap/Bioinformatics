__author__ = 'ilap'


import sys
sys.setrecursionlimit(1000000000)

import threading
def composition (read, k):

    result = []
    for i in range (0, len (read) - k + 1):
        kmer = read[i:i+k]
        result.append (kmer)

    return sorted (result)

'''
String Spelled by a Genome Path Problem. Reconstruct a string from its genome path.
     Input: A sequence of k-mers Pattern1, ... ,Patternn such that the last k - 1 symbols of Patterni are
            equal to the first k-1 symbols of Patterni + 1 for 1 <= i <= n-1.
     Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Patterni  (for 1 <= i <= n).

'''
def genomePathProblem (reads):
    genome = ""
    klen = len (reads[0])
    rlen = len (reads)

    tarr = reads[:]
    for i in range (0, rlen-1):

        suffix = reads[i][1:klen]
        prefix = reads[i+1][0:klen-1]

        if (prefix == suffix):
            genome += reads[i][0]
    genome += reads[i+1]

    return genome

'''
Overlap Graph Problem: Construct the overlap graph of a collection of k-mers.
     Input: A collection Patterns of k-mers.
     Output: The overlap graph Overlap(Patterns).
'''
import collections
def overlapGraphProblem (reads):

    rlen = len (reads)
    klen = len (reads[0])
    idx = 0
    graph = { }
    for kmer in reads:
        graph[kmer] = []
        suffix = kmer[1:klen]
        for  i in range (0, rlen):
            if (i == idx): continue # skip the actual one.
            prefix = reads[i][0:klen-1]
            if ( prefix == suffix):
                graph[kmer].append (reads[i])

        idx += 1
    return graph

'''
Hamiltonian Path Problem: Construct a Hamiltonian path in a graph.
     Input: A directed graph.
     Output: A path visiting every node in the graph exactly once (if such a path exists).
'''

def hamiltonianPathProblem (reads):
    genome=""
    return genome

import numpy as np

def pathGraph (read, k):

    mlen = len (read) - (k-2)
    matrix = np.array ([[0]*mlen]*mlen, np.uint32)
    nodes = []

    for i in range (0, mlen - 1):

        kmer = read[i:i+k-1]

        nodes.append (kmer)
        matrix[i][i+1] = 1

    nodes.append(read[mlen-1:mlen-1+k])

    return (nodes, matrix)


def deBrujin (read, k):

    (nodes, matrix) = pathGraph(read, k)
    nlen = len (nodes)


    for node in nodes:

        tnodes = nodes[:]
        idx = nodes.index(node)

        for j in range (len (nodes)-1, idx, -1):
            if (node == nodes[j]): # merge

                matrix[:, idx] = matrix[:, idx] + matrix[:, j]
                matrix[idx, :] = matrix[idx, :] + matrix[j, :]

                del tnodes[j]
                matrix = np.delete (matrix, j, 0)
                matrix = np.delete (matrix, j, 1)


        nodes = tnodes[:]

    return (nodes, matrix)


def printDeBurjin (nodes, matrix):

    result = []
    for i in range (0, len (nodes)):

        indexes = [k for k,x in enumerate(matrix[i,:]) if x != 0]
        res =  nodes[i] + " -> "
        ra = []
        str = ""
        for j in indexes:
            num = matrix[i,j]
            for k in range (0, num):
                ra.append(nodes[j])
            str = res +  ','.join (ra)
        result.append (str)

    print '\n'.join(sorted (result))


def deBrujinPatterns (reads):

    reads = sorted (reads)
    plen = len (reads)
    klen = len (reads[0])
    graph = {}
    for kmer in reads:

        prefix = kmer[0:klen-1]
        suffix = kmer[1:klen]
        #print "KMER", kmer, prefix, suffix

        if prefix not in graph:
            graph[prefix] = []
        #print "APP", prefix, suffix
        graph[prefix].append (suffix)

    #print graph
    return graph

def printGraph (graph):
    for key in sorted(graph):
        to = graph[key]
        #print "TO", to
        if len (to) != 0:
            print "%s -> %s" % (key, ','.join(to))

def findPath(graph, start, end, path=[]):
        print "START", start
        path = path + [start]
        if start == end:
            return path
        if not graph.has_key(start):
            return None
        for node in graph[start]:
            print "NODE", node
            if node not in path:
                newpath = findPath (graph, node, end, path)
                if newpath: return newpath
        return None

def deletePath (graph, from_node, to_node=-1):
    lock = threading.Lock()

    lock.acquire() # will block if lock is already held
    if from_node in graph:
        nodes = graph[from_node]
        if len (nodes) > 1:
            del graph[from_node][0]
        else:
            del graph[from_node]
    lock.release()
    return graph

def walkPath (graph, start, unused_nodes, path=[]):

    if path:
        end = path[0]
    else:
        end = -1

    if start == end:
        return path

    path = path + [start]

    if not graph.has_key(start):
        return None

    for node in graph[start]:
        graph = deletePath(graph,start)

        lock = threading.Lock ()
        lock.acquire ()
        if graph.has_key(start):
            if not start in unused_nodes:
                unused_nodes.append (start)
        elif start in unused_nodes:
            unused_nodes.remove(start)
        lock.release ()

        newpath = walkPath (graph, node, unused_nodes, path)

        if newpath:
            return newpath

    return None

import random
from random import randint

def eucladianCycle (graph, start_node=-1):

    import copy
    unused_edges = copy.deepcopy(graph)
    unused_nodes = []

    #[start] = graph[start_node]

    cycle = []

    if start_node == -1:
        start = randint (0, len (unused_edges)-1)
    else:
        start = start_node

    while unused_edges:
        temp_cycle = walkPath(unused_edges, start, unused_nodes)

        if unused_nodes:
            start = unused_nodes[0]

        temp_cycle += cycle
        idx = temp_cycle.index(start)
        temp_cycle = temp_cycle[idx:]+temp_cycle[0:idx]

        cycle = temp_cycle
    cycle.append(cycle[0])
    #print "CY", cycle


    return cycle

import random
from random import randint

def eucladianPath (graph, start_node=-1):

    import copy
    unused_edges = copy.deepcopy(graph)
    unused_nodes = []

    #[start] = graph[start_node]

    cycle = []

    if start_node == -1:
        start = randint (0, len (unused_edges)-1)
    else:
        start = start_node

    while unused_edges:
        temp_cycle = walkPath(unused_edges, start, unused_nodes)

        if unused_nodes:
            start = unused_nodes[0]

        temp_cycle += cycle
        idx = temp_cycle.index(start)
        temp_cycle = temp_cycle[idx:]+temp_cycle[0:idx]

        cycle = temp_cycle
    #cycle.append(cycle[0])
    #print "CY", cycle


    return cycle

def eucladianPathPatterns (graph, start_node):

    import copy
    unused_edges = copy.deepcopy(graph)
    unused_nodes = []

    #[start] = graph[start_node]

    cycle = []


    while unused_edges:
        temp_cycle = walkPath(unused_edges, start, unused_nodes)

        if unused_nodes:
            start = unused_nodes[0]

        temp_cycle += cycle
        idx = temp_cycle.index(start)
        temp_cycle = temp_cycle[idx:]+temp_cycle[0:idx]

        cycle = temp_cycle
    #cycle.append(cycle[0])
    #print "CY", cycle


    return cycle

import random
from random import randint

import numpy as np


def getUnbalancedNode (graph):

    from_node = -1
    to_node = -1

    for key in graph.keys ():
        for node in graph[key]:
            if not graph.has_key(node):
                from_node = node
                break

    if from_node != -1:
        print "from is -1"
        # increase the graph
        graph[from_node] = []

    mlen = len (graph)
    #print mlen, from_node
    matrix = np.array ([[0]*mlen]*mlen, np.uint32)


    for key in graph.keys():
        nodes = graph[key]
        key_idx = graph.keys().index(key)
        for node in nodes:
            if graph.has_key(node):
                node_idx = graph.keys().index(node)
                #print key_idx, node_idx, node
                matrix[key_idx][node_idx] = 1

    for i in range (0, mlen-1):
        matrix[i,i] = 0
        col_count = np.sum(matrix[:,i])
        row_count = np.sum(matrix[i,:])
        k = graph.keys()[i]
        if row_count > col_count:
            to_node = k
        elif row_count < col_count:
            from_node = k

    graph[from_node].append (to_node)

    #print graph


    print "FROM TO:", from_node, to_node


    return (from_node, to_node, graph)


def deletePath2 (graph, from_node, to_node=-1):
    lock = threading.Lock()

    lock.acquire() # will block if lock is already held
    if from_node in graph:
        nodes = graph[from_node]
        if len (nodes) > 1:
            del graph[from_node][nodes.index (to_node)]
        else:
            del graph[from_node]
    lock.release()
    return graph

def walkThrough (graph, node_arr):

    print "AAAAAAAAA", node_arr
    unused_edges = graph.copy()
    arr = node_arr[:]
    from_node = arr[0]
    arr = arr[1:len(arr)]

    print "ARR", arr
    for node in arr:
        to_node = node
        graph = deletePath2(unused_edges, from_node, to_node)
        print "REMOVE", node, unused_edges, from_node, to_node

        from_node = to_node

def generateAllBinKmers (pattern, k):

    if k == 0: return pattern

    result = ""

    for N in ["0", "1"]:
        tval = generateAllBinKmers(pattern + N, k-1)
        if k == 1:
            result += tval + " "
        else:
            result += tval

    return result


###################################################################################
#### MAIN
###################################################################################
READS=''.join ('''
000
001
010
011
100
101
110
111
''').split ()
READS= ''.join (generateAllBinKmers("", 9)).split ()
graph = deBrujinPatterns(READS)


print "GRAPH", graph
#(from_path, to_path, g) =  getUnbalancedNode(graph)
g = graph
print g
cy = eucladianCycle(g, "00001010")

'''cy = eucladianPath(g, to_path)
idx = cy.index(to_path)
cy = cy[idx:len(cy)] + cy[0:idx]'''

rest = ""
for i in range (0, len(cy)-1):
    rest += cy[i][0]
#rest += cy[len(cy)-1]
print "REST", rest

#### 4.9
EUCL_GRAPH={
    0: [2],
    1: [3],
    2: [1],
    3: [0,4],
    6: [3,7],
    7: [8],
    8: [9],
    9: [6],
}

(from_path, to_path, g) =  getUnbalancedNode(EUCL_GRAPH)

print "AAAAA"

cy = eucladianPath(g, from_path)
idx = cy.index(to_path)
cy = cy[idx:len(cy)] + cy[0:idx]

print '->'.join ( str (v) for v in cy)
'''


### 4.7 Eucladian....
EUCL_GRAPH={
    0: [3],
    1: [0],
    2: [6,1],
    3: [2],
    4: [2],
    5: [4],
    6: [5,8],
    7: [9],
    8: [7],
    9: [6],
}

cy = eucladianCycle(EUCL_GRAPH)
print '->'.join ( str (v) for v in cy)
#printEucledianGraph('''

## 4.5

#READS=''.join ('''GAGG
#CAGG
#GGGG
#GGGA
#CAGG
#AGGG
#GGAG''').split ()

#printGraph(deBrujinPatterns(READS))

#READ = "TAATGCCATGGGATGTT"

#print composition(READ, 3)
#(N, M) = deBrujin(READ, 3)
#printDeBurjin(N, M)


###################################################################################
## 4.4

#print "FI", x,  x.filled(1)

#READS=pathGraph("TAATGCCATGGGATGTT", 3) #TAATGCCATGGGATGTT
#print READS

#(N, M) = deBrujin("TAATGCCATGGGATGTT", 4)
#printDeBurjin(N, M)
#(N, M) = deBrujin("TAATGCCATGGGATGTT", 3)
#printDeBurjin(N, M)
#exit()

###################################################################################
### 4.3
#READS= ''.join ('''ATGCG
#GCATG
#CATGC
#AGGCA
#GGCAT''').split ()
##### 4.


#graph = overlapGraphProblem(READS)
'''graph = overlapGraphProblem(READS)

for key in sorted(graph):
    to = graph[key]
    if len (to) != 0:
        print "%s -> %s" % (key, to[0])

mask = np.ones(10, dtype=bool)
print mask
mask[[0,2,4]] = False
print mask
#http://docs.scipy.org/doc/numpy/reference/generated/numpy.delete.html
#result = arr[mask,...]
'''
###################################################################################
### 4.3.3


#print READS
#print genomePathProblem (READS)


###################################################################################
### 4.1
#read = "TGCTGTAGCTTCGCATACTTAGCTAAGTTATACGCAATGACTGGAAAAAGTGCAACTCTCGACTGAGCACATGGTTCCCGTCACTACAAATTTATCCTGTGACCAGACTCCGGCTTAATTTACCTTGTACCACTTGCCCAGTTTGAGTAAGCCAAATACACTGGAACGTACGTGCACTCCCCGTCTTCTGGAAAAGGTGGTTCGTCTTCGATTACGGTACCCAGTTATATTGATTGAATGCCGCGCGTGTGCCAACCACCTTGATTGAATAGCCTACGGTAATTGCTTAATGCAGTGCGAAGTGAGTACCCGTATTCCCGTCCGGGAACCGGCGGTTCGCAGTTCTCGGGGGCATGTCATAAATCGGATCTCAGGTCTAATTTAATAATAGTCAATGGAGGGCGACCTTATGGCGTTGAGCCACCACCAGGATACTTGACATCTCCGTGTATACTACAGATGAGAACCTAAGGCCGAAGCGCGCTACCGCCAGCCCAAGGTAAAACTTGATTATAGAGTAGCGTTCAGATTTTACTACTTTAACATATGGATGAGAGACCCGCGGGATTTGTTTGGTAATCAGCGTATACGTACGAAACTGGAGCCTACGGCGATCTTAATACCCCGGACAGAATGGGTTCAGTCGTCTCATTTCGGACTACGGAAAGATGTAAGTGTGCCCGACCGTACTTGGACAAGACTTTGGCTGCTTATCAGACAAGGACCTTTCCGGAGTAATAGCGACGCCAGGGTATGATCATATGCTATGAGACTATCAATGGGACTAGCATACTTGTACGCGGTACCCCGGGATAGGCACCTGGGGTGCGGGAAATGACATGACTAGCATCCCGGTGTCCACGGAGTGGGCCGCAGAAAACTCCGCTTAGCTTACGTCCAGGACCAACTGTAAGATATGCGTAGGACCTGCCAAGCACAGTAATGCTGCGTGCAGTTTGTTCTTGACCCAAAATCGACCCGGAAAAGCTGCTGCAGCATCATCATATCCACCTCGTTGGGTACTGCTGGCCCTCTAGTTGGCTCCGGTATGAGAGCCTTATACTGGTTTTTCTCCTTGTCAAGTCCCTGGTCTCATCGCTGAAAAAGACACTTAGGCCTCCTTATAGGTCTCTAGACAACCCTCTTTCTTCTGTAGACTGGAATCTCCGCATCAGCCCGACTTGCGCCTAGATTCTAAGATGAAGAAGGCGATACTAATGACTTCAATATATCATTAATCGCATAGAAACTTATGATATGTACGCTATTAAGGCGCCGTAATCCACCAGGCTTGGGGTCCGAGTAGTCTATTCCGTCTTGAAATTGGGTGCGGTCTTGATCGTGACTATTAATAGCTTCCTTGCGTTCCCCTGGCCTAATACGATCTTGAGCCACAGTGTCTTTGGCTGAGTCAATACTGGAACTGCGCCGAGGTGTACAACATAACAATTATGCTGCGAAGCGATTACAGTCCCGAGATTCAAGGAGCCCATAAAAAGCGGACGGGTAAATAGTTCTGGAGAGATGCCAAACGTACGGTCCGAATTGACTTGCGGCCACTAAGGTTCGGTAGCACCTGTCAATTAGGGATCACCTAAAGTGCGTTGCTCCACTTTAGGAACGGCACTTTGACATTCGCCCCCGGGCGTCTGAAACGGTTTTGCCACATCTAACGAGGCGTACGTAAATACGGCGTTCTGATAACCTATGTAGAATTACAAGTGCCCACTACATCACGAATTCGGGTCCCTTGTGCATAAAGTAAGCTTTTTTAGCTGGGGTTATCGGCGAAAGTAAAAGGCAATCTGAGAGACGATTGTCTGTCCCTACGTGACCGGGTCGGACATACTAGGAAGAATGGACCTACAGACCACCCGTCGTCGCTTATATTGCTAGGGGCAAAATGATGGACGAGGTGGTTAATAAACGTACGTGCGAGTCGCGTTCAAACCTTTTTACCCGCTTTCACAAATATATTTGACTACGGCGTACCAACCTAACTATCTTTACTCCCTCATCGCAAGGGGAGTAGGACTTTTGTACACGGCGCCCTACTTATTCACCCGAACTACTCGTGATATCGTGCGAAGGAGTTCCGAACCATTCTATTGCAGGGGGGCTCTTCTCCTCGCAATCGCAGACCTTGGGCTACTCAATAGGCACCAAGTTTACTTGGACATTTTTGCAGTGTCTCCTACACATGTCGCTTTGGGGGCTCTGGTAAATCAGGTGTTGGTCGGGGGCACATATCGGGGGGAAAGTACATTCCTGCGTAGGTCGACCCTTGTCATGGGGGGGGTTCCACACTGTATCGCTCGGTGAACACAGAGTTGGGAACGAGGGTCACCGCTCGTAAAGACTTTTGCTTTACTAGCGCGGTTCGCGCTTGATGATACATAGATAAGCAGCGCAGCGATTCGGCAGGTATTCTAGCGTGCCTAGAACCGGCCGCCAGTTAGACATTTAACACAAGTGTCCTGCAATCACGAAGAGCTGAGGTGATCCGTCCTTCAGGCCTATCTAGGAAAGCATCCGGCTAGGAAAACTAATCGACACTATTATATACACATAACCCCTCGTGTCTCGCCCTATAGAAAACTGATAGAATAAATGGCTCAAAATAGTCACGTGTTAGCCCTAGTAGTAGGCGGCCCTATCGGCACAATGCTTCTGTGCGACGACTTGTGGTGAAACCCTATTAACTCAACCGGGGATCGTAATAACGTGTGAGCCCCTAGGAGGATGGAACGCTTAAGACTGACAGCGGATGCCCAGGAATCCATCCAATACACGATGGATTTAGGGCCTGATCTAAGCTGCTGCTAGGACTGCGGCTGGTGTCTCGAGTCGCTAATCCGCCCTCACTTACGATACTGAAATGCGGCGCCCCAATTTGAAGCCCCCCATTTGTTGTGGACGAAAGTGTCCAATTCGAGGAACCGCGCTCGTATCGTGGCTATGACATTGTACGTCAGGTTCGGGGGGGCCTGGCGAAACAGAACCCGATCGCTCTCACCACGGAATGTACTTCGTATGAAGTAAGTCGTTATATAACCTCCCGCTCTGCAAGGAACGGCAGCATCTGGTTATCGGAACCTTATATCGAGCAACCGACAGAAGGCTACTGCAGGCCGCGACGGTACCCATCGTTGCCTCTGACTAATGGCGTCGGGTGTAACTCACTGTCAAAGACTGCATTGATAGAAATAAGCGGCGTGCATAAGCACGTCCTGTATGCCCTGCGCGAGTTCTGAACGCTTGCACTTTGACAGGGTAGGAGCAACGCTATGGAAATCGGGGTCTGAGGTACATTGACTGCCGATGAGCTCTTAGAGGGATAAAAGCGCTGGTGAGGATCCTTAATCCGAAGGCCACCACCAGAATGCTGCCCGAAAGGATATCCTCGCGGTAAGCAGGTGGACTGGCGCAGGAGGAGTAGAGATGACGCCTTGTATTATACCTTCTAACGCAAATTAGTAGCCCGAAGAATGGCCGACCATAGTCGTCGGATACGATGGATTCCCGTTTTCCTACCGCATCGCCGTAGCCCCGCAAGACGCTCGGTTCTATGGAGGCCTCCCTGGTAACCTGCTGTTTGATAAGAGACGGACTGACATACAGATTCACTGCAATGGTTGCTTGCATGGTGCCATCCTGCGCGGGATGTCGTATCAGAACCTCGACGGACGGGAATAACTTGCGGTCGCGAACCAGGCACATTCAAATCTATCCCTGATCTGTTAGGCGTACCTCTATTCTCCCCTTTTTACATTTCACGGAAAAACCAACCCAAGATGCATGAGTTGTTAACTGCCTCCTGTTTTCAAGTAGTGGGTAAGATACATACCCGATCTCGAAAGCACATGTAGATACACGGGCCTATTGTCGCTTTTCGTCAGCCCTCTCACGATATTGGCCTACGAAACAGGTCTGTAGATGCACTGATGAAAGCCCAGAGAGTCGACCGCCCCGCAGAGCGCAATCGAGACTTTACAAGATATCTTGGTTGTTAAGTAAATAAGGTTCCCTTCCCGGAGCGACACGGTCTAGAGGCGGGCCCCATACATTCGGATATCCATTCTGGCGCTCTCCACAATTTTCCTGTTTATTCCCATTCTGCCCCTTTACGGGCGGCCGTCGTTCCGGGGCTAAGTAAGATTCCTACGAAGCCCAGCGAGGAGAGACTGTCTCGTCGCTCGACTGGTAAACTGATGCGGGGCTAGAAATGATGCCTGCCTGGCCGGTCTGAATTGCCCACTAGCCGCAGGCATGTAGATGGAGGGTGTTGCCCGCAATGCCCCGAGCTTAATCCAACTGCTCTTGCAGCGGCGTCTCTAAAAGATTGCAGAGTTAAATCAGCCATAGAGATTGTTATCCGAGCCGCCAGTGGGTAGGTGATCTGCAAAACATAGGTCAATCAACCACAGCTATACTTATATTGAAACAGATGGCGCTCGATGCTCAATTACAGGGGGTGTTTCGCGCTTGCTGTCCCCTAAAGAACTTTCCCTCGCTTGCCTACCCCGTCCCTAGACTAACTATTATCTGAGATGTGACTGGTTGATGGCTAGCCCTTCCGAAGACACCCAGATAAGATCCGGACGTAATGAGCGACACCTTACCTTTCCATGTGAACCCGCTACGGCCGGTTTTTCGCCAACCTACGGCCGTCAGGCTCGTGCATCACGGTTCTAGTGTCTCAAAAATGCACGTAGGCGTTGTGGTAGTGGATGACGGCTTGTCGTTGAGGTGTCCATTTCTACTCTATCTGACCCACGTGGTATTCATGATGCAAAGGCATGCGCTGTAATATCTACTCCGGAATTCGAGTGCTTAATTGCTAATATCGCGCTTAAACCCTGTTGTGGGTAAACGAACACACAGTAGCTGGCTAAGCACAGCGGTGGTTAAACGGCGCATTTCCGAGCCAGACCCGCTTTCAAAGGGCAGCAGAGATGGGGTGGCAATAGGAGGATGCGGAGCCTAGGCCAGTCTTCGGT"

#read = "CAATCCAAC"
#print '\n'.join (composition(read, 100))


