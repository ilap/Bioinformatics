__author__ = 'ilap'

import threading

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

def deBrujinReadPairs (read, k, d):

    (nodes, matrix) = pathGraphReadPairs(read, k, d)
    nlen = len (nodes)


    for node in nodes:

        tnodes = nodes[:]
        idx = nodes.index(node)

        for j in range (len (nodes)-1, idx, -1):
            if (node == nodes[j]): # merge
                print "GLUING", node, nodes[j]
                matrix[:, idx] = matrix[:, idx] + matrix[:, j]
                matrix[idx, :] = matrix[idx, :] + matrix[j, :]

                del tnodes[j]
                matrix = np.delete (matrix, j, 0)
                matrix = np.delete (matrix, j, 1)


        nodes = tnodes[:]

    return (nodes, matrix)


def deBrujinReadPairs (read, k, d):

    (nodes, matrix) = pathGraphReadPairs(read, k, d)
    nlen = len (nodes)


    for node in nodes:

        tnodes = nodes[:]
        idx = nodes.index(node)

        for j in range (len (nodes)-1, idx, -1):
            if (node == nodes[j]): # merge
                print "GLUING", node, nodes[j]
                matrix[:, idx] = matrix[:, idx] + matrix[:, j]
                matrix[idx, :] = matrix[idx, :] + matrix[j, :]

                del tnodes[j]
                matrix = np.delete (matrix, j, 0)
                matrix = np.delete (matrix, j, 1)


        nodes = tnodes[:]

    return (nodes, matrix)


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

def pairedComposition (pattern, k, d):

    plen = len (pattern)

    process_len = plen - (2*k+d)

    result = []
    print "P", process_len
    for i in range (0, process_len+1):
        #val = "("+pattern[i:i+k] + "|" + pattern[i+k+d:i+2*k+d] +")"
        val = [pattern[i:i+k],pattern[i+k+d:i+2*k+d]]
        result.append(val)

    #result  = ' '.join (sorted (result))

    #return sorted (result)
    return result

def printReadPairs (read_pairs):

    result = []
    for read_pair in read_pairs:
        result.append("("+read_pair[0] +"|"+read_pair[1]+")")

    return ' '.join (result)

def prefixReadPair (read_pair):

    klen = len (read_pair[0])

    result = [read_pair[0][0:klen-1], read_pair[1][0:klen-1]]
    return result

def suffixReadPair (read_pair):

    klen = len (read_pair[0])

    result = [read_pair[0][1:klen], read_pair[1][1:klen]]
    return result

import numpy as np
def pathGraphReadPairs (text, k, d):

    read_pairs = pairedComposition(text, k, d)
    mlen = len (read_pairs)+1
    matrix = np.array ([[0]*mlen]*mlen, np.uint32)
    nodes = []

    for i in range (0, mlen - 1):

        kmer = prefixReadPair(read_pairs[i])
        nodes.append (kmer)
        matrix[i][i+1] = 1

    nodes.append(suffixReadPair(read_pairs[len (read_pairs)-1]))
    return (nodes, matrix)

def composition (read, k):

    result = []
    for i in range (0, len (read) - k + 1):
        kmer = read[i:i+k]
        result.append (kmer)

    return sorted (result)

def compositionReadPairsGraph (text, k, d):
    graph = {}


    (nodes, matrix) = pathGraphReadPairs(text, k, d)
    nlen = len (nodes)

    for i in range(0, nlen):
        values = matrix[i,:]
        arr = values.tolist()
        if 1 in arr:
            idx = arr.index (1)
            graph[''.join(nodes[i])] = ''.join(nodes[idx])
        #idx = values.index


    print graph
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




def getUnbalancedNode (graph):

    from_node = -1
    to_node = -1

    for key in graph.keys ():
        for node in graph[key]:
            if not graph.has_key(node):
                from_node = node
                graph[from_node] = []
                break

    mlen = len (graph)
    print mlen, from_node
    matrix = np.array ([[0]*mlen]*mlen, np.uint32)


    for key in graph.keys():
        nodes = graph[key]
        key_idx = graph.keys().index(key)
        for node in nodes:
            if graph.has_key(node):
                node_idx = graph.keys().index(node)
                print key_idx, node_idx, node
                matrix[key_idx][node_idx] += 1

    print graph.keys ()
    print matrix
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

'''
String Reconstruction from Read-Pairs Problem: Reconstruct a string from its paired composition.
     Input: A collection of paired k-mers PairedReads and an integer d.
     Output: A string Text with (k,d)-mer composition equal to PairedReads (if such a
     string exists).
'''
def reconstructFromReadPairs (paired_kmers, d):

    plen = len (paired_kmers)
    k = (len(paired_kmers[0])-1)/2
    result = [None]*(plen+2*k+d)
    rlen = len (result)

    print "result", k, d, rlen, result
    for i in range (0, plen-1):
        f= paired_kmers[i]
        t = paired_kmers[i+1]

        pair1= f[0:k] + t[0:k]
        pair2= f[k+1:2*k+1] + t[k+1:2*k+1]
        #print "PAURS:", list(pair1), pair2
        result[i:2*k] = list(pair1)
        result[i+2*k+d:i+4*k+d] = list(pair2)
        print "RES", pair1, pair2, result




    print result




    return result

def printGraph (graph):
    for key in sorted(graph):
        to = graph[key]
        #print "TO", to
        if len (to) != 0:
            print "%s -> %s" % (key, ','.join(to))


## Main

READ_PAIRS= ''.join ('''
GAGA|TTGA
TCGT|GATG
CGTG|ATGT
TGGT|TGAG
GTGA|TGTT
GTGG|GTGA
TGAG|GTTG
GGTC|GAGA
GTCG|AGAT
''').split ()

RP_GRAPHS= {
"A|A": ["G|G"],
"A|T": ["G|G"],
"G|G": ["C|C","C|C","C|C"],
"C|C": ["T|A", "A|T","T|T"],
"T|T": ["G|G"],
}

RP_GRAPHS2={
    0: [2],
    1: [3],
    2: [1],
    3: [0,4],
    6: [3,7],
    7: [8],
    8: [9],
    9: [6],
}

########
print READ_PAIRS
READS=''.join ('''
GAGA|TTGA
TCGT|GATG
CGTG|ATGT
TGGT|TGAG
GTGA|TGTT
GTGG|GTGA
TGAG|GTTG
GGTC|GAGA
GTCG|AGAT
''').split()
k=4
d=2

print READS
READS=   [['GAGA','TTGA'], ['TCGT','GATG'], ['CGTG','ATGT'], ['TGGT','TGAG'], ['GTGA','TGTT'], ['GTGG','GTGA'], ['TGAG','GTTG'], ['GGTC','GAGA'], ['GTCG','AGAT']]


RP_GRAPHS = deBrujinReadPairs(READS, k, d)

exit ()#######

(from_path, to_path, g) =  getUnbalancedNode(RP_GRAPHS)

print "G", from_path, to_path, g

cy = eucladianPath(g, from_path)

idx = cy.index(to_path)
cy = cy[idx:len(cy)] + cy[0:idx]

reconstructFromReadPairs(cy, d)

print '->'.join ( str (v) for v in cy)



'''STR = "TAATGCCATGGGATGTT"

read_pairs =  pairedComposition(STR, 3, 1)

print printReadPairs(read_pairs)

print pathGraphReadPairs("TAATGCCATGGGATGTT", 3, 1)

print prefixReadPair(["TAA", "GCC"])
print suffixReadPair(["TAA", "GCC"])
print "DEB"
print deBrujinReadPairs("TAATGCCATGGGATGTT", 3, 1)
#print compositionReadPairsGraph("TAATGCCATGGGATGTT", 3, 1)
'''