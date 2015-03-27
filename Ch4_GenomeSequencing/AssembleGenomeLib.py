__author__ = 'ilap'

## Imports
import random
import numpy as np
import threading
from random import randint
import os
#print "GETCWD", os.getcwd()

import sys
sys.setrecursionlimit(1000000000)

### Functions
def prefixReadPair (read_pair):

    klen = (len (read_pair)-1)/2

    result = read_pair[0:klen-1] +"|"+ read_pair[klen+1:len(read_pair)-1]
    return result

def suffixReadPair (read_pair):

    klen = (len (read_pair)-1)/2

    result = read_pair[1:klen] +"|"+ read_pair[klen+2:len(read_pair)]
    return result

def compositionReadPairsGraph (nodes, matrix):
    #XXXprint matrix
    graph = {}

    nlen = len (nodes)

    for i in range(0, nlen):
        values = matrix[i,:]
        arr = values.tolist()
        if 1L in arr:
            graph[nodes[i]] = []
        for j in range (0, nlen):
            if arr[j] > 0:
                #XXXprint "JJJJJJJJJ", j, arr[j]
                vnodes = []
                for k in range (0,arr[j]):
                    graph[nodes[i]].append (nodes[j])

    #XXXprint graph
    return graph


def pathGraphReadPairs (read_pairs, k, d):

    mlen = len (read_pairs)+1
    matrix = np.array ([[0]*mlen]*mlen, np.uint32)
    nodes = []

    #for i in range (0, mlen - 1):
    for node in read_pairs:
        prefix = prefixReadPair(node)
        suffix = suffixReadPair(node)

        if prefix not in nodes:
            nodes.append (prefix)

        if suffix not in nodes:
            nodes.append (suffix)

        pidx = nodes.index(prefix)
        sidx = nodes.index(suffix)

        #print prefix, suffix, pidx, sidx, nodes
        matrix[pidx][sidx] += 1

    #print "AAAA", len (nodes), nodes
    return (nodes, matrix)


def deBrujinReadPairs (read_pairs, k, d):

    (nodes, matrix) = pathGraphReadPairs(read_pairs, k, d)
    nlen = len (nodes)

    return compositionReadPairsGraph (nodes, matrix)

def getUnbalancedNode (graph):

    from_node = -1
    to_node = -1
    can_delete = False
    #print "MLEN", len (graph)
    for key in graph.keys ():
        for node in graph[key]:
            if not graph.has_key(node):
                from_node = node
                graph[from_node] = []
                can_delete = True
                break

    #print "MLEN", len (graph)
    mlen = len (graph)
    print mlen, from_node
    matrix = np.array ([[0]*mlen]*mlen, np.uint32)


    for key in graph.keys():
        nodes = graph[key]
        key_idx = graph.keys().index(key)
        for node in nodes:
            if graph.has_key(node):
                node_idx = graph.keys().index(node)
                #XXX print key_idx, node_idx, node
                matrix[key_idx][node_idx] += 1

    #print graph.keys ()
    #print matrix
    #print "MLEN", len (graph)
    for i in range (0, mlen):
        matrix[i,i] = 0
        col_count = np.sum(matrix[:,i])
        row_count = np.sum(matrix[i,:])
        k = graph.keys()[i]
        #print "GRAPH", i, k
        if row_count > col_count:
            #print "HEUREKA"
            to_node = k
        elif row_count < col_count:
            from_node = k

    graph[from_node].append (to_node)


    #print graph


    print "FROM TO:", from_node, to_node


    return (from_node, to_node, graph)


def deletePath (graph, from_node, to_node=-1):
    lock = threading.Lock()

    lock.acquire() # will block if lock is already held
    if from_node in graph:
        nodes = graph[from_node]
        if len (nodes) > 1:
            if to_node == -1:
                idx = 0
            else:
                #print "NOEDS", nodes
                idx = nodes.index (to_node)
            del graph[from_node][idx]
        else:
            del graph[from_node]
    lock.release()
    return graph

def deletePath2 (graph, from_node, to_node=-1):
    #print "GRAPH",  graph
    lock = threading.Lock()

    lock.acquire() # will block if lock is already held
    if from_node in graph:
        nodes = graph[from_node]
        #print "NODE", nodes
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


def eucladianPath (graph, start_node=-1, end_node=-1):

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
        #XXXprint "@##################"
        #print "NEW", printGraph(unused_edges)
        temp_cycle = walkPath(unused_edges, start, unused_nodes)

        if unused_nodes:
            start = unused_nodes[0]
        #XXXprint "$$$$$$$$$$$$$$$$$$$$$"
        #XXXprint "NEXT", printGraph(unused_edges)
        #XXXprint "UNISED", unused_nodes, start, unused_edges, "CYCLE", cycle
        temp_cycle += cycle
        idx = temp_cycle.index(start)
        temp_cycle = temp_cycle[idx:]+temp_cycle[0:idx]

        cycle = temp_cycle
    #cycle.append(cycle[0])
    #print "CY", cycle

    idx = cycle.index(end_node)
    cycle = cycle[idx:len(cycle)] + cycle[0:idx]

    return cycle

def prefixReadPair (read_pair):

    klen = (len (read_pair)-1)/2

    result = read_pair[0:klen-1] +"|"+ read_pair[klen+1:len(read_pair)-1]
    return result

def suffixReadPair (read_pair):

    klen = (len (read_pair)-1)/2

    result = read_pair[1:klen] +"|"+ read_pair[klen+2:len(read_pair)]
    return result

def compositionReadPairsGraph (nodes, matrix):
    #XXXprint matrix
    graph = {}

    nlen = len (nodes)

    for i in range(0, nlen):
        values = matrix[i,:]
        arr = values.tolist()
        if 1L in arr:
            graph[nodes[i]] = []
        for j in range (0, nlen):
            if arr[j] > 0:
                #XXXprint "JJJJJJJJJ", j, arr[j]
                vnodes = []
                for k in range (0,arr[j]):
                    graph[nodes[i]].append (nodes[j])

    #XXXprint graph
    return graph


def pathGraphReadPairs (read_pairs, k, d):

    mlen = len (read_pairs)+1
    matrix = np.array ([[0]*mlen]*mlen, np.uint32)
    nodes = []

    #for i in range (0, mlen - 1):
    for node in read_pairs:
        prefix = prefixReadPair(node)
        suffix = suffixReadPair(node)

        if prefix not in nodes:
            nodes.append (prefix)

        if suffix not in nodes:
            nodes.append (suffix)

        pidx = nodes.index(prefix)
        sidx = nodes.index(suffix)

        #print prefix, suffix, pidx, sidx, nodes
        matrix[pidx][sidx] += 1

    #print "AAAA", len (nodes), nodes
    return (nodes, matrix)


def deBrujinReadPairs (read_pairs, k, d):

    (nodes, matrix) = pathGraphReadPairs(read_pairs, k, d)
    nlen = len (nodes)

    return compositionReadPairsGraph (nodes, matrix)

def printGraph (graph):
    #XXXprint "KEYS", graph.keys()
    print graph
    for key in graph:
        to = graph[key]
        #print "TO", to
        if len (to) != 0:
            print "%s -> %s" % (key, ','.join(to))

'''
String Reconstruction from Read-Pairs Problem: Reconstruct a string from its paired composition.
     Input: A collection of paired k-mers PairedReads and an integer d.
     Output: A string Text with (k,d)-mer composition equal to PairedReads (if such a
     string exists).
'''
def reconstructFromReadPairs (paired_kmers, k, d):

    plen = len (paired_kmers)-1
    res_len = plen+2*k+d-1

    #XXXprint plen, res_len

    klen = (len(paired_kmers[0])-1)/2

    result = ['']*res_len
    for i in range (0, res_len):
        if i < plen:
            #result[i] = paired_kmers[i][0:klen]
            val = paired_kmers[i][0]
            #XXXprint "VAL", val[0]
            result[i] =val[0]
        else:
            #XXXprint "VALT", 2*klen+d
            val = paired_kmers[i-(2*k+d-2)][2*klen]
            #XXXprint "IIII", val, i, i-(2*k+d-2)
            #print paired_kmers[i]
            result[i] = val
    #XXXprint "RRREEES", result

    return ''.join (result)

def stringSpelledByGappedPatterns(read_pairs, d):

    #XXXprint read_pairs
    klen = (len(read_pairs[0])-1)/2
    result = ""

    prefix_string = ""
    suffix_string = ""
    for pat in read_pairs:
        #XXXprint "PAT", pat
        prefix_string  += pat[0]
        suffix_string += pat[klen+1]

    prefix_string += read_pairs[len(read_pairs)-1][1:klen]
    suffix_string += read_pairs[len(read_pairs)-1][klen+2:]

    #print "Preffix and suffix", prefix_string, suffix_string

    plen = len (prefix_string)
    for i in range(klen+d, plen):
        #XXXprint "III", i, prefix_string[i], suffix_string[i-klen-d]
        if prefix_string[i] != suffix_string[i-klen-d]:
            #print "ERROR", i, read_pairs[i]
            break

    #return
    result = prefix_string[0:klen+d+1]+ suffix_string
    return result

    #return prefix_string + suffix_string[plen-klen-d:]

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

