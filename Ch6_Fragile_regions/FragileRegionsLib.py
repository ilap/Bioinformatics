__author__ = 'ilap'

from math import *
import numpy as np

def getArrayFromTuppleString (string):

    return [int (i) for i in string.split ()]

def kSortinReversal (P, start, end):

    A=P[start:end+1]
    #print "A", A

    A=np.array(A)
    A=A*-1

    A = A[::-1]
    A = A.tolist ()
    #print "AA", P[0:start], A, P[end+1:]

    P = P[0:start] + A + P[end+1:]
    return P

global distance
distance = 0
def greedySort (P):

    global distance
    print "P", P
    plen = len (P)
    result = ""
    for k in range(plen):
        if k+1 not in P:
            idx = P.index(-(k+1))
        else:
            idx = P.index(k+1)
        #print "IDX", k, idx
        if (k+1) != P[k]:
            P = kSortinReversal(P,k, idx)
            distance += 1
            result = result +  "("+' '.join ('{0:+}'.format(x) for x in P)+")\n"
        if -(k+1) == P[k]:
            P[k] = -P[k]
            distance += 1
            result = result +  "("+' '.join ('{0:+}'.format(x) for x in P)+")\n"
    print "Distance", distance
    return result


def adjency (Pi, Pi1):
    return (Pi1 - Pi) == 1

global breakpoints
breakpoints = 0
def greedySortWithBreakPoint (P):

    global breakpoints
    plen = len (P)
    P = [0] + P + [plen+1]
    print "P", P

    result = ""
    for k in range(plen+1):

        if not adjency(P[k], P[k+1]):
            breakpoints += 1
            #print "AAA", k, P[k], P[k+1], breakpoints


    return breakpoints

def chromosomeToCycle (chrom):
    nodes = []
    #print chrom
    for j in range (0, len (chrom)):
        i = chrom [j]

        #print "J", j,i
        if i > 0:
            nodes.append(2*i-1)
            nodes.append(2*i)
        else:
            nodes.append(-2*i)
            nodes.append(-2*i-1)
    return nodes

def cycleToChromosome (nodes):
    #print "NODES", nodes
    chrom = []
    for j in range (len (nodes)/2):
        #print "N", nodes[2*j], nodes[2*j+1]
        if nodes[2*j] < nodes[2*j+1]:
            chrom.append ('{0:+}'.format(nodes[2*j+1]/2))
        else:
            chrom.append ('{0:+}'.format(-(nodes[2*j]/2)))
    return "("+' '.join (chrom)+")"

def coloredEdges (perms):
    cycles = []

    for chrom in perms:
        nodes = chromosomeToCycle(chrom)
        nodes.append (nodes[0])
        edges = []
        print "xNODES", chrom, nodes
        for j in range (1, len (chrom)+1):
            edges.append ((nodes[2*j-1],nodes[2*j]))
        cycles.append(edges)

    return cycles

def getPermsFromTuppleStrings (string):
    result = []
    arrays=string.split ("(")
    arrays=arrays[1:]

    #print "ARR", arrays
    for arr in arrays:
    # return [int (i) for i in string.split ()]
        arr=arr.split(")")
        arr=arr[0:-1]
        arr = arr[0].split()
        arr = [int (x) for x in arr]
        result.append (arr)
    return result

def graphToGenome (graph):
    P = {}
    for nodes in graph:
        chrom = cycleToChromosome(nodes)
        print chrom