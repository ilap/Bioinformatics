__author__ = 'ilap'
import sys
sys.setrecursionlimit(1000000000)
'''
Ch 5.2
Longest Common Subsequence Problem: Find a longest common subsequence of two strings.
     Input: Two strings.
     Output: A longest common subsequence of these strings.
     Method one:
     Any time we remove a symbol on the left, we add it to a growing alignment of
     ATGTTATA and ATCGTCC on the right.
     When only one symbol is removed in a turn, we align it with a space symbol.
'''
def longestCommonSubsequene (str1, str2):

    reps = len (str2)
    tstr1 = str2
    tstr2 = str1

    if len (str1) >= len (str2):
        reps = len (str1)
        tstr1 = str2
        tstr2 = str1
    gstr1 = ""
    gstr2 = ""
    score = 0

    for i in range (0, reps):
        a = i

'''
Manhattan Tourist Problem: Find a longest path in a rectangular city.
     Input: A weighted n x m rectangular grid with n + 1 rows and m + 1 columns.
     Output: A longest path from source (0,0) to sink (n, m) in the grid.
'''

def manhattamPath ():
    result = ""
    return result

'''
Longest Path in a Directed Graph Problem: Find a longest path between two nodes in an edge-weighted directed graph.
     Input: An edge-weighted directed graph with source and sink nodes.
     Output: A longest path from source to sink in the directed graph.
'''
def maximumPathInDAG (weighted_dag):
    result = ""

    return result

'''
Greedy change:
It does not always return with the minumum number of coins e.g. 48
it returns 40, 5, 1, 1, 1 instead of 24, 24
'''

import bisect
def greedyChange (money):

    denominations = [120, 40, 30, 24, 20, 10, 5, 4, 1]
    denominations.sort()

    change = []

    while money:

        index = bisect.bisect(denominations, money)
        coin = denominations[index-1]
        change.append(coin)
        money -= coin

    return change

'''
Change Problem: Find the minimum number of coins needed to make change.
     Input: An integer money and an array Coins of d positive integers.
     Output: The minimum number of coins with denominations Coins that changes money.
'''
import sys
def recursiveChange (money, coins):

    if money == 0:
        return 0

    min_num_coins = sys.maxint
    for i in range (0, len (coins)):
        coin = coins[i]
        if money >= coin:
            num_coins = recursiveChange(money - coin, coins)
            if num_coins + 1  < min_num_coins:
                min_num_coins = num_coins + 1
                return min_num_coins

'''
The following dynamic programming algorithm calculates MinNumCoins(money) with
runtime O(money * |Coins|).
'''
def dynamicChange (money, coins):

    min_num_coins = [0]
    clen = len (coins)
    max_coin = coins[0]

    for m in range (1, money+1):
        min_num_coins.append(sys.maxint)

        midx = m
        if m > max_coin:
            min_num_coins = min_num_coins[1:]
            midx = max_coin

        for i in range (0, clen):
            coin = coins[i]
            if m >= coin:
                if min_num_coins [midx-coin] + 1 < min_num_coins[midx]:
                    min_num_coins[midx] = min_num_coins[midx-coin] + 1

    return min_num_coins[midx]

'''
Manhattam revisited
Dynamic programming...
    n = column
    m = row

    downs is n*(m+1)
    right is (n+1)*m
'''
import numpy as np
def manhattamTouristDP (n, m, down, right):

    n += 1
    m += 1
    profile = [0]*(n*m)
    s = np.array (profile)
    s = s.reshape(n,m)

    print s

    for i in range (1, n):
        s[i,0] = s[i-1,0] + down [i-1,0]

    for j in range (1, m):
        s[0,j] = s[0,j-1]+ right[0,j-1]


    for i in range (1, n):
        for j in range (1, m):
            v1 = s[i-1,j] + down[i-1, j]

            v2 = s[i,j-1] + right[i, j-1]

            print i, j, v1, v2
            s[i,j] = max (v1, v2)

    print s


    return s[n-1,m-1]

def generateMatrix (graph):
    result = []

    glen = len (graph)
    matrix = np.zeros (shape=(glen,glen))

    rows = graph.keys ()
    cols = rows[:]
    for i in range (0, glen):
        node = cols[i]
        out_nodes = graph[node]

        olen = len (out_nodes)
        for j in range (0, olen):
        #for onode in out_nodes:
            onode = out_nodes[j]
            if type (onode) is tuple:
                (onode, dummy) = onode

            if not onode in cols:
                cols.append (onode)
                col = np.zeros ((glen, 1))
                matrix = np.append(matrix, col, axis=1)

            idx = cols.index (onode)
            print i, idx, node, onode
            matrix[i, idx] += 1


    return (rows, cols, matrix)

import copy
def topologicalOrder (graph):
    result = []
    graph = copy.deepcopy(graph)

    (rows, cols, matrix) = generateMatrix(graph)

    print matrix
    candidates = []
    print"KEY", graph.keys ()
    for node in graph.keys ():
        idx = cols.index (node)
        if not sum(matrix[:,idx]):
            candidates.append(node)

    #print "Candidates", candidates
    while candidates:
        #print "LEN CAND", len (candidates)
        node = candidates[0] #len(candidates)-1]
        candidates.remove(node)
        result.append(node)

        if graph.has_key (node):
            out_nodes = graph[node]
            del (graph[node])
        else:
            out_nodes = []

        for onode in out_nodes:
            #print "COLS", onode, cols
            if type (onode) is tuple:
                (onode, dummy) = onode
            idx_col = cols.index (onode)
            idx_row = cols.index (node)
            if matrix [idx_row, idx_col] >= 1:
                matrix [idx_row, idx_col] -= 1
            if not sum(matrix[:,idx_col]):
                candidates.append(onode)

    if graph:
        print "Thi input graph is not a DAG"
        return []

    return result

import re
def convertToTupleDAG (str):
    result = {}
    for line in str.splitlines ():
        if not line:
            continue

        key = int (re.findall(r'(\d+).*->', line)[0])
        #NOT TUPLES nodes = re.findall(r'(\d+):\d+', line)
        nodes = [(int (x),int(y)) for (x,y) in re.findall(r'(\d+):(\d+)', line)]
        if not nodes:

            nodes = [int (x) for x in re.findall(r'\d+', line)]
            nodes = nodes[1:]
        print "NODES", nodes
        if result.has_key (key):
            result[key] += nodes
        else:
            result[key] = nodes

    return result

def LCSBacktrack (v, w):

    v_len = len (v)+1
    w_len = len (w)+1

    s = np.zeros (shape=(v_len, w_len))


    for i in range (1, v_len):
        for j in range (1, w_len):
            down = s[i-1,j]
            left = s[i, j-1]
            downleft = s[i-1, j-1]

            if v[i-1] == w[j-1]:
                s[i,j] = downleft + 1
            else:
                s[i,j] = max (down, left)


    return s


def outputLCS (backtrack, v, w, i, j):
    if i == 0 or j == 0:
            return ""

    print backtrack[i,j]
    if v[i-1] == w[j-1]:
        return outputLCS(backtrack, v, w, i - 1, j - 1) + v[i-1]
    else:
        if backtrack[i, j-1] > backtrack[i-1, j]:
            return outputLCS(backtrack, v, w, i, j-1)
        else:
            return outputLCS(backtrack, v, w, i - 1, j)


def find_path(graph, start, end, path=[]):
    print "FFP", start, end, path
    path = path + [start]
    if start == end:
        return path
    if not graph.has_key(start):
        return None
    for node in graph[start]:
        print "NODE:", start, " -->", node, path
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath: return newpath
    return None





def backtrackDAG (source, sink, graph, s=-1, e=-1):
    print "GGRAPH", graph

    arr =  topologicalOrder(graph)

    DAG = arr[arr.index(source):arr.index(sink)+1]

    print "DAG", DAG, graph
    maxes =  findLongestPath (DAG, source, sink, graph)

    print "MAXES", maxes
    max_val = -1
    idx = 0
    for max in maxes:
        tval = max[0]
        if max_val < tval:
            max_val = tval
            max_idx = idx

        idx += 1

    if max_val > 0:
        print "MAXES", maxes[max_idx]
        result = maxes[max_idx][1:]
        return '->'.join (str (x) for x in result)
    return "ERROR"



def findLongestPath (arr, start, end, graph, path=[0], weight=0):

        print "PATH", start, weight, path
        path = path + [start]

        path[0] = path[0] + weight



        if start == end:
            print path[0]
            return [path]

        if not graph.has_key(start):
            return []

        paths = []

        for node in graph[start]:
            #print "NODE", node
            if type (node) is tuple:
                (node, dummy) = node



            if node not in path:
                newpaths = findLongestPath(arr, node, end, graph, path, dummy)
                for newpath in newpaths:
                    paths.append(newpath)

            weight = 0

        return paths

import numpy as np
def globalAlignmentGraph (v, w, ro=0, mu=0):

    inc=1
    v_len = len (v)+1
    w_len = len (w)+1



    s = np.zeros (shape=(v_len, w_len))
    for i in range (v_len):
        s[i,0] = i*-ro
    for i in range (w_len):
        s[0,i] = i*-ro


    (row, col, matrix) = generateScoreMatrixFromFile("BLOSUM62.txt")


    for i in range (1, v_len):
        for j in range (1, w_len):
            v_char = v[i-1]
            w_char = w[j-1]
            v_idx = row.index(v_char)
            w_idx = col.index(w_char)

            down = s[i-1,j] - ro
            left = s[i, j-1] - ro
            downleft = s[i-1, j-1] + matrix[v_idx, w_idx]


            s[i,j] = max (down, left, downleft)

            #print i, j, v[i-1], w[j-1], down, left, downleft, matrix[v_idx, w_idx]


    return (s, row, col, matrix)

def generateScoreMatrixFromFile (filename):
    lines = open(filename,'rU').read().split("\n")

    col = ''.join (lines[0].split())
    row = ""
    lines = lines[1:]

    mlen = len (col)

    matrix = np.zeros (shape=(mlen,mlen))

    for i in range (mlen):
        for j in range (mlen+1):
            line = ' '.join(lines[i].split()).split()

            if not j:
                row += line[j]
            else:
                matrix[i,j-1] = line[j]

    return (row, col, matrix)
####
global __v_alignment
__v_alignment = str ()

def backTrackScoreMatrix (backtrack, v, w, i, j, ro, row, col, matrix):

    global __v_alignment
    v_len = len (v)
    w_len = len (w)

    print "IJ", i, j
    if v_len >= w_len:
        if i == 0:
            print __v_alignment
            return ""
    else:
        if j == 0:
            print __v_alignment
            return ""

    v_idx= row.index(v[i-1])
    w_idx= col.index(w[j-1])

    diag = backtrack[i-1,j-1]+ matrix[v_idx, w_idx]
    down = backtrack[i-1,j] - ro
    right = backtrack[i,j-1] - ro

    max_v = max (down, right, diag)

    #print "(", i,",",j,") ",max_v, " DRG", down, right, diag,

    if max_v == down:
        #print "down"
        __v_alignment = v[i-1] + __v_alignment
        return backTrackScoreMatrix(backtrack, v, w, i-1, j, ro, row, col, matrix) + "-"
    elif max_v == right:
        #print "right"
        __v_alignment = "-" + __v_alignment
        return backTrackScoreMatrix(backtrack, v, w, i, j-1, ro, row, col, matrix) + w[j-1]
    else:
        #print "diag"
        __v_alignment = v[i-1] + __v_alignment
        return backTrackScoreMatrix(backtrack, v, w, i-1, j-1, ro, row, col, matrix) + w[j-1]

import numpy as np
def levenshteinDistance(v, w):

    v_len = len (v)+1
    w_len = len (w)+1

    s = np.zeros (shape=(v_len, w_len))
    for i in range (v_len):
        s[i,0] = i
    for i in range (w_len):
        s[0,i] = i


    for i in range (1, v_len):
        for j in range (1, w_len):
            v_char = v[i-1]
            w_char = w[j-1]

            down = s[i-1,j]+1
            left = s[i, j-1]+1
            downleft = s[i-1, j-1]

            if v_char == w_char:
               s[i,j] = downleft
            else:
                s[i,j] = min (down, left, downleft+1)
    #print s
    return int (s[v_len-1, w_len-1])
