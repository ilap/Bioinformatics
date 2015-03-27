__author__ = 'ilap'

import math
import time

def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0)
        return ret
    return wrap

global MASS_TABLE
global PEPTID
global MASS

PEPTID = 0
MASS = 1
MASS_TABLE = \
[["G", "A", "S", "P", "V", "T", "C", "I", "L", "N", "D", "K", "Q", "E", "M", "H", "F", "R", "Y", "W"],
 [57,   71,  87,  97,  99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]]

SMASS_TABLE = \
[["G", "A", "S", "P", "V", "T", "C", "I", "N", "D", "K", "E", "M", "H", "F", "R", "Y", "W"],
 [57,   71,  87,  97,  99, 101, 103, 113, 114, 115,  128, 129, 131, 137, 147, 156, 163, 186]]

'''
Input: Length of peptide e.g. NQEL => 4
Output: Number of linear subpeptides e.g.
N Q E L NQ QE EL LN NQE QEL ELN LNQ NQEL => 13
'''
def getNumOfCyclicSubPeptides (len):
    return len * (len - 1) + 1

'''
Input: Animo acid e.g. "G"
Output: mass of animo acid e.g 57
'''
def getPeptidMassByName (amino_acid):
    return MASS_TABLE[1][MASS_TABLE[0].index(amino_acid)]


'''
Input: spectrum of an peptid e.g.
Output:
'''
def getCyclicalPeptideLengthFromSpectrum (spectrum):

    result = 0
    idx = 1
    while idx**2 - idx < spectrum:
        idx += 1


    print "IDX:" , idx
    if idx**2 - idx == spectrum:
        result = idx

    return result

'''
Input: Length of peptide e.g. NQEL => 4
Output: Number of linear subpeptides e.g.
N Q E L NQ QE EL LN NQE QEL ELN NQEL => 10
'''
def getNumOfLinearPeptides (len):
    return len * (len + 1 ) / 2

def getLinearPeptideSizeFromSpectrum (spectrum):

    result = 0
    idx = 0
    tval = 0
    while result <= spectrum - idx:

        idx += 1
        result += idx

    return idx

def peptidAsMassFast (peptide):
    #print "PEPTIDE", peptide
    mass = 0
    plen = len (peptide)

    for i in range (0,plen):
        mass += getPeptidMassByName(peptide[i])

    return mass

def peptidAsMass (peptide):
    #print "PEPTIDE", peptide
    mass = 0
    plen = len (peptide)
    #print "MASS0", plen, peptide
    if plen == 1:
        mass += getPeptidMassByName(peptide)
        #print "MASS1", mass, peptide
        return mass
    elif plen == 0:
        return mass
    mass = getPeptidMassByName(peptide[plen-1]) + peptidAsMass (peptide[0:plen-1])
    #print "MASS2", mass, peptide
    return mass

def peptidAsMasses (peptide, masses = None):

    plen = len (peptide)
    if masses == None: masses = []

    if plen == 1:
        mass = getPeptidMassByName(peptide)
        return masses.append(mass)

    mass = getPeptidMassByName(peptide[0])
    masses.append(mass)
    peptidAsMasses (peptide[1:plen], masses)
    return masses

def getPeptidFrommasses (masses):
    result = ""
    for i in masses:
        result += MASS_TABLE[0][MASS_TABLE[1].index(i)]
    return result
'''
Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.
     Input: An amino acid string Peptide.
     Output: cycloSpectrum(Peptide).
'''
def cycloSpectrum (peptid):

    plen = len (peptid)

    #if plen == 1:
    #    return getPeptidMassByName(peptid)

    new_peptid = peptid + peptid
    arr_len = getNumOfCyclicSubPeptides(plen) + 1 # add the extra ("", 0) as first element in arr

    #print arr_len

    cyclo_spec = []

    arr = [0] * arr_len
    cyclo_spec.append (arr)

    arr = [""] * arr_len
    cyclo_spec.append (arr)

    idx = 0
    for i in range (0, plen - 1):
        for j in range (1, plen + 1):
            idx = i * plen + j
            chr = new_peptid[(j-1):(j+i)]
            mass = peptidAsMass(chr)
            cyclo_spec[0][idx] = mass
            cyclo_spec[1][idx] = chr

    cyclo_spec [0][arr_len - 1] = peptidAsMass(peptid)
    cyclo_spec [1][arr_len - 1] = peptid
            #print "IDX:", idx-1, "CHR:", arr_len, chr
    #print cyclo_spec
    return sorted (cyclo_spec[0])

'''
Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.
     Input: An amino acid string Peptide.
     Output: cycloSpectrum(Peptide).
'''
def linearSpectrum (peptid):

    plen = len (peptid)

    #print "PEPTID:", peptid
    #if plen == 1:
    #    return getPeptidMassByName(peptid)

    new_peptid = peptid + peptid
    arr_len = getNumOfLinearPeptides(plen) + 1

    #print arr_len

    linear_spec = []

    arr = [0] * arr_len
    linear_spec.append (arr)

    arr = [""] * arr_len
    linear_spec.append (arr)

    idx = 1
    for i in range (0, plen - 1):
        for j in range (1, plen + 1 - i):
            chr = new_peptid[(j-1):(j+i)]
            mass = peptidAsMass(chr)
            linear_spec[0][idx] = mass
            linear_spec[1][idx] = chr
            #print chr, mass


            idx += 1


    linear_spec [0][arr_len - 1] = peptidAsMass(peptid)
    linear_spec [1][arr_len - 1] = peptid

    #print linear_spec

    return sorted (linear_spec[0])

def linearSpectrumFast (peptid):

    plen = len (peptid)

    #print "PEPTID:", peptid
    #if plen == 1:
    #    return getPeptidMassByName(peptid)

    new_peptid = peptid + peptid
    arr_len = getNumOfLinearPeptides(plen) + 1

    #print arr_len

    linear_spec = []

    arr = [0] * arr_len
    linear_spec.append (arr)

    arr = [""] * arr_len
    linear_spec.append (arr)
    #print linear_spec


    idx = 1
    for i in range (0, plen - 1):
        for j in range (1, plen + 1 - i):
            chr = new_peptid[(j-1):(j+i)]
            mass = peptidAsMass(chr)
            linear_spec[0][idx] = mass
            linear_spec[1][idx] = chr
            #print chr, mass


            idx += 1


    linear_spec [0][arr_len - 1] = peptidAsMass(peptid)
    linear_spec [1][arr_len - 1] = peptid

    #print linear_spec

    return sorted (linear_spec[0])

'''
Parent Mass is the know mass of the unknown peptide in the experimental spectrum
'''
def parentMass (spectrum):
    return max (spectrum)

'''
Cyclopeptide Sequencing Problem: Given an ideal experimental spectrum, find a cyclic peptide whose theoretical spectrum matches the experimental spectrum.
     Input: A collection of (possibly repeated) integers Spectrum corresponding to an ideal
     experimental spectrum.
     Output: An amino acid string Peptide such that Cyclospectrum(Peptide) = Spectrum (if such a
     string exists).
'''
def BFcycloPeptideSequencing (spectrum):
    #TODO: output Peptide

    return ""

'''

'''
def expandPeptides (peptides):
    result = []

    if peptides == []:
        for mass in SMASS_TABLE[MASS]:
                a = [mass]
                result.append (a)
    else:
        for i in peptides:
            summa = sum (i)

            #print"SUM:", summa
            for mass in SMASS_TABLE[MASS]:
                #print "LENI", i[len (i)-1], mass
                tmp = i[:] # or list (i)
                #if mass not in i and mass + summa in spectrum:
                    #print "MASS+SUM", mass + summa

                tmp.append (mass)
                    #tmp.append (mass)
                    #print "TMP:", tmp, i
                result.append (tmp)
    print "expandPeptides", result
    return result

def expandStringPeptides (peptides):
    result = []

    if peptides == []:
        for peptide in SMASS_TABLE[PEPTID]:
                result.append (peptide)
    else:
        for i in peptides:
            #summa = sum (i)

            #print"SUM:", summa

            for peptide in SMASS_TABLE[PEPTID]:
                #print "LENI", i[len (i)-1], mass
                 # or list (i)
                #if mass not in i and mass + summa in spectrum:
                    #print "MASS+SUM", mass + summa
                tmp = i + peptide
                #tmp.append (peptide)
                    #tmp.append (mass)
                    #print "TMP:", tmp, i
                result.append (tmp)
    #print "expandPeptides", result
    return result

def convertMassesToPeptid (masses):
    result = ""
    #print "MASSES:", masses
    for mass in masses:
        if mass in MASS_TABLE[1]:
            idx = MASS_TABLE[1].index (mass)
            result += MASS_TABLE[0][idx]
    return result

def isPeptideConsistent (peptide, spectrum):
    #print "PEPTIDE", peptide
    if sum (peptide) not in spectrum:
        #print "FFALSE". sum (peptide), peptide
        return False
    tmp = []
    for amino_acid in peptide:
        tmp.append(amino_acid)
        summa = sum (tmp)

     #   print "ISCON", amino_acid
        if amino_acid in spectrum:
            pcount = peptide.count (amino_acid)
            scount = spectrum.count (amino_acid)
            #print "PS:", pcount, scount
            if pcount > scount:
                #print "FALSE"
                return False
        else:
            #print "FALSE"
            return False
    #print "TRUE"
    return True

def scorePeptid (peptid, experimental):
    score = 0

    #print "scorePeptid", peptid
    theoretical = cycloSpectrum(peptid)
    unique = set (theoretical)

    for mass in unique:

        if mass in experimental:
            tcount = theoretical.count (mass)
            multiplicity = experimental.count (mass)
            #print "PS:", tcount, multiplicity
            if tcount > multiplicity:
                score += multiplicity
            else:
                score += tcount
    #print "Score", len (experimental), score
    return score

def linearPeptidScorePeptid (peptid, experimental):
    score = 0
    theoretical = linearSpectrum(peptid)
    unique = set (theoretical)

    for mass in unique:

        if mass in experimental:
            tcount = theoretical.count (mass)
            multiplicity = experimental.count (mass)
            #print "PS:", tcount, multiplicity
            if tcount > multiplicity:
                score += multiplicity
            else:
                score += tcount
    #print "TLEN", len (experimental)
    return score

def linearPeptidScorePeptidFast2 (peptid, experimental):
    score = 0
    theoretical = linearSpectrumFast(peptid)
    unique = set (theoretical)

    for mass in unique:

        if mass in experimental:
            tcount = theoretical.count (mass)

            if tcount > 1:
                multiplicity = experimental.count (mass)
                #print "PS:", tcount, multiplicity
                if tcount > multiplicity:
                    score += multiplicity
                else:
                    score += tcount
            score += tcount
    #print "TLEN", len (experimental)
    return score
'''
Banch-and-bound algorythm for peptide sequencing.
'''
def cycloPeptideSequencing (spectrum):
    #TODO: cycloPeptidSequencing
    print "SPectrum:", spectrum

    peptides = expandPeptides ([])
    #print "peptides:", peptides
    idx=0
    result = []
    while len (peptides) != 0:
        idx += 1
        #print "PEPTIDES", peptides
        newpeptids = peptides [:]
        for peptide in peptides:
            if  sum (peptide) == parentMass (spectrum):
                peptide_str = convertMassesToPeptid (peptide)
                #print "STR:", peptide_str, peptide
                #print ' '.join (str (x) for x in  peptide)
                #print '-'.join (map (str,peptide))

                if cycloSpectrum (peptide_str) == spectrum:
                    print '-'.join (map (str,peptide))
                    result.append (peptide_str)
                newpeptids.remove(peptide)
            elif not isPeptideConsistent (peptide, spectrum):
                newpeptids.remove(peptide)

        #print "peptides:", len (newpeptids), newpeptids
        if len (newpeptids) == 0:
            break
        peptides = expandPeptides(newpeptids)
        #if idx == 1: break
    return result

@timing
def trimLeaderBoard (leader_board, spectrum, N, amino_acid = 0, amino_acid_mass = 0):

    plen = len (leader_board)
    linear_scores = []

    idx = 0
    print "##### TRIM for begin"
    for peptide in leader_board:

        score = linearPeptidScorePeptid(peptide, spectrum)
        linear_scores.append([score])
        linear_scores[idx].append (peptide)
        idx += 1
    print "##### TRIM for end"
    #print "TRIM sort"
    linear_scores = sorted (linear_scores,None,None,True)

    #print "TRIM for 2"
    for j in range (N, len (linear_scores)):
        if linear_scores[j][0] < linear_scores[N-1][0]:
            leader_board = [lb[1] for lb in linear_scores[:j]]
            #leader_board[j:] = []
            #print len (new), new
            #print len (leader_board), leader_board
            return leader_board

    #"TRRIM create array"
    leader_board = [lb[1] for lb in linear_scores]
    return leader_board

@timing
def trimLeaderBoardFast (leader_board, spectrum, N, amino_acid = 0, amino_acid_mass = 0):

    plen = len (leader_board)
    linear_scores = [[0, '']]*plen

    print "######### FOR"
    for idx in range (0, plen):
        peptide = leader_board[idx]
        score = linearPeptidScorePeptid(peptide, spectrum)
        linear_scores[idx] = [score, peptide]
    print "######### SORT"
    linear_scores = sorted (linear_scores,None,None,True)

    print "######### FOR 2"
    for j in range (N, len (linear_scores)):
        if linear_scores[j][0] < linear_scores[N-1][0]:
            leader_board = [lb[1] for lb in linear_scores[:j]]
            #leader_board[j:] = []
            #print len (new), new
            #print len (leader_board), leader_board
            return leader_board

    leader_board = [lb[1] for lb in linear_scores]
    return leader_board

'''
Banch-and-bound algorythm for peptide sequencing.
'''
def leaderBoardCycloPeptideSequencing3 (spectrum, N):
    #TODO: cycloPeptidSequencing
    print "SPectrum:", spectrum

    leader_board = expandPeptides ([])
    leader_peptide = [""]
    idx=0
    result = []

    while len (leader_board) != 0:
        idx += 1
        newpeptids = leader_board[:]

        for peptide in newpeptids:
            if  sum (peptide) == parentMass (spectrum):
                p =  convertMassesToPeptid (peptide)
                lp = convertMassesToPeptid (leader_peptide)
                ps = linearPeptidScorePeptid(p, spectrum)
                ls = linearPeptidScorePeptid(lp, spectrum)
                print "SCORE:", p, lp, ps, ls
                if  ps > ls :
                    print '-'.join (map (str,peptide))
                    result.append (p)
                    leader_peptide.append (p)
                newpeptids.remove(peptide)
            elif sum (peptide) > parentMass (spectrum):
                newpeptids.remove(peptide)

        leader_board = trimLeaderBoard (leader_peptide, spectrum, N)
        print "leader_board:", len (leader_board), leader_board
        if len (newpeptids) == 0:
            break
        leader_board = expandPeptides(leader_board)
        if idx == 2: break
    return result


@timing
def leaderBoardCycloPeptideSequencing (spectrum, N):

    leader_board = ['']
    leader_peptide = ""
    result = []

    idx = 0
    while len (leader_board) != 0:
        print "@@@ expand"
        leader_board = expandStringPeptides (leader_board)
        #print leader_board
        #newpeptids = leader_board[:]
        scores = []
        print "### nlb"
        nlb = leader_board[:]
        print "@@@ FOR"
        for peptide in leader_board:

            #print "P", peptide[0], leader_board
            peptide_mass = peptidAsMass(peptide)
            parent_mass = parentMass(spectrum)
            #print "Mass", peptide_mass, parent_mass, len (leader_board)
            if  peptide_mass == parent_mass:
                ps = scorePeptid(peptide, spectrum)
                ls = scorePeptid(leader_peptide, spectrum)

                #print "SCORE:", ps, ls
                if  ps > ls :
                    #print "SCOREPEP:", ps, "SCORE2", ls

                    #print '-'.join (map (str,peptide)), ps, ls
                    print '-'.join ( str (x) for x in peptidAsMasses(peptide))
                    #result.append (peptide)
                    leader_peptide = peptide
                    #return result
            elif peptide_mass > parent_mass:
                #print "NLP, remove"
                nlb.remove(peptide)
                #print "remove:, size", peptide, len (leader_board)
        #print "SCORES:", scores
        print "LEN1:", len (leader_board)
        leader_board = trimLeaderBoardFast (nlb, spectrum, N)
        print "LEN2:", len (leader_board)
        #print "LEN3:", len (leader_board)
        #idx += 1
        #if idx == 23: break

    return result

'''
Spectral Convolution Problem: Compute the convolution of a spectrum.
     Input: A collection of integers Spectrum.
     Output: The list of elements in the convolution of Spectrum. If an element has
     multiplicity k, it should appear exactly k times; you may return the elements in any order.
'''
def convolution (spectrum):
    slen = len (spectrum)
    result = []

    sorted_spectrum = sorted (spectrum)
    print "SORTED", sorted_spectrum[1]
    for idx in range (0, slen - 1):
        tmp = []
        for j in range (slen - 1, idx, -1):
            #print idx, j
            val = sorted_spectrum[j] - sorted_spectrum [idx]
            if val == 0: continue
            tmp.append (val)
            result.append(val)
        print tmp
    print [str (x) for x in result]
    return  result
### MAIN
T = "441 99 597 409 684 131 0 230 787 156 625 113 312 540 128 574 213 827 313 310 471 475 200 415 785 571 771 216 526 469 287 753 728 358 671 57 572 156 259 253 413 528 728 97 643 512 654 103 756 668 344 631 884 356 241 443 781 372"
#T = "0 99 113 114 128 227 257 299 355 356 370 371 484"

#T = "0 113 114 128 129 227 242 242 257 355 356 370 371 484"
#T="0 136 137 323"
#T =""
T = "0 57 118 179 236 240 301"
SPECTRUM =  [int (x) for x in T.split (' ')]
#print "LENS", len (T), len (SPECTRUM)
c = sorted (convolution(SPECTRUM))
print ' '.join (str (x) for x in c)
from itertools import groupby as g
print "MOST", max(g(sorted(c)), key=lambda(x, v):(len(list(v)),-c.index(x)))[0]
print "COUNT", c.count (99)

import collections
print "MOST COMMONS", collections.Counter(c).most_common()
#print len (c), len (SPECTRUM), 74*73/2
## Solve cyclicSpectrum
#print getNumOfLinearPeptides (4)
#print getNumOfCyclicSubPeptides (4)
SPECTRUM = [0, 97, 97, 99, 101, 103, 196, 198, 198, 200, 202, 295, 297, 299, 299, 301, 394, 396, 398, 400, 400, 497]
#SPECTRUM = [ 0, 113, 128, 186, 241, 299, 314, 427]
#T = "0 87 99 99 113 114 128 147 186 186 186 227 227 242 246 285 299 314 333 333 341 355 413 413 428 432 432 446 454 519 527 541 541 545 560 560 618 632 640 640 659 674 688 727 731 746 746 787 787 787 826 845 859 860 874 874 886 973"
#SPECTRUM =  [int (x) for x in T.split (' ')]

#p = expandPeptides([])
#print p

#p = expandPeptides(p)
#print p

#AA = cycloPeptideSequencing(SPECTRUM)
#print AA
#
#a = cycloSpectrum("NQEL")
#print a

#T = "0 71 71 87 97 97 97 99 99 99 99 99 101 101 103 103 113 113 113 113 113 113 114 115 128 128 129 129 129 131 131 137 137 147 147 147 156 163 163 163 163 170 170 184 186 198 200 204 204 210 211 212 212 214 214 218 226 226 226 227 227 228 234 236 241 242 242 243 244 246 250 251 258 262 265 275 281 283 289 294 294 301 301 301 310 317 317 317 319 324 325 326 326 327 329 339 340 341 345 348 350 351 354 355 356 357 359 362 363 364 371 374 378 380 388 390 390 395 398 400 404 406 407 416 425 425 430 430 430 442 447 448 450 451 452 456 457 461 461 462 465 469 470 470 472 473 477 477 485 487 489 491 494 501 503 503 503 505 508 518 527 529 532 537 537 543 543 545 553 553 555 559 560 562 565 569 570 581 588 588 590 590 590 590 598 598 598 600 600 604 606 609 611 616 617 618 620 631 631 632 636 640 640 642 651 656 656 656 668 671 674 675 680 681 684 689 697 697 697 702 703 708 712 713 713 717 719 719 730 731 733 735 737 741 744 745 745 746 746 748 751 751 753 754 757 759 767 774 779 780 787 801 802 808 810 811 812 815 825 826 830 831 832 832 834 834 836 837 842 844 845 848 849 850 851 854 854 855 860 860 861 864 871 875 898 900 907 907 907 909 911 914 916 923 924 931 931 932 937 943 945 947 948 948 949 950 952 957 957 959 961 962 963 964 965 967 968 972 978 978 988 994 995 1000 1001 1006 1020 1031 1034 1035 1037 1038 1040 1042 1046 1047 1056 1059 1060 1060 1060 1060 1061 1061 1063 1063 1065 1065 1068 1070 1075 1077 1078 1081 1087 1091 1096 1096 1102 1104 1106 1109 1113 1119 1134 1136 1141 1143 1149 1151 1157 1158 1159 1159 1160 1161 1162 1163 1164 1165 1169 1171 1174 1174 1178 1178 1188 1188 1190 1190 1193 1196 1199 1205 1205 1207 1210 1215 1217 1222 1226 1235 1248 1249 1249 1250 1258 1259 1260 1267 1269 1273 1275 1276 1278 1282 1287 1287 1290 1291 1294 1296 1296 1300 1302 1303 1304 1304 1306 1306 1307 1307 1309 1311 1314 1319 1327 1334 1335 1348 1348 1352 1352 1357 1362 1362 1368 1377 1377 1378 1381 1388 1389 1397 1400 1401 1401 1404 1405 1407 1410 1410 1413 1415 1419 1420 1422 1425 1430 1432 1432 1434 1434 1439 1448 1451 1451 1453 1453 1459 1461 1461 1463 1465 1476 1490 1491 1494 1497 1502 1504 1505 1507 1509 1511 1512 1514 1515 1517 1518 1520 1521 1525 1529 1531 1532 1538 1538 1538 1547 1548 1552 1561 1564 1566 1568 1574 1576 1578 1579 1581 1588 1588 1590 1593 1604 1605 1608 1608 1615 1616 1620 1621 1626 1628 1637 1640 1641 1643 1646 1649 1651 1651 1652 1653 1658 1660 1660 1661 1665 1665 1667 1668 1675 1675 1676 1677 1677 1678 1679 1687 1691 1701 1702 1705 1711 1721 1724 1731 1739 1739 1741 1744 1750 1751 1752 1752 1752 1755 1757 1758 1759 1772 1775 1777 1777 1778 1778 1778 1779 1780 1788 1789 1790 1792 1803 1804 1807 1808 1814 1814 1818 1822 1824 1830 1831 1838 1838 1839 1840 1849 1852 1855 1857 1861 1865 1871 1872 1876 1879 1880 1885 1889 1890 1891 1891 1891 1894 1895 1902 1904 1905 1907 1909 1909 1910 1915 1915 1918 1935 1936 1938 1940 1942 1944 1962 1965 1966 1967 1970 1971 1977 1978 1979 1987 1988 1993 1994 1994 1999 2001 2003 2004 2004 2005 2008 2008 2008 2015 2018 2020 2022 2025 2028 2032 2033 2036 2041 2041 2043 2049 2051 2053 2057 2058 2059 2065 2071 2076 2081 2091 2102 2103 2106 2107 2107 2108 2109 2117 2118 2119 2123 2124 2128 2129 2130 2133 2134 2135 2136 2141 2142 2146 2146 2148 2154 2156 2156 2157 2158 2166 2174 2178 2178 2180 2189 2199 2204 2204 2205 2205 2206 2206 2208 2216 2219 2220 2222 2228 2230 2233 2235 2237 2244 2245 2247 2249 2255 2255 2255 2261 2261 2265 2266 2269 2270 2275 2276 2283 2283 2288 2291 2298 2305 2305 2309 2318 2319 2320 2321 2329 2332 2333 2334 2334 2341 2341 2346 2350 2350 2353 2359 2360 2362 2367 2368 2368 2369 2370 2373 2375 2379 2382 2384 2389 2392 2398 2408 2408 2412 2412 2413 2418 2418 2419 2428 2431 2433 2434 2438 2447 2447 2447 2449 2452 2454 2461 2463 2469 2470 2472 2472 2474 2475 2481 2483 2490 2492 2495 2497 2497 2497 2506 2510 2511 2513 2516 2525 2526 2526 2531 2532 2532 2536 2536 2546 2552 2555 2560 2562 2565 2569 2571 2574 2575 2576 2576 2577 2582 2584 2585 2591 2594 2594 2598 2603 2603 2610 2610 2611 2612 2615 2623 2624 2625 2626 2635 2639 2639 2646 2653 2656 2661 2661 2668 2669 2674 2675 2678 2679 2683 2683 2689 2689 2689 2695 2697 2699 2700 2707 2709 2711 2714 2716 2722 2724 2725 2728 2736 2738 2738 2739 2739 2740 2740 2745 2755 2764 2766 2766 2770 2778 2786 2787 2788 2788 2790 2796 2798 2798 2802 2803 2808 2809 2810 2811 2814 2815 2816 2820 2821 2825 2826 2827 2835 2836 2837 2837 2838 2841 2842 2853 2863 2868 2873 2879 2885 2886 2887 2891 2893 2895 2901 2903 2903 2908 2911 2912 2916 2919 2922 2924 2926 2929 2936 2936 2936 2939 2940 2940 2941 2943 2945 2950 2950 2951 2956 2957 2965 2966 2967 2973 2974 2977 2978 2979 2982 3000 3002 3004 3006 3008 3009 3026 3029 3029 3034 3035 3035 3037 3039 3040 3042 3049 3050 3053 3053 3053 3054 3055 3059 3064 3065 3068 3072 3073 3079 3083 3087 3089 3092 3095 3104 3105 3106 3106 3113 3114 3120 3122 3126 3130 3130 3136 3137 3140 3141 3152 3154 3155 3156 3164 3165 3166 3166 3166 3167 3167 3169 3172 3185 3186 3187 3189 3192 3192 3192 3193 3194 3200 3203 3205 3205 3213 3220 3223 3233 3239 3242 3243 3253 3257 3265 3266 3267 3267 3268 3269 3269 3276 3277 3279 3279 3283 3284 3284 3286 3291 3292 3293 3293 3295 3298 3301 3303 3304 3307 3316 3318 3323 3324 3328 3329 3336 3336 3339 3340 3351 3354 3356 3356 3363 3365 3366 3368 3370 3376 3378 3380 3383 3392 3396 3397 3406 3406 3406 3412 3413 3415 3419 3423 3424 3426 3427 3429 3430 3432 3433 3435 3437 3439 3440 3442 3447 3450 3453 3454 3468 3479 3481 3483 3483 3485 3491 3491 3493 3493 3496 3505 3510 3510 3512 3512 3514 3519 3522 3524 3525 3529 3531 3534 3534 3537 3539 3540 3543 3543 3544 3547 3555 3556 3563 3566 3567 3567 3576 3582 3582 3587 3592 3592 3596 3596 3609 3610 3617 3625 3630 3633 3635 3637 3637 3638 3638 3640 3640 3641 3642 3644 3648 3648 3650 3653 3654 3657 3657 3662 3666 3668 3669 3671 3675 3677 3684 3685 3686 3694 3695 3695 3696 3709 3718 3722 3727 3729 3734 3737 3739 3739 3745 3748 3751 3754 3754 3756 3756 3766 3766 3770 3770 3773 3775 3779 3780 3781 3782 3783 3784 3785 3785 3786 3787 3793 3795 3801 3803 3808 3810 3825 3831 3835 3838 3840 3842 3848 3848 3853 3857 3863 3866 3867 3869 3874 3876 3879 3879 3881 3881 3883 3883 3884 3884 3884 3884 3885 3888 3897 3898 3902 3904 3906 3907 3909 3910 3913 3924 3938 3943 3944 3949 3950 3956 3966 3966 3972 3976 3977 3979 3980 3981 3982 3983 3985 3987 3987 3992 3994 3995 3996 3996 3997 3999 4001 4007 4012 4013 4013 4020 4021 4028 4030 4033 4035 4037 4037 4037 4044 4046 4069 4073 4080 4083 4084 4084 4089 4090 4090 4093 4094 4095 4096 4099 4100 4102 4107 4108 4110 4110 4112 4112 4113 4114 4118 4119 4129 4132 4133 4134 4136 4142 4143 4157 4164 4165 4170 4177 4185 4187 4190 4191 4193 4193 4196 4198 4198 4199 4199 4200 4203 4207 4209 4211 4213 4214 4225 4225 4227 4231 4231 4232 4236 4241 4242 4247 4247 4247 4255 4260 4263 4264 4269 4270 4273 4276 4288 4288 4288 4293 4302 4304 4304 4308 4312 4313 4313 4324 4326 4327 4328 4333 4335 4338 4340 4344 4344 4346 4346 4346 4354 4354 4354 4354 4356 4356 4363 4374 4375 4379 4382 4384 4385 4389 4391 4391 4399 4401 4401 4407 4407 4412 4415 4417 4426 4436 4439 4441 4441 4441 4443 4450 4453 4455 4457 4459 4467 4467 4471 4472 4474 4474 4475 4479 4482 4483 4483 4487 4488 4492 4493 4494 4496 4497 4502 4514 4514 4514 4519 4519 4528 4537 4538 4540 4544 4546 4549 4554 4554 4556 4564 4566 4570 4573 4580 4581 4582 4585 4587 4588 4589 4590 4593 4594 4596 4599 4603 4604 4605 4615 4617 4618 4618 4619 4620 4625 4627 4627 4627 4634 4643 4643 4643 4650 4650 4655 4661 4663 4669 4679 4682 4686 4693 4694 4698 4700 4701 4702 4702 4703 4708 4710 4716 4717 4717 4718 4718 4718 4726 4730 4730 4732 4732 4733 4734 4740 4740 4744 4746 4758 4760 4774 4774 4781 4781 4781 4781 4788 4797 4797 4797 4807 4807 4813 4813 4815 4815 4815 4816 4816 4829 4830 4831 4831 4831 4831 4831 4831 4841 4841 4843 4843 4845 4845 4845 4845 4845 4847 4847 4847 4857 4873 4873 4944"
#SPECTRUM =  [int (x) for x in T.split (' ')]

p =[]
#for i in range (1, 4):
#    p =  expandStringPeptides(p)
#    print p


#T = "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"
#SPECTRUM =  [int (x) for x in T.split (' ')]
#leaderBoardCycloPeptideSequencing(SPECTRUM, 1000)
#print '-'.join (str (x) for x in peptidAsMasses(""))

T ="0 71 87 97 101 103 103 113 113 113 114 115 115 115 128 128 128 128 129 131 163 163 163 174 184 186 186 198 211 215 216 217 228 231 232 241 242 243 243 246 256 273 276 278 287 291 294 299 302 312 313 314 330 345 346 349 349 357 359 369 370 371 371 379 394 401 404 404 409 415 415 427 436 443 458 459 460 462 464 465 473 476 480 484 488 499 507 516 522 530 532 532 533 544 556 564 565 567 574 577 578 578 586 590 591 595 601 612 635 636 644 647 648 651 657 659 660 662 671 678 679 682 691 693 695 699 705 706 708 749 751 754 763 764 764 772 772 772 775 775 775 776 795 802 807 810 810 819 822 834 837 846 852 862 868 877 878 879 879 885 887 889 903 903 923 924 925 927 935 935 938 940 947 949 949 960 965 980 990 992 1005 1006 1008 1018 1027 1038 1040 1042 1048 1048 1050 1052 1053 1054 1062 1063 1063 1066 1066 1077 1077 1093 1121 1121 1133 1141 1143 1151 1155 1155 1164 1166 1166 1167 1168 1176 1179 1181 1181 1181 1190 1191 1192 1211 1234 1236 1238 1240 1248 1252 1256 1269 1270 1277 1279 1283 1294 1294 1294 1295 1295 1296 1305 1339 1339 1342 1344 1349 1350 1353 1354 1366 1367 1367 1376 1384 1384 1398 1407 1408 1411 1411 1420 1423 1423 1442 1454 1457 1457 1463 1463 1467 1467 1470 1480 1481 1481 1512 1512 1513 1523 1526 1526 1530 1530 1536 1536 1539 1551 1570 1570 1573 1582 1582 1585 1586 1595 1609 1609 1617 1626 1626 1627 1639 1640 1643 1644 1649 1651 1654 1654 1688 1697 1698 1698 1699 1699 1699 1710 1714 1716 1723 1724 1737 1741 1745 1753 1755 1757 1759 1782 1801 1802 1803 1812 1812 1812 1814 1817 1825 1826 1827 1827 1829 1838 1838 1842 1850 1852 1860 1872 1872 1900 1916 1916 1927 1927 1930 1930 1931 1939 1940 1941 1943 1945 1945 1951 1953 1955 1966 1975 1985 1987 1988 2001 2003 2013 2028 2033 2044 2044 2046 2053 2055 2058 2058 2066 2068 2069 2070 2090 2090 2104 2106 2114 2114 2115 2116 2125 2131 2141 2147 2156 2159 2171 2174 2183 2183 2186 2191 2198 2217 2218 2218 2218 2221 2221 2221 2229 2229 2230 2239 2242 2244 2285 2287 2288 2294 2298 2300 2302 2311 2314 2315 2322 2331 2333 2334 2336 2342 2345 2346 2349 2357 2358 2381 2392 2398 2402 2403 2407 2415 2415 2416 2419 2426 2428 2429 2437 2449 2460 2461 2461 2463 2471 2477 2486 2494 2505 2509 2513 2517 2520 2528 2529 2531 2533 2534 2535 2550 2557 2566 2578 2578 2584 2589 2589 2592 2599 2614 2622 2622 2623 2624 2634 2636 2644 2644 2647 2648 2663 2679 2680 2681 2691 2694 2699 2702 2706 2715 2717 2720 2737 2747 2750 2750 2751 2752 2761 2762 2765 2776 2777 2778 2782 2795 2807 2807 2809 2819 2830 2830 2830 2862 2864 2865 2865 2865 2865 2878 2878 2878 2879 2880 2880 2880 2890 2890 2892 2896 2906 2922 2993"
N = 443

T = "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"

N = 1000
SPECTRUM =  [int (x) for x in T.split (' ')]
leaderBoardCycloPeptideSequencing(SPECTRUM, N)

exit()
#T="97-129-97-147-99-71-186-71-113-163-115-71-113-128-103-87-128-101-137-163-114"
#masses =  [int (x) for x in T.split ('-')]
#print "PEPTID", getPeptidFrommasses(masses)
T="0 97 97 129 194 196 226 226 244 258 323 323 452"
SPECTRUM =  [int (x) for x in T.split (' ')]

print "LPS PEEP", linearPeptidScorePeptid("PEEP", SPECTRUM)

T ="0 71 71 71 131 131 131 156 198 199 199 202 202 202 333 333 333 404 404"
SPECTRUM =  [int (x) for x in T.split (' ')]
print "SCR MAMA", scorePeptid("MAMA", SPECTRUM)
#print "LPS", linearPeptidScorePeptidFast("PEEP", SPECTRUM)

T = "GFAQHVMEGIGLDVKFTNIISCFFDHEWSTCHCKHHNSINHTMSMVF LIGDDDEADNCMMMVQSIKWKTLLRYGAFFTFPFYSYAILHVFYVLW KPMWWAFIFGFCDMKNCFDAPFWMHNSVQWEQHYRCNDVKMMSQLCW MAPRDIRMYFDKYHETAALDSQWIIQQIYHLMNVRKLNRTNRFTSVG FEKYHQQQILIDAQRVRLVHTVARAGPGWVQTGGWQQTCPRYKPYAW NVNPCERSSPPNFSWFMSFWADNSDYGDVIFCCPSVLRTMEMQSKKG WDTDTFFQKAMLKKDETADQIFNLRPYSLTCHNENILGNDNQEKQAG TLGSGENDKGHTVGAGHKGHPEREFEAPIERHEHPRVMMTKVGCYWI VCGHHHEQTVIMKAFDAWKVGFLGPIVAWVIFPAVYLWGKSLCPWTN YDSPTTYLSTHCHRLTNRMVHENPVICPPQDFAKYLIQSGWEFPLVA KDPINQTGDTNVRNFNVGCFCGCYFQWERHDGTPMHFWFSQKLSLTW HMKKLFWGIMKHHILFDFVNQPAFTNKAKGPTPHKAEELIRNLGQEK FNDRQRLVCHTNQCCAYKNKVVCSGGGSEISTNAHTYHFLALGHQVG MYYSAWTEPYYPPTLQIWWWYWKYGCTACQTGPHTMVFVQPTCKCVH YYGYRQCSWCQRWTVRRMLCWIDVLHKALHWHVCLLFHQALYGFSHE WASIGAIMRSAKDMYESLEFHKTHCTYFVYMVCKEARPGWTFFIEWV"


#T ="0 71 87 101 113 158 184 188 259 271 372"
#N = 2
#SPECTRUM =  [int (x) for x in T.split (' ')]

T = "GFAQHVMEGIGLDVKFTNIISCFFDHEWSTCHCKHHNSINHTMSMVF LIGDDDEADNCMMMVQSIKWKTLLRYGAFFTFPFYSYAILHVFYVLW KPMWWAFIFGFCDMKNCFDAPFWMHNSVQWEQHYRCNDVKMMSQLCW MAPRDIRMYFDKYHETAALDSQWIIQQIYHLMNVRKLNRTNRFTSVG FEKYHQQQILIDAQRVRLVHTVARAGPGWVQTGGWQQTCPRYKPYAW NVNPCERSSPPNFSWFMSFWADNSDYGDVIFCCPSVLRTMEMQSKKG WDTDTFFQKAMLKKDETADQIFNLRPYSLTCHNENILGNDNQEKQAG TLGSGENDKGHTVGAGHKGHPEREFEAPIERHEHPRVMMTKVGCYWI VCGHHHEQTVIMKAFDAWKVGFLGPIVAWVIFPAVYLWGKSLCPWTN YDSPTTYLSTHCHRLTNRMVHENPVICPPQDFAKYLIQSGWEFPLVA KDPINQTGDTNVRNFNVGCFCGCYFQWERHDGTPMHFWFSQKLSLTW HMKKLFWGIMKHHILFDFVNQPAFTNKAKGPTPHKAEELIRNLGQEK FNDRQRLVCHTNQCCAYKNKVVCSGGGSEISTNAHTYHFLALGHQVG MYYSAWTEPYYPPTLQIWWWYWKYGCTACQTGPHTMVFVQPTCKCVH YYGYRQCSWCQRWTVRRMLCWIDVLHKALHWHVCLLFHQALYGFSHE WASIGAIMRSAKDMYESLEFHKTHCTYFVYMVCKEARPGWTFFIEWV"
LEADERBOARD  =  [str (x) for x in T.split (' ')]

T ="0 71 87 97 97 99 101 101 101 101 101 103 103 113 113 113 113 113 113 113 113 114 114 115 115 128 128 129 129 129 129 129 129 131 131 131 137 137 147 147 156 156 156 163 163 163 186 186 198 200 200 202 202 204 214 214 216 216 216 226 228 230 234 242 242 242 242 242 242 244 245 253 257 257 257 259 260 266 266 268 269 271 276 276 276 276 278 283 284 285 287 293 294 299 301 303 313 317 317 317 327 329 331 343 343 347 354 355 355 356 359 363 363 370 370 371 371 372 379 382 384 388 389 397 397 400 405 407 408 408 413 413 415 415 416 418 418 420 422 428 430 430 432 434 439 444 455 456 458 459 473 476 484 484 485 485 487 499 500 501 506 507 510 510 511 513 515 518 522 526 527 528 529 529 531 533 533 537 540 541 544 544 545 547 558 559 562 568 569 571 572 574 585 586 588 597 598 607 610 612 616 619 620 624 624 625 626 631 636 641 644 646 650 657 657 660 662 663 669 669 670 671 672 674 675 678 681 684 685 689 691 692 696 698 699 700 700 701 711 733 733 735 738 738 739 741 753 771 772 772 773 775 775 775 778 779 782 783 783 786 788 789 794 794 797 798 798 800 801 802 804 805 806 808 810 813 815 828 837 840 846 854 862 864 864 866 869 882 884 885 886 888 889 899 901 902 902 903 904 907 908 911 911 914 914 924 925 926 928 935 935 937 937 941 941 945 952 955 961 961 975 975 977 984 987 988 992 995 999 1002 1013 1015 1016 1017 1017 1017 1022 1025 1032 1032 1038 1039 1039 1040 1044 1051 1055 1058 1058 1059 1062 1065 1066 1070 1070 1072 1074 1084 1088 1097 1099 1100 1101 1103 1104 1105 1106 1118 1118 1121 1130 1133 1135 1142 1150 1151 1153 1153 1154 1154 1156 1165 1171 1172 1186 1187 1194 1196 1198 1200 1201 1201 1202 1207 1212 1213 1214 1215 1216 1217 1217 1218 1218 1219 1231 1233 1234 1234 1236 1248 1255 1259 1264 1267 1272 1279 1281 1285 1298 1300 1301 1309 1311 1315 1315 1315 1316 1318 1319 1319 1321 1325 1330 1331 1334 1335 1338 1341 1343 1344 1346 1347 1348 1352 1363 1363 1372 1379 1395 1396 1398 1402 1410 1414 1414 1422 1422 1428 1429 1429 1430 1431 1432 1433 1434 1435 1435 1438 1447 1448 1450 1452 1460 1461 1471 1472 1472 1476 1476 1478 1481 1494 1499 1509 1515 1517 1521 1525 1532 1534 1535 1535 1535 1543 1546 1548 1551 1557 1557 1561 1561 1561 1562 1566 1566 1573 1577 1581 1585 1585 1589 1603 1605 1606 1608 1609 1609 1612 1618 1628 1634 1635 1636 1638 1647 1663 1664 1664 1664 1672 1672 1674 1674 1677 1680 1682 1686 1690 1694 1699 1702 1708 1713 1715 1716 1718 1718 1722 1723 1724 1740 1741 1743 1748 1749 1756 1763 1777 1777 1778 1781 1785 1792 1801 1803 1803 1805 1810 1811 1811 1815 1819 1819 1822 1823 1827 1830 1837 1837 1842 1844 1849 1850 1853 1869 1872 1876 1878 1878 1879 1887 1905 1906 1906 1912 1916 1916 1916 1924 1925 1929 1934 1936 1939 1940 1940 1940 1948 1948 1950 1965 1972 1979 1979 1991 1993 1996 2000 2000 2005 2007 2007 2009 2017 2019 2019 2026 2031 2034 2035 2039 2042 2053 2053 2058 2061 2064 2076 2079 2079 2087 2087 2094 2096 2097 2103 2113 2118 2120 2120 2123 2125 2132 2132 2135 2135 2136 2138 2140 2147 2150 2167 2171 2171 2192 2193 2200 2200 2200 2205 2209 2216 2216 2219 2224 2226 2226 2232 2233 2233 2233 2235 2237 2250 2252 2262 2266 2266 2268 2276 2284 2287 2296 2300 2303 2313 2318 2334 2334 2334 2336 2336 2337 2339 2347 2349 2355 2355 2356 2363 2363 2365 2366 2379 2380 2389 2397 2397 2413 2413 2416 2418 2434 2435 2437 2446 2447 2449 2449 2452 2466 2468 2468 2473 2476 2490 2492 2493 2494 2494 2494 2502 2510 2511 2519 2526 2533 2536 2542 2547 2549 2550 2553 2560 2578 2579 2579 2586 2597 2605 2605 2615 2620 2622 2623 2624 2624 2625 2627 2631 2632 2634 2648 2650 2650 2650 2651 2655 2689 2691 2706 2710 2711 2715 2733 2733 2733 2735 2738 2742 2742 2744 2747 2751 2753 2756 2761 2763 2768 2771 2779 2802 2807 2813 2820 2824 2825 2828 2836 2836 2836 2843 2846 2862 2862 2866 2871 2873 2876 2876 2881 2882 2900 2903 2907 2924 2933 2933 2937 2938 2944 2944 2949 2953 2957 2965 2972 2975 2975 2975 2987 2995 2999 3004 3018 3034 3034 3037 3037 3038 3039 3058 3062 3066 3066 3070 3073 3078 3078 3082 3088 3101 3104 3112 3119 3130 3135 3138 3149 3151 3152 3163 3163 3171 3172 3179 3190 3191 3193 3195 3217 3217 3220 3229 3238 3241 3241 3243 3244 3248 3250 3251 3276 3286 3291 3292 3300 3308 3318 3320 3321 3330 3335 3342 3351 3351 3354 3354 3357 3358 3358 3372 3387 3405 3407 3420 3422 3429 3431 3433 3433 3434 3439 3448 3455 3471 3483 3485 3486 3486 3488 3510 3514 3521 3533 3534 3534 3535 3546 3551 3552 3568 3585 3596 3599 3599 3600 3611 3611 3614 3615 3622 3634 3635 3647 3649 3664 3664 3681 3682 3696 3697 3708 3713 3727 3728 3728 3728 3735 3748 3750 3751 3771 3777 3797 3809 3812 3827 3828 3837 3841 3841 3842 3851 3857 3864 3864 3868 3868 3884 3898 3913 3940 3940 3942 3943 3955 3965 3969 3970 3970 3977 3981 4013 4014 4027 4027 4044 4053 4054 4056 4057 4083 4096 4098 4099 4110 4126 4140 4140 4140 4145 4155 4158 4167 4171 4184 4209 4211 4212 4223 4253 4255 4259 4268 4272 4284 4296 4296 4299 4303 4313 4352 4368 4373 4374 4397 4397 4397 4400 4409 4409 4416 4428 4465 4469 4487 4501 4510 4510 4526 4529 4538 4560 4566 4572 4584 4630 4639 4639 4639 4643 4651 4673 4673 4681 4685 4752 4752 4752 4768 4782 4786 4786 4802 4829 4853 4867 4881 4881 4883 4915 4915 4942 4968 4968 4982 4994 5028 5044 5069 5069 5071 5095 5097 5157 5157 5170 5184 5198 5210 5258 5270 5299 5311 5313 5371 5373 5412 5426 5474 5486 5527 5575 5587 5642 5688 5743 5844"
SPECTRUM =  [int (x) for x in T.split (' ')]
N = 5
#print "TRIM:", ' '.join (trimLeaderBoard(LEADERBOARD, SPECTRUM, N, 0, 0))

#print "TRIM1:", ' '.join (trimLeaderBoardFast(LEADERBOARD, SPECTRUM, N, 0, 0))

#T ="0 71 113 129 147 200 218 260 313 331 347 389 460"
#SPECTRUM =  [int (x) for x in T.split (' ')]



#print "SCORE", linearPeptidScorePeptid("YNPEPFVAWAIYDAIKCSKTH",  SPECTRUM)
#print "SCORE", linearPeptidScorePeptid("PEPFVAWAIYDAIKCSKTHYN",  SPECTRUM)

#print "SCOREN", scorePeptid("YNPEPFVAWAIYDAIKCSKTH",  SPECTRUM)
#print "SCOREN", scorePeptid("PEPFVAWAIYDAIKCSKTHYN",  SPECTRUM)

#print peptidAsMass("YNPEPFVAWAIYDAIKCSKTH")


#print linearPeptidScorePeptid( "NQEL", [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484])
#T="71-131-114-113-113-115-99-97-103-137"
#A =  [int (x) for x in T.split ('-')]

#print     convertMassesToPeptid(A)
#p = expand(SPECTRUM, p)
#print p
#print "LS:", linearSpectrum("NQEL")
#p = expand(SPECTRUM, p)
#print p
T = "0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333"

T ="0 97 97 99 101 103 196 198 198 200 202 295 297 299 299 301 394 396 398 400 400 497"
SPECTRUM =  [int (x) for x in T.split (' ')]

print "CONS", isPeptideConsistent([97, 103, 97], SPECTRUM)

pp = ["TVQ", "VAQ", "CTV", "ETC", "QCV", "TCE"]
for p in pp:
    #m =convertMassesToPeptid(p)
    m = peptidAsMasses(p)
    print "MMMM", m
    print isPeptideConsistent( m,SPECTRUM)

for i in ["TLAM", "MAIT", "MTAL", "MTAI", "TAIM", "ALTM"]:
    print "M2P", cycloSpectrum(i)

print pep
#print getNumOfCyclicPeptides(31315)
#print cycloSpectrum("LEQN")
#PEPTID = "NQEL"
#print peptidAsMasses (PEPTID)
#print peptidAsMass(PEPTID)
#print MASS_TABLE[MASS][19]
#print getNumOfPetids(1024)

#print getPeptidMass("VVV")
#res = cycloSpectrum(PEPTID)
#masses = res[0]
#print ' '.join (res[1])
#print parentMass(masses)
#print "SUM:" , sum (masses)
#print ' '.join ( str (v) for v in sorted (masses))

#print cycloPeptideSequencing(masses)
#print getLinearPeptideSizeFromSpectrum(41181)

#T = "GFAQHVMEGIGLDVKFTNIISCFFDHEWSTCHCKHHNSINHTMSMVF LIGDDDEADNCMMMVQSIKWKTLLRYGAFFTFPFYSYAILHVFYVLW KPMWWAFIFGFCDMKNCFDAPFWMHNSVQWEQHYRCNDVKMMSQLCW MAPRDIRMYFDKYHETAALDSQWIIQQIYHLMNVRKLNRTNRFTSVG FEKYHQQQILIDAQRVRLVHTVARAGPGWVQTGGWQQTCPRYKPYAW NVNPCERSSPPNFSWFMSFWADNSDYGDVIFCCPSVLRTMEMQSKKG WDTDTFFQKAMLKKDETADQIFNLRPYSLTCHNENILGNDNQEKQAG TLGSGENDKGHTVGAGHKGHPEREFEAPIERHEHPRVMMTKVGCYWI VCGHHHEQTVIMKAFDAWKVGFLGPIVAWVIFPAVYLWGKSLCPWTN YDSPTTYLSTHCHRLTNRMVHENPVICPPQDFAKYLIQSGWEFPLVA KDPINQTGDTNVRNFNVGCFCGCYFQWERHDGTPMHFWFSQKLSLTW HMKKLFWGIMKHHILFDFVNQPAFTNKAKGPTPHKAEELIRNLGQEK FNDRQRLVCHTNQCCAYKNKVVCSGGGSEISTNAHTYHFLALGHQVG MYYSAWTEPYYPPTLQIWWWYWKYGCTACQTGPHTMVFVQPTCKCVH YYGYRQCSWCQRWTVRRMLCWIDVLHKALHWHVCLLFHQALYGFSHE WASIGAIMRSAKDMYESLEFHKTHCTYFVYMVCKEARPGWTFFIEWV"
#LEADERBOARD  =  [str (x) for x in T.split (' ')]

#T ="0 71 87 97 97 99 101 101 101 101 101 103 103 113 113 113 113 113 113 113 113 114 114 115 115 128 128 129 129 129 129 129 129 131 131 131 137 137 147 147 156 156 156 163 163 163 186 186 198 200 200 202 202 204 214 214 216 216 216 226 228 230 234 242 242 242 242 242 242 244 245 253 257 257 257 259 260 266 266 268 269 271 276 276 276 276 278 283 284 285 287 293 294 299 301 303 313 317 317 317 327 329 331 343 343 347 354 355 355 356 359 363 363 370 370 371 371 372 379 382 384 388 389 397 397 400 405 407 408 408 413 413 415 415 416 418 418 420 422 428 430 430 432 434 439 444 455 456 458 459 473 476 484 484 485 485 487 499 500 501 506 507 510 510 511 513 515 518 522 526 527 528 529 529 531 533 533 537 540 541 544 544 545 547 558 559 562 568 569 571 572 574 585 586 588 597 598 607 610 612 616 619 620 624 624 625 626 631 636 641 644 646 650 657 657 660 662 663 669 669 670 671 672 674 675 678 681 684 685 689 691 692 696 698 699 700 700 701 711 733 733 735 738 738 739 741 753 771 772 772 773 775 775 775 778 779 782 783 783 786 788 789 794 794 797 798 798 800 801 802 804 805 806 808 810 813 815 828 837 840 846 854 862 864 864 866 869 882 884 885 886 888 889 899 901 902 902 903 904 907 908 911 911 914 914 924 925 926 928 935 935 937 937 941 941 945 952 955 961 961 975 975 977 984 987 988 992 995 999 1002 1013 1015 1016 1017 1017 1017 1022 1025 1032 1032 1038 1039 1039 1040 1044 1051 1055 1058 1058 1059 1062 1065 1066 1070 1070 1072 1074 1084 1088 1097 1099 1100 1101 1103 1104 1105 1106 1118 1118 1121 1130 1133 1135 1142 1150 1151 1153 1153 1154 1154 1156 1165 1171 1172 1186 1187 1194 1196 1198 1200 1201 1201 1202 1207 1212 1213 1214 1215 1216 1217 1217 1218 1218 1219 1231 1233 1234 1234 1236 1248 1255 1259 1264 1267 1272 1279 1281 1285 1298 1300 1301 1309 1311 1315 1315 1315 1316 1318 1319 1319 1321 1325 1330 1331 1334 1335 1338 1341 1343 1344 1346 1347 1348 1352 1363 1363 1372 1379 1395 1396 1398 1402 1410 1414 1414 1422 1422 1428 1429 1429 1430 1431 1432 1433 1434 1435 1435 1438 1447 1448 1450 1452 1460 1461 1471 1472 1472 1476 1476 1478 1481 1494 1499 1509 1515 1517 1521 1525 1532 1534 1535 1535 1535 1543 1546 1548 1551 1557 1557 1561 1561 1561 1562 1566 1566 1573 1577 1581 1585 1585 1589 1603 1605 1606 1608 1609 1609 1612 1618 1628 1634 1635 1636 1638 1647 1663 1664 1664 1664 1672 1672 1674 1674 1677 1680 1682 1686 1690 1694 1699 1702 1708 1713 1715 1716 1718 1718 1722 1723 1724 1740 1741 1743 1748 1749 1756 1763 1777 1777 1778 1781 1785 1792 1801 1803 1803 1805 1810 1811 1811 1815 1819 1819 1822 1823 1827 1830 1837 1837 1842 1844 1849 1850 1853 1869 1872 1876 1878 1878 1879 1887 1905 1906 1906 1912 1916 1916 1916 1924 1925 1929 1934 1936 1939 1940 1940 1940 1948 1948 1950 1965 1972 1979 1979 1991 1993 1996 2000 2000 2005 2007 2007 2009 2017 2019 2019 2026 2031 2034 2035 2039 2042 2053 2053 2058 2061 2064 2076 2079 2079 2087 2087 2094 2096 2097 2103 2113 2118 2120 2120 2123 2125 2132 2132 2135 2135 2136 2138 2140 2147 2150 2167 2171 2171 2192 2193 2200 2200 2200 2205 2209 2216 2216 2219 2224 2226 2226 2232 2233 2233 2233 2235 2237 2250 2252 2262 2266 2266 2268 2276 2284 2287 2296 2300 2303 2313 2318 2334 2334 2334 2336 2336 2337 2339 2347 2349 2355 2355 2356 2363 2363 2365 2366 2379 2380 2389 2397 2397 2413 2413 2416 2418 2434 2435 2437 2446 2447 2449 2449 2452 2466 2468 2468 2473 2476 2490 2492 2493 2494 2494 2494 2502 2510 2511 2519 2526 2533 2536 2542 2547 2549 2550 2553 2560 2578 2579 2579 2586 2597 2605 2605 2615 2620 2622 2623 2624 2624 2625 2627 2631 2632 2634 2648 2650 2650 2650 2651 2655 2689 2691 2706 2710 2711 2715 2733 2733 2733 2735 2738 2742 2742 2744 2747 2751 2753 2756 2761 2763 2768 2771 2779 2802 2807 2813 2820 2824 2825 2828 2836 2836 2836 2843 2846 2862 2862 2866 2871 2873 2876 2876 2881 2882 2900 2903 2907 2924 2933 2933 2937 2938 2944 2944 2949 2953 2957 2965 2972 2975 2975 2975 2987 2995 2999 3004 3018 3034 3034 3037 3037 3038 3039 3058 3062 3066 3066 3070 3073 3078 3078 3082 3088 3101 3104 3112 3119 3130 3135 3138 3149 3151 3152 3163 3163 3171 3172 3179 3190 3191 3193 3195 3217 3217 3220 3229 3238 3241 3241 3243 3244 3248 3250 3251 3276 3286 3291 3292 3300 3308 3318 3320 3321 3330 3335 3342 3351 3351 3354 3354 3357 3358 3358 3372 3387 3405 3407 3420 3422 3429 3431 3433 3433 3434 3439 3448 3455 3471 3483 3485 3486 3486 3488 3510 3514 3521 3533 3534 3534 3535 3546 3551 3552 3568 3585 3596 3599 3599 3600 3611 3611 3614 3615 3622 3634 3635 3647 3649 3664 3664 3681 3682 3696 3697 3708 3713 3727 3728 3728 3728 3735 3748 3750 3751 3771 3777 3797 3809 3812 3827 3828 3837 3841 3841 3842 3851 3857 3864 3864 3868 3868 3884 3898 3913 3940 3940 3942 3943 3955 3965 3969 3970 3970 3977 3981 4013 4014 4027 4027 4044 4053 4054 4056 4057 4083 4096 4098 4099 4110 4126 4140 4140 4140 4145 4155 4158 4167 4171 4184 4209 4211 4212 4223 4253 4255 4259 4268 4272 4284 4296 4296 4299 4303 4313 4352 4368 4373 4374 4397 4397 4397 4400 4409 4409 4416 4428 4465 4469 4487 4501 4510 4510 4526 4529 4538 4560 4566 4572 4584 4630 4639 4639 4639 4643 4651 4673 4673 4681 4685 4752 4752 4752 4768 4782 4786 4786 4802 4829 4853 4867 4881 4881 4883 4915 4915 4942 4968 4968 4982 4994 5028 5044 5069 5069 5071 5095 5097 5157 5157 5170 5184 5198 5210 5258 5270 5299 5311 5313 5371 5373 5412 5426 5474 5486 5527 5575 5587 5642 5688 5743 5844"
#SPECTRUM =  [int (x) for x in T.split (' ')]
#N = 5
#print ' '.join (trimLeaderBoard(LEADERBOARD, SPECTRUM, N, 0, 0))