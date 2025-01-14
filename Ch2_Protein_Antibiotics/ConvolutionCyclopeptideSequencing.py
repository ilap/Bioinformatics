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

OMASS_TABLE = \
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

def getPeptidFromMasses (masses):
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
        for mass in MASS_TABLE[MASS]:
                a = [mass]
                result.append (a)
    else:
        for i in peptides:
            summa = sum (i)

            #print"SUM:", summa
            for mass in MASS_TABLE[MASS]:
                #print "LENI", i[len (i)-1], mass
                tmp = i[:] # or list (i)
                #if mass not in i and mass + summa in spectrum:
                    #print "MASS+SUM", mass + summa

                tmp.append (mass)
                    #tmp.append (mass)
                    #print "TMP:", tmp, i
                result.append (tmp)
    #print "expandPeptides", result
    return result

def expandStringPeptidesAlphabe (peptides, alphabet):
    result = []

    if peptides == []:
        for peptide in alphabet[PEPTID]:
                result.append (peptide)
    else:
        for i in peptides:
            #summa = sum (i)

            #print"SUM:", summa

            for peptide in alphabet[PEPTID]:
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

def expandStringPeptides (peptides):
    result = []

    if peptides == []:
        for peptide in MASS_TABLE[PEPTID]:
                result.append (peptide)
    else:
        for i in peptides:
            #summa = sum (i)

            #print"SUM:", summa

            for peptide in MASS_TABLE[PEPTID]:
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


        if len (newpeptids) == 0:
            break
        peptides = expandPeptides(newpeptids)


    return result

def trimLeaderBoard (leader_board, spectrum, N, amino_acid = 0, amino_acid_mass = 0):

    plen = len (leader_board)
    linear_scores = []

    idx = 0
    for peptide in leader_board:

        score = linearPeptidScorePeptid(peptide, spectrum)
        linear_scores.append([score])
        linear_scores[idx].append (peptide)
        idx += 1

    linear_scores = sorted (linear_scores,None,None,True)

    for j in range (N, len (linear_scores)):
        if linear_scores[j][0] < linear_scores[N-1][0]:
            leader_board = [lb[1] for lb in linear_scores[:j]]

            return leader_board

    leader_board = [lb[1] for lb in linear_scores]

    return leader_board

'''
Banch-and-bound algorythm for peptide sequencing.
'''
@timing
def leaderBoardCycloPeptideSequencing (spectrum, N):

    leader_board = ['']
    leader_peptide = ""
    result = []
    print "MASS:", MASS_TABLE
    while len (leader_board) != 0:
        leader_board = expandStringPeptides(leader_board)
        new_leader_board = leader_board[:]

        for peptide in leader_board:

            peptide_mass = peptidAsMass(peptide)
            parent_mass = parentMass(spectrum)

            if  peptide_mass == parent_mass:
                ps = scorePeptid(peptide, spectrum)
                ls = scorePeptid(leader_peptide, spectrum)

                #print "SCORE:", ps, ls
                # REMOVE = if you want the highest....
                if  ps > ls :
                    print "SCORE:", ps, ls
                    print '-'.join ( str (x) for x in peptidAsMasses(peptide))
                    leader_peptide = peptide
            elif peptide_mass > parent_mass:
                new_leader_board.remove(peptide)

        leader_board = trimLeaderBoard (new_leader_board, spectrum, N)


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
    #print "SORTED", sorted_spectrum[1]
    for idx in range (0, slen - 1):
        tmp = []
        for j in range (slen - 1, idx, -1):

            val = sorted_spectrum[j] - sorted_spectrum [idx]
            if val == 0: continue
            tmp.append (val)
            result.append(val)
        #print tmp
    #print [str (x) for x in result]
    return  result

def buildAlphaBet (arr):
    base = 33
    alen = len (arr)
    global MASS_TABLE
    MASS_TABLE = [['']*alen,[0]*alen]
    for i in range (0, alen):
        MASS_TABLE[0][i] = chr (base + i)
        MASS_TABLE[1][i] = arr[i]


def convolutionCyclycPeptideSequencing (spectrum, M, N):

    c = [x for x in convolution(spectrum) if x <= 200 and x >= 57]

    from collections import Counter
    x =  Counter (c)

    conv =  x.most_common(M)

    import operator
    ss =  sorted(x.items(), key=operator.itemgetter(1), reverse=True)

    if len (ss) <= M:
        return dict (ss).keys ()

    alphabet = []
    needed = -1
    for i in range (0, len (ss)):
        (key, value) = ss[i]
        #print "KV", i, needed, key, value
        if i > M and value < needed:
            break
        elif i == M:
            needed = value
        else:
            alphabet.append(key)

    buildAlphaBet(alphabet)

    print spectrum, N
    print leaderBoardCycloPeptideSequencing(spectrum, N)


### MAIN
# T="0 71 71 87 87 97 97 99 101 101 101 101 103 113 113 113 113 113 115 115 115 115 128 128 129 129 131 131 137 147 156 156 156 163 163 163 163 172 184 186 186 199 204 210 212 212 214 214 214 216 218 226 227 228 230 231 232 243 244 244 257 259 260 260 262 273 276 278 278 281 285 287 292 293 293 301 302 303 312 314 319 325 326 327 327 327 327 328 331 332 341 358 360 364 373 374 374 375 377 389 390 391 393 394 398 402 403 407 407 409 415 416 418 425 428 429 440 440 441 441 444 445 449 456 464 465 472 475 488 490 492 494 497 499 502 503 503 504 504 505 505 512 516 517 520 521 521 531 540 541 544 550 554 556 559 560 569 570 572 573 577 578 587 592 592 601 603 605 612 613 617 618 618 619 621 621 622 625 628 634 634 653 655 657 659 661 669 671 672 675 678 680 684 688 691 696 701 704 705 707 714 715 716 716 718 721 722 723 724 732 734 735 743 748 749 749 756 768 769 775 776 781 781 785 791 792 799 800 804 806 808 809 816 819 822 822 825 829 829 830 835 836 836 843 846 847 850 852 852 870 872 877 879 887 890 895 900 904 905 905 906 912 913 919 919 921 923 928 931 932 935 937 938 942 942 944 947 948 948 953 959 961 965 976 985 985 985 992 999 1005 1008 1013 1015 1016 1020 1024 1024 1028 1032 1035 1036 1036 1036 1043 1046 1049 1050 1051 1061 1061 1062 1062 1066 1068 1069 1073 1075 1076 1089 1090 1098 1100 1114 1117 1122 1123 1131 1132 1133 1133 1133 1137 1137 1137 1139 1141 1148 1149 1152 1155 1160 1162 1162 1163 1174 1174 1176 1179 1183 1191 1198 1199 1201 1201 1204 1219 1220 1224 1225 1229 1232 1232 1234 1244 1245 1248 1250 1250 1252 1252 1253 1261 1261 1262 1263 1264 1265 1273 1280 1288 1296 1302 1302 1305 1312 1313 1316 1316 1318 1325 1332 1335 1335 1339 1345 1347 1349 1350 1351 1354 1354 1360 1362 1363 1363 1367 1369 1376 1377 1387 1388 1391 1392 1392 1393 1397 1403 1410 1415 1417 1425 1425 1428 1433 1434 1440 1444 1451 1460 1463 1464 1464 1464 1464 1464 1465 1469 1475 1477 1479 1488 1491 1492 1497 1501 1502 1505 1510 1520 1525 1525 1526 1526 1530 1531 1534 1540 1540 1543 1547 1548 1556 1561 1564 1564 1566 1577 1578 1579 1588 1590 1590 1592 1593 1597 1603 1606 1619 1623 1626 1627 1627 1627 1628 1634 1635 1637 1637 1638 1639 1641 1648 1648 1653 1662 1665 1671 1674 1676 1677 1681 1689 1693 1693 1695 1703 1706 1708 1719 1724 1724 1727 1729 1734 1735 1738 1740 1742 1742 1750 1752 1752 1753 1754 1754 1756 1760 1763 1765 1766 1778 1784 1785 1790 1792 1794 1804 1804 1806 1808 1821 1821 1823 1832 1837 1837 1837 1851 1853 1853 1853 1855 1857 1857 1862 1865 1866 1866 1867 1871 1871 1875 1879 1879 1883 1885 1889 1890 1891 1897 1905 1908 1916 1919 1923 1928 1933 1936 1941 1941 1949 1952 1964 1966 1966 1966 1966 1967 1968 1968 1977 1980 1984 1984 1985 1988 1992 2000 2004 2004 2009 2010 2012 2020 2020 2022 2028 2031 2034 2034 2036 2037 2041 2046 2051 2052 2053 2056 2065 2079 2079 2081 2081 2087 2091 2093 2097 2097 2105 2111 2113 2115 2116 2122 2123 2124 2129 2131 2133 2135 2138 2139 2140 2147 2148 2152 2153 2159 2160 2165 2168 2168 2169 2178 2180 2184 2187 2192 2197 2198 2206 2206 2209 2215 2224 2226 2237 2237 2240 2244 2244 2244 2247 2250 2253 2255 2260 2263 2264 2266 2268 2269 2277 2280 2286 2287 2293 2293 2294 2296 2297 2300 2301 2307 2309 2310 2311 2315 2315 2324 2334 2337 2339 2341 2352 2354 2356 2359 2369 2378 2378 2382 2383 2384 2386 2392 2393 2396 2397 2399 2400 2400 2406 2407 2413 2416 2424 2425 2427 2429 2430 2433 2438 2440 2443 2446 2449 2449 2449 2453 2456 2456 2467 2469 2478 2484 2487 2487 2495 2496 2501 2506 2509 2513 2515 2524 2525 2525 2528 2533 2534 2540 2541 2545 2546 2553 2554 2555 2558 2560 2562 2564 2569 2570 2571 2577 2578 2580 2582 2588 2596 2596 2600 2602 2606 2612 2612 2614 2614 2628 2637 2640 2641 2642 2647 2652 2656 2657 2659 2659 2662 2665 2671 2673 2673 2681 2683 2684 2689 2689 2693 2701 2705 2708 2709 2709 2713 2716 2725 2725 2726 2727 2727 2727 2727 2729 2741 2744 2752 2752 2757 2760 2765 2770 2774 2777 2785 2788 2796 2802 2803 2804 2808 2810 2814 2814 2818 2822 2822 2826 2827 2827 2828 2831 2836 2836 2838 2840 2840 2840 2842 2856 2856 2856 2861 2870 2872 2872 2885 2887 2889 2889 2899 2901 2903 2908 2909 2915 2927 2928 2930 2933 2937 2939 2939 2940 2941 2941 2943 2951 2951 2953 2955 2958 2959 2964 2966 2969 2969 2974 2985 2987 2990 2998 3000 3000 3004 3012 3016 3017 3019 3022 3028 3031 3040 3045 3045 3052 3054 3055 3056 3056 3058 3059 3065 3066 3066 3066 3067 3070 3074 3087 3090 3096 3100 3101 3103 3103 3105 3114 3115 3116 3127 3129 3129 3132 3137 3145 3146 3150 3153 3153 3159 3162 3163 3167 3167 3168 3168 3173 3183 3188 3191 3192 3196 3201 3202 3205 3214 3216 3218 3224 3228 3229 3229 3229 3229 3229 3230 3233 3242 3249 3253 3259 3260 3265 3268 3268 3276 3278 3283 3290 3296 3300 3301 3301 3302 3305 3306 3316 3317 3324 3326 3330 3330 3331 3333 3339 3339 3342 3343 3344 3346 3348 3354 3358 3358 3361 3368 3375 3377 3377 3380 3381 3388 3391 3391 3397 3405 3413 3420 3428 3429 3430 3431 3432 3432 3440 3441 3441 3443 3443 3445 3448 3449 3459 3461 3461 3464 3468 3469 3473 3474 3489 3492 3492 3494 3495 3502 3510 3514 3517 3519 3519 3530 3531 3531 3533 3538 3541 3544 3545 3552 3554 3556 3556 3556 3560 3560 3560 3561 3562 3570 3571 3576 3579 3593 3595 3603 3604 3617 3618 3620 3624 3625 3627 3631 3631 3632 3632 3642 3643 3644 3647 3650 3657 3657 3657 3658 3661 3665 3669 3669 3673 3677 3678 3680 3685 3688 3694 3701 3708 3708 3708 3717 3728 3732 3734 3740 3745 3745 3746 3749 3751 3751 3755 3756 3758 3761 3762 3765 3770 3772 3774 3774 3780 3781 3787 3788 3788 3789 3793 3798 3803 3806 3814 3816 3821 3823 3841 3841 3843 3846 3847 3850 3857 3857 3858 3863 3864 3864 3868 3871 3871 3874 3877 3884 3885 3887 3889 3893 3894 3901 3902 3908 3912 3912 3917 3918 3924 3925 3937 3944 3944 3945 3950 3958 3959 3961 3969 3970 3971 3972 3975 3977 3977 3978 3979 3986 3988 3989 3992 3997 4002 4005 4009 4013 4015 4018 4021 4022 4024 4032 4034 4036 4038 4040 4059 4059 4065 4068 4071 4072 4072 4074 4075 4075 4076 4080 4081 4088 4090 4092 4101 4101 4106 4115 4116 4120 4121 4123 4124 4133 4134 4137 4139 4143 4149 4152 4153 4162 4172 4172 4173 4176 4177 4181 4188 4188 4189 4189 4190 4190 4191 4194 4196 4199 4201 4203 4205 4218 4221 4228 4229 4237 4244 4248 4249 4252 4252 4253 4253 4264 4265 4268 4275 4277 4278 4284 4286 4286 4290 4291 4295 4299 4300 4302 4303 4304 4316 4318 4319 4319 4320 4329 4333 4335 4352 4361 4362 4365 4366 4366 4366 4366 4367 4368 4374 4379 4381 4390 4391 4392 4400 4400 4401 4406 4408 4412 4415 4415 4417 4420 4431 4433 4433 4434 4436 4449 4449 4450 4461 4462 4463 4465 4466 4467 4475 4477 4479 4479 4479 4481 4481 4483 4489 4494 4507 4507 4509 4521 4530 4530 4530 4530 4537 4537 4537 4546 4556 4562 4562 4564 4564 4565 4565 4578 4578 4578 4578 4580 4580 4580 4580 4580 4590 4592 4592 4592 4592 4594 4596 4596 4606 4606 4622 4622 4693"
#SPECTRUM = [int (x) for x in T.split (' ')]
#PEP = "VYYEVDWTMGRQIDPDEYPIAQCTRHRATILTLPDWQM"
# GOOD
# GOOD print scorePeptid(PEP, SPECTRUM)
#import datetimeprint datetime.datetime.now ().time ()
#T="0 71 87 97 97 99 101 101 101 101 101 103 103 113 113 113 113 113 113 113 113 114 114 115 115 128 128 129 129 129 129 129 129 131 131 131 137 137 147 147 156 156 156 163 163 163 186 186 198 200 200 202 202 204 214 214 216 216 216 226 228 230 234 242 242 242 242 242 242 244 245 253 257 257 257 259 260 266 266 268 269 271 276 276 276 276 278 283 284 285 287 293 294 299 301 303 313 317 317 317 327 329 331 343 343 347 354 355 355 356 359 363 363 370 370 371 371 372 379 382 384 388 389 397 397 400 405 407 408 408 413 413 415 415 416 418 418 420 422 428 430 430 432 434 439 444 455 456 458 459 473 476 484 484 485 485 487 499 500 501 506 507 510 510 511 513 515 518 522 526 527 528 529 529 531 533 533 537 540 541 544 544 545 547 558 559 562 568 569 571 572 574 585 586 588 597 598 607 610 612 616 619 620 624 624 625 626 631 636 641 644 646 650 657 657 660 662 663 669 669 670 671 672 674 675 678 681 684 685 689 691 692 696 698 699 700 700 701 711 733 733 735 738 738 739 741 753 771 772 772 773 775 775 775 778 779 782 783 783 786 788 789 794 794 797 798 798 800 801 802 804 805 806 808 810 813 815 828 837 840 846 854 862 864 864 866 869 882 884 885 886 888 889 899 901 902 902 903 904 907 908 911 911 914 914 924 925 926 928 935 935 937 937 941 941 945 952 955 961 961 975 975 977 984 987 988 992 995 999 1002 1013 1015 1016 1017 1017 1017 1022 1025 1032 1032 1038 1039 1039 1040 1044 1051 1055 1058 1058 1059 1062 1065 1066 1070 1070 1072 1074 1084 1088 1097 1099 1100 1101 1103 1104 1105 1106 1118 1118 1121 1130 1133 1135 1142 1150 1151 1153 1153 1154 1154 1156 1165 1171 1172 1186 1187 1194 1196 1198 1200 1201 1201 1202 1207 1212 1213 1214 1215 1216 1217 1217 1218 1218 1219 1231 1233 1234 1234 1236 1248 1255 1259 1264 1267 1272 1279 1281 1285 1298 1300 1301 1309 1311 1315 1315 1315 1316 1318 1319 1319 1321 1325 1330 1331 1334 1335 1338 1341 1343 1344 1346 1347 1348 1352 1363 1363 1372 1379 1395 1396 1398 1402 1410 1414 1414 1422 1422 1428 1429 1429 1430 1431 1432 1433 1434 1435 1435 1438 1447 1448 1450 1452 1460 1461 1471 1472 1472 1476 1476 1478 1481 1494 1499 1509 1515 1517 1521 1525 1532 1534 1535 1535 1535 1543 1546 1548 1551 1557 1557 1561 1561 1561 1562 1566 1566 1573 1577 1581 1585 1585 1589 1603 1605 1606 1608 1609 1609 1612 1618 1628 1634 1635 1636 1638 1647 1663 1664 1664 1664 1672 1672 1674 1674 1677 1680 1682 1686 1690 1694 1699 1702 1708 1713 1715 1716 1718 1718 1722 1723 1724 1740 1741 1743 1748 1749 1756 1763 1777 1777 1778 1781 1785 1792 1801 1803 1803 1805 1810 1811 1811 1815 1819 1819 1822 1823 1827 1830 1837 1837 1842 1844 1849 1850 1853 1869 1872 1876 1878 1878 1879 1887 1905 1906 1906 1912 1916 1916 1916 1924 1925 1929 1934 1936 1939 1940 1940 1940 1948 1948 1950 1965 1972 1979 1979 1991 1993 1996 2000 2000 2005 2007 2007 2009 2017 2019 2019 2026 2031 2034 2035 2039 2042 2053 2053 2058 2061 2064 2076 2079 2079 2087 2087 2094 2096 2097 2103 2113 2118 2120 2120 2123 2125 2132 2132 2135 2135 2136 2138 2140 2147 2150 2167 2171 2171 2192 2193 2200 2200 2200 2205 2209 2216 2216 2219 2224 2226 2226 2232 2233 2233 2233 2235 2237 2250 2252 2262 2266 2266 2268 2276 2284 2287 2296 2300 2303 2313 2318 2334 2334 2334 2336 2336 2337 2339 2347 2349 2355 2355 2356 2363 2363 2365 2366 2379 2380 2389 2397 2397 2413 2413 2416 2418 2434 2435 2437 2446 2447 2449 2449 2452 2466 2468 2468 2473 2476 2490 2492 2493 2494 2494 2494 2502 2510 2511 2519 2526 2533 2536 2542 2547 2549 2550 2553 2560 2578 2579 2579 2586 2597 2605 2605 2615 2620 2622 2623 2624 2624 2625 2627 2631 2632 2634 2648 2650 2650 2650 2651 2655 2689 2691 2706 2710 2711 2715 2733 2733 2733 2735 2738 2742 2742 2744 2747 2751 2753 2756 2761 2763 2768 2771 2779 2802 2807 2813 2820 2824 2825 2828 2836 2836 2836 2843 2846 2862 2862 2866 2871 2873 2876 2876 2881 2882 2900 2903 2907 2924 2933 2933 2937 2938 2944 2944 2949 2953 2957 2965 2972 2975 2975 2975 2987 2995 2999 3004 3018 3034 3034 3037 3037 3038 3039 3058 3062 3066 3066 3070 3073 3078 3078 3082 3088 3101 3104 3112 3119 3130 3135 3138 3149 3151 3152 3163 3163 3171 3172 3179 3190 3191 3193 3195 3217 3217 3220 3229 3238 3241 3241 3243 3244 3248 3250 3251 3276 3286 3291 3292 3300 3308 3318 3320 3321 3330 3335 3342 3351 3351 3354 3354 3357 3358 3358 3372 3387 3405 3407 3420 3422 3429 3431 3433 3433 3434 3439 3448 3455 3471 3483 3485 3486 3486 3488 3510 3514 3521 3533 3534 3534 3535 3546 3551 3552 3568 3585 3596 3599 3599 3600 3611 3611 3614 3615 3622 3634 3635 3647 3649 3664 3664 3681 3682 3696 3697 3708 3713 3727 3728 3728 3728 3735 3748 3750 3751 3771 3777 3797 3809 3812 3827 3828 3837 3841 3841 3842 3851 3857 3864 3864 3868 3868 3884 3898 3913 3940 3940 3942 3943 3955 3965 3969 3970 3970 3977 3981 4013 4014 4027 4027 4044 4053 4054 4056 4057 4083 4096 4098 4099 4110 4126 4140 4140 4140 4145 4155 4158 4167 4171 4184 4209 4211 4212 4223 4253 4255 4259 4268 4272 4284 4296 4296 4299 4303 4313 4352 4368 4373 4374 4397 4397 4397 4400 4409 4409 4416 4428 4465 4469 4487 4501 4510 4510 4526 4529 4538 4560 4566 4572 4584 4630 4639 4639 4639 4643 4651 4673 4673 4681 4685 4752 4752 4752 4768 4782 4786 4786 4802 4829 4853 4867 4881 4881 4883 4915 4915 4942 4968 4968 4982 4994 5028 5044 5069 5069 5071 5095 5097 5157 5157 5170 5184 5198 5210 5258 5270 5299 5311 5313 5371 5373 5412 5426 5474 5486 5527 5575 5587 5642 5688 5743 5844"

#SPECTRUM = [int (x) for x in T.split (' ')]
#P = "GFAQHVMEGIGLDVKFTNIISCFFDHEWSTCHCKHHNSINHTMSMVF LIGDDDEADNCMMMVQSIKWKTLLRYGAFFTFPFYSYAILHVFYVLW KPMWWAFIFGFCDMKNCFDAPFWMHNSVQWEQHYRCNDVKMMSQLCW MAPRDIRMYFDKYHETAALDSQWIIQQIYHLMNVRKLNRTNRFTSVG FEKYHQQQILIDAQRVRLVHTVARAGPGWVQTGGWQQTCPRYKPYAW NVNPCERSSPPNFSWFMSFWADNSDYGDVIFCCPSVLRTMEMQSKKG WDTDTFFQKAMLKKDETADQIFNLRPYSLTCHNENILGNDNQEKQAG TLGSGENDKGHTVGAGHKGHPEREFEAPIERHEHPRVMMTKVGCYWI VCGHHHEQTVIMKAFDAWKVGFLGPIVAWVIFPAVYLWGKSLCPWTN YDSPTTYLSTHCHRLTNRMVHENPVICPPQDFAKYLIQSGWEFPLVA KDPINQTGDTNVRNFNVGCFCGCYFQWERHDGTPMHFWFSQKLSLTW HMKKLFWGIMKHHILFDFVNQPAFTNKAKGPTPHKAEELIRNLGQEK FNDRQRLVCHTNQCCAYKNKVVCSGGGSEISTNAHTYHFLALGHQVG MYYSAWTEPYYPPTLQIWWWYWKYGCTACQTGPHTMVFVQPTCKCVH YYGYRQCSWCQRWTVRRMLCWIDVLHKALHWHVCLLFHQALYGFSHE WASIGAIMRSAKDMYESLEFHKTHCTYFVYMVCKEARPGWTFFIEWV"
#PEP = [str (x) for x in P.split (' ')]
# GOOD print ' '.join (trimLeaderBoard(PEP, SPECTRUM ,5))
# GOOD print ' '.join (str (x) for x in linearSpectrum("WARLYHPARDKAEWDTTMKHYYNIDNPLSNYENWFWMHNRELDCYD"))


T="0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"
M=20
N=1000
SPECTRUM = [int (x) for x in T.split (' ')]
print convolutionCyclycPeptideSequencing (SPECTRUM, M, N)
