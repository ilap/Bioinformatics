__author__ = 'ilap'

from SequenceCompareLib import *

import numpy as np
def localAlignmentGraph (v, w, ro=0, mu=0):

    inc=1
    v_len = len (v)+1
    w_len = len (w)+1



    s = np.zeros (shape=(v_len, w_len))
    #for i in range (v_len):
    #    s[i,0] = i*-5
    #for i in range (w_len):
    #    s[0,i] = i*-5


    (row, col, matrix) = generateScoreMatrixFromFile("PAM250.txt")

    max_i = -sys.maxint
    max_j = -sys.maxint
    max_v = -sys.maxint
    for i in range (1, v_len):
        for j in range (1, w_len):
            v_char = v[i-1]
            w_char = w[j-1]
            v_idx = row.index(v_char)
            w_idx = col.index(w_char)

            down = s[i-1,j] - ro
            left = s[i, j-1] - ro
            downleft = s[i-1, j-1] + matrix[v_idx, w_idx]


            max_t = max (0, down, left, downleft)
            if max_t > max_v:
                max_v = max_t
                max_i = i
                max_j = j
            s[i,j] = max_t



            #print i, j, v[i-1], w[j-1], down, left, downleft, matrix[v_idx, w_idx]


    return (s, row, col, matrix, max_i, max_j)


global _v_alignment
_v_alignment = str ()

def backTrackLocalScoreMatrix (backtrack, v, w, i, j, ro, row, col, matrix):

    global _v_alignment
    v_len = len (v)
    w_len = len (w)
    if v_len >= w_len:
        if i == 0:
            print _v_alignment
            return ""
    else:
        if j == 0:
            print _v_alignment
            return ""

    v_idx= row.index(v[i-1])
    w_idx= col.index(w[j-1])

    diag = backtrack[i-1,j-1]+ matrix[v_idx, w_idx]
    down = backtrack[i-1,j] - ro
    right = backtrack[i,j-1] - ro

    max_v = max (0, down, right, diag)

    print "(", i,",",j,") ",max_v, " DRG", down, right, diag,

    if not max_v:
        print _v_alignment
        return ""
    elif max_v == down:
        print "down"
        _v_alignment = v[i-1] + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, i-1, j, ro, row, col, matrix) + "-"
    elif max_v == right:
        print "right"
        _v_alignment = "-" + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, i, j-1, ro, row, col, matrix) + w[j-1]
    else:
        print "diag"
        _v_alignment = v[i-1] + _v_alignment
        return backTrackLocalScoreMatrix(backtrack, v, w, i-1, j-1, ro, row, col, matrix) + w[j-1]




v="HIWWFDNQGKEPNTIWHTETGMLGNGEGISNFMYRFPLWQAGDKVWNFCEMPTELKYTKNEDKHCDIITTADVCCHYPWDFDLQFQNGDPCPIASITYYHTPRNFFFEMPMYPALSWKCTKLMYKVQPLYKTKVGWDDAGLNSLHVGSVKMHRLTASRCWLLSDYMDIQSTVIVITCFAWFNVTCQHFWFPGHMFLQFAMKCDCHKSYSWNITMMESYFCVIHNIGGDFFQRPWYYQLHLCCCNAAEAHYLVLCETFILMNDLSHMVLQFKLPSVWAELHEHIGLQLNLCWFDRSPQRRRIMWENSLETLPYKIYHIPMLMEDKRIEHNSFSSGPFYWPCSDNHCGAKMSFFLQAFNPKMQVRWNFRIWKLMASGLRKHCALLDGYGAHCSDFGDLCQQYLPPCTLRKAFMMKDGPQKGEMQYINIWWERFQWVGAPIQKVWPSTCSEKHTEVFGVSGDPWPLQKRSGAPAPINTNTHIHKWTESTDQGQTAMWAPWMKPKTCFQVYFVACGIFQRWPKCCKYCPDCSGTLVMHGERRSHPFNKPWMRSIGFTYPSIVREDHLFTCWWETVTFMQYMEDFRGQRRFPEEEDGMWNMKKSYGYYLEYNFPVMKMDTFMTTAQSRTTGCEISDHTVVCAGGRNKQAEFQVLTSYGEFMPIHRGQQAYFTGYKGMALTTPWTYQPRYPVLHDAKEEDRVCFDVKMRINWHCLFEKNRTWAHAPCGGKTVQMGMWLVVFRNRFITIVQNKWKWDCLHFLSGYSCDSHHVDIAMCTAPHFYMNISRSPQFEEPWHCVPECDTYIGHYPPKISGVNTPTWMNGVIAMVQTQRASSGAHAEMVRDTSSHRRNSANCQYYRNWCPNTAMNWPDYETFLDGTHKRKARPGQRESFDRNMWMCMCAASLQ"
w="CCVTATAMYSLQYDSRSEIAMCHWKHFVHMGTMMIGNQSNSEGALTEKLFYNAVPALYMDSGSMWIYHIWVEKYYEWAGIQGPTRGWSLYQGYFMFGMIEMKRSRWFFEPEVRMSPHEEQELHLQQVLKMDACKEHFMRLVQRQINMTENPTMVADHWQHEECNIRNLDPLYPDKVYQCGDKFNEVEAGFKINRPIMCESPEYMLAVWMAESYEYDLWCQVFWVSYYDEAKFLFVALQFNCQYAHYNKEIYMAHQHPRSPQIKKLRGHQDFWENSHETSSGLKGMYHYPMLMHNSSTSGPPCSDNHCGGFTVFNQCYQMSFTLQAFNPKMTMVRWNFRIWFQTLIWVELMASGLRKHCALLGARDPGEDNGDLCQQYLAFKAMCTLRWAFMMKDGPQTYETAGKGEMQYENWVGQKVFCSTCSEKHTEVFGFMKGDPWPLQKCGAPFPINTNTFCVCHIEIHKWTESTDMDTYDASTQGQTAMWAPWMKPKTCFQACSHNPLEYFEWWNMGYAAFGIFQRWPYCCKYCPDCSGTKVMEGERNNFKNYNPFNTPWMRSIGFTYPSIVTEDSLFTCWWETRGMRRCWSTPVNEEEDFMWKKSYGYYLDTPDGCCPDTTGCEISTVNSAGGRNKQAEFQVKHRGQQAYFFSSCYTVLGKDMMMQGGCPMIRTVWVDIHFFLEVHDFAPRHKHGCLPPRKVAAIMDAPYMFTDENKPNRGPQHGAYNNKDGECAMPVCPYEAMMPEDDYKVFSLSDHMWQKHAKRKEIALRAPNYLNGWADMNCIEMPLGERKMQHVIKTMSIAPRAISAVIQGRFGFWYTRVHYQGLYSHVWKFQAQFDAPTIKWQCNCCWQCKFVYTGWTLIFQIFPPWWEYMHMYNVAAGCGPLAHEDLEYMTDSDLA"
#generateScoreMatrixFromFile("BLOSUM62.txt")

v="MEANLY"
w="PENALTY"

ro = 5
(b, row, col, matrix, max_i, max_j) =  localAlignmentGraph(v, w, ro)
print b
print int (b[max_i,max_j]), max_i, max_j

print ''.join (backTrackLocalScoreMatrix(b, v, w, max_i, max_j, ro, row, col, matrix))


