__author__ = 'ilap'

import time

def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0)
        return ret
    return wrap

global CODONS
global RNA_TABLE

CODONS = open('RNA_Codon_Table.txt','rU').read()

def getRNATable ():
    rna_table = []
    rna_table.append ([])
    rna_table.append ([])

    #print rna_table

    for i in range (0, 64):
        idx= i*6
        codon=CODONS[idx:idx+3]
        short_codon = CODONS[idx+4]

        #print "CODON: ", codon, short_codon

        rna_table[0].append (codon)
        if short_codon == "*": short_codon = ""
        rna_table[1].append (short_codon)

    return rna_table

RNA_TABLE = getRNATable()

'''
Input: DNA Sequence assume 5 prime to 3 prime
Output: transcirbed DNA (complement).
Example:
In: GAAACT
OUT: AGTTC
'''
def reverseComplement (dna):

    arr=[['A', 'C', 'G', 'T'],['T', 'G', 'C', 'A']]
    dlen = len(dna)

    res = ""
    #print "DNA LEN: ", dlen
    for i in range (dlen-1, -1, -1):
        #print "III: ", i
        chr = dna[i:i+1]
        idx = arr[0].index (chr)
        comp_chr = arr[1][idx]
        res += comp_chr
    return res

def transcribeRNA (dna):

    rna = ""
    for i in range (0, len (dna)):
        if dna[i] == "T":
            char = "U"
        else:
            char = dna[i]
        rna += char

    return rna

def rnaToDNA (rna):
    dna = ""

    for i in range (0, len (rna)):
        if rna[i] == "U":
            char = "T"
        else:
            char = rna[i]
        dna += char

    return dna

def translateCodons (rna):

    rna_len = len (rna)

    protein = ""
    protein = ""
    for i in range (0, rna_len, 3):
        codon = rna[i:i+3]

        idx = RNA_TABLE[0].index (codon)
        protein += RNA_TABLE[1][idx]

    return protein

def peptideToCodons (peptide):

    #print "PEPTIDE: " + peptide
    codons = []

    try:
        indexes = [index for index, value in enumerate (RNA_TABLE[1]) if value == peptide]
        for i in indexes:
            #print "CODONS1: " + RNA_TABLE[0][i]
            codons.append(RNA_TABLE[0][i])
    except StopIteration:
        print "Index issue"

    return codons

'''
Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
     Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
     Output: All substrings of Text encoding Peptide (if any such substrings exist).
'''
@timing
def encodePeptide (rna, peptides):

    plen = len (peptides)
    first_codons = peptideToCodons(peptides[0])

    result = []
    #print "First COdons:" , first_codons
    rna_len = len (rna) - (plen - 1) * 3
    #print rna_len
    for i in range (0, rna_len):
        #print "AAAAAAA", rna_len
        codon = rna[i:i+3]
        if codon in first_codons:
            # start checking
            tmp_codons = codon
            for j in range (1, plen):
                #print "JJJJ:" , j
                tmp_codon = rna[i+j*3:i+j*3+3]
                codons = peptideToCodons(peptides[j])
                #print "TMP IN CODS", tmp_codon, codons
                if tmp_codon in codons:
                    #print "TMPCOD IN CODONS:", j, tmp_codon, codons
                    tmp_codons += tmp_codon
                    if j == plen - 1 and len (tmp_codons) == plen * 3:
                        #print plen, j,  " IHVE FOUND ONE:" + tmp_codons, len (tmp_codons)
                        result.append (tmp_codons)
                        tmp_codons = ""

                #print i, j, codon, "CODONSAA:" + rna[i+j*3:i+j*3+3]

    # codons = []
    # for i in range (0, len (peptides)):
    #    codons += peptideToCodons(peptides[i])

    #print "CODONS:", codons

    return result

########### MAIN

DNA = open('Bacillus_Brevis_Genome.txt','rU').read().split("\n")
DNA = ''.join (DNA)
peptides = "VKLFPWFNGY"

#print DNA
print translateCodons("CCCCGUACGGAGAUGAAA")
print translateCodons("CCUCGUACAGAAAUCAAC")
print translateCodons("CCCAGUACCGAAAUUAAC")
print translateCodons("CCGAGGACCGAAAUCAAC")

#print translateCodons("CYCLIC")
RNA = transcribeRNA(DNA)
print "TRANS DONE"
rDNA = reverseComplement(DNA)
rRNA = transcribeRNA(rDNA)

print "STart Encoding"
enc_RNA = encodePeptide(RNA ,peptides)
enc_rRNA = encodePeptide(rRNA ,peptides)

print "PRINT",  enc_RNA, enc_rRNA
for rna in enc_RNA:
    print rnaToDNA(rna)

for rrna in enc_rRNA:
    print reverseComplement (rnaToDNA (rrna))

#print "P2C: ", peptideToCodons("S")
print

