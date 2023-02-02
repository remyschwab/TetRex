#! /usr/bin/python3
"""
A small python program for reverse translating an amino sequence into nuclelic acid
codon regular expressions
- Remy Schwab, November 18th 2022
"""

import sys

fourWobble = "AC|GT||."

codon_postfix_lut = {
    "M": "AT.G.",       # Methionine
    "W": "TG.G.",       # Tryptophan
    "K": "AA.AG|.",     # Lysine
    "D": "GA.TC|.",     # Aspartic Acid
    "E": "GA.AG|.",     # Glutamic Acid
    "H": "CA.TC|.",     # Histadine
    "N": "AA.CT|.",     # Asparagine
    "Q": "CA.GA|.",     # Glutamine
    "Y": "TA.TC|.",     # Tyrosine
    "F": "TT.CT|.",     # Phenylalanine
    "C": "TG.TC|.",     # Cysteine
    "I": "AT.AC|T|.",   # Isoleucine
    "V": "GT.AC|GT||.", # Valine
    "T": "AC.AC|GT||.", # Threonine
    "P": "CC.AC|GT||.", # Proline
    "G": "GG.AC|GT||.", # Glycine
    "A": "GC.AC|GT||.", # Alanine
    "R": "CG.AC|GT||.", # Arginine NOT COMPLETE
    "L": "CT.AC|GT||.", # Leucine NOT COMPLETE
    "S": "TC.AC|GT||."  # Serine NOT COMPLETE
}

aa_codon_lut = {
    "M": "ATG",                   # Methionine
    "W": "TGG",                   # Tryptophan
    "K": "AA(A|G)",               # Lysine
    "D": "GA(T|C)",               # Aspartic Acid
    "E": "GA(A|G)",               # Glutamic Acid
    "H": "CA(T|C)",               # Histadine
    "N": "AA(C|T)",               # Asparagine
    "Q": "CA(G|A)",               # Glutamine
    "Y": "TA(T|C)",               # Tyrosine
    "F": "TT(C|T)",               # Phenylalanine
    "C": "TG(T|C)",               # Cysteine
    "I": "AT(A|C|T|)",            # Isoleucine
    "V": "GT(A|C|G|T)",           # Valine
    "T": "AC(A|C|G|T)",           # Threonine
    "P": "CC(A|C|G|T)",           # Proline
    "G": "GG(A|C|G|T)",           # Glycine
    "A": "GC(A|C|G|T)",           # Alanine
    "R": "(AG(G|A)|CG(A|C|G|T))", # Arginine
    "L": "(CT(A|C|G|T)|TT(A|G))", # Leucine
    "S": "(TC(A|C|G|T)|AG(C|T))"  # Serine
}

codon_count_table = {
    "M": 1, # Methionine
    "W": 1, # Tryptophan
    "K": 2, # Lysine
    "D": 2, # Aspartic Acid
    "E": 2, # Glutamic Acid
    "H": 2, # Histadine
    "N": 2, # Asparagine
    "Q": 2, # Glutamine
    "Y": 2, # Tyrosine
    "F": 2, # Phenylalanine
    "C": 2, # Cysteine
    "I": 4, # Isoleucine
    "V": 4, # Valine
    "T": 4, # Threonine
    "P": 4, # Proline
    "G": 4, # Glycine
    "A": 4, # Alanine
    "R": 6, # Arginine
    "L": 6, # Leucine
    "S": 6  # Serine
}

def computeCombos(aa_seq):
    total = 1
    for aa in aa_seq:
        total = total*codon_count_table[aa]
    return total

def computeRegEx(aa_sequence):
    return "".join([aa_codon_lut[aa] for aa in aa_sequence])
            

def main():
    na_regex = computeRegEx(sys.argv[1])
    print(na_regex)
    print(computeCombos(sys.argv[1]))


if __name__ == "__main__":
    main()
