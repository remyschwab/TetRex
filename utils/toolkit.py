#! /usr/bin/python3
"""
A small python program for reverse translating an amino sequence into nuclelic acid
codon regular expressions in reverse polish postfix notation.
- Remy Schwab, November 18th 2022
"""

import sys

fourWobble = "AC|GT||."

aa_codon_lut = {
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

def computePostfix(aa_sequence):
    postfix = ""
    for i in range(len(aa_sequence)):
        postfix += aa_codon_lut[aa_sequence[i]]
        if i != 0:
            postfix += '.'
    return postfix
            

def main():
    postfix_expression = computePostfix(sys.argv[1])
    print(postfix_expression)

if __name__ == "__main__":
    main()
