#! /opt/homebrew/bin/python3.10

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
    "I": "AT(A|C|T)",            # Isoleucine
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


AA_SET = set(codon_count_table.keys())
AA_UNION = "(" + "|".join([_ for _ in AA_SET]) + ")"


def computeCombos(aa_seq):
    total = 1
    for aa in aa_seq:
        total = total*codon_count_table[aa]
    return total


def computeRegEx(aa_sequence):
    return "".join([aa_codon_lut[aa] for aa in aa_sequence])


def make_union_iter(sub_disjunction):
    return "(" + "|".join([_ for _ in sub_disjunction if _ not in '[]']) + ")"


def convert_prosite_pattern(pattern):
    """
    Scanning prosite patterns one character at a time is limited
    Break the pattern into tokens at the pattern delimiter '-'
    """
    posix_pattern = []
    tkn_lst = pattern.split("-")
    for tkn in tkn_lst:
        match tkn:
            case 'x': ## Wildcard
                posix_pattern.append(".")
            case 'x(0,1)': ## Optional Wildcard
                posix_pattern.append(".?")
            case dis if dis.startswith('['): ## Disjunction
                posix_pattern.append(make_union_iter(dis))
            case dis if dis.startswith('{'): ## Negative Disjunction
                posix_pattern.append(make_union_iter(AA_SET-set(tkn)))
            case _: ## Just a plain old AA
                posix_pattern.append(tkn)
    output_query = "".join(posix_pattern)
    output_query = output_query.replace(">","$")
    output_query = output_query.replace("<", "^")
    return output_query


def main():
    PROSITE_PATTERN = sys.argv[1]
    print(convert_prosite_pattern(PROSITE_PATTERN))


if __name__ == "__main__":
    main()
