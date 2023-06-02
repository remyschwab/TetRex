import string

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.pyplot import figure
from pysam import FastaFile


def init_amino_acids():
    amino_acids = [alpha for alpha in string.ascii_uppercase]
    amino_acids.remove('B')
    amino_acids.remove('J')
    amino_acids.remove('O')
    amino_acids.remove('U')
    amino_acids.remove('X')
    amino_acids.remove('Z')
    return amino_acids


def init_tracker(amino_acids):
    tracker = {residue: -1 for residue in amino_acids}
    return tracker


def track(track_map, seq):
    record_dists = []
    for i, residue in enumerate(seq):
        for aa in track_map.keys():
            if track_map[aa] != -1:
                record_dists.append((residue, aa, i-track_map[aa]))
        track_map[residue] = i
    return record_dists


def track_fasta(fa_path):
    fa = FastaFile(fa_path)
    amino_acids = init_amino_acids()
    ref_lengths_map = dict(zip(fa.references, fa.lengths))
    dist_lst = []
    for ref in ref_lengths_map.keys():
        record_tracker = init_tracker(amino_acids)
        record_seq = fa.fetch(ref, 0, ref_lengths_map[ref])
        dists = track(record_tracker, record_seq)
        dist_lst.append(len(dists))
    return dist_lst
