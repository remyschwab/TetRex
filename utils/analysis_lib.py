import string

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import os

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


def compare_to_prosite(tetrex_output):
    prosite_path = '/Users/rschwab/repos/TetRex/data/ps_scan/ps_scan.pl'
    dat_path = '/Users/rschwab/repos/TetRex/data/REFERENCES/prosite.dat'
    true_positives = []
    tetrex_results = []
    for res in os.listdir(tetrex_output):
        handle_path = tetrex_output+res
        tetrex_count = 0
        prosite_count = 0
        pattern = res.split('.')[0]
        for line in open(handle_path):
            tetrex_count += 1
            output = subprocess.run([prosite_path, '-p', pattern, '-d', dat_path, line], capture_output=True, text=True)
            if output.stdout != '':
                prosite_count += 1
        tetrex_results.append(tetrex_count)
        true_positives.append(prosite_count)
    return (tetrex_results, true_positives)

