import pathlib
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction.Restriction import RestrictionBatch
from itertools import combinations
import math

import matplotlib.pyplot as plt
import numpy as np

project_folder = pathlib.Path().cwd()
data_folder = project_folder.joinpath('data')
file = [file for file in data_folder.glob('*.fna')][0]

def _soft_round(x, base=100): return int(math.floor(x/base)) * base

def main():
    
    total_fragments: list = []
    for chrom in SeqIO.parse(file, 'fasta'):
        resbat = RestrictionBatch(['EcoRI', 'BamHI', 'ApaI'])
        result = resbat.search(chrom.seq)
        slice_positions_list: list = sorted([slice_pos for _, enzyme_slice_list in result.items() for slice_pos in enzyme_slice_list])
        
        fragments_list: list = []
        chrom_len = len(chrom.seq)
        for i_pos, pos in enumerate(slice_positions_list):
            start_pos = slice_positions_list[i_pos] if i_pos > 0 else 0
            end_pos = slice_positions_list[i_pos + 1] if i_pos < len(slice_positions_list) - 1 else chrom_len
            fragments_list.append(end_pos - start_pos)
        total_fragments += fragments_list
    

    bins = np.arange(0, 2200, 100)

    fig, ax = plt.subplots(figsize=(9, 5))
    _, bins, patches = plt.hist([np.clip(total_fragments, bins[0], bins[-1])], bins=bins,)

    xlabels = bins[1:].astype(str)
    xlabels[-1] += '+'

    N_labels = len(xlabels)
    plt.xlim([0, 2000])
    plt.xticks(100 * np.arange(N_labels) + 50)
    ax.set_xticklabels(xlabels)
    fig.savefig(f'temp_{chrom}.png', dpi=fig.dpi)
if __name__=='__main__':
    main()
    #print(project_folder)
