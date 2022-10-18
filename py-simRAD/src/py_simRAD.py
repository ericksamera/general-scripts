#!/usr/bin/env python3
"""
Purpose: Perform double restriction digest on a given genome.
"""
__author__ = "Erick Samera; adapted from Joon Lee"
__version__ = "1.1.0"
__comments__ = "stable enough; more biologically accurate?"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    ArgumentDefaultsHelpFormatter)
from pathlib import Path
# --------------------------------------------------
from itertools import combinations
from math import ceil, floor
from Bio import SeqIO
from Bio.Restriction.Restriction import RestrictionBatch
import numpy as np
import pandas as pd
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        #usage='%(prog)s',
        description="Perform double restriction digest on a given genome.",
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    parser.add_argument(
        '-o',
        '--out',
        dest='output_path',
        metavar='DIR',
        type=Path,
        required=True,
        help="path of output file (.csv)")
    
    group_binning_parser = parser.add_argument_group('fragment size binning options')
    group_binning_parser.add_argument(
        '-m',
        '--min',
        dest='bin_min',
        metavar='INT',
        type=int,
        default=100,
        help="min value for bins")
    group_binning_parser.add_argument(
        '-s',
        '--step',
        dest='bin_step',
        metavar='INT',
        type=int,
        default=100,
        help="step for bins")
    group_binning_parser.add_argument(
        '-M',
        '--max',
        dest='bin_max',
        metavar='INT',
        type=int,
        default=2100,
        help="max value for bins")

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if args.input_path: args.input_path = args.input_path.resolve()
    if args.output_path: args.output_path = args.output_path.resolve()

    return args
# --------------------------------------------------
def _ceil_round(x, base=100)-> int: return ceil(int(x)/base)*base
def _generate_fragment_counts_dict(rst_enz_combos_arg: list, bins_arg: range) -> dict:
    """
    Generates a fragment counts dict:
    
    Parameters:
        rst_enz_combos_arg: list
            the generated list of restriction fragment combinations
        bins_arg: range
            a range of bins to be generated
    
    Returns:
        (dict)
            dictionary of each combination and its count of fragments in each bin
    """
    fragment_combination_counts = {str('-'.join(list(combination))): {} for combination in rst_enz_combos_arg}
    for combination in rst_enz_combos_arg:
        combination_str = str('-'.join(list(combination)))
        enzyme_site_str = [f"{str(enzyme)}: {enzyme.elucidate()}" for enzyme in RestrictionBatch(list(combination))]
        fragment_combination_counts[combination_str]['restriction_sites'] = ', '.join(enzyme_site_str)
        fragment_combination_counts[combination_str]['total'] = 0
        fragment_combination_counts[combination_str][f'< {bins_arg[0]}'] = 0
        fragment_combination_counts[combination_str].update({bin_num: 0 for bin_num in bins_arg})
        fragment_combination_counts[combination_str][f'> {bins_arg[-1]}'] = 0
    return fragment_combination_counts
def _generate_restriction_fragments(slice_positions_arg: list, end_pos_arg: int) -> list:
    """
    From a list of cutting positions, do the cutting and generate a list of fragment sizes.

    Parameters:
        slice_positions_arg: list
            list of slice positions from the restriction enzymes
        end_pos_arg: int
            the last position in a given chromosome
    
    Returns:
        (list)
            list of restriction restriction fragments of given lengths
    """
    
    restriction_fragments: list = []

    for i_pos, _ in enumerate(slice_positions_arg):
        start_pos = slice_positions_arg[i_pos] if i_pos > 0 else 0
        end_pos = slice_positions_arg[i_pos + 1] if i_pos < len(slice_positions_arg) - 1 else end_pos_arg
        restriction_fragments.append(end_pos - start_pos)
    return restriction_fragments
def _generate_restriction_fragments2(seq_arg, restriction_enzymes_arg: RestrictionBatch) -> list:
    """
    From a list of cutting positions, do the cutting and generate a list of fragment sizes, \
    but this way is more biologically accurate, technically. It cuts with the first enzyme, \
    then cuts again with the second enzyme (if applicable).

    Parameters:
        seq_arg: Seq
            chromosomal sequence to cut
        restriction_enzymes_arg: RestrictionBatch
            a set of enzymes to cut it with
    
    Returns:
        (list)
            list of restriction restriction fragments of given lengths
    """
    enzymes = [i for i in restriction_enzymes_arg]
    first_restriction = enzymes[0]
    second_restriction = enzymes[1] if len(enzymes)>1 else None
    first_pass = first_restriction.catalyse(seq_arg)
    if second_restriction:
        second_pass: list = []
        for fragment in first_pass:
            second_pass += second_restriction.catalyse(fragment)
    else:
        second_pass = first_pass
    return second_pass
def _subtract_8_fix(rst_enz_combos_arg: list, fragment_combination_counts_arg: dict) -> dict:
    """
    In at least the total counts, the R implementation of simRAD seems to be off by 8, for some reason. \
    This function deals with this.

    Parameters:
        rst_enz_combos_arg: list
            the restriction fragment combinations list
        fragment_combination_counts_arg: dict
            the fragments by combination dictionary
        
    Returns:
        (dict)
            a fixed dictionary maybe with the 8-off error dealt with
    """    
    for combination in rst_enz_combos_arg:
        combination_str = '-'.join(list(combination))
        for key in fragment_combination_counts_arg[combination_str]:
            if all([
                    isinstance(fragment_combination_counts_arg[combination_str][key], int),
                    fragment_combination_counts_arg[combination_str][key] > 8
                    ]):
                fragment_combination_counts_arg[combination_str][key] -= 8
    return fragment_combination_counts_arg
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    restriction_enzymes_list: list = ['SbfI', 'EcoRI', 'SphI', 'PstI', 'MspI', 'MseI']

    args = get_args()

    # generate bins for csv output
    bins: range = range(args.bin_min, args.bin_max, args.bin_step)

    # list of restriction enzyme pairs, including individual restrictions
    rst_enz_combinations: list = \
        sorted([enzyme_pair for enzyme_pair in combinations(restriction_enzymes_list, 2)] +\
        [(enzyme, ) for enzyme in restriction_enzymes_list])

    # dictionary of fragment counts, iteratively generate bins for counts, and include oversize bin 
    fragment_combination_counts = _generate_fragment_counts_dict(rst_enz_combinations, bins)
    
    for combination in rst_enz_combinations:
        combination_str = '-'.join(list(combination))
        
        for chr in SeqIO.parse(args.input_path, 'fasta'):
            
            # switch which implementation
            biologically_accurate = True

            # ignore the mitochrondrial genome, only do nuclear
            if 'mitochondrion' in chr.description: continue

            if not biologically_accurate:            
                # get all unique slice positions
                restriction_batch = RestrictionBatch(list(combination))
                restriction_result = restriction_batch.search(chr.seq.upper())
                # using (end - beginning), add 0 position so that first fragment is the entire length from start to that position
                slice_positions: list = sorted(set([0] + [slice_pos for _, enzyme_slice_list in restriction_result.items() for slice_pos in enzyme_slice_list]))

                # generate "fragments" by cutting between slice positions
                restriction_fragments = _generate_restriction_fragments(slice_positions, len(chr.seq))

                # iteratively add fragments to bins and total
                for fragment_len in restriction_fragments:
                    bin_key = _ceil_round(fragment_len, args.bin_step)
                    bin_keys = [int(key) for key in fragment_combination_counts[combination_str] if str(key).isnumeric()]
                    
                    # if the bin_key is 0, this is probably a mistake
                    if bin_key == 0: continue
                    
                    # count fragments in the correct bin
                    if bin_key in bin_keys:
                        fragment_combination_counts[combination_str][bin_key] += 1
                    elif bin_key < min(bin_keys):
                        fragment_combination_counts[combination_str][f'< {bins[0]}'] += 1
                    elif bin_key > min(bin_keys):
                        fragment_combination_counts[combination_str][f'> {bins[-1]}'] += 1
                    
                    # add to the total
                    fragment_combination_counts[combination_str]['total'] += 1

            elif biologically_accurate:
                restriction_enzymes = RestrictionBatch(list(combination))
                restriction_fragments = _generate_restriction_fragments2(chr.seq, restriction_enzymes)

                for fragment_len in restriction_fragments:
                    bin_key = _ceil_round(fragment_len, args.bin_step)
                    bin_keys = [int(key) for key in fragment_combination_counts[combination_str] if str(key).isnumeric()]
                    
                    # if the bin_key is 0, this is probably a mistake
                    if bin_key == 0: continue
                    
                    # count fragments in the correct bin
                    if bin_key in bin_keys:
                        fragment_combination_counts[combination_str][bin_key] += 1
                    elif bin_key < min(bin_keys):
                        fragment_combination_counts[combination_str][f'< {bins[0]}'] += 1
                    elif bin_key > min(bin_keys):
                        fragment_combination_counts[combination_str][f'> {bins[-1]}'] += 1

    pd.DataFrame.from_dict(fragment_combination_counts, orient='index').to_csv(args.output_path)
# # --------------------------------------------------
if __name__ == '__main__':
    main()
