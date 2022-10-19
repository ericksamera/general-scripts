#!/usr/bin/env python3
"""
Purpose: Perform double restriction digest on a given genome.
"""
__author__ = "Erick Samera; Michael Ke; ACK: Joon Lee"
__version__ = "2.0.0"
__comments__ = "more biologically accurate"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
from itertools import combinations
from math import ceil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import Restriction
import pandas as pd
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        #usage='%(prog)s',
        description="Perform double restriction digest on a given genome.",
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        'input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    parser.add_argument(
        '-o',
        '--out',
        dest='output_path',
        metavar='FILE',
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
        help="min value for bins (inclusive), will contain 0 to first bin size (default=100)")
    group_binning_parser.add_argument(
        '-s',
        '--step',
        dest='bin_step',
        metavar='INT',
        type=int,
        default=100,
        help="step for bins (default=100)")
    group_binning_parser.add_argument(
        '-M',
        '--max',
        dest='bin_max',
        metavar='INT',
        type=int,
        default=2000,
        help="max value for bins (inclusive) (default=2000)")

    group_rst_enz_parser = parser.add_argument_group('restriction enzyme options')
    group_rst_enz_parser.add_argument(
        '--enzymes',
        dest='enzymes',
        metavar='STR',
        type=str,
        default="SbfI;EcoRI;SphI;PstI;MspI;MseI",
        help='enzymes to use, ex: "EcoRI;BamHI" (default: "SbfI;EcoRI;SphI;PstI;MspI;MseI")')
    group_rst_enz_parser.add_argument(
        '-b',
        '--buffer',
        dest='buffer',
        metavar='INT',
        type=int,
        default=None,
        help="if given, min length on fragment ends to find restriction site (default: None)")
    group_rst_enz_parser.add_argument(
        '--fast',
        dest='use_fast',
        action='store_true',
        help="use the fast algorithm (but maybe inaccurate)")
    group_rst_enz_parser.add_argument(
        '--output_fastas',
        dest='output_fasta_dir',
        metavar='DIR',
        type=Path,
        required=False,
        help="if given, dir of output for fasta fragments (SLOW!) (default: None)")

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if args.input_path: args.input_path = args.input_path.resolve()
    if args.output_path: args.output_path = args.output_path.resolve()
    # --------------------------------------------------
    invalid_enzymes = [enzyme for enzyme in args.enzymes.split(';') if enzyme not in Restriction.AllEnzymes.elements()]
    if invalid_enzymes: 
        parser.error(f"Couldn't process the following enzymes: {' '.join(invalid_enzymes)}")

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
        enzyme_site_str = [f"{str(enzyme)}: {enzyme.elucidate()}" for enzyme in Restriction.RestrictionBatch(list(combination))]
        fragment_combination_counts[combination_str]['restriction_sites'] = ', '.join(enzyme_site_str)
        fragment_combination_counts[combination_str]['total'] = 0
        fragment_combination_counts[combination_str].update({bin_num: 0 for bin_num in bins_arg})
        fragment_combination_counts[combination_str][f'> {bins_arg[-1]}'] = 0
    return fragment_combination_counts
def _generate_restriction_fragments_fast(slice_positions_arg: list, end_pos_arg: int) -> list:
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
def _generate_restriction_fragments(seq_arg: SeqRecord, restriction_enzymes_arg: Restriction.RestrictionBatch, buffer_arg: int, fasta_output_path_arg: Path) -> list:
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
    enzymes_str = [str(i) for i in restriction_enzymes_arg]
    first_restriction = enzymes[0]
    second_restriction = enzymes[1] if len(enzymes)>1 else None

    # generate restriction fragments using the first enzyme,
    first_pass = first_restriction.catalyse(seq_arg) if not buffer_arg else [fragment for fragment in first_restriction.catalyse(seq_arg) if len(fragment) > buffer_arg]


    if second_restriction:

        # generate a list of restriction fragments
        second_pass: list = []
        
        for fragment in first_pass:

            # if the buffer arg is passed and the length is greater than the buffer, if not greater dont bother catalyzing
            filtered_fragment: SeqRecord = fragment
            if buffer_arg:
                if len(fragment) > buffer_arg: filtered_fragment = fragment[buffer_arg:-buffer_arg]
                else: second_restriction_results = [fragment]

            second_restriction_results = second_restriction.catalyse(filtered_fragment)
            second_pass += second_restriction_results

            if fasta_output_path_arg:
                with open(fasta_output_path_arg.joinpath(f'{"-".join(enzymes_str)}.fasta'), 'a') as fasta_output:
                    for test_fragment in second_restriction_results:
                        enzyme_sides: dict = {
                            '5': str(first_restriction) if test_fragment[:10] == fragment[:10] else str(second_restriction),
                            '3': str(first_restriction) if test_fragment[-10:] == fragment[-10:] else str(second_restriction),
                            'seq': test_fragment
                        }
                        SeqIO.write(SeqRecord(enzyme_sides['seq'], id=f"{enzyme_sides['5']}-{enzyme_sides['3']}", description=''), fasta_output, 'fasta')
    else:
        second_pass = first_restriction.catalyse(seq_arg)
    return second_pass
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    restriction_enzymes_list: list = ['SbfI', 'EcoRI', 'SphI', 'PstI', 'MspI', 'MseI']

    args = get_args()

    # generate bins for csv output
    bins: range = range(args.bin_min, args.bin_max + 1, args.bin_step)

    # list of restriction enzyme pairs, including individual restrictions
    rst_enz_combinations: list = \
        sorted([enzyme_pair for enzyme_pair in combinations(restriction_enzymes_list, 2)] +\
        [(enzyme, ) for enzyme in restriction_enzymes_list])

    # dictionary of fragment counts, iteratively generate bins for counts, and include oversize bin 
    fragment_combination_counts = _generate_fragment_counts_dict(rst_enz_combinations, bins)

    for combination in rst_enz_combinations:
        combination_str = '-'.join(list(combination))

        restriction_fragments_per_combination: list = []
        for chr in SeqIO.parse(args.input_path, 'fasta'):
            # ignore the mitochrondrial genome, only do nuclear
            if 'mitochondrion' in chr.description: continue

            if args.use_fast:
                # get all unique slice positions
                restriction_batch = Restriction.RestrictionBatch(list(combination))
                restriction_result = restriction_batch.search(chr.seq.upper())
                # using (end - beginning), add 0 position so that first fragment is the entire length from start to that position
                # and remove duplicate positions from isoschizomers or similar cut sites
                slice_positions: list = sorted(set([0] + [slice_pos for _, enzyme_slice_list in restriction_result.items() for slice_pos in enzyme_slice_list]))

                # generate "fragments" by cutting between slice positions
                restriction_fragments = _generate_restriction_fragments_fast(
                    slice_positions_arg=slice_positions, 
                    end_pos_arg=len(chr.seq))
                # iteratively add fragments to bins and total
                for fragment_len in restriction_fragments:
                    bin_key = _ceil_round(fragment_len, args.bin_step)
                    bin_keys = [int(key) for key in fragment_combination_counts[combination_str] if str(key).isnumeric()]
                    
                    # if the bin_key is 0, this is probably a mistake
                    if bin_key == 0: continue
                    
                    # count fragments in the correct bin
                    if bin_key in bin_keys:
                        fragment_combination_counts[combination_str][bin_key] += 1
                    elif bin_key > min(bin_keys):
                        fragment_combination_counts[combination_str][f'> {bins[-1]}'] += 1
                    
                    # add to the total
                    fragment_combination_counts[combination_str]['total'] += 1
                restriction_fragments_per_combination += restriction_fragments

            elif not args.use_fast:
                restriction_enzymes = Restriction.RestrictionBatch(list(combination))
                restriction_fragments = _generate_restriction_fragments(
                    seq_arg=chr.seq, 
                    restriction_enzymes_arg=restriction_enzymes, 
                    buffer_arg=args.buffer, 
                    fasta_output_path_arg=args.output_fasta_dir)

                for fragment in restriction_fragments:
                    bin_key = _ceil_round(len(fragment), args.bin_step)
                    bin_keys = [int(key) for key in fragment_combination_counts[combination_str] if str(key).isnumeric()]
                    
                    # if the bin_key is 0, this is probably a mistake
                    if bin_key == 0: continue

                    # count fragments in the correct bin
                    if bin_key in bin_keys:
                        fragment_combination_counts[combination_str][bin_key] += 1
                    elif bin_key > min(bin_keys):
                        fragment_combination_counts[combination_str][f'> {bins[-1]}'] += 1

                    # add to the total
                    fragment_combination_counts[combination_str]['total'] += 1
                restriction_fragments_per_combination += [len(fragment) for fragment in restriction_fragments]

    # generate DataFrame output
    pd.DataFrame.from_dict(fragment_combination_counts, orient='index').to_csv(args.output_path)
# # --------------------------------------------------
if __name__ == '__main__':
    main()
