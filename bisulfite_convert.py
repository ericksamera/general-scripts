#!/usr/bin/env python3
"""
Purpose: Performs in-silico bisulfite conversion on a given sequence file.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "stable enough"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    ArgumentDefaultsHelpFormatter)
from pathlib import Path
# --------------------------------------------------
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        #usage='%(prog)s',
        description="Performs in-silico bisulfite conversion on a given sequence file.",
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'input_path',
        type=Path,
        help="path of input file (.fasta/.gb) or dir of input files")
    parser.add_argument(
        '-o',
        '--out',
        dest='output_path',
        metavar='DIR',
        type=Path,
        help="path of output dir")
    parser.add_argument(
        '--output-type',
        dest='output_type',
        metavar='TYPE',
        type=str,
        choices=['gb', 'fasta'],
        default='gb',
        help="output type")
    # --------------------------------------------------

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    args.input_path = args.input_path.resolve()
    args.output_path = args.output_path.resolve()

    return args
# --------------------------------------------------
def _bs_convert_seq(sequence_arg: str) -> str:
    """
    Function returns a bisulfite converted sequence.

    Parameters:
        sequence_arg: str
            string of sequence
    
    Returns:
        str
            bisulfite-converted sequence
    """
    for nuc in 'CATCAT':
        sequence_arg = sequence_arg.replace(f'C{nuc}', f'T{nuc}')
    return sequence_arg
def _generate_CpG_feature(sequence_arg: str) -> SeqFeature.SeqFeature:
    """
    Function generates a SeqFeature object containing CpG positions.

    Parameters:
        sequence_arg: str
            string of sequence, doesn't matter if bisulfite converted.
    
    Returns:
        SeqFeature object
    """
    CpG_locations: list = []
    for nuc_i, _ in enumerate(sequence_arg):
            if nuc_i < len(sequence_arg)-1:
                if (sequence_arg[nuc_i] == 'C') and (sequence_arg[nuc_i+1] == 'G'):
                    CpG_locations.append(SeqFeature.FeatureLocation(nuc_i, nuc_i+2))
    
    CpG_SeqFeature = SeqFeature.SeqFeature(
        SeqFeature.CompoundLocation(CpG_locations),
        strand=1,
        type='modified_base',
        qualifiers={'mod_base':'5mc'}
        )

    return CpG_SeqFeature
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()
    
    inferred_filetypes = {
        '.fasta': 'fasta',
        '.fa': 'fasta',
        '.gb': 'gb'
    }

    output_dir = args.output_path.parent if args.output_path.is_file() else args.output_path
    output_dir.mkdir(parents=True, exist_ok=True)

    files_to_write: list = []

    if args.input_path.is_dir():
        for file in args.input_path.glob(r'*[.gb|.fasta|.fa]'):
            SeqObject = SeqIO.read(file, inferred_filetypes[file.suffix])
            file_to_write = SeqRecord(_bs_convert_seq(SeqObject.seq), id=SeqObject.id, description='', annotations={"molecule_type": "DNA"})

            if args.output_type == 'gb':
                generated_CpG_features = _generate_CpG_feature(SeqObject.seq)
                if generated_CpG_features:
                    file_to_write.features.append(generated_CpG_features) 
                   
            files_to_write.append((f'bs_{file.stem}.{args.output_type}', file_to_write))
            
    
    elif args.input_path.is_file():
        file = args.input_path
        SeqObject = SeqIO.read(file, inferred_filetypes[file.suffix])
        file_to_write = SeqRecord(_bs_convert_seq(SeqObject.seq), id=SeqObject.id, description='', annotations={"molecule_type": "DNA"})
        
        if args.output_type == 'gb':
            generated_CpG_features = _generate_CpG_feature(SeqObject.seq)
            if generated_CpG_features:
                file_to_write.features.append(generated_CpG_features)

        files_to_write.append((f'bs_{file.stem}.{args.output_type}', file_to_write))
    
    for filename, SeqObject in files_to_write:
        SeqIO.write(SeqObject, output_dir.joinpath(filename), f'{args.output_type}')
# --------------------------------------------------
if __name__ == '__main__':
    main()
