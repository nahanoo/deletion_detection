import argparse
from .deletion_detection import Deletion
from os.path import join


def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect deletions in evolved bacterial strains.')
    parser.add_argument(
        'mutant', help='fasta file of the mutated strain.')
    parser.add_argument(
        'ancestral', help='genbank file of the ancestor'
    )
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument('--plot', help='plot alignments and annotations', action='store_true')

    return parser.parse_args()

def main():
    args = parse_args()
    # Class innit
    d = Deletion(args)
    # Chunk assembly into smaller sequences using a sliding window
    d.chunk_assembly()
    # Map chunked mutant assembly to reference
    d.map_chunks()
    # Gets all areas with no coverage
    d.get_deletions()
    # Gets deletions location in mutant
    d.get_location()
    # Annotates location
    d.annotate()
    # Dumps to file
    d.dump()
    if args.plot:
        # Plot alignment of no coverage region in reference
        d.plot_deletions()
        # Plots annotations
        d.plot_annotation()
    # Cleans temporary files
    d.clean()