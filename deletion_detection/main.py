import argparse
from deletion_detection import Deletion
from os.path import join


def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect insertions in evolved bacterial strains. Fast and simpple.')
    parser.add_argument(
        'anceteral', help='fasta file of the ancesteral strain'
    )
    parser.add_argument(
        'mutant', help='genbank file of the mutated strain.')
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument('--plot', help='if this flag is added the alignment of every insertion\
        is plotted.', action='store_true')

    return parser.parse_args()


args = parse_args()
d = Deletion(args)
d.chunk_assembly()
d.mapper(d.reference_fasta, join(d.out_dir, "chunked_sequences.fasta"),
         join(d.out_dir, "aligned.sam"))
d.get_deletions()
d.get_mut_origin()
d.dump()
if args.plot:
    d.plot_deletions()
    d.plot_annotation()
# i.clean()
