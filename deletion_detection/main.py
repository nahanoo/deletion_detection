import argparse
from deletion_detection import Deletion
from deletion_detection import Annotation
from deletion_detection import NoAlignment

def parse_args():
    """Parsing required and optional arguemtns."""
    parser = argparse.ArgumentParser(description='Deletion detection based on sorted BAM files.\
        By default it only detects deletions which are located in reads.\
        Adding the --output_no_alignment_regions will also output the regions\
        in the reference where no reads aligned. This package requires SAMtools>=1.11 in your PATH.\
        ')
    parser.add_argument('bam_file',help='path to sorted BAM file.')
    parser.add_argument('output_dir',help='output direcotry.')
    parser.add_argument('--genbank',help='if specified deletions will be annotated based on genbank files')
    parser.add_argument('--no_alignment_only',action='store_true',\
        help='if specified only no alignment regions will be outputted')
    parser.add_argument('--min_counts',type=int,help='minimal observed \
        reads with the deletion [default is 5]')
    parser.add_argument('--min_frequency',type=float,help='minimal observed \
        frequency of the deletion (observed reads with\
        the deletion divided by the coverage at this position) [default is 0.8]')
    parser.add_argument('--min_mapping_quality',type=int,help='minimal sequence mapping quality\
        considered for coverage calculations [default is 60]')
    return parser.parse_args()

def main():
    args = parse_args()
    #Setting default values for min_counts and min_frequency if not specified
    if args.min_counts is None:
        args.min_counts = 5
    if args.min_frequency is None:
        args.min_frequency = 0.8
    if args.min_mapping_quality is None:
        args.min_mapping_quality = 60
    
    n = NoAlignment()
    #Get entire coverage
    n.get_coverage(args.bam_file,args.min_mapping_quality)
    n.get_no_alignments()
    #Get regions with 0 coverage.
    #Samtools depth -J is used so in read deletions are counted as coverage.
    n.get_no_alignments_length()
    n.create_df()
    if args.genbank:
        a = Annotation()
        a.parse_genbank(args.genbank)
        a.annotate(n.output)
        n.write_no_alignments(a.annotation,args.output_dir)
    else:
        n.write_no_alignments(n.output,args.output_dir)

    if not args.no_alignment_only:
        #Output in-read deletions
        d = Deletion()
        #Get all in-read deletions
        d.get_deletions(args.bam_file)
        #Apply filter
        d.filter_deletions(args.min_counts,args.min_frequency)
        d.create_df()
        #Write to tsv
        if args.genbank:
            a = Annotation()
            a.parse_genbank(args.genbank)
            a.annotate(d.output)
            d.write_deletions(a.annotation,args.output_dir)
        else:
            d.write_deletions(d.output,args.output_dir)