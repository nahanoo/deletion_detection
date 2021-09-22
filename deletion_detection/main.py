import argparse
from deletion_detection import Deletion
from deletion_detection import NoAlignment

def parse_args():
    """Parsing required and optional arguemtns."""
    parser = argparse.ArgumentParser(description='Deletion detection based on sorted BAM files.\
        By default it only detects deletions which are located in reads.\
        Adding the --output_no_alignment_regions will also output the regions\
        in the reference where no reads aligned. This feature requires SAMtools>=1.11 in your PATH.\
        ')
    parser.add_argument('bam_file',help='path to sorted BAM file.')
    parser.add_argument('output_dir',help='output direcotry.')
    parser.add_argument('--output_no_alignment_regions',required=False,action='store_true',help='outputs the regions\
        in the reference where no reads aligned. This feature requires SAMtools>=1.11 in your PATH.')
    parser.add_argument('--min_counts',type=int,help='minimal observed reads with the deletion [default is 5]')
    parser.add_argument('--min_frequency',type=float,help='minimal observed frequency of the deletion (observed reads with\
        the deletion divided by the coverage at this position) [default is 0.8]')
    return parser.parse_args()

def main():
    args = parse_args()
    #Setting default values for min_counts and min_frequency if not specified
    if args.min_counts is None:
        args.min_counts = 5
    if args.min_frequency is None:
        args.min_frequency = 0.8
    #Output in-read deletions
    if not args.output_no_alignment_regions:
        d = Deletion()
        #Get all in-read deletions
        d.get_deletions(args.bam_file)
        #Apply filter
        d.filter_deletions(args.min_counts,args.min_frequency)
        #Write to tsv
        d.write_deletions(args.output_dir)
    #Additionally to in-read deletions output also no alignment regions (if requested with flag)
    if args.output_no_alignment_regions:
        d = Deletion()
        #Get all in-read deletions
        d.get_deletions(args.bam_file)
        #Apply filter
        d.filter_deletions(args.min_counts,args.min_frequency)
        #Write in-read deletions to tsv
        d.write_deletions(args.output_dir)
        n = NoAlignment()
        #Get entire coverage
        n.get_coverage(args.bam_file)
        n.get_no_alignments()
        #Get regions with 0 coverage.
        #Samtools depth -J is used so in read deletions are counted as coverage.
        n.get_no_alignments_length()
        #Write to tsv
        n.write_no_alignments(args.output_dir)