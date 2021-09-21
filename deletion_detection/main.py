import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Deletion detection based on sorted BAM files.\
        By default it only detects deletions which are located in reads.\
        Adding the --output_no_alignment_regions will also output the regions\
        in the reference where no reads aligned. This feature requires SAMtools>=1.11 in your PATH.\
        ')
    parser.add_argument('bam_file',help='path to sorted BAM file.')
    parser.add_argument('output_dir',help='ouput direcotry.')
    parser.add_argument('--output_no_alignment_regions',required=False,action='store_true',help='outputs the regions\
        in the reference where no reads aligned. This feature requires SAMtools>=1.11 in your PATH.')
    return parser.parse_args()

def main():
    args = parse_args()

args = parse_args()