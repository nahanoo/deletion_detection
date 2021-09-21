# Deletion detection for haploid genomes

## Introduction

This is a very simple deletion detector for haploid genomes.  
As an input it uses sorted BAM files. Using `pysam` this tool iterates over every read in the Alignment file. The CIGAR code of every read is checked for deletions. Every position and how many times this deletion was found is stored. Additionally the coverage at given position is stored.  
All deletions are then filtered for `--min_counts` and `--min_frequency`. If no values are specified using those flags default min_count is 5 and default min_frequency is 0.8.  
If specified with `--output_no_alignment_regions` all regions with no read alignments are outputted additionally to in-read deletions.
Note that this feature requires SAMtools>=1.11 in your PATH. 
In-read deletion detection for haploid genomes.

## Installation

This package can be used using pip. Since it is still in testing it is not available on PyPi yet.
```
git clone git@github.com:nahanoo/deletion_detection.git
cd deletion_detection
pip install .
```

## Usage

The tool can then be called using  
`deletion_detection BAMFILE OUTPUTDIR`

Help page called with `deletion_detection -h`:
```
usage: detect_deletions [-h] [--output_no_alignment_regions] [--min_counts MIN_COUNTS] [--min_frequency MIN_FREQUENCY]
                        bam_file output_dir

Deletion detection based on sorted BAM files. By default it only detects deletions which are located in reads. Adding the
--output_no_alignment_regions will also output the regions in the reference where no reads aligned. This feature requires
SAMtools>=1.11 in your PATH.

positional arguments:
  bam_file              path to sorted BAM file.
  output_dir            ouput direcotry.

optional arguments:
  -h, --help            show this help message and exit
  --output_no_alignment_regions
                        outputs the regions in the reference where no reads aligned. This feature requires SAMtools>=1.11 in your PATH.
  --min_counts MIN_COUNTS
                        minimal observed reads with the deletion [default is 5]
  --min_frequency MIN_FREQUENCY
                        minimal observed frequency of deletion (observed reads with deletion divided by coverage) [default is 5]
```
