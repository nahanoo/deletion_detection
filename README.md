# Deletion detection for haploid genomes

## Introduction

This is a very simple deletion detector for haploid genomes. Due to it's lightweight design it's pretty fast. Using 16 CPUs it takes about 2 minutes to run on a BAM file with the size of 1GB. 
As an input it uses sorted BAM files. Using `pysam` this tool iterates over every read in the Alignment file. The CIGAR code of every read is checked for deletions. Every position and how many times this deletion was found is stored as well as the coverage at given position.  
All deletions are then filtered for `--min_counts` and `--min_frequency`. If no values are specified using those flags default min_counts is 5 and default min_frequency is 0.8.  
If specified with `--output_no_alignment_regions` all regions with no read alignments are outputted additionally to in-read deletions.
This package was tested with PacBio HiFi data using [minimap2](https://github.com/lh3/minimap2) for alignment. It works best by using corrected reads. For detecting short deletions this package should also work well with Illumina data. Uncorrected nanopore reads could be a bit more tricky because of the indel issues that come with nanopore data. It should still work pretty well, however the outputted positions will be affected by the flags `--min_counts` and `--min_frequency`.

## Installation

This package requires SAMtools>=1.11 in your PATH.  
This package can be installed using pip. Since it is still in testing it is not available on PyPi yet.
```
git clone git@github.com:nahanoo/deletion_detection.git
cd deletion_detection
pip install .
```

## Usage

Help page called with `deletion_detection -h`:
```
usage: detect_deletions [-h] [--output_no_alignment_regions] [--min_counts MIN_COUNTS]
                        [--min_frequency MIN_FREQUENCY] [--min_mapping_quality MIN_MAPPING_QUALITY]
                        bam_file output_dir

Deletion detection based on sorted BAM files. By default it only detects deletions which are located in
reads. Adding the --output_no_alignment_regions will also output the regions in the reference where no reads
aligned. This feature requires SAMtools>=1.11 in your PATH.

positional arguments:
  bam_file              path to sorted BAM file.
  output_dir            output direcotry.

optional arguments:
  -h, --help            show this help message and exit
  --output_no_alignment_regions
                        outputs the regions in the reference where no reads aligned. This feature requires
                        SAMtools>=1.11 in your PATH.
  --min_counts MIN_COUNTS
                        minimal observed reads with the deletion [default is 5]
  --min_frequency MIN_FREQUENCY
                        minimal observed frequency of the deletion (observed reads with the deletion divided
                        by the coverage at this position) [default is 0.8]
  --min_mapping_quality MIN_MAPPING_QUALITY
                        minimal sequence mapping quality considered for coverage calculations [default is 60]
```
