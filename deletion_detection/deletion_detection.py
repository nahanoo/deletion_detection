import pysam
import pandas as pd
import os
import subprocess
from io import StringIO

class Deletion():
    def get_coverage(self,path):
        cmd = ['samtools','depth','-a','-J',path]
        process = subprocess.run(cmd,capture_output=True)
        coverage = pd.read_csv(StringIO(process.stdout.decode()),sep='\t')
        coverage.columns=['chromosome','position','coverage']
        self.coverage = dict()
        for chromosome,position,coverage in zip(coverage['chromosome'],\
            coverage['position'],coverage['coverage']):
            self.coverage[(chromosome,position)] = coverage

    def get_deletions(self,path):
        """This functions iterates over all primary reads in the alignment file.
        The CIGAR code of every read is checked for deletions.
        If deletions are found, the deletion count and the coverage at given position is stored.
        """
        self.deletions = dict()
        alignment = pysam.AlignmentFile(path,'rb')
        reads = []
        for read in alignment:
            if not read.is_secondary:
                reads.append(read)
        total_reads = len(reads)
        self.get_coverage(path)
        for counter,read in enumerate(reads):
            #We set the block position to zero which stands for the position in the reference.
            block_pos = 0
            #We set the cigar position to zero which stands for the position within the cigar code.
            cigar_pos = 0
            #Iterating over the cigar code:
            for (block_type, block_len) in read.cigar:
                """Block_type contains information about variant type. Depending on the block type
                we need to move forward on the reference side (block_pos) or the query side (cigar_pos).
                There is a very useful table about this on https://samtools.github.io/hts-specs/SAMv1.pdf page 8"""
                #0 is match or mismatch
                if block_type == 0:
                    block_pos += block_len
                    cigar_pos += block_len
                #1 is insertion to the reference
                elif block_type == 1:
                    cigar_pos += block_len
                #2 deletion from the reference
                elif block_type == 2:
                    #Deletion is stored in dictionary.
                    #Values are ([count,coverage],length)
                    #pysam is 0 index based, to make it compareable with samtools I switched to 1 indexing for storing.
                    if self.deletions.get((read.reference_name,read.pos+block_pos+1)) != None:
                        self.deletions[(read.reference_name,read.pos+block_pos+1)][0][0] += 1
                        coverage = self.coverage[(read.reference_name,read.pos+block_pos+1)]
                        self.deletions[(read.reference_name,read.pos+block_pos+1)][0][1] = coverage
                    else:
                        self.deletions[(read.reference_name,read.pos+block_pos+1)] = ([1,None],block_len)
                    block_pos += block_len
                #3 is skipped region from the reference
                elif block_type == 3:
                    block_pos += block_len
                #4 is soft clipping
                elif block_type == 4:
                    cigar_pos += block_len

    def filter_deletions(self,min_count,min_freq):
        """This function filters deletions. Min_count stands for the minimal
        observed count of deletions and min_freq minimal observed deletion frequency (count/coverage)."""
        self.filtered_deletions = dict()
        for k,((count,coverage),length) in self.deletions.items():
            if (count >= min_count) and (count/coverage >= min_freq):
                self.filtered_deletions[k] = ([count,coverage],length)

    def create_df(self):
        self.output = pd.DataFrame(columns=['chromosome','position','length','count','coverage','type'],index=range(len(self.filtered_deletions.values())))
        for counter,((chromosome,position),(count_coverage,length)) in enumerate(self.filtered_deletions.items()):
            self.output.at[counter,'chromosome'] = chromosome
            self.output.at[counter,'position'] = position
            self.output.at[counter,'length'] = length
            self.output.at[counter,'count'] = count_coverage[0]
            self.output.at[counter,'coverage'] = count_coverage[1]
            self.output.at[counter,'type'] = 'in-read deletion'

    def write_deletions(self,df,out,prefix):
        """This function writes detected deletions to tsv."""
        df.to_csv(os.path.join(out,prefix+'.reads.tsv'),index=False,sep='\t')
