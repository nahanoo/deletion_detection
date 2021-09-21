import pysam
import pandas as pd

class Deletion():
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
        for read in reads:
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
                    if (read.reference_name,read.pos+block_pos+1) in self.deletions.keys():
                        self.deletions[(read.reference_name,read.pos+block_pos+1)][0][0] += 1
                    else:
                        coverage = alignment.count(read.reference_name,read.pos+block_pos,read.pos+block_pos+1)
                        self.deletions[(read.reference_name,read.pos+block_pos+1)] = ([1,coverage],block_len)
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

    def write_deletions(self,out):
        """This function writes detected deletions to tsv."""
        df = pd.DataFrame(columns=['chromosome','position','length','count','coverage','type'],index=range(len(self.filtered_deletions.values())))
        for counter,((chromosome,position),(count_coverage,length)) in enumerate(self.filtered_deletions.items()):
            df.at[counter,'chromosome'] = chromosome
            df.at[counter,'position'] = position
            df.at[counter,'length'] = length
            df.at[counter,'count'] = count_coverage[0]
            df.at[counter,'coverage'] = count_coverage[1]
            df.at[counter,'type'] = 'In-read deletion'
        df.to_csv(out,index=False,sep='\t')