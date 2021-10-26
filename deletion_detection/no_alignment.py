import pysam
from io import StringIO
import pandas as pd
import subprocess
import os

class NoAlignment():
    """This class is used to detect areas with zero coverage.
    This is useful if you want to find plasmids or large 
    genomic regions which were potentially deleted, something
    you would not find by investigating aligned reads.
    Remember that contigs at the start and the end often have
    poor assembly quality which can lead to zero coverage. 
    """

    def get_coverage(self,path,min_mapping_quality):
        """This function calculates coverage using samtools depth.
        Depending on the size of the bam file this can take a while
        and you may want to run it on a cluster. The flag "-J" is used
        which means that deletions in reads count as coverage.
        Therefore coverage is only zero in areas where no reads align.
        This function requires that samtools>=1.11 is in your PATH variable.
        Depending on the sequencing data it's worth experimenting with the
        --min_mapping_quality flag to get the desired results.
        """
        self.bamfile_path = path
        cmd = ['samtools','depth','-aa','-J','-Q',str(min_mapping_quality),self.bamfile_path]
        process = subprocess.run(cmd,capture_output=True)
        self.coverage = pd.read_csv(StringIO(process.stdout.decode()),sep='\t')
        self.coverage.columns=['chromosome','position','coverage']

    def get_no_alignments(self):
        """Hero we get all positions with zero reads aligned.
        """
        self.no_alignments_df = self.coverage[self.coverage['coverage'] == 0]

    def get_no_alignments_length(self):
        """This function detects the length of the unaligned region. 
        It does it by walking up the unaligned region based on a starting position
        and increasing the count of the length if the position is increased by one.
        """
        #For convenience we group the no alignment df as dictionary.
        no_alignment = {chromosome:list() for chromosome in set(self.no_alignments_df['chromosome'])}
        for chromosome,position in zip(self.no_alignments_df['chromosome'],self.no_alignments_df['position']):
            no_alignment[chromosome].append(position)

        #All starting positions are stored in a dictionary.
        #Values are the length of the continuous no alignment region.
        start_positions = dict()
        for chromosome,positions in no_alignment.items():
            #For every chromosome the first base is always the first
            #starting position.
            start_position = positions[0]
            start_positions[(chromosome,start_position)] = 1
            for position in positions[1:]:
                #If the no alignment region is continuous the length counter in increased.
                if start_position + start_positions[(chromosome,start_position)] == position:
                    start_positions[(chromosome,start_position)] += 1
                else:
                    #If not continuous new starting position is defined.
                    start_position = position
                    start_positions[(chromosome,start_position)] = 1
        self.no_alignments = start_positions

    def create_df(self):
        """Creates df for outputting or for annotating."""
        self.output = pd.DataFrame(columns=['chromosome','position','length','type'],index=range(len(self.no_alignments.values())))
        for counter,((chromosome,position),length) in enumerate(self.no_alignments.items()):
            self.output.at[counter,'chromosome'] = chromosome
            self.output.at[counter,'position'] = position
            self.output.at[counter,'length'] = length
            self.output.at[counter,'type'] = 'no alignment region'

    def write_no_alignments(self,df,out):
        """This function writes detected no alignment regions to tsv."""
        df.to_csv(os.path.join(out,'no_alignment_regions.tsv'),index=False,sep='\t')
