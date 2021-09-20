import pysam
from io import StringIO
import pandas as pd
import subprocess

class Coverage():
    """This class is used to detect areas with zero coverage.
    This is useful if you want to find plasmids or large 
    genomic regions which were potentially deleted, something
    you would not find by investigating aligned reads.
    Remember that contigs at the start and the end often have
    poor assembly quality which can lead to zero coverage. 
    """

    def get_coverage(self,path):
        """This function calculates coverage using samtools depth.
        Depending on the size of the bam file this can take a while
        and you may want to run it on a cluster. The flag "-J" is used
        which means that deletions in reads count as coverage.
        Therefore coverage is only zero in areas where no reads align.
        """
        self.bamfile_path = path
        cmd = ['samtools','depth','-a','-J',self.bamfile_path]
        process = subprocess.run(cmd,capture_output=True)
        self.coverage = pd.read_csv(StringIO(process.stdout.decode()),sep='\t')
        self.coverage.columns=['chromosome','position','coverage']

    def get_length(self):
        #Splitting df in dictionary with contigs as keys
        zero_coverage = {chromosome:list() for chromosome in set(c.zero_coverage['chromosome'])}
        for chromosome,pos in zip(self.zero_coverage['chromosome'],self.zero_coverage['position']):
            zero_coverage[chromosome].append(pos)
        start_positions = dict()
        for chromosome,positions in zero_coverage.items():
            start_position = positions[0]
            start_positions[start_position] = 1
            previous_position = start_position
            for pos in positions[1:]:
                if previous_position + 1 == pos:
                    start_positions[start_position] += 1
                    previous_position = pos
                else:
                    start_position = pos
                    start_positions[start_position] = 1
                    previous_position = pos
        return zero_coverage

    def get_zero_coverage(self):
        """Hero we get all positions with zero reads aligned.
        """
        self.zero_coverage = self.coverage[self.coverage['coverage'] == 0]


c = Coverage()
c.get_coverage('testdata/at/mapped_reads.sorted.bam')
c.get_zero_coverage()
z = c.get_length()


