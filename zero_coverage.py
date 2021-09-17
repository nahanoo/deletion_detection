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
        """This function calucaltes coverage using samtools depth.
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

    def get_zero_coverage_length(self,zero_coverage):
        pass

    def get_zero_coverage(self):
        df = self.coverage[self.coverage['coverage'] == 0]



c = Coverage()
c.get_coverage('testdata/mapped_reads.sorted.bam')

