from setuptools import setup, find_packages

setup(name='deltec',
      version='1.0',
      description='Detects deletion and regions with no alignments in BAM file derived from haploid genomes.',
      author='Eric Ulrich',
      url='https://github.com/nahanoo/deletion_detection',
      packages=find_packages(),
      install_requires=['pandas',
                        'pysam']
     )