from setuptools import setup, find_packages

setup(name='deletion_detection',
      version='1.0',
      description='Detects deletions and regions with no alignments in BAM files derived from haploid genomes.',
      author='Eric Ulrich',
      url='https://github.com/nahanoo/deletion_detection',
      packages=['deletion_detection'],
      install_requires=['pandas',
                        'pysam'],
      entry_points={
          'console_scripts': [
              'detect_deletions = deletion_detection.main:main'
          ]
      }
     )
