from setuptools import setup

setup(name='deletion_detection',
      version='1.0',
      description='Detects deletions from assemblies of evolved bacteria.',
      author='Eric Ulrich',
      url='https://github.com/nahanoo/deletion_detection',
      packages=['deletion_detection'],
      install_requires=['pandas',
                        'pysam',
                        'Bio',
                        'dna_features_viewer'],
      entry_points={
          'console_scripts': [
              'detect_deletions = deletion_detection.main:main'
          ]
      }
     )
