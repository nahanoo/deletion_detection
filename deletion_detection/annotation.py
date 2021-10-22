import pandas as pd
from deletion_detection import Deletion
from no_alignment import NoAlignment
from Bio import SeqIO

class Annotation():
    def __init__(self,input):
        self.input = input

    def normalize_input(self):
        raw = self.input
        self.input = {key:None for key in raw.keys()}
        for key,value in raw.items():
            self.input[key] = value[1]

    def parse_genbank(self,genbank):
        contigs = [contig for contig in SeqIO.parse(genbank,'genbank')]
        return contigs

def parse_genbank():
    features = dict()
    contigs = [contig for contig in SeqIO.parse('testdata/at.gbk','genbank')]
    for contig in contigs:
        for f in contig.features:
            feature = dict()
            feature['start'] = f.location.start
            feature['end'] = f.location.end
            feature['strand'] = f.location.strand
            feature['type'] = f.type
            q = f.qualifiers
            try:
                feature['gene'] = q['gene']
            except KeyError:
                pass
            try:
                feature['inference'] = q['inference']
            except KeyError:
                pass
            try:
                feature['product'] = q['product']
            except KeyError:
                pass
            features[(contig.id,f.location.start)] = feature
    return features
f = parse_genbank()
    
"""
d = Deletion()
d.get_deletions('testdata/mapped_corrected_reads.sorted.bam')
d.filter_deletions(10,0.8)

n = NoAlignment()
n.get_coverage('testdata/mapped_corrected_reads.sorted.bam',60)
n.get_no_alignments()
n.get_no_alignments_length()
a = Annotation(d.filtered_deletions)
a.normalize_input()
c = a.parse_genbank('testdata/at.gbk')"""
