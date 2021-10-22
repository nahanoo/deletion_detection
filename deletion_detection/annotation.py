import pandas as pd
from deletion_detection import Deletion
from no_alignment import NoAlignment
from Bio import SeqIO
import numpy as np

class Annotation():
    def parse_genbank(self,genbank):
        contigs = [contig for contig in SeqIO.parse(genbank,'genbank')]
        return contigs

def annotate(features,df):
    out = pd.DataFrame(columns=df.columns.to_list()+\
        ['coding_type','nt_pos','aa_pos','gene','product'])
    loc = 0
    for i,row in df.iterrows():
        c,p,l = row['chromosome'],row['position'],row['length']
        for (start,end),feature in features[c].items():
            for nucleotide in range(start,end):
                if nucleotide in range(p,p+l):
                    try:
                        gene = feature['gene']
                    except KeyError:
                        gene = np.nan
                    row['position'] 
                    entries = [feature['type'],np.nan,np.nan,gene,feature['product']]
                    print(row.to_list()+entries)
                    out.loc[loc] = row.to_list()+entries
                    loc += 1
                    break
    return out

def parse_genbank():
    contigs = [contig for contig in SeqIO.parse('testdata/at.gbk','genbank')]
    features = {key:dict() for key in [contig.id for contig in contigs]}
    for contig in contigs:
        for f in contig.features:
            feature = dict()
            feature['type'] = f.type
            feature['strand'] = f.location.strand
            q = f.qualifiers
            if not f.type == 'source':
                try:
                    feature['gene'] = q['gene'][0]
                except KeyError:
                    pass
                try:
                    feature['product'] = q['product'][0]
                except KeyError:
                    pass
                features[contig.id][(int(f.location.start),int(f.location.end))] = feature
    return features

f = parse_genbank()
df = pd.concat([pd.read_csv('testdata/in_read_deletions.tsv',sep='\t'),
    pd.read_csv('testdata/no_alignment_regions.tsv',sep='\t')])
df.index = range(len(df))
deletions = {(chromosome,position):length for chromosome,position,length in zip(df['chromosome'],df['position'],df['length'])}
df = annotate(f,df)
    
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
