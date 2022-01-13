import pandas as pd
from Bio import SeqIO
import numpy as np

class Annotation():
    def parse_genbank(self,genbank):
        """Creating a dictionary with all the features from the genbank file."""
        self.contigs = [contig for contig in SeqIO.parse(genbank,'genbank')]
        self.features = {key:dict() for key in [contig.id for contig in self.contigs]}
        for contig in self.contigs:
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
                    #Storing the feature dictionary with start and end position of feature as key
                    self.features[contig.id][(int(f.location.start),int(f.location.end))] = feature
    
    def get_gc_content(self,chromosome,position,length):
        contigs = [contig.id for contig in self.contigs]
        sequence = self.contigs[contigs.index(chromosome)][position:position+length]
        if len(sequence.seq) > 0:
            return 100.0*len([base for base in sequence if base in "GC"])/len(sequence)
        if len(sequence.seq) == 0:
            return 'end of contig'

    def annotate(self,df):
        """Iterates over all deletions and no aignment regions and annotates
        deletions."""
        self.annotation = pd.DataFrame(columns=df.columns.to_list()+\
            ['gc_content','coding_type','strand','nt_pos','gene','product'])
        loc = 0
        for i,row in df.iterrows():
            #Iterating over all deletions
            c,p,l = row['chromosome'],row['position'],row['length']
            #We need to store if deletion is annotated in order that it is still outputted
            #in case no feature is found.
            #Iterating over all features. Follwing code is inefficient
            #but necessary to annotate long no alignment regions
            gc_content = self.get_gc_content(c,p,l)
            annotated = False
            for (start,end),feature in self.features[c].items():
                for nucleotide in range(start,end):
                    #Checks if any position in a feature is in a deletion
                    if nucleotide in range(p,p+l):
                        annotated = True
                        try:
                            gene = feature['gene']
                        except KeyError:
                            gene = np.nan
                        nt = [start,end]
                        nt = [str(n) for n in nt]
                        nt = nt[0]+'-'+nt[1]
                        entries = [gc_content,feature['type'],feature['strand'],nt,gene,feature['product']]
                        self.annotation.loc[loc] = row.to_list()+entries
                        loc += 1
                        break
            if not annotated:
                entries = [gc_content,np.nan,np.nan,np.nan,np.nan,np.nan]
                self.annotation.loc[loc] = row.to_list()+entries
                loc += 1