import pandas as pd
import plotly.express as px
from get_samples import Samples
import os

s = Samples()

class Plot():
    """Plotting class for black queen hypothesis analysis."""
    def boxplot(self,df):
        """This boxplot function takes a df as argument which typically has
        treatments as columns and samples as rows."""
        fig = px.box(df,points='all')
        return fig


def parse_deletions(name_dir_tuples,fname):
    base_dir = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'
    treatments = list(set([s.treatments[name] for name,folder in name_dir_tuples]))
    out = pd.DataFrame(columns=treatments,index=[name for name,folder in name_dir_tuples])
    for name,folder in name_dir_tuples:
        treatment = s.treatments[name]
        df = pd.read_csv(os.path.join(base_dir,folder,fname),sep='\t')
        out.at[name,treatment] = df['length'].sum()
    out.fillna(0)
    return out

def plot_no_alignments():
    pass 

p = Plot()
df = parse_deletions(s.ct_name_dir,'no_alignment_regions.tsv')
f = p.boxplot(df)
f.write_image('ct.png')
df = parse_deletions(s.ct_name_dir,'in_read_deletions.tsv')
f = p.boxplot(df)
f.write_image('in_read_ct.png')

df = parse_deletions(s.at_name_dir,'no_alignment_regions.tsv')
f = p.boxplot(df)
f.write_image('at.png')
df = parse_deletions(s.at_name_dir,'in_read_deletions.tsv')
f = p.boxplot(df)
f.write_image('in_read_at.png')

df = parse_deletions(s.oa_name_dir,'no_alignment_regions.tsv')
f = p.boxplot(df)
f.write_image('oa.png')
df = parse_deletions(s.oa_name_dir,'in_read_deletions.tsv')
f = p.boxplot(df)
f.write_image('in_read_oa.png')

df = parse_deletions(s.ms_name_dir,'no_alignment_regions.tsv')
f = p.boxplot(df)
f.write_image('ms.png')
df = parse_deletions(s.ms_name_dir,'in_read_deletions.tsv')
f = p.boxplot(df)
f.write_image('in_read_ms.png')
