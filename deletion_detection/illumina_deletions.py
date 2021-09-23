import pandas as pd
from get_samples import Samples
import glob

s = Samples()
base_dir = '/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/'

def split_to_files():
    df = pd.read_csv('deletions.tsv',\
        names=['sample','chromosome','position','t1','t2','t3','t4'],sep='\t')
    df = df.loc[((df['t2']==1) | (df['t2']==0)) & (df['t3']==1)]
    df.insert(0,'strain',None)
    for i,chromosome,sample in zip(df.index,df['chromosome'],df['sample']):
        if chromosome[0] == 'C':
            strain = 'Ct'+str(sample).replace('.','')
            df.at[i,'strain'] = strain
        if chromosome[0] == 'M':
            df.at[i,'strain'] = 'Ms'+str(sample).replace('.','')
        if chromosome[0] == 'O':
            df.at[i,'strain'] = 'Oa'+str(sample).replace('.','')
        if chromosome[0] == 'A':
            df.at[i,'strain'] = 'At'+str(sample).replace('.','')
    name_dirs = s.ct_name_dir+s.at_name_dir+s.oa_name_dir+s.ms_name_dir
    for strain in set(df['strain']):
        out = df[df['strain'] == strain]
        keys = [key for key,folder in name_dirs if key[0:4] == strain]
        for key in keys:
            out.to_csv(base_dir+dict(name_dirs)[key]+'/all_deletions.tsv',\
                columns=['strain','chromosome','position'],sep='\t',index=False)

def regroup_deletions():
    for f in glob.glob(base_dir+'*/all_deletions.tsv'):
        df = pd.read_csv(f,sep='\t')
        deletions = {chromosome:list() for chromosome in set(df['chromosome'])}
        for chromosome,position in zip(df['chromosome'],df['position']):
            deletions[chromosome].append(position)

        #All starting positions are stored in a dictionary.
        #Values are the length of the continuous no alignment region.
        start_positions = dict()
        for chromosome,positions in deletions.items():
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
        out = pd.DataFrame(columns=['chromosome','position','length'],\
            index=range(len(start_positions)))
        for counter,((chromosome,position),length) in enumerate(start_positions.items()):
            out.at[counter,'chromosome'] = chromosome
            out.at[counter,'position'] = position
            out.at[counter,'length'] = length
        d = '/'.join(f.split('/')[:-1])
        out.to_csv(d+'/grouped_deletions.tsv',sep='\t',index=False)