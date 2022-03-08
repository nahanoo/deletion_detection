from os.path import join, exists
from os import mkdir, remove
from io import StringIO
from subprocess import call, run as r, DEVNULL, STDOUT
from Bio import SeqIO
import pandas as pd
from plotting import plot_alignment, plot_genbank
import pysam


class Deletion():
    def __init__(self, args):
        self.out_dir = args.out_dir
        self.reference_gbk = args.anceteral
        self.mutant_fasta = args.mutant

        self.step = 10000
        self.window = 50000

        self.mutant_contigs = [contig for contig in SeqIO.parse(
            self.mutant_fasta, 'fasta')]

        self.reference_contigs = [contig for contig in SeqIO.parse(
            self.reference_gbk, 'genbank')]
        self.reference_fasta = join(self.out_dir, 'reference.fasta')
        with open(join(self.out_dir, 'reference.fasta'), 'w') as handle:
            SeqIO.write(self.reference_contigs, handle, 'fasta')
        self.genbank = self.parse_genbank()

        self.bam = join(self.out_dir, "aligned.sorted.bam")
        self.deletions = None

        self.annotated = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'start', 'end', 'protein'])
        self.deleted_plasmids = pd.DataFrame(
            columns=['chromosome', 'position', 'length'])
        self.deleted_plasmids_annotated = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'start', 'end', 'protein'])

        if not exists(join(self.out_dir, 'plots')):
            mkdir(join(self.out_dir, 'plots'))

    def parse_genbank(self):
        genbank = {contig.id: {} for contig in self.reference_contigs}
        for contig in self.reference_contigs:
            for feature in contig.features:
                try:
                    start = feature.location.start
                    end = feature.location.end
                    product = feature.qualifiers['product']
                    genbank[contig.id][(start, end)] = product[0]
                except KeyError:
                    pass
        return genbank

    def chunker(self, seq, window_size, step):
        """Creates chunks of a sequence. window_size defines
        chunk size and step the amount of basepairs the window
        is moved forward."""
        # List which stores all chunks
        seqs = []
        seqlen = len(seq)
        self.step = step
        for counter, q in enumerate(range(0, seqlen, step)):
            # Returns ether entire sequence or window depending on sequence length
            j = seqlen if q + window_size > seqlen else q + window_size
            chunk = seq[q:j]
            # Add chunk id to sequence id
            chunk.id = chunk.id + "." + str(counter)
            seqs.append(chunk)
            if j == seqlen:
                break
        return seqs

    def chunk_assembly(self):
        """Chunks an assembly of multiple contigs into different 
        chunks using a sliding window algorightm (see chunker function)."""
        assembly_chunks = []
        for contig in self.mutant_contigs:
            # Creates chunks of every contig
            assembly_chunks += self.chunker(contig, self.window, self.step)
        target = join(self.out_dir,
                      "chunked_sequences.fasta")
        # Dumps chunks to fasta
        with open(target, "w") as handle:
            SeqIO.write(assembly_chunks, handle, "fasta")

    def mapper(self):
        """Maps chunked sequence to all ancesteral genomes
        experiment with minimap2.
        Minimap2 settings are set to accureate PacBio reads."""
        reads = join(self.out_dir, "chunked_sequences.fasta")
        sam = join(self.out_dir, "aligned.sam")
        cmd = [
            "minimap2",
            "-ax",
            "asm5",
            self.reference_fasta,
            reads,
            ">",
            sam,
        ]
        # Calling minimap and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'sort', '-o', self.bam, sam]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'index', self.bam]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)

    def get_deletions(self):
        cmd = ['samtools', 'depth', '-aa', '-Q', '0', self.bam]
        process = r(cmd, capture_output=True)
        df = pd.read_csv(StringIO(process.stdout.decode()), sep='\t')
        df.columns = ['chromosome', 'position', 'coverage']
        df = df[df['coverage'] == 0]
        self.deletions = self.concat_deletions(df)
        self.deletions['position'] = self.deletions['position'] - 1
        mask = []
        a = pysam.AlignmentFile(self.bam)
        for c, p, l in zip(self.deletions['chromosome'], self.deletions['position'], self.deletions['length']):
            if l == a.get_reference_length(c):
                self.deleted_plasmids.loc[len(self.deleted_plasmids)] = [c,p,l]
                mask.append(False)
            else:
                reads = []
                for read in a.fetch(c, p-100, p+l+100):
                    reads.append(read.qname)
                deletion_break = self.check_break(reads)
                if l < 500:
                    deletion_break = False
                mask.append(deletion_break)
        self.deletions = self.deletions[mask]
        self.deletions.to_csv(
            join(self.out_dir, 'deletions.tsv'), sep='\t', index=False)

    def concat_deletions(self, df):
        out = pd.DataFrame(columns=['chromosome', 'position', 'length'])
        i = -1
        prev_pos = 0
        prev_contig = None
        for contig, pos in zip(df['chromosome'], df['position']):
            if (prev_contig == contig) & (pos - 1 == prev_pos):
                out.at[i, 'length'] += 1
            else:
                i += 1
                out.loc[i] = [contig, pos, 1]
            prev_pos = pos
            prev_contig = contig
        return out

    def check_break(self, reads):
        contigs = ['.'.join(read.split('.')[:-1]) for read in reads]
        if len(set(contigs)) > 1:
            return False

        index = sorted([int(read.split('.')[-1]) for read in reads])
        prev = index[0]
        for i in index[1:]:
            if i - 1 == prev:
                pass
            else:
                return False
            prev += 1

        return True

    def annotate(self):
        i = 0
        for counter, row in self.deletions.iterrows():
            c = row['chromosome']
            p = row['position']
            l = row['length']
            a = pysam.AlignmentFile(self.bam)
            for read in a.fetch(c, p-1, p):
                c_mut = '.'.join(read.qname.split('.')[:-1])
                p_mut = int(read.qname.split('.')[-1]) * self.step + read.qend
                break
        
            for (start, end), product in self.genbank[c].items():
                if (start >= p) & (end <= p + l):
                    self.annotated.loc[i] = [
                        c_mut, p_mut, l, start, end, product]
                    i += 1
        self.annotated = self.annotated.drop_duplicates()
        self.annotated.to_csv(
            join(self.out_dir, 'deletions.annotated.tsv'), sep='\t', index=False)


        print(self.deleted_plasmids)
        i = 0
        for counter, row in self.deleted_plasmids.iterrows():
            c = row['chromosome']
            p = row['position']
            l = row['length']
            for (start, end), product in self.genbank[c].items():
                if (start >= p) & (end <= p + l):
                    self.deleted_plasmids_annotated.loc[i] = [
                        c_mut, p_mut, l, start, end, product]
                    i += 1
        self.deleted_plasmids_annotated = self.deleted_plasmids_annotated.drop_duplicates()
        self.deleted_plasmids_annotated.to_csv(
            join(self.out_dir, 'deleted_plasmids.annotated.tsv'), sep='\t', index=False)

    def plot_deletions(self):
        out = join(self.out_dir, 'plots', 'alignments')
        if not exists(out):
            mkdir(out)

        for chromosome, position in zip(self.deletions['chromosome'], self.deletions['position']):
            plot_alignment(self.bam, chromosome, position, out)

    def plot_annotation(self):
        out = join(self.out_dir, 'plots', 'annotations')
        if not exists(out):
            mkdir(out)
        for i, row in self.deletions.iterrows():
            plot_genbank(
                self.reference_contigs, row['chromosome'], row['position'],
                row['position']+row['length'], out)

    def clean(self):
        remove(self.mutant_fasta)
