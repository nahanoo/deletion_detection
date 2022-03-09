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
        self.ref_deletions = None
        self.no_coverage = None
        self.mut_deletions = pd.DataFrame(
            columns=['chromosome', 'position', 'length'])
        self.mut_deletions_annotated = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'products'])

        self.plasmids = pd.DataFrame(
            columns=['chromosome', 'position', 'length'])
        self.plasmids_annotated = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'products'])

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

    def mapper(self, reference, reads, out):
        """Maps chunked sequence to all ancesteral genomes
        experiment with minimap2.
        Minimap2 settings are set to accureate PacBio reads."""
        cmd = [
            "minimap2",
            "-ax",
            "asm5",
            reference,
            reads,
            ">",
            out,
        ]
        bam = out.replace('.sam', '.sorted.bam')
        # Calling minimap and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'sort', '-o', bam, out]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'index', bam]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)

    def get_deletions(self):
        cmd = ['samtools', 'depth', '-aa', '-Q', '0', self.bam]
        process = r(cmd, capture_output=True)
        df = pd.read_csv(StringIO(process.stdout.decode()), sep='\t')
        df.columns = ['chromosome', 'position', 'coverage']
        df = df[df['coverage'] == 0]
        self.ref_deletions = self.concat_deletions(df)
        self.ref_deletions['position'] = self.ref_deletions['position'] - 1
        self.no_coverage = self.ref_deletions
        a = pysam.AlignmentFile(self.bam)
        mask = []
        for c, p, l in zip(self.ref_deletions['chromosome'], self.ref_deletions['position'], self.ref_deletions['length']):
            if l == a.get_reference_length(c):
                self.plasmids.loc[len(self.plasmids)] = [c, p, l]
                mask.append(False)
            else:
                mask.append(True)
        self.ref_deletions = self.ref_deletions[mask]
        for c, p, l in zip(self.plasmids['chromosome'], self.plasmids['position'], self.plasmids['length']):
            i = len(self.plasmids_annotated)
            products = self.annotate_position(c, p, l)
            self.plasmids_annotated.loc[i] = [c, p, l, products]

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

    def get_mut_origin(self):
        reference_contigs = {
            contig.id: contig for contig in self.reference_contigs}
        for c, p, l in zip(self.ref_deletions['chromosome'], self.ref_deletions['position'], self.ref_deletions['length']):
            a = pysam.AlignmentFile(self.bam)
            padding = 2000
            if p - 2000 < 0:
                start = 0
            else:
                start = p - 2000
            if p + 2000 > a.get_reference_length(c):
                end = a.get_reference_length(c)
            else:
                end = p + 2000
            seq = reference_contigs[c][start:p]
            seq += reference_contigs[c][p+l:end]
            seq_id = c + '.' + str(p)
            seq.id = seq_id
            rel_pos = p - start
            tmp_seq = join(self.out_dir, 'tmp_seq.fasta')
            tmp_sam = join(self.out_dir, 'tmp_seq.sam')
            tmp_bam = tmp_sam.replace('.sam', '.sorted.bam')
            with open(tmp_seq, 'w') as handle:
                SeqIO.write(seq, handle, 'fasta')
            self.mapper(self.mutant_fasta, tmp_seq, tmp_sam)
            a = pysam.AlignmentFile(tmp_bam)
            for read in a:
                if not (read.is_unmapped):
                    i = len(self.mut_deletions)
                    self.mut_deletions.loc[i] = [
                        read.reference_name, read.reference_start + rel_pos, l]
                products = self.annotate_position(c, p, l)
                self.mut_deletions_annotated.loc[i] = [
                    read.reference_name, read.reference_start + rel_pos, l, products]

    def annotate_position(self, c, p, l):
        products = []
        for (start, end), product in self.genbank[c].items():
            if not set(range(start, end)).isdisjoint(range(p, p+l)):
                products.append(product)
        return products

    def format_out(self, df):
        out = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'product'])
        i = 0
        for c, p, l, products in zip(df['chromosome'], df['position'], df['length'], df['products']):
            for product in products:
                out.loc[i] = [c, p, l, product]
                i += 1
        return out

    def dump(self):
        self.no_coverage.to_csv(
            join(self.out_dir, 'no_coverage.tsv'), sep='\t', index=False)
        self.mut_deletions.to_csv(
            join(self.out_dir, 'deletions.tsv'), sep='\t', index=False)
        self.format_out(self.mut_deletions_annotated).to_csv(
            join(self.out_dir, 'deletions.annotated.tsv'), sep='\t', index=False)
        self.plasmids.to_csv(
            join(self.out_dir, 'plasmids.tsv'), sep='\t', index=False)
        self.format_out(self.plasmids_annotated).to_csv(
            join(self.out_dir, 'plasmids.annotated.tsv'), sep='\t', index=False)

    def plot_deletions(self):
        out = join(self.out_dir, 'plots', 'alignments')
        if not exists(out):
            mkdir(out)

        for chromosome, position in zip(self.ref_deletions['chromosome'], self.ref_deletions['position']):
            plot_alignment(self.bam, chromosome, position, out)

    def plot_annotation(self):
        out = join(self.out_dir, 'plots', 'annotations')
        if not exists(out):
            mkdir(out)
        for i, row in self.ref_deletions.iterrows():
            plot_genbank(
                self.reference_contigs, row['chromosome'], row['position'],
                row['position']+row['length'], out)

    def clean(self):
        remove(self.mutant_fasta)
