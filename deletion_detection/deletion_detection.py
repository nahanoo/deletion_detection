from os.path import join, exists
from os import mkdir, remove
from io import StringIO
from subprocess import call, run as r, DEVNULL, STDOUT
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from .plotting import plot_alignment, plot_genbank
import pysam


class Deletion():
    def __init__(self, args):
        # Input file paths
        self.out_dir = args.out_dir
        self.reference_gbk = args.ancestral
        self.mutant_fasta = args.mutant
        # Bam file path
        self.bam = join(self.out_dir, "aligned.sorted.bam")

        # Chunk parameters
        self.step = 100
        self.window = 500

        # List contigs
        self.mutant_contigs = [contig for contig in SeqIO.parse(
            self.mutant_fasta, 'fasta')]
        self.reference_contigs = [contig for contig in SeqIO.parse(
            self.reference_gbk, 'genbank')]
        self.reference_features = self.parse_genbank()
        # Dumping input reference as fasta for mapping
        self.reference_fasta = join(self.out_dir, 'reference.fasta')
        with open(join(self.out_dir, 'reference.fasta'), 'w') as handle:
            SeqIO.write(self.reference_contigs, handle, 'fasta')

        # Output dataframes
        self.no_coverage = pd.DataFrame(
            columns=['chromosome', 'position', 'length'])
        self.deletions = pd.DataFrame(
            columns=['chromosome', 'position',
                     'length', 'chromosome_origin', 'position_origin'])
        self.deletions_annotated = pd.DataFrame(
            columns=['chromosome', 'position',
                     'length', 'chromosome_origin', 'position_origin', 'products'])
        self.plasmids = pd.DataFrame(
            columns=['chromosome', 'position', 'length'])
        self.plasmids_annotated = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'products'])

        # List storing files paths to trash
        self.trash = []
        self.trash.append(self.reference_fasta)

        # Create plot directory
        if not exists(join(self.out_dir, 'plots')):
            mkdir(join(self.out_dir, 'plots'))

    def parse_genbank(self):
        """Parses all features of a genbanka and stores locations
        and products in dictionary"""
        genbank = {contig.id: {} for contig in self.reference_contigs}
        for contig in self.reference_contigs:
            for feature in contig.features:
                # Some features don't have all desired keys
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
            chunk_id = seq.id
            chunk_seq = seq.seq[q:j]
            # Add chunk id to sequence id
            chunk = SeqRecord(seq=chunk_seq,id=chunk_id + "." + str(counter))
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
        self.chunks = join(self.out_dir,
                      "chunked_sequences.fasta")
        # Dumps chunks to fasta
        with open(self.chunks, "w") as handle:
            SeqIO.write(assembly_chunks, handle, "fasta")
        # Delete chunked sequence
        self.trash.append(self.chunks)

    def mapper(self, reference, reads, out):
        """Maps long accurate sequences to references with minimap2."""
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
        # Files to trash
        if out not in self.trash:
            self.trash.append(out)
            self.trash.append(bam)
            self.trash.append(bam+'.bai')

    def map_chunks(self):
        """Maps chunked sequences to reference"""
        self.mapper(self.reference_fasta, join(self.out_dir, "chunked_sequences.fasta"),
                    join(self.out_dir, "aligned.sam"))

    def get_deletions(self):
        """Gets all regions with 0 coverage and seperates between
        in strand deletion or deleted plasmids."""
        cmd = ['samtools', 'depth', '-aa', '-Q', '0', self.bam]
        process = r(cmd, capture_output=True)
        df = pd.read_csv(StringIO(process.stdout.decode()), sep='\t')
        df.columns = ['chromosome', 'position', 'coverage']
        # Masks for regions with no coverage
        df = df[df['coverage'] == 0]
        # Concats positions into start position + length
        self.concat_deletions(df)
        # Switching to 0 based index
        self.no_coverage['position'] = self.no_coverage['position'] - 1
        a = pysam.AlignmentFile(self.bam)
        for i, row in self.no_coverage.iterrows():
            c, p, l = row
            if l == a.get_reference_length(c):
                # Entire plasmid is deleted
                self.plasmids.loc[len(self.plasmids)] = [c, p, l]
            else:
                self.deletions.loc[len(self.deletions)] = [None, None, l, c, p]

    def concat_deletions(self, df):
        """Concats following positions into deletions
        with length equals n following positions."""
        i = -1
        prev_pos = 0
        prev_contig = None
        for contig, pos in zip(df['chromosome'], df['position']):
            if (prev_contig == contig) & (pos - 1 == prev_pos):
                self.no_coverage.at[i, 'length'] += 1
            else:
                i += 1
                self.no_coverage.loc[i] = [contig, pos, 1]
            prev_pos = pos
            prev_contig = contig

    def get_location(self):
        """Extracts sequence from reference around detected deletion.
        Sequence is aligned to mutant to find position in the mutant."""
        # Switching to dict style
        reference_contigs = {
            contig.id: contig for contig in self.reference_contigs}
        i = 0
        deletions = pd.DataFrame(columns=self.deletions.columns)
        for c, p, l in zip(self.deletions['chromosome_origin'],
                           self.deletions['position_origin'], self.deletions['length']):
            a = pysam.AlignmentFile(self.bam)
            # Getting correct paddins considering pos possiblye
            # being at the start or the end of a contig
            padding = 2000
            if p - padding < 0:
                start = 0
            else:
                start = p - padding
            if p + padding > a.get_reference_length(c):
                end = a.get_reference_length(c)
            else:
                end = p + padding

            # Getting sequence around deletion in reference
            seq = reference_contigs[c][start:p]
            seq += reference_contigs[c][p+l:end]
            seq_id = c + '.' + str(p)
            seq.id = seq_id

            # Location of deletion in sequence
            rel_pos = p - start

            # Mapping files
            tmp_seq = join(self.out_dir, 'tmp_seq.fasta')
            tmp_sam = join(self.out_dir, 'tmp_seq.sam')
            tmp_bam = tmp_sam.replace('.sam', '.sorted.bam')
            with open(tmp_seq, 'w') as handle:
                SeqIO.write(seq, handle, 'fasta')

            # Mapping to mutant
            self.mapper(self.mutant_fasta, tmp_seq, tmp_sam)

            a = pysam.AlignmentFile(tmp_bam)
            for read in a:
                if not (read.is_unmapped):
                    # Idetifying positions
                    deletions.loc[i] = [
                        read.reference_name, read.reference_start + rel_pos, l, c, p]
                    i += 1
        self.deletions = deletions
        # Deleting tmp_seq
        self.trash.append(tmp_seq)

    def annotate(self):
        """Annotates the dataframes of deleted plasmids and
        in strand deletions."""
        i = 0
        for c, p, l in zip(self.plasmids['chromosome'], self.plasmids['position'], self.plasmids['length']):
            products = self.annotate_position(c, p, l)
            for product in products:
                self.plasmids_annotated.loc[i] = [c, p, l, product]
                i += 1

        i = 0
        for counter, row in self.deletions.iterrows():
            c, p, l, c_o, p_o = row
            products = self.annotate_position(c_o, p_o, l)
            for product in products:
                self.deletions_annotated.loc[i] = [c, p, l, c_o, p_o, product]
                i += 1

    def annotate_position(self, c, p, l):
        """Returns products in a region in a genbank."""
        products = []
        for (start, end), product in self.reference_features[c].items():
            if not set(range(start, end)).isdisjoint(range(p, p+l)):
                products.append(product)
        return products

    def dump(self):
        """Dumping all dataframes."""
        self.no_coverage.to_csv(
            join(self.out_dir, 'no_coverage.tsv'), sep='\t', index=False)

        self.deletions.to_csv(
            join(self.out_dir, 'deletions.tsv'), sep='\t', index=False)
        self.deletions_annotated.to_csv(
            join(self.out_dir, 'deletions.annotated.tsv'), sep='\t', index=False)

        self.plasmids.to_csv(
            join(self.out_dir, 'plasmids.tsv'), sep='\t', index=False)
        self.plasmids_annotated.to_csv(
            join(self.out_dir, 'plasmids.annotated.tsv'), sep='\t', index=False)

    def plot_deletions(self):
        """Plots alignments."""
        out = join(self.out_dir, 'plots', 'alignments')
        if not exists(out):
            mkdir(out)

        for chromosome, position, length in zip(self.deletions['chromosome_origin'], 
        self.deletions['position_origin'], self.deletions['length']):
            plot_alignment(self.bam, chromosome, position, length, self.window, out)

    def plot_annotation(self):
        """Plots annotation"""
        out = join(self.out_dir, 'plots', 'annotations')
        if not exists(out):
            mkdir(out)
        for i, row in self.deletions.iterrows():
            plot_genbank(
                self.reference_contigs, row['chromosome_origin'], row['position_origin'],
                row['position_origin']+row['length'], out)

    def clean(self):
        """Deletes temporary files."""
        for item in self.trash:
            remove(item)
