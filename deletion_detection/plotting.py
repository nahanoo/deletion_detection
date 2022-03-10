from dna_features_viewer import GraphicFeature, GraphicRecord
import pysam
from os.path import join

def plot_alignment(bam, chromosome, position, length, window, out):
    """Plots alignment of regions with no coverage in
    reference."""

    # Gets correct padding
    a = pysam.AlignmentFile(bam)
    padding = window
    if position - padding < 0:
        start = 0
    else:
        start = position - padding

    ref_length = a.get_reference_length(chromosome)
    if position + length + padding > ref_length:
        end = ref_length
    else:
        end = position +length + padding

    # Gets reads of no coverage regions
    reads = []
    for read in a.fetch(chromosome, start, end):
        if (not read.is_unmapped) & (not read.is_secondary):
            reads.append(read)

    # Creates features for plotting
    features = []
    for read in reads:
        if read.is_reverse:
            strand = 1
        else:
            strand = -1
        gf = GraphicFeature(start=read.reference_start, end=read.reference_end, strand=strand,
                            color="#673AB7", label=read.query_name)
        features.append(gf)
    record = GraphicRecord(sequence_length=end, features=features)
    record = record.crop((start, end))
    f_name = '.'.join([chromosome, str(position), 'pdf'])
    record.plot_on_multiple_pages(join(out,f_name),
                                  nucl_per_line=end-start,
                                  lines_per_page=10,
                                  plot_sequence=False
                                  )

def plot_genbank(genbank_list, chromosome, start,end, out):
    """Plots genbank annotation of deleted sequence"""
    genbank = {contig.id: contig for contig in genbank_list}
    features = []
    for feature in genbank[chromosome].features:
        f_start = feature.location.start
        f_end = feature.location.end
        if (f_start >= start) & (f_end <= end):
            try:
                gf = GraphicFeature(start=f_start, end=f_end, strand=feature.location.strand,
                                    color="#673AB7", label=feature.qualifiers['product'][0])
                features.append(gf)
            except KeyError:
                pass
    record = GraphicRecord(sequence_length=end, features=features)
    record = record.crop((start, end))
    f_name = '.'.join([chromosome, str(start), 'pdf'])
    record.plot_on_multiple_pages(join(out, f_name),
                                  nucl_per_line=(end-start)/3,
                                  lines_per_page=10,
                                  plot_sequence=False
                                  )
