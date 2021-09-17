import pysam

class Deletion():
    def __init__(self):
        self.deletions = []

class Unmapped():
    def __init__(self):
        self.contig = None
        self.position = None
        self.length = None
        self.reference_sequence = None



class Detector():

    def __init__(self):
        self.deletions = []

    def read_sam_file(self,path):
        return pysam.AlignmentFile(path,'rb')

    def get_deletions(self,alignment):
        # we iterate over every read of the alignment
        for read in alignment:
            # we set the block position to zero which stands
            # for the position in the reference
            block_pos = 0
            # we set the cigar position to zero which stands
            # for the position within the cigar code
            cigar_pos = 0
            # we iterate over the cigar code of every read
            # cigar code contains the information about alignment
            # (match)
            for (block_type, block_len) in read.cigar:
                # here we grab whether the alignment is match
                # deletions insertion and so on
                # depending on the block type we need to move forward
                # on the reference side (block_pos) or the query
                # side (cigar_pos). there is a very useful table about this
                # on https://samtools.github.io/hts-specs/SAMv1.pdf page 8
                # 0 is match or mismatch
                if block_type == 0:
                    block_pos += block_len
                    cigar_pos += block_len
                # 1 is insertion to the reference
                elif block_type == 1:
                    cigar_pos += block_len
                # 2 deletion from the reference
                elif block_type == 2:
                    deletion = dict()
                    deletion['contig'] = read.rname
                    deletion['position'] = read.pos+block_pos
                    deletion['length'] = block_len
                    self.deletions.append(deletion)
                    block_pos += block_len
                # 3 is skipped region from the reference
                elif block_type == 3:
                    block_pos += block_len
                # 4 is soft clipping
                elif block_type == 4:
                    cigar_pos += block_len
                
    def sort_deletions(self,deletions):
        return sorted(deletions, key=lambda k: k['length'],reverse=True)

def test():
    d = Deltec()
    alignment = d.read_sam_file('../testdata/alignments/aligned.sam')
    d.get_deletions(alignment)
    return d.sort_deletions(d.deletions)