from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.motifs import create
from Bio.Phylo.TreeConstruction import DistanceCalculator

class GenomeAnalyzer:
    def __init__(self, fasta_path):
        self.records = list(SeqIO.parse(fasta_path, "fasta"))
        
    def gc_content(self):
        return {rec.id: GC(rec.seq) for rec in self.records}
    
    def find_motifs(self, motif_length=6):
        sequences = [str(rec.seq) for rec in self.records]
        motifs = create(sequences, 'ACGT')
        return motifs.consensus
    
    def build_phylogeny(self):
        aligner = DistanceCalculator('identity')
        matrix = aligner.get_distance(self.records)
        return matrix

    def orf_finder(self, min_length=300):
        for record in self.records:
            for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                for frame in range(3):
                    trans = nuc[frame:].translate()
                    trans_len = len(trans)
                    aa_start = 0
                    while aa_start < trans_len:
                        aa_end = trans.find("*", aa_start)
                        if aa_end == -1: aa_end = trans_len
                        if aa_end - aa_start >= min_length/3:
                            yield {
                                'id': record.id,
                                'start': frame + aa_start*3,
                                'end': frame + aa_end*3,
                                'strand': strand
                            }
                        aa_start = aa_end + 1
