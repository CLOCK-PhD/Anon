#!/usr/bin/python3

import json
import pysam
from difflib import SequenceMatcher
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

class Kmer:
    def __init__(self, sequence):
        self.sequence = sequence

    @staticmethod
    def reverse_complement(seq):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement[base] for base in reversed(seq))

    def __repr__(self):
        return f"Kmer(sequence='{self.sequence}')"

class UnalignedKmer(Kmer):
    def __init__(self, sequence):
        super().__init__(sequence)

    def __repr__(self):
        return f"Unaligned{super().__repr__()}"


class SingleAlignedKmer(Kmer):
    def __init__(self, sequence, chromosome, start, end, strand):
        super().__init__(sequence)
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        if strand == '-':
            self.sequence = self.reverse_complement(self.sequence)

    def __repr__(self):
        return (f"SingleAligned{super().__repr__()}, chromosome='{self.chromosome}', "
                f"start={self.start}, end={self.end}, strand='{self.strand}')")


class MultiAlignedKmer(Kmer):
    def __init__(self, sequence, alignments):
        super().__init__(sequence)
        self.alignments = alignments  # alignments should be a list of tuple or objects containing the alignment details

    def __repr__(self):
        return f"MultiAligned{super().__repr__()}, alignments={len(self.alignments)}"



def adjust_positions(ref_seq, kmer, start, end):
    matcher = SequenceMatcher(None, kmer, ref_seq)
    match = matcher.find_longest_match(0, len(kmer), 0, len(ref_seq))
    
    # Calculate how much to extend the sequence fetch at both ends
    extend_start = match.a
    extend_end = len(kmer) - (match.a + match.size)

    # Adjust start and end ensuring they are within bounds
    new_start = max(0, start - extend_start)
    new_end = end + extend_end
    
    return new_start, new_end

# Conversion des noms des chromosomes pour correspondre à RefSeq
chromosome_map = {
    'chr1': 'NC_000001.11',
    'chr2': 'NC_000002.12',
    'chr3': 'NC_000003.12',
    'chr4': 'NC_000004.12',
    'chr5': 'NC_000005.10',
    'chr6': 'NC_000006.12',
    'chr7': 'NC_000007.14',
    'chr8': 'NC_000008.11',
    'chr9': 'NC_000009.12',
    'chr10': 'NC_000010.11',
    'chr11': 'NC_000011.10',
    'chr12': 'NC_000012.12',
    'chr13': 'NC_000013.11',
    'chr14': 'NC_000014.9',
    'chr15': 'NC_000015.10',
    'chr16': 'NC_000016.10',
    'chr17': 'NC_000017.11',
    'chr18': 'NC_000018.10',
    'chr19': 'NC_000019.10',
    'chr20': 'NC_000020.11',
    'chr21': 'NC_000021.9',
    'chr22': 'NC_000022.11',
    'chrX': 'NC_000023.11',
    'chrY': 'NC_000024.10'
}

# En fait c'est pareil que p14
"""chromosome_map_p13 = {
    'chr1': 'NC_000001.11',
    'chr2': 'NC_000002.12',
    'chr3': 'NC_000003.12',
    'chr4': 'NC_000004.12',
    'chr5': 'NC_000005.10',
    'chr6': 'NC_000006.12',
    'chr7': 'NC_000007.14',
    'chr8': 'NC_000008.11',
    'chr9': 'NC_000009.12',
    'chr10': 'NC_000010.11',
    'chr11': 'NC_000011.10',
    'chr12': 'NC_000012.12',
    'chr13': 'NC_000013.11',
    'chr14': 'NC_000014.9',
    'chr15': 'NC_000015.10',
    'chr16': 'NC_000016.10',
    'chr17': 'NC_000017.11',
    'chr18': 'NC_000018.10',
    'chr19': 'NC_000019.10',
    'chr20': 'NC_000020.11',
    'chr21': 'NC_000021.9',
    'chr22': 'NC_000022.11',
    'chrX': 'NC_000023.11',
    'chrY': 'NC_000024.10'
}"""


##############################################################################
# Load JSON data
with open('TCGA_BC_aggregated.json', 'r') as file:
#with open('TCGA_OV_aggregated.json', 'r') as file:
#with open('BEAUTY_aggregated.json', 'r') as file:
#with open('DLBCL_aggregated.json', 'r') as file:
    data = json.load(file)

# Access k-mer information
kmers_data = data['kmers']

count = 0
align_count = 0
no_align_count = 0
multimap_count = 0

# Load data and create Kmer objects
kmer_objects = []

for kmer_data in kmers_data:
    count += 1
    if not kmer_data['alignments']:
        kmer = UnalignedKmer(kmer_data['kmer'])
        no_align_count += 1
    elif len(kmer_data['alignments']) == 1:
        alignment = kmer_data['alignments'][0]
        kmer = SingleAlignedKmer(kmer_data['kmer'], alignment['chromosome'], alignment['start'], alignment['end'], alignment['strand'])
        align_count += 1
    else:
        kmer = MultiAlignedKmer(kmer_data['kmer'], kmer_data['alignments'])
        multimap_count += 1

    kmer_objects.append(kmer)

# Load the reference genome
# GRCh38p14
#reference_genome = pysam.FastaFile('/home/remycosta/phd/Anon/data/grch38p14/GCF_000001405.40_GRCh38.p14_genomic.fna')
# GRCh38p13
reference_genome = pysam.FastaFile('/home/remycosta/phd/Anon/data/GRCh38p13/GCF_000001405.39_GRCh38.p13_genomic.fna')
# FAIRE DES TRUCS

match_count = 0
no_match_count = 0
size_pb_count = 0
smaller_size = 0
not_found_count = 0

for kmer in kmer_objects:
    if isinstance(kmer, SingleAlignedKmer):
        # Map the chromosome name to the corresponding RefSeq ID
        refseq_chromosome = chromosome_map.get(kmer.chromosome)
        if refseq_chromosome:
            # Fetch the sequence from the reference genome
            ref_sequence = (reference_genome.fetch(refseq_chromosome, kmer.start, kmer.end)).upper()
            #print(f"Ref:\t{ref_sequence}")
            #print(f"k-mer:\t{kmer.sequence}")
            # Compare the sequences
            if kmer.sequence == ref_sequence:
                match_count += 1
            else:
                if len(kmer.sequence) > len(ref_sequence):
                    size_pb_count += 1
                    
                    print(f"Ref:\t{ref_sequence}")
                    print(f"k-mer:\t{kmer.sequence}")
                    # MAKE ALIGNMENT HERE
                    # Perform alignment considering potential insertions
                    # Custom scoring: match = 2, mismatch = -1, gap opening in ref = -2, gap extension in ref = -1, gap opening in k-mer = -0.5, gap extension in k-mer = -0.1
                    alignments = pairwise2.align.globalms(kmer.sequence, ref_sequence, 
                                                          2,    # score for a match
                                                          -1,   # score for a mismatch
                                                          -0.5, # penalty for opening a gap in k-mer (less penalizing)
                                                          -0.1) # penalty for extending a gap in k-mer (less penalizing)
                    
                    # Display the best alignment
                    if alignments:
                        print("Best alignment:")
                        print(format_alignment(*alignments[0]))
                    print("-----")

                    print("-----")
                else :
                    no_match_count += 1
        else:
            not_found_count += 1
            print(f"Chromosome {kmer.chromosome} not found in reference.")
            print("-----")


print("\n----- Recap -----")
print(f"Nombre de k-mers: {count}")
print(f"Nombre d'alignement uniques: {align_count}")
print(f"Nombre de k-mers multimap {multimap_count}")
print(f"Nombre de k-mers sans alignement : {no_align_count}")
print(f"Total : {align_count + multimap_count + no_align_count}")


print("----- Alignements uniques -----")
print(f"Match:\t\t{match_count}")
print(f"Size pb:\t{size_pb_count}")
print(f"No match:\t{no_match_count}")
print(f"Not found:\t{not_found_count}")
print(f"total:\t\t{match_count + no_match_count + size_pb_count}")
print("--------------------")