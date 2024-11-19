#!/usr/bin/python3

import random
import argparse

def parse_vcf(vcf_file):
    """
    Parse the VCF file to extract REF alleles and their positions by chromosome.
    """
    chrom_refs = {}
    print(f"Opening VCF file: {vcf_file}")
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                print(f"Skipping malformed line: {repr(line)}")
                continue  # Skip lines that do not have enough columns
            chrom, pos, ref = parts[0], int(parts[1]), parts[3]
            if chrom not in chrom_refs:
                chrom_refs[chrom] = {}
            chrom_refs[chrom][pos] = ref
            print(f"Parsed line: {chrom}, {pos}, {ref}")
    print("Finished parsing VCF file.")
    return chrom_refs

def generate_random_sequence_with_refs(chrom, positions, length):
    """
    Generate a random DNA sequence with REF alleles at the specified positions for a given chromosome.
    """
    nucleotides = ['a', 'c', 'g', 't']
    sequence = []
    current_pos = 1

    for pos in sorted(positions):
        # Add random nucleotides until we reach the position
        while current_pos < pos:
            sequence.append(random.choice(nucleotides))
            current_pos += 1
        # Add the reference allele at the specified position
        ref = positions[pos]
        for nucleotide in ref:
            if current_pos <= length:
                sequence.append(nucleotide)
                current_pos += 1
            else:
                break

    # Fill the rest of the sequence with random nucleotides if necessary
    while len(sequence) < length:
        sequence.append(random.choice(nucleotides))

    return ''.join(sequence[:length])

def print_fasta(chrom, sequence):
    """
    Print the sequence in FASTA format.
    """
    print(f">{chrom}")
    for i in range(0, len(sequence), 60):
        print(sequence[i:i+60])

def main():
    parser = argparse.ArgumentParser(description="Generate random DNA sequences with REF alleles at specified positions from a VCF file.")
    parser.add_argument('vcf_file', help="Path to the VCF file")
    
    args = parser.parse_args()
    vcf_file = args.vcf_file
    
    chrom_refs = parse_vcf(vcf_file)
    
    # Display a summary of the information found for each chromosome
    for chrom, positions in chrom_refs.items():
        print(f"Chromosome: {chrom}")
        print(f"Number of positions: {len(positions)}")
        for pos, ref in positions.items():
            print(f"  Position: {pos}, REF: {ref}")
        print()
    
    # Ask for sequence length and generate sequence for each chromosome
    for chrom, positions in chrom_refs.items():
        length = int(input(f"Enter the length of the DNA sequence for chromosome {chrom}: "))
        sequence = generate_random_sequence_with_refs(chrom, positions, length)
        print_fasta(chrom, sequence)

if __name__ == "__main__":
    main()
