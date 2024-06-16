#!/usr/bin/python3

import random

def generate_random_sequence(length):
    """
    Generate a random DNA sequence of a given length.
    """
    nucleotides = ['a', 't', 'c', 'g']
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    return sequence

def print_fasta(sequence, header=">random_sequence"):
    """
    Print the sequence in FASTA format.
    """
    print(header)
    for i in range(0, len(sequence), 60):
        print(sequence[i:i+60])

def main():
    # Define the length of the random sequence
    length = int(input("Enter the length of the DNA sequence: "))
    sequence = generate_random_sequence(length)
    print_fasta(sequence)

if __name__ == "__main__":
    main()
