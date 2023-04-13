#!/usr/bin/python3

"""
Genome K-mer Generator

Programme de test pour générer les k-mers uniques de tout le génome
dans un dictionnaire.
"""

from Bio import SeqIO

def main():

# Lecture du fichier fasta
    seq = []
    with open(fastaFile) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq.upper())

if __name__ == '__main__':
    print("le chat")
    main()