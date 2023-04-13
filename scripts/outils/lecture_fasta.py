#!/usr/bin/python3

"""
Test pour vérifier la lecture de séquence
"""

import argparse
from Bio import SeqIO

# Gestion des arguments
parser = argparse.ArgumentParser()

# Création des arguments
parser.add_argument("-i", "--input", dest="fasta_file", help = "fasta input file")
parser.add_argument("-p", "--position", dest="pos", type=int, help="position qu'on veut récupérer")
parser.add_argument("-k", "--kmer-size", dest="ksize", default=31, help="kmer size")

# Récupération des valeurs des arguments et attribution
args = parser.parse_args()
input_file = args.fasta_file
pos = args.pos
kmerSize = args.ksize

# 2. Récupérer la séquence du chromosome en mémoire
def main() :
    seq = []
    with open(input_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = record.seq

    grab_seq = (kmerSize-1)
    snp = seq[pos-1]
    print(f"SNP à la position {pos} : {snp}")
    print(f"kmer : {seq[(pos-1)-grab_seq:(pos-1)+grab_seq]}")



if __name__ == '__main__':
    main()