#!/usr/bin/python3

# Extraire des kmers depuis la séquence de référence du chromosome directement, à partir des infos vcf ref snp

import re
import sys
import argparse
import os
from Bio import SeqIO
from numpy import equal
from matplotlib.pyplot import close
from pprint import pprint

# Gestion des arguments
parser = argparse.ArgumentParser()

# Création des arguments
parser.add_argument("-i", "--input", dest="fasta_file", help = "fasta input file")
parser.add_argument("-k", "--kmer_size", dest="kmer_size", default=21, help="Select k-mer size")
parser.add_argument("-nt","-- nucleotide_postition", default=1, dest="nt_pos", help="Position du nucléotide à partir duquel générer le kmer")

# Récupération des valeurs des arguments et attribution
args = parser.parse_args()
input_file = args.fasta_file
kmer_size = int(args.kmer_size)
nt_pos = int(args.nt_pos) - 1

#print(kmer_size//2)

# Récupérer la séquence
seq = []
with open(input_file) as handle :
    for record in SeqIO.parse(handle, "fasta"):
        seq = record.seq
        print(record.description)
        #print(record.seq[nt_pos])
        #print(record.seq[nt_pos - kmer_size//2 : nt_pos + kmer_size//2])

snp_pos_list = [10001, 10002, 10003, 10007, 10008, 10009, 10013, 10028, 10040, 10041, 10042]

#print(f"Position du SNP : {seq[nt_pos]}")
#print(f"kmer : {seq[nt_pos - kmer_size//2 : nt_pos + kmer_size//2]}")
print(f"Longueur de la séquence : {len(seq)}")

# Récupère le kmer à la position indiquée dans la séquence
for i in snp_pos_list:
    snp_pos = i -1
    print(f"Position du snp : {i}")
    print(f"SNP : {seq[snp_pos]}")
    print(f"kmer : {seq[snp_pos - kmer_size//2 : snp_pos + kmer_size//2 + 1]}")