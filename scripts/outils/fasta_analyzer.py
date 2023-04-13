#!/usr/bin/python3

# Extraire les Primary assemblies des chromosomes
# Données d'entrées : grch38p7

from Bio import SeqIO
from pprint import pprint

selected_sequences = []

with open("grch38p7/grch38p7.fasta") as handle :
    for record in SeqIO.parse(handle, "fasta"):
        #print(record.description)
        seq_desc = record.description
        x = seq_desc.split()
        #pprint(x)
        if x[0].startswith("NC") and len(record) > 30000:
            print(record.description)
            #print(len(record))
            selected_sequences.append(record)

SeqIO.write(selected_sequences, "grch38p7_prim_assembly.fasta", "fasta")