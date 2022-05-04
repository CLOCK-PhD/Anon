#!/usr/bin/python3

# Programme pour récupérer les séquences fasta et les mettre dans des fichiers séparés

from hashlib import new
import argparse
from Bio import SeqIO

# Gestion des arguments
parser = argparse.ArgumentParser()

# Création des arguments
parser.add_argument("-i", "--input", dest="fasta_file", help = "fasta input file")

# Récupération des valeurs des arguments et attribution
args = parser.parse_args()
input_file = args.fasta_file

with open(input_file) as handle :
    for record in SeqIO.parse(handle, "fasta"):
        if "NC" in record.description: # Pour récupérer les séquences qui commencent par NC
            output_name = record.description.split(",")[0]
            output_name = output_name.replace(" ", "_")
            print(output_name)
            with open(f"{output_name}.fasta", "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")