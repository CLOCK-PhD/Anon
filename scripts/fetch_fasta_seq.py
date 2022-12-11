#!/usr/bin/python3

# Programme pour récupérer les séquences fasta et les mettre dans des fichiers séparés

from hashlib import new
import argparse
from Bio import SeqIO
import re

# Conversion du nom du chromosome
def convertChromName(c:str)->str:
    res = re.search("NC_00+([0-9]{1,2}).*$", c)
    if res :
        chrom = res.group(1)
        if res.group(1) == "23":
            return "X"
        elif res.group(1) == "24":
            return "Y"
        else :
            return res.group(1)
    else :
        return c

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
            output_name = convertChromName(record.description)
            #output_name = output_name.replace(" ", "_")
            print(output_name)
            with open(f"{output_name}.fasta", "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")