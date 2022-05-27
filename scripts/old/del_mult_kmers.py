#!/usr/bin/python3

# Supprimer les kmers multiples dans les fichiers kmer.tsv produits par le générateur de kmers

# TROP GOURMAND EN MEMOIRE

from genericpath import isfile
import argparse
import os
import csv
from heapq import merge
from os.path import isfile, join
from pprint import pprint

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", dest="kmer_lists_dir", help = "Directory with the kmers list")
parser.add_argument("-o", "--ouput", dest="output_dir", default="kmer_snp.tsv", help="Output file name")

args = parser.parse_args()
kmer_dir = args.kmer_lists_dir
output_file = args.output_dir

# lister tous les fichiers du dossier
dir_files = [f for f in os.listdir(kmer_dir) if isfile(join(kmer_dir, f))]

kmers_trie = []
file1 = []
kmer_dict = {}

for f in dir_files:
    file_name = f"{kmer_dir}/{f}"
    with open(file_name, newline="") as csvfile:
        tsv_file = csv.reader(csvfile, delimiter="\t")
        # Mettre le fichier courant en mémoire
        for line in tsv_file :
            #file1.append(line)
            kmer_dict[line[0]]=line[1]
    os.remove(file_name)

# Exporter en fichier tsv
output_file_name = f"{kmer_dir}/{output_file}"
with open(output_file_name, "w", encoding="utf-8") as f:
    tsv_header = "k-mer seq\trs_id\n"
    f.write(tsv_header)
    for kmer, id in kmer_dict.items() :
        line_output = f"{kmer}\t{id}\n"
        f.write(line_output)