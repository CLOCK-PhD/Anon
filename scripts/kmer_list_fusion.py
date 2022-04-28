#!/usr/bin/python3

# Programme pour lire et trier les fichiers outputs de gen_kmer_snp.py

from genericpath import isfile
import argparse
import os
import csv
import pandas as pd
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

# VERSION FACILE
kmers_trie = []
file1 = []

for f in dir_files:
    file_name = f"{kmer_dir}/{f}"
    with open(file_name, newline="") as csvfile:
        tsv_file = csv.reader(csvfile, delimiter="\t")
        for line in tsv_file :
            file1.append(line)
        kmers_trie = list(merge(kmers_trie, file1))
        kmers_trie.sort()
    os.remove(file_name)

# Exporter en fichier tsv
output_file_name = f"{kmer_dir}/{output_file}"
with open(output_file_name, "w", encoding="utf-8") as f:
    tsv_header = "k-mer seq\trs_id\n"
    f.write(tsv_header)
    for kmer in kmers_trie :
        line_output = f"{kmer[0]}\t{kmer[1]}\n"
        f.write(line_output)