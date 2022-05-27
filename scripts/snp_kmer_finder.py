#!/usr/bin/python3

"""
snp_kmer_finder : un programme pour rechercher les kmers des snp présents dans une séquence.

EN COURS DE DEV

Objectif :
1. Ouvrir la séquence de référence
2. La parcourir pour trouver chaque k-mer
3. Pour chaque k-mer, rechercher s'il est présent dans l'index
4. Relever les k-mers présents et les compter. (vérifier la position dans la séquence)
5. Faire un fichier output de comptage de k-mers

A FAIRE : Intégrer une rechercher dichotomique (binary search) pour la recherche de suffixe

IDEE : Actualiser l'indexe :
    1. Faire une recherche de k-mers sur la séquence de référence pour supprimer les kmers qui
    apparaissent plusieurs fois.
    2. Réutiliser le programme avec la séquence d'une autre personne.
"""

import ahocorasick
import re
import sys
import argparse
import os
import heapq
import linecache
from Bio import SeqIO
from typing import OrderedDict
from pprint import pprint
from os.path import isfile, join

def main():
    # Variables
    ksize = 21
    prefix_size = 5
    index_dir = "../data/snp_chr_name_test"
    found_kmer_dict = {}

    # Charger la séquence
    seq_file = "../data/grch37p13/NC_000024.9_Homo_sapiens_chromosome_Y.fasta"
    with open(seq_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
    
    # Parcourir les kmers de la séquence
    n_kmers = len(seq) - ksize + 1
    count = 0
    for i in range(n_kmers):
        count += 1
        kmer = str(seq[i:i+ksize])
        prefix = kmer[:prefix_size]
        suffix = kmer[prefix_size:]
        try :
            with open(f"{index_dir}/{prefix}", "r")as f:
                for line in f: # On changera ça pour une recherche dichotomique après
                    if suffix == line.split("\t")[0] :
                        try:
                            found_kmer_dict[prefix+line.split("\t")[0]].append(line.split("\t")[1])
                        except KeyError:
                            found_kmer_dict[prefix+line.split("\t")[0]] = [line.split("\t")[1]]
        except FileNotFoundError:
            print(f"{prefix} is not a prefix file")
        if count == 100000:
            break

    #pprint(found_kmer_dict)
    for key, value in found_kmer_dict.items() :
        print(f"{key} : {len(value)}")


if __name__ == "__main__":
    main()