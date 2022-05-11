#!/usr/bin/python3

"""
Programme pour merge les tsv avec suppression des doublons
"""

import heapq
import argparse
import sys

# A FAIRE : préciser un fichier de sortie pour éviter de tout supprimer
# A FAIRE : ajouter une option pour supprimer tous les fichiers après le tri

def main():
    parser = argparse.ArgumentParser(description='Merge multiple sorted files')
    parser.add_argument('files', metavar='files', type=str, nargs='+', help='input files')

    args = parser.parse_args()

    files = [open(filename, 'r') for filename in args.files]

    """
    # Sans suppression
    merged = heapq.merge(*files)
    for line in merged:
        sys.stdout.write(line)"""

    # Suppression des doublons
    header = "#Kmer_seq\trs_id\tchromosome\tsnp_position\tkmer_position\n"
    sys.stdout.write(header)
    merged = heapq.merge(*files)
    prev_line = ""
    for line in merged:
        kmer = line.split("\t")[0]
        if kmer == prev_line :
            prev_line = kmer
        else:
            sys.stdout.write(line)
            prev_line = kmer

    # Suppression des fichiers des kmers :


if __name__ == '__main__':
    main()