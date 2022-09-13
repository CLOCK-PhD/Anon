#!/usr/bin/python3

"""
Un programme pour convertir le fichier de sortie "mathching_kmers.tsv" de create_kmer_prefix_dict.py
en fichier fasta afin de pouvoir faire des tests avec

REMARQUES :
Transipedia en ligne ne prend que des fichiers fasta d'environ 1Mo, ce qui représente avec le format actuel environ 20000 kmers.
On peut envisager de faire une sortie d'un dossier contenant les kmers par fichiers de 20k kmers.

"""

import argparse

def main():

    # Gestion des arguments
    parser = argparse.ArgumentParser()

    # Création des arguments
    parser = argparse.ArgumentParser(description='Extract kmers at SNP positions from a reference VCF and fasta file')
    parser.add_argument("-i", "--input_file", dest="matching_kmers_file", help="Path to the file containing the matching kmers")

    args = parser.parse_args()
    matching_kmers_file = args.matching_kmers_file          # fichier matching_kmers.tsv

    #matching_kmers_file = "../data/test_matching_kmers.tsv"

    with open(matching_kmers_file) as f :
        for line in f:
            the_line = line.split("\t")
            kmer_seq = the_line[0]
            rs_id = the_line[1]
            chrom = the_line[2]
            snp_pos = the_line[3]
            kmer_pos = the_line[4].split("\n")[0]
            print(f">{rs_id};chrom={chrom};snp_pos={snp_pos};kmer_pos={kmer_pos};ref=GRCh37p13")
            print(f"{the_line[0]}")


if __name__ == '__main__':
    main()