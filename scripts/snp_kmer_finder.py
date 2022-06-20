#!/usr/bin/python3

"""
snp_kmer_finder : rechercher les kmers des snp présents dans une séquence à partir d'un index.

Prend en entrée :
    une séquence au format fasta (option -f)
    l'index de k-mer fournit par le programme kmer_snp_gen_index.py (option -f)
Argument optionnels :
    taille du k-mer (-k)
    taille du préfixe (-p)
    --no-prefix : fichier de sortie sous forme de fichier unique

Fournit en sortie :
    - un fichier .tsv trié contenant:
        - les k-mers
        - le rs_id de leur snp
        - le nombre de fois qu'ils ont été trouvés
    - Optionnel : un fichier contenant les préfixes qui ne sont pas dans l'indexe

Le fichier de sortie .tsv est sous la forme :
kmer_seq    rs_id   comptage
"""

"""
AMELIORATION : Algo de recherche plus efficace que la simple recherche par ligne actuelle
    -Intégrer une rechercher dichotomique (binary search) pour la recherche de suffixe
        Testé, pas approuvé. À approfondir.
"""

import argparse
import os
from tqdm import tqdm
from Bio import SeqIO
from typing import OrderedDict
from pprint import pprint

def uniquify(path:str) -> str:
    """
    Génère un nom de dossier unique pour FileExistsError

    Parameters :
        path (str): Nom du chemin du dossier à créer

    Returns :
        path(str):  Nouveau nom du chemin du dossier
    """
    filename, extension = os.path.splitext(path)
    counter = 1

    while os.path.exists(path):
        path = filename + "_" + str(counter) + extension
        counter += 1

    return path

def main():

    # Gestion des arguments
    parser = argparse.ArgumentParser()

    # Création des arguments
    parser = argparse.ArgumentParser(description='Extract kmers at SNP positions from a reference VCF and fasta file')
    parser.add_argument("-i", "--index", dest="index", help="Path to the directory containing the index")
    parser.add_argument("-f", "--fasta", dest="seq", help="Path to the sequence fasta file")
    parser.add_argument("-k", "--kmer_size", dest="kmer_size", default=21, type=int, help="Select k-mer size")
    parser.add_argument("-p", "--prefix_size", dest="prefix_size", default=5, type=int, help="Prefix size")
    parser.add_argument('--show-prefix', dest="no_pref" ,action='store_true', help="Add this option to create a simple output file in lexicographic order instead of a prefix index")

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    index_dir = args.index
    seq_file = args.seq
    ksize = args.kmer_size
    prefix_size = args.prefix_size
    output_unknown_prefix = args.no_pref
    output_file_name = "kmer_finder_results.tsv"
    uniquify(output_file_name)

    found_kmer_dict = {}
    not_in_index = []

    print("Loading sequence")
    with open(seq_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
    print(f"Sequence length : {len(seq)}")
    
    # Parcourir les kmers de la séquence
    n_kmers = len(seq) - ksize + 1
    print(f"Number of k-mers to analyse : {n_kmers}")
    print("Searching sequence k-mers in the index...")
    pbar = tqdm(total=n_kmers)
    for i in range(n_kmers):
        pbar.update(1)
        kmer = str(seq[i:i+ksize])
        prefix = kmer[:prefix_size]
        suffix = kmer[prefix_size:]
        try :
            with open(f"{index_dir}/{prefix}", "r")as f:
                for line in f:
                    if suffix == line.split("\t")[0] :
                        try:
                            found_kmer_dict[(prefix+line.split("\t")[0], line.split("\t")[1])].append(i)
                        except KeyError:
                            found_kmer_dict[(prefix+line.split("\t")[0], line.split("\t")[1])] = [i]
                        break
        except FileNotFoundError:
            not_in_index.append(prefix)

    pbar.close()
    
    # Tri du dictionnaire
    found_kmer_dict = OrderedDict(sorted(found_kmer_dict.items()))

    # Afficher le k-mer et le nombre de fois qu'il apparait dans la séquence
    with open(f"{index_dir}/{output_file_name}", "w") as f:
        for key, value in found_kmer_dict.items() :
            f.write(f"{key[0]}\t{key[1]}\t{len(value)}\n")

    if output_unknown_prefix:
        pprint(not_in_index)

if __name__ == "__main__":
    main()