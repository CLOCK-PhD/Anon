#!/usr/bin/python3

"""
snp_kmer_finder : rechercher les kmers des snp présents dans une séquence à partir d'un index.

EN COURS DE DEV

Objectif :
1. Ouvrir la séquence de référence
2. La parcourir pour trouver chaque k-mer
3. Pour chaque k-mer, rechercher s'il est présent dans l'index
4. Relever les k-mers présents et les compter. (vérifier la position dans la séquence)
5. Faire un fichier output de comptage de k-mers

A FAIRE : Intégrer une rechercher dichotomique (binary search) pour la recherche de suffixe
A FAIRE : Retirer les conditions de test

IDEE : Actualiser l'indexe :
    1. Faire une recherche de k-mers sur la séquence de référence pour supprimer les kmers qui
    apparaissent plusieurs fois.
    2. Réutiliser le programme avec la séquence d'une autre personne.
"""

from tqdm import tqdm
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
    not_in_index = []

    # Charger la séquence
    print("Loading sequence")
    seq_file = "../data/grch37p13/NC_000024.9_Homo_sapiens_chromosome_Y.fasta"
    with open(seq_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
    
    print("Analysing sequence k-mers...")
    # Parcourir les kmers de la séquence
    n_kmers = len(seq) - ksize + 1
    count = 0
    pbar = tqdm(total=100000)
    for i in range(n_kmers):
        pbar.update(1)
        count += 1
        kmer = str(seq[i:i+ksize])
        prefix = kmer[:prefix_size]
        suffix = kmer[prefix_size:]
        try :
            with open(f"{index_dir}/{prefix}", "r")as f:
                for line in f: # On changera ça pour une recherche dichotomique après
                    if suffix == line.split("\t")[0] :
                        # Key = k-mer, Value = [rs_id]
                        """try:
                            found_kmer_dict[prefix+line.split("\t")[0]].append(line.split("\t")[1])
                        except KeyError:
                            found_kmer_dict[prefix+line.split("\t")[0]] = [line.split("\t")[1]]"""
                        # Key = (k-mer, rs_id), Value = comptage
                        """try:
                            found_kmer_dict[(prefix+line.split("\t")[0], line.split("\t")[1])] += 1
                        except KeyError:
                            found_kmer_dict[(prefix+line.split("\t")[0], line.split("\t")[1])] = 1"""
                        # Key = (k-mer, rs_id), Value = [positions]
                        try:
                            found_kmer_dict[(prefix+line.split("\t")[0], line.split("\t")[1])].append(i)
                        except KeyError:
                            found_kmer_dict[(prefix+line.split("\t")[0], line.split("\t")[1])] = [i]
        except FileNotFoundError:
            not_in_index.append(prefix)
        if count == 100000:
            break

    pbar.close()

    pprint(not_in_index)
    # Afficher le k-mer et le nombre de fois qu'il apparait dans la séquence
    """for key, value in found_kmer_dict.items() :
        print(f"{key} : {value}")"""


if __name__ == "__main__":
    main()