#!/usr/bin/python3

"""
Créer un dictionnaire qui sert d'index préfixe pour les k-mers d'une séquence :
    clés : préfixe(str)
    valeur : liste qui contient tous les suffixes associés au préfixe (str)

Étapes :
    1. Charger la séquence en mémoire
    2. Parcourir chaque k-mer de la séquance
    3. Pour chaque k-mer : découper en préfixe et suffixe
    4. Remplissage du dictionnaire :
        Si le préfixe n'existe pas, l'ajouter nouvelle clé du dictionnaire
        Ajouter le suffixe dans une liste en tant que valeur (en liste)
        Si le préfixe existe, ajouter le suffixe dans la valeur (sans doublon)
    5. Ordonner le dictionnaire
    6. Ordonner chaque liste
    7. Faire des recherches dans l'index 
"""
import sys
from typing import OrderedDict
from Bio import SeqIO
from tqdm import tqdm
from pprint import pprint
from os import listdir
from os.path import isfile, join

def encode(pref:str, prefix_list:list)->int:
    try :
        return prefix_list.index(pref)
    except ValueError :
        return "nope"
    print("le chaton")

def main():

    # 0. Création des variables
    kmers_dict = {}
    ksize = 21
    prefix_size = 5

    # 1. Charger la séquence en mémoire
    seq_file = "../data/grch37p13/NC_000024.9_Homo_sapiens_chromosome_Y.fasta"
    print("Loading sequence...")
    with open(seq_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
    print(f"Sequence length : {len(seq)}")

    # 2. Parcourir les k-mers de la séquence
    # Parcourir les kmers de la séquence
    n_kmers = len(seq) - ksize + 1
    print(f"Number of k-mers to create : {n_kmers}")
    print("Creating prefix dictionnary...")
    pbar = tqdm(total=n_kmers)
    
    # 3. Pour chaque k-mer, découper en préfixe et suffixe (range(n_kmers))
    for i in range(n_kmers):
        pbar.update(1)
        kmer = str(seq[i:i+ksize])
        prefix = kmer[:prefix_size]
        suffix = kmer[prefix_size:]

        # 4. Remplissage du dictionnaire
        try :
            if suffix not in kmers_dict[prefix]:
                kmers_dict[prefix].add(suffix)
        except KeyError:
            kmers_dict[prefix] = {suffix}

    pbar.close()

    # 5. Ordonner les clés du dictionnaire :
    print("Ordering dictionnary...")
    kmers_dict = OrderedDict(sorted(kmers_dict.items()))
    print("\tDone.")
    
    # 6. Ordonner les valeurs du dictionnaire pour chaque clé :
    print("Ordering values...")
    """for value in kmers_dict.values():
        value = list(value)
        value.sort()
        #value = sorted(value)"""
    for key in kmers_dict :
        kmers_dict[key] = list(kmers_dict[key])
        kmers_dict[key].sort()
    print("\tDone.")

    print(len(kmers_dict))

    # Parcourir chaque élément de la liste en valeur pour chaque clé
    kmer_count = 0
    for key in kmers_dict:
        for value in kmers_dict[key]:
            kmer_count += 1

    print(f"Nombre de kmers : {kmer_count}")

    # 7. Faire les recherches dans l'index

    # Dossier contenant les fichiers :
    output_dir = "../data/snp_chr_name_test"
    # Lister tous les fichiers :
    kmer_files = [f for f in listdir(output_dir) if isfile(join(output_dir, f))]
    kmer_files.sort()

    # Ouvrir tous les fichiers
    f_desc = [open(f"{output_dir}/{filename}", "r") for filename in kmer_files]

    for key in kmers_dict:
        id_pref = encode(key, kmer_files)
        if id_pref != "nope":
            for value in kmers_dict[key]:
                suffix = value
                for line in f_desc[id_pref]:
                    if line.split("\t")[0] == suffix:
                        print(f"Trouvé : {key}{value}")



    # Fermer tous les fichiers
    for f in f_desc:
        f.close()

    print("le chat")

if __name__ == '__main__':
    main()