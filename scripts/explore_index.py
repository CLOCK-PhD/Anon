#!/usr/bin/python3

"""
Effectue une recherche des k-mers d'une séquence donnée dans son index correspondant.
Fournit en sortie les k-mers communs (susceptibles de servir à l'identification).
CETTE VERSION CHARGE TOUT LE FICHIER PREFIXE EN MEMOIRE

Créer un dictionnaire qui sert d'index préfixe pour les k-mers d'une séquence :
    clés : préfixe(str)
    valeur : liste qui contient tous les suffixes associés au préfixe (str)

Effectue la recherche à partir du dictionnaire sur l'index.

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
    7. Recherche des k-mers dans l'index
    8. Écrire le fichier de sortie

    URGENT : Régler le problème de la recherche sur l'index purifié où la recherche est très lente
            Voir si on peut renvoyer le curseur de lecture à un endroit adapté

    ? A FAIRE : Avec le changement de méthode, on perd le comptage ; l'intégrer à nouveau.
        Voir pour faire un dictionnaire suffixe : comptage plutôt qu'un set
        Comparer si on gagne ou perd du temps
"""
import os
import argparse
from typing import OrderedDict
from Bio import SeqIO
from tqdm import tqdm
from pprint import pprint
from os import listdir
from os.path import isfile, join

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

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    index_dir = args.index          # Dossier contenant l'index
    seq_file = args.seq             # Fichier fasta de la séquence
    ksize = args.kmer_size          # Taille du k-mer
    prefix_size = args.prefix_size  # Taille du préfixe

    # 0. Création des variables
    kmers_dict = {}
    matching_kmers_dict = {}
    output_file_name = "matching_kmers_v2.tsv" # Nom du fichier de sortie
    output_file_name = uniquify(output_file_name)

    # À supprimer à la fin des tests
    index_dir = "../data/kmer_snp/chrY_index" # A supprimer à la fin des tests
    seq_file = "../data/grch37p13/NC_000024.9_Homo_sapiens_chromosome_Y.fasta"

    # 1. Charger la séquence en mémoire
    print("Loading sequence...")
    with open(seq_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
    print(f"Sequence length : {len(seq)}")

    # 2. Parcourir les k-mers de la séquence
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
    for key in kmers_dict :
        kmers_dict[key] = list(kmers_dict[key])
        kmers_dict[key].sort()
    print("\tDone.")

    # Parcourir chaque élément de la liste en valeur pour chaque clé
    kmer_count = 0
    for key in kmers_dict:
        for value in kmers_dict[key]:
            kmer_count += 1

    # 7. Recherche des k-mers dans l'index
    kmer_files = [f for f in listdir(index_dir) if isfile(join(index_dir, f))] # Liste des les fichiers
    kmer_files.sort()

    # DEUXIEME TEST - On parcourt chaque ligne et on vérifie si son suffixe est dans la liste du dico
    print("Recherche dans l'index...")
    count = 0
    not_in_index = []
    pbar2 = tqdm(total=len(kmers_dict))
    for key in kmers_dict:
        pbar2.update(1)
        current_file = f"{index_dir}/{key}"
        try :
            with open(current_file, "r") as f:
                lines = f.readlines()
                for l in lines:
                    suff = l.split("\t")[0]
                    if suff in kmers_dict[key]:
                        count += 1
                        matching_kmers_dict[key+suff] = l.split("\t")[1:]           
        except :
            not_in_index.append(key) # Préfixe absent de l'index
    pbar2.close()




    # RECHERCHE AVEC CHARGEMENT DU FICHIER PREFIXE COMPLET
    """print("Recherche dans l'index...")
    not_in_index = []   # Liste des préfixes qui ne sont pas dans l'index
    count = 0           # Comptage du nombre de k-mers contenus dans le dictionnaire
    pbar2 = tqdm(total=len(kmers_dict))
    for key in kmers_dict:
        pbar2.update(1)
        current_file = f"{index_dir}/{key}"
        try :
            with open(current_file, "r") as f:
                for value in kmers_dict[key]:
                    f.seek(reading_position)
                    curr_suffix = value
                    lines = f.readlines()
                    
        except :
            not_in_index.append(key) # Préfixe absent de l'index
    pbar2.close()"""

    # RECHERCHE AVEC LECTURE PAR LIGNE
    """print("Recherche dans l'index...")
    not_in_index = []   # Liste des préfixes qui ne sont pas dans l'index
    count = 0           # Comptage du nombre de k-mers contenus dans le dictionnaire
    pbar2 = tqdm(total=len(kmers_dict))
    for key in kmers_dict:
        pbar2.update(1)
        current_file = f"{index_dir}/{key}"
        reading_position = 0 # Position du curseur pour la lecture du fichier
        try :
            with open(current_file, "r") as f:
                #f.seek(reading_position)
                for value in kmers_dict[key]:
                    f.seek(reading_position)
                    curr_suffix = value
                    line = f.readline()
                    if line.split("\t")[0] == curr_suffix :
                            reading_position = f.tell() # Enregistrer la position et y retourner directement
                            count += 1
                            # Ajouter la ligne dans le dico des k-mers trouvés
                            matching_kmers_dict[key+curr_suffix] = line.split("\t")[1:]
                    while line:
                        line = f.readline()
                        if line.split("\t")[0] == curr_suffix :
                            count += 1
                            reading_position = f.tell()
                            matching_kmers_dict[key+curr_suffix] = line.split("\t")[1:]
                            break

                        elif line.split("\t")[0] > curr_suffix:
                            break # Arrêt de la recherche
        except :
            not_in_index.append(key) # Préfixe absent de l'index
    pbar2.close()"""
    
    print()
    print(f"Taille de l'index :                         {len(kmer_files)}")
    print(f"Nombre de préfixes trouvés :                {len(kmers_dict)}")
    print(f"Nombre de kmers dans le dictionnaire :      {kmer_count}")
    print(f"Nombre de kmers en commun avec l'index :    {count}")
    print(f"Nombre de préfixes absents de l'index :     {len(not_in_index)}")

    # 8. Écrire le dictionnare des k-mers retrouvés dans le fichier ouput :
    print(f"Creating output file : {output_file_name}")
    matching_kmers_dict = OrderedDict(sorted(matching_kmers_dict.items()))
    with open(f"{index_dir}/{output_file_name}", "w") as f:
        for key, value in matching_kmers_dict.items() :
            f.write(f"{key}\t{value[0]}\t{value[1]}\t{value[2]}\t{value[3]}")

if __name__ == '__main__':
    main()