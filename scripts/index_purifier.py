#!/usr/bin/python3

"""
index_purifier.py

Programme pour éliminer les k-mers des snps qu'on retrouve dans la séquence de référence
grâce au programme snp_kmer_finder.py.

Produit un nouvel index qui ne contient plus les k-mers contenus dans l'output de snp_kmer_finder.py.

Objectif :
    1. Créer un dossier pour le nouvel index (OK)
    2. Lire le fichier output kmer_snp_finder.py (OK)
    3. Charger le fichier en mémoire dans un dictionnaire ordonné : Clé = préfixe, Valeur = [suffixes] (OK)
        - Ordonner les listes en Valeur (OK)
    4. Pour chaque clé, parcourir le fichier préfixe de l'index qui correspond (OK)
    5. Réécrire chaque ligne de l'index dont le suffixe n'est pas contenu dans le dictionnaire (OK)
    6. Supprimer le fichier d'index correspondant (en attente de la fin des tests)
    7. Vérifier les fichiers créés dans le nouvel index (OK)
    8. Copier les fichiers qui manquent qui n'ont pas été affectés par la purification. (OK)
"""

import argparse
from collections import OrderedDict
import os
from pprint import pprint
from os import listdir
from os.path import isfile, join
import shutil


# Générer un nom de dossier unique
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
    parser = argparse.ArgumentParser(description='Purify the index by making a new one without k-mers found on the reference sequence')
    parser.add_argument("-d", "--directory", dest="directory", help="Path to the directory containing the index")
    parser.add_argument("-kf", "--kmer-file", dest="kmer_file", help="Path to the sequence fasta file")
    parser.add_argument("-p", "--prefix_size", dest="prefix_size", default=5, type=int, help="Prefix size")

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    dir = args.directory
    kmer_to_del_file = args.kmer_file
    prefix_size = args.prefix_size

    output_dir = f"{dir}/purified_index"
    kmer_dict = {}
    new_index_files_list = []
    index_files_list = [f for f in listdir(dir) if isfile(join(dir, f))]

    # Variables utilisées lors des tests
    #dir = "../data/snp_chr_name_test"
    #kmer_to_del_file = "../data/snp_chr_name_test/07_kmer_finder_results.tsv"
    #prefix_size = 5
    
    # 1. Création du fichier d'index - On enlève les commentaires quand le reste sera OK
    try :
        os.makedirs(output_dir)
    except FileExistsError:
        output_dir = uniquify(output_dir)
        os.makedirs(output_dir)
    print(f"Output directory name : {output_dir}")

    # 2. Ouverture du fichier kmer_found
    with open(kmer_to_del_file, "r") as f:
        for line in f:
            kmer = line.split("\t")[0]
            prefix = kmer[:prefix_size]
            suffix = kmer[prefix_size:]

            # 3. Remplir le dictionnaire
            try:
                kmer_dict[prefix].append(suffix)
            except KeyError:
                kmer_dict[prefix] = [suffix]
    
    # Tri du dico (pour être sûr)
    kmer_dict = OrderedDict(sorted(kmer_dict.items()))

    # Tri des listes en valeur du dict
    for value in kmer_dict.values():
        value.sort()

    # 4. Parcours du fichier préfixe à la recherche des k-mers à supprimer
    for key, value in kmer_dict.items():
        index_file = f"{dir}/{key}"
        new_index_file = f"{output_dir}/{key}"
        with open(index_file, "r") as f:
            # 5. Réécriture du fichier de l'index sans les k-mers multiples
            with open(new_index_file, "w") as nif:
                #print(f"{index_file} opened !")
                for line in f:
                    file_suffix = line.split("\t")[0]
                    if file_suffix not in value:
                        nif.write(line)

    # 7. Vérifier les fichiers crées dans le nouvel index
    new_index_files_list = [f for f in listdir(output_dir) if isfile(join(output_dir, f))]
    new_index_files_list.sort()
    index_files_list.sort()

    # 8a. Lister les fichiers qui manquent
    for i in index_files_list: # A virer quand on aura un dir spécifique pour l'index original
        if len(i) != prefix_size :
            print(f"{i} removed from list (prefix size)")
            index_files_list.remove(i)
    # Liste des fichiers à copier
    s = set(new_index_files_list)
    files_to_copy = [x for x in index_files_list if x not in s]

    # 8b. Copier les fichiers quand manquent
    for i in files_to_copy:
        #os.replace(f"{dir}/{i}", f"{output_dir}/{i}") # Pour déplacer, à faire une fois les tests terminés
        shutil.copy(f"{dir}/{i}", f"{output_dir}/{i}") # Pour copier, pendant les tests

    print("le chat")


if __name__ == "__main__" :
    main()