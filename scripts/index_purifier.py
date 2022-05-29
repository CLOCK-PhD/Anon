#!/usr/bin/python3

"""
index_purifier.py

Un programme pour éliminer les k-mers des snps qu'on retrouve dans la séquence de référence
grâce au programme snp_kmer_finder.py.

Produit un nouvel index qui ne contient plus les k-mers contenus dans l'output de snp_kmer_finder.py.

EN DEVELOPPEMENT

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

A FAIRE : Gestion des arguments
"""
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

    dir = "../data/snp_chr_name_test"
    output_dir = "../data/snp_chr_name_test/000purified_index" # A modifier avec argparse
    kmer_to_del_file = "../data/snp_chr_name_test/07_kmer_finder_results.tsv"
    kmer_dict = {}
    prefix_size = 5
    new_index_files_list = []
    index_files_list = [f for f in listdir(dir) if isfile(join(dir, f))]
    
    # 1. Création du fichier d'index - On enlève les commentaires quand le reste sera OK
    try :
        os.makedirs(output_dir)
    except FileExistsError:
        output_dir = uniquify(output_dir)
        os.makedirs(output_dir)


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
    print(len(kmer_dict))
    #pprint(kmer_dict)

    # 4. Parcours du fichier préfixe à la recherche des k-mers à supprimer
    for key, value in kmer_dict.items():
        index_file = f"{dir}/{key}"
        new_index_file = f"{output_dir}/{key}"
        with open(index_file, "r") as f:
            # 5. Réécriture du fichier de l'index sans les k-mers multiples
            with open(new_index_file, "w") as nif:
                print(f"{index_file} opened !")
                for line in f:
                    file_suffix = line.split("\t")[0]
                    if file_suffix not in value:
                        nif.write(line)
                # 7. Vérifier les fichiers créés dans le nouvel index
                new_index_files_list.append(key)

    # 8. Lister les fichiers qui manquent
    for i in index_files_list:
        if len(i) != prefix_size or i in new_index_files_list:
            index_files_list.remove(i)

    # 8. Copier les fichiers quand manquent
    for i in index_files_list:
        #os.replace(f"{dir}/{i}", f"{output_dir}/{i}") # Pour déplacer, à faire une fois les tests terminés
        shutil.copy(f"{dir}/{i}", f"{output_dir}/{i}") # Pour copier, pendant les tests


    print("le chat")

if __name__ == "__main__" :
    main()