#!/usr/bin/python3

"""
index_merger.py

Programme pour merger tous les index crées par le script ksg.py

Objectif : merger les noms de chaque fichier préfix des différents chromosomes
    en un seul, afin de réunir tous les k-mers porteurs de SNP dans un seul fichier,
    et les contenir tous dans un seul dossier.

Exemple : réunion des fichiers AAAAA du chromosome 1 et 2 en un seul fichier AAAAA
    qui les contient, et pareil pour AAAAC, AAAAT, etc...

Étapes :
    0. Gestion des arguments
    1. Création de la liste des préfixes
    2. Création du dossier de sortie contenant l'index final
    3. Ouverture des fichiers de préfixes homonymes
    . Merge des fichiers de préfixes homonymes
    . (Optionnel) Suppression des fichiers d'index initiaux
    . Extraction des k-mers identiques retrouvés dans les fichiers préfixes homonymes
    . Enjoie
"""

# Premier essai : Approche heap merge qui était efficace lors des tests précédents
# pour réunir les fichiers.

# A FAIRE : On va retrouver des k-mers identiques dans des chromosomes différents.
# Il sera nécessaire de les supprimer les deux et de les mettre dans un dictionnaire
# comme lors de la création des indexes chromosomiques.
#   Idée : soit à faire au vol, ce qui risque d'être un peu chiant
#           Soit à faire une fois le merge effectué (hmmmm...)

# A FAIRE : Gestion des arguments une fois les tests réussis

# OK : PREMIERS TESTS : sur les fichiers output des chromosomes Y, 5 et 11
# A FAIRE :  Créer les fichiers de préfixe de sortie contenant les kmers uniques
# A FAIRE : Gérer plusieurs préfixes à la suite


import os
import glob
import heapq
import argparse
import sys

import pandas as pd

from itertools import product
from pprint import pprint

# Générer un nom de dossier unique
def uniquify(path:str) -> str:
    """
    Génère un nom de dossier unique pour FileExistsError

    Parameters :
        path    (str):  Nom du chemin du dossier à créer

    Returns :
        path    (str):  Nouveau nom du chemin du dossier
    """
    filename, extension = os.path.splitext(path)
    counter = 1

    while os.path.exists(path):
        path = filename + "_" + str(counter) + extension
        counter += 1

    return path

def main():

    # On le déplace parce qu'il faut en faire un par préfixe
    #indexDict = {}  # Dictionnaire de tous les fichiers réunis
    #dupDict = {}    # Dictionnaire sans doublon
    prefixSize = 5
    inputDir = "../data/ksg_test/"  # Dossier contenant les indexes pour chaque chromosome
    input_dir_Y = "../data/ksg_test/kmer_snp_index_chrY/"
    input_dir_11 = "../data/ksg_test/kmer_snp_index_chr11/"
    input_dir_5 = "../data/ksg_test/kmer_snp_index_chr5/"


    # 1. Créer les suffixes possibles:
    prefixes = product("ACGT", repeat=prefixSize)
    prefList = []
    for pref in prefixes :
        prefList.append(''.join(pref))

    print(f"Nombre de prédixes créés : {len(prefList)}")

    # 2. Création du dossier de sortie contenant l'index final - OK
    """outputDirectory = "../data/ksg_test/index_full_genome"
    try :
        os.makedirs(outputDirectory)
    except FileExistsError:
        outputDirectory = uniquify(outputDirectory)
        os.makedirs(outputDirectory)
    print(f"Dossier {outputDirectory} créé.")"""


    # 3. Ouverture des fichiers préfixes homonymes

    # Afficher le contenu des dossiers
    indexFileNames = [input_dir_Y, input_dir_11, input_dir_5]
    #pprint(indexFileNames)
    my_files = []

    # TEST EN COURS : Récupérer tous les noms des dossiers des indexes de chromosomes


    # OK - TEST : Ouvrir tous les préfixes à la suite pour faire les opérations
    for p in prefList :
        print(f"\nPréfixe en cours : {p}")
        indexDict = {}  # Dictionnaire de tous les fichiers réunis
        dupDict = {}    # Dictionnaire sans doublon
        les_fichiers = []
        for f in indexFileNames :
            les_fichiers.append(f + p)
    
        #pprint(les_fichiers)
        
        # Liste avec tous les noms de fichiers d'indexe homonymes pour tous les chrom
        # La suite des opérations se passe donc ici

        for f in les_fichiers :
            with open(f, "r") as le_index:
                print(f"Reading file {f}")
                # test en virant le range()
                #for n in range(1000):
                for l in le_index:
                    line = le_index.readline().strip("\n").split("\t")
                    # line[0] = suffixe; line[1:] : tout le reste des infos
                    #print(line[0])
                    #print(f"\t{line[1:]}")

                    # Remplir le dictionnaire - premier test
                    try :
                        indexDict[line[0]].append(line[1:])
                    except KeyError :
                        indexDict[line[0]] = [line[1:]]

        # Afficher le dictionnaire crée
        print(f"\tDictionnaire initial : {len(indexDict)}")
        #pprint(indexDict)

        # Parcours du dictionnaire pour liste les éléments à supprimer :
        kmerToDel = []
        for suffix, val in indexDict.items() :
            # un k-mer est unique ssi il n'a qu'un variant dans sa liste
            # donc la liste de ses valeur ne contient qu'un élément (qui est une liste)
            if len(val) > 1 :
                #print(f"{suffix}\n\t{val}")
                kmerToDel.append(suffix)

        # Afficher les k-mers à supprimer
        #pprint(kmerToDel)

        # Déplacer les k-mers à supprimer dans l'autre dictionnaire
        for suffix in kmerToDel:
            dupDict[suffix] = indexDict.pop(suffix)
            
        # Afficher le dictionnaire des dups
        print(f"\tDictionnaire des dups : {len(dupDict)}")
        #pprint(dupDict)
        # Afficher le dictionnaire index final
        print(f"\tDictionnaire final : {len(indexDict)}")
        #pprint(indexDict)





    # Premiers tests - ici ça marche pour juste un nom de fichier donné
    # Ajouter le nom d'un préfixe aux dossiers pour ouvrir les bons fichiers
    """for f in indexFileNames :
        my_files.append(f + "AAAAA")"""
    
    # 4. 
    ### TEST : Plus simple, méthode d'origine avec un dictionnaire

    #pprint(my_files)   # Afficher tous les noms de fichiers créés

    """for f in my_files :
        with open(f, "r") as le_index :
            print(f"Reading file {f}")
            for n in range(5):
                line = le_index.readline().strip("\n").split("\t")
                # line[0] = suffixe; line[1:] : tout le reste des infos
                #print(line[0])
                #print(f"\t{line[1:]}")

                # Remplir le dictionnaire - premier test
                try :
                    indexDict[line[0]].append(line[1:])
                except KeyError :
                    indexDict[line[0]] = [line[1:]]"""

    # Afficher le dictionnaire crée
    """print("\nDictionnaire initial :")
    #pprint(indexDict)
    print(len(indexDict))"""
    
    # Parcours du dictionnaire pour liste les éléments à supprimer :
    """kmerToDel = []
    for suffix, val in indexDict.items() :
        # un k-mer est unique ssi il n'a qu'un variant dans sa liste
        if len(val) > 1 :
            #print(f"{suffix}\n\t{val}")
            kmerToDel.append(suffix)"""

    # Afficher les k-mers à supprimer
    #pprint(kmerToDel)

    # Déplacer les k-mers à supprimer dans l'autre dictionnaire
    """for suffix in kmerToDel:
        dupDict[suffix] = indexDict.pop(suffix)"""
    

if __name__ == '__main__':
    main()