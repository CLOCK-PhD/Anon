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
    4. Merge des fichiers de préfixes homonymes
    . (Optionnel) Suppression des fichiers d'index initiaux
    5. Extraction des k-mers identiques retrouvés dans les fichiers préfixes homonymes
    6. Écriture du nouvel index dans un dossier
    7. Enjoie
"""




####################################################################
# URGENT : TRIER LES DICTIONNAIRES AVANT ECRITURE DES FICHIERS !!! #
####################################################################




# A FAIRE : Gestion des arguments une fois les tests réussis

import os
from glob import glob

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

    prefixSize = 5  # Taille du préfixe
    inputDir = "../data/ksg_test/"  # Dossier contenant les indexes pour chaque chromosome
    # Attention de ne pas mettre le nouvel index dedans.

    # 1. Créer les suffixes possibles
    prefixes = product("ACGT", repeat=prefixSize)
    prefList = [] # Liste qui va contenir tous les préfixes possibles
    for pref in prefixes :
        prefList.append(''.join(pref))
    #print(f"Nombre de prédixes créés : {len(prefList)}")

    # 2. Création du dossier de sortie contenant l'index final
    outputDirectory = "../data/index_full_genome"
    try :
        os.makedirs(outputDirectory)
    except FileExistsError:
        outputDirectory = uniquify(outputDirectory)
        os.makedirs(outputDirectory)
    print(f"Dossier {outputDirectory} créé.")


    # 3. Ouverture des fichiers préfixes homonymes

    # Liste contenant tous les noms de dossiers d'indexes
    indexDirList = glob(inputDir + "*/")
    #pprint(indexDirList)

    # Ouvrir tous les préfixes à la suite pour faire les opérations
    for p in prefList :
        print(f"\nPréfixe en cours : {p}")
        indexDict = {}  # Dictionnaire de tous les fichiers réunis
        dupDict = {}    # Dictionnaire sans doublon
        les_fichiers = [] # Contient les fichiers homonymes des différents indexes
        for f in indexDirList :
            les_fichiers.append(f + p)
        # Afficher les fichiers utilisés :
        #pprint(les_fichiers)
        
        # Liste avec tous les noms de fichiers d'indexe homonymes pour tous les chrom
        # La suite des opérations se passe donc ici

        for f in les_fichiers :
            with open(f, "r") as currentIndex:
                #print(f"Reading file {f}")
                # 4. 
                for l in currentIndex:
                    line = currentIndex.readline().strip("\n").split("\t")
                    # Remplir le dictionnaire - premier test
                    try :
                        indexDict[line[0]].append(line[1:])
                    except KeyError :
                        indexDict[line[0]] = [line[1:]]

        # Afficher le dictionnaire crée
        print(f"\tDictionnaire initial : {len(indexDict)}")
        #pprint(indexDict)

        # 5. Extraction des k-mers identiques retrouvés dans les fichiers préfixes homonymes 
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

        # 6. Écriture du nouvel index dans un dossier
        with open(f"{outputDirectory}/{p}", "a") as f:
            for suff, var in indexDict.items() :
                the_var = "\t".join(var[0])
                line = f"{suff}\t{the_var}\n"
                f.write(line)

#7. Enjoie
if __name__ == '__main__':
    main()