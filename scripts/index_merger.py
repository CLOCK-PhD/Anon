#!/usr/bin/python3

"""
index_merger.py

Programme pour merger tous les index crées par le script ksg.py

Objectif : merger les noms de chaque fichier préfix des différents chromosomes
    en un seul, afin de réunir tous les k-mers porteurs de SNP dans un seul fichier,
    et les contenir tous dans un seul dossier.

Exemple : réunion des fichiers AAAAA du chromosome 1 et 2 en un seul fichier AAAAA
    qui les contient, et pareil pour AAAAC, AAAAT, etc...
"""

# Premier essai : Approche heap merge qui était efficace lors des tests précédents
# pour réunir les fichiers.

# A FAIRE : On va retrouver des k-mers identiques dans des chromosomes différents.
# Il sera nécessaire de les supprimer les deux et de les mettre dans un dictionnaire
# comme lors de la création des indexes chromosomiques.

# PREMIERS TESTS : sur les fichiers output des chromosomes Y et 11

from itertools import product
from pprint import pprint

def main():

    prefixSize = 5
    input_dir_11 = "../data/ksg_test/"


    # Créer les suffixes possibles:
    prefixes = product("ACGT", repeat=prefixSize)
    prefList = []
    for pref in prefixes :
        prefList.append(''.join(pref))

    print(len(prefList))

if __name__ == '__main__':
    main()