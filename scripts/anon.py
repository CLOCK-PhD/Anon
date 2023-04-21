#!/usr/bin/python3

"""
anon.py

Programme qui va chercher les reads dans l'index des k-mers porteurs de SNP.

Retourne des informations quant aux k-mers contenus dans le read, relatifs aux
risque de ré-identification.

Objectifs :
    . À définir
"""

# Pour l'instant, tests à la con pour les bases :
# OK : Ouvrir un fichier fastq
# OK : Lire le fichier fastq
# OK : Récupérer la séquence du fichier fastq
# OK : Découper la séquence fastq en k-mers
# A FAIRE : rechercher les fastq dans l'index
#       On va essayer deux méthodes :
#       - Charger TOUS les fichiers en mémoire pour la recherche (main_fuckMyRAM())
#       - Ouvrir les fichiers à la demande (main())
# A FAIRE : Compter les occurences
# A FAIRE : Créer un output comme dans le mail d'Alban

# IDEE :  Est-ce que c'est mieux de découper nos kmers dans un dictionnaire
# puis de faire une recherche ? Ou alors rechercher directement ?

import sys
import pandas as pd

from Bio import SeqIO
from pprint import pprint


def main():
    
    fastqFile = "../data/raw_seq_data/ERR020236_1.fastq"    # Fichier fastq
    kmerSize = 31                                           # Taille du k-mer
    indexDir = "../data/index_full_genome/"

    # TEST Lecture du fichier
    """fichier = "../data/index_full_genome/AAACG"
    df = pd.read_table(fichier, header=None)

    print(df)"""


    # TEST PRINCIPAL
    # Lire le fichier fastq
    for record in SeqIO.parse(fastqFile, 'fastq'):
        seq = str(record.seq.upper())
        readId = record.id
    
        # On continue la boucle si la séquence contient un N
        if "N" in seq :
            continue

        print(readId)
        print(seq)

        kmerReadDict = {} # ?

        # Division du read en k-mers
        n_kmers = len(seq) - kmerSize + 1
        for i in range(n_kmers):
            rKmer = seq[i:i+kmerSize]
            #print(rKmer)
            
            # Division du k-mer en préfixe/suffixe
            pref = rKmer[:5]    # Fichier de l'index à ouvrir
            suff = rKmer[5:]    # Suffixe à chercher dans la première colone de l'index
            
            # TEST AVEC UN DICO POUR PAS OUVRIR UN FICHIER A CHAQUE KMER
            # Remplissage du dictionnaire
            try :
                kmerReadDict[pref].append(suff)
            except KeyError :
                kmerReadDict[pref] = [suff]

        # Attention à l'indetation, un \t plus loin et c'est le drame
        #pprint(kmerReadDict)

        # Parcours du dictionnaire
        for pref, suff in kmerReadDict.items():
            indexFile = indexDir + pref
            df = pd.read_table(indexFile, header=None)
            for s in suff :
                """result = df[df[0].str.contains(s)]
                print(result)"""
                result = df[0].str.contains(s)
                print(result)
                if result :
                    print("YAY")


def main_fuckMyRAM():
    print("soon")

    # TEST OUVERTURE TOUS LES FICHIERS
    # mais d'abord un seul
    fichier = "../data/index_full_genome/AAACG"
    fichierDico = {}
    with open(fichier) as f :
        bla = f.readlines()
        pprint(bla)

if __name__ == '__main__':
    main()



