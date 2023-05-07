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
from sys import getsizeof

# TEST AVEC DICO DE DICO POUR OUVRIR CHAQUE FICHIER DINDEXE UNE FOIS - TROP DE RAM
# A voir si on fait pas une découpe tous les 100k reads pour ensuite faire une recherche,
# remplir le dico des résultats, et boucler là dessus jusqu'à la fin
def main2():
    
    fastqFile = "../data/raw_seq_data/ERR020236_1.fastq"    # Fichier fastq
    kmerSize = 31                                           # Taille du k-mer
    indexDir = "../data/index_full_genome/"

    # TEST Lecture du fichier
    """fichier = "../data/index_full_genome/AAACG"
    df = pd.read_table(fichier, header=None)

    print(df)"""

    count = 0
    # TEST PRINCIPAL
    bigDict = {}    # Dictionnaire de dictionnaires
                    # C'est lui qui contient les préfixes
                    # Les sous dico seront : key = readID, value = [suffixes]

    # Lire le fichier fastq
    for record in SeqIO.parse(fastqFile, 'fastq'):
        seq = str(record.seq.upper())
        readId = record.id
    
        # On continue la boucle si la séquence contient un N
        if "N" in seq :
            #print("Il y a un N dans la séqunce")
            continue # Si N, on repart sur une nouvelle boucle sans faire ce qui suit

        #print(readId)
        #print(seq)
        count += 1 # Pour les tests
        print(count)

        #kmerReadDict = {} # Dictionnaire qui contient tous les suffixes

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

            #bigDict[pref] = {}

            if pref in bigDict.keys() :
                try :
                    bigDict[pref][readId].append(suff)
                except KeyError :
                    bigDict[pref][readId] = [suff]
            else :
                bigDict[pref] = {readId : [suff]}
                
            
            """try :
                bigDict[pref][readId].append(suff)
            except KeyError :
                    bigDict[pref] = {readId : [suff]}"""

            #pprint(kmerReadDict)

        #if count == 1000000 : break
    
    print(len(bigDict))
    print(f"Taille : {getsizeof(bigDict)}")
    #pprint(bigDict)

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
            #print("Il y a un N dans la séqunce")
            continue

        print(readId)
        print(seq)

        kmerReadDict = {} # Dictionnaire qui contient tous les suffixes

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

        # Parcours du dictionnaire - TEST AVEC PANDAS
        for pref, suff in kmerReadDict.items():
            indexFile = indexDir + pref
            print(f"Ouverture de l'indexe {indexFile}")
            df = pd.read_table(indexFile, header=None)

            # Test recherche du suffixe
            
            """for s in mes_suff:
                print(f"\tRecherche du suffixe : {s}")
                if s in df[0] :
                    print("yay")
                else : print("booo")"""


            # Recherche du préfixe
            """for s in suff :
                print(f"\tRecherche du suffixe : {s}")
                if s in df[0] :
                    print("yay")
                else : print("booo")"""

            # tests pandas
            #print(pref)

            #my_suff = ["AAAAAAAAAAAAAAAAAAGAAGCTCA", "AAAAAAAAAAAAAAAAGAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAGGTTT", "AAAAAAAAAAAAAAAAAAGGAAGGCG"]

            #print(df[0])

            """for s in suff :
                print(s)
                #print((df[0].eq(s)).any())
                #Autre truc
                #print(df[0].str.contains(s).any())

                # Fonctionne : affiche la ligne trouvée
                result = df[df[0].str.contains(s)]
                #print(result)
                if not result.empty :
                    print(result)
                #print(type(result))
                
                # Marche aussi mais affiche tout
                #result = df[0].str.contains(s)
                #print(result)"""

        # Parcours du dictionnaire des rKmers - méthode dictionnaire
        """for pref, suff in kmerReadDict.items():
            indexFile = indexDir + pref    
            # TEST EN METTANT LE FICHIER EN MEMOIRE EN DICTIONNAIRE
            indexDict = {}
            with open(indexFile, "r") as f:
                for line in f:
                    line = f.readline().strip("\n")
                    print(line)"""


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
    main2()



