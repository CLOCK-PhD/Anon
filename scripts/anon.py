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
# OK : A FAIRE : rechercher les fastq dans l'index
#       On va essayer deux méthodes :
#       - Charger TOUS les fichiers en mémoire pour la recherche (main3())
#       - Ouvrir les fichiers à la demande (main())
# A FAIRE : Compter les occurences
# A FAIRE : Créer un output comme dans le mail d'Alban

# IDEE :  Est-ce que c'est mieux de découper nos kmers dans un dictionnaire
# puis de faire une recherche ? Ou alors rechercher directement ?

import sys
import glob
import os
import pandas as pd
import json
import bisect

from Bio import SeqIO
from pprint import pprint
from sys import getsizeof
from tqdm import tqdm


def dicho(wordList, word):

    # Test pour adapter à mon input pourri

    if len(wordList) == 0:
        return False
    
    i = len(wordList) // 2
    if wordList[i].split("\t")[0] == word:
        #return True
        return wordList[i].split("\t")[1]
    
    if wordList[i].split("\t")[0] > word:
        # recherche dans la première moitié
        return dicho(wordList[:i], word)
    else :
        # recherche dans la seconde moité
        return dicho(wordList[i+1:], word)
    
def dicho2(wordList, word):
    i = bisect.bisect_left(wordList, word)
    if i == len(wordList):
        return False
    return wordList[i] == word

"""# TEST AVEC DICO DE DICO POUR OUVRIR CHAQUE FICHIER DINDEXE UNE FOIS - TROP DE RAM
# A voir si on fait pas une découpe tous les 100k reads pour ensuite faire une recherche,
# remplir le dico des résultats, et boucler là dessus jusqu'à la fin"""
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

    resDict = {}    # Dictionnaire contenant les résultats
                    # clé = read
                    # valeur = les SNP contenus

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
        #print(count)

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

        if count == 100 :
            break
    
    # BON ICI CA PEUT MARCHER
    # MAIS PANDAS
    # ON VA FAIRE LE TEST AVEC TOUT L'INDEX EN MEMOIRE
    kmer_count = 0
    for prefixe, read in bigDict.items() :
        #print(prefixe)
        indexFile = indexDir + prefixe
        print(f"Ouverture de l'indexe {indexFile}")
        df = pd.read_table(indexFile, header=None)
        for readName, suffixe in read.items() :
            print(f"\t{readName}")
            for e in suffixe :
                print(f"\t\t{e}")
                kmer_count += 1


    print(len(bigDict))
    print(f"Taille : {getsizeof(bigDict)}")
    print(f"Nombre de kmers générés : {kmer_count}")
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


def main3():
    
    # A FAIRE : tqdm

    fastqFile = "../data/raw_seq_data/ERR020236_1.fastq"    # Fichier fastq
    kmerSize = 31                                           # Taille du k-mer
    indexDir = "../data/index_full_genome/"

    inMemoryFilesDict = {} # Dictionnaire contenant l'indexe en mémoire
    resultDict = {}         # Dictionnaire des reads ayant des k-mers porteurs de SNP

    # Liste de tous les fichiers
    indexFiles = sorted(os.listdir(indexDir))
    #print(indexFiles)

    count = 0 # Pour les tests
    pbar = tqdm(total=len(indexFiles))
    for p in indexFiles :
        count += 1
        pbar.update(1)
        inMemoryFilesDict[p] = []
        indexFilePath = indexDir + p
        with open(indexFilePath, "r") as f:
            #print(f"Opening {indexFilePath}")
            for line in f :
                line = f.readline().strip("\n")
                #print(line) # Jusque là, on ouvre et on lit
                # Mise en mémoire du fichier dans un dictionnaire
                inMemoryFilesDict[p].append(line)
        #if count == 2 : break
    pbar.close()
    
    #print(inMemoryFilesDict.keys())
    #pprint(inMemoryFilesDict) # jusque là on dirait qu'on a tous ce qu'il nous faut

    # Rappel :
    # Dictionnaire avec clé = préfixe, valeur = liste de chaque ligne
    # ligne 0 : suffixe, 1:rsid, 2:?, 3: SNP, 4:snppos, 5:kmerpos, 6 et 7 osef

    # En vrai ce serait mieux avec pandas hein...
    # ----------------------------------------------
    # Lire le fichier fastq

    # Comptage de lignes
    print("Comptage du nombre de reads...")
    with open(fastqFile, "rb") as f:
        num_lines = int(sum(1 for _ in f)/4)

    print(f"Il y a {num_lines} reads dans le fichier fastq")

    #searchCount = 0 # Pour les tests
    pbar2 = tqdm(total=num_lines)
    for record in SeqIO.parse(fastqFile, 'fastq'):
        pbar2.update(1)
        seq = str(record.seq.upper())
        readId = record.id
    
        # On continue la boucle si la séquence contient un N
        if "N" in seq :
            #print("Il y a un N dans la séqunce")
            continue
        
        #searchCount += 1 # Pour les tests
        #print(f"{searchCount} - {readId}")
        #print(seq)

        kmerReadDict = {} # Dictionnaire qui contient tous les suffixes du read en cours

        # Division du read en k-mers
        n_kmers = len(seq) - kmerSize + 1
        for i in range(n_kmers):
            rKmer = seq[i:i+kmerSize]
            #print(rKmer)
            
            # Division du k-mer en préfixe/suffixe
            pref = rKmer[:5]    # Fichier de l'index à ouvrir
            suff = rKmer[5:]    # Suffixe à chercher dans la première colone de l'index

            # LENTEUR EXTREME - Tentative sans passer par un dico
            """for e in inMemoryFilesDict[pref]:
                if suff == e.split("\t")[0] :
                    #print("yay")
                    try :
                        resultDict[readId].append(e.split("\t")[1])
                    except KeyError :
                        resultDict[readId] = [e.split("\t")[1]]"""
            
            # ESSAI - RECHERCHE DICHOTOMIQUE
            res = dicho(inMemoryFilesDict[pref], suff)
            if res :
                try :
                    resultDict[readId].append(res)
                except KeyError :
                    resultDict[readId] = [res]

        #if searchCount == 100:break
    
    #pprint(resultDict)

    # Affichage des résultats :
    
    for read, snps, in resultDict.items() :
        print(f"Le read {read} est associé aux SNPs : ")
        for s in snps :
            print(f"\t{s}")
    
    outputFile = "../data/anon.json"
    with open(outputFile, "w") as of :
        json.dump(resultDict, of, indent=4)


if __name__ == '__main__':
    main3()