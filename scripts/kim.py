#!/usr/bin/python3

"""
Parce que j'en ai marre.

DEROULEMENT

A - Création de l'index
    1. On ouvre le vcf
    2. On va lire la ligne et récuprer les infos selon ce qui nous intéresse
    3. On va récupérer la séquence correspondante dans le fasta
    4. On découpe la séquence en k-mer avec les fonctions dédiées aux VC
    5. On balance dans un dataframe
    5bis. On supprime les doublons
    6. On génère le dataframe en tsv
    7. On génère un index du dataframe

B - Vérification des kmers présents dans le génome de référence
    1. On se démerde pour trouver une méthode moins conne qu'avant

C- Requête
    1. Maintenant on a un putain d'index grâce à Pandas
    2. On lit le fastq (pysam again)
    3. On requête l'index
    4. Merveilleux.
"""

# NOTES
"""
Infos qui m'intéressent (ça peut être modifié):
 - rsid
 - position de la variation dans le k-mer (snpPos - relPos)
 - nombre de kmers créés
 - un bool
"""

import re
import os
import gzip
import glob
import pysam
import argparse
import pandas as pd

from pysam import VariantFile
from itertools import product
from pprint import pprint
from Bio import SeqIO
from typing import OrderedDict
from tqdm import tqdm
from sys import getsizeof


#########################
# MES FONCTIONS DE MORT ################################################################################
#########################


##########################################
# Fonctions converties depuis c++
##########################################
def isDegenerate(base):
    # Define a function to check if a base is degenerate
    return base in ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']

def kmer_generator(umer, k):
    # Generate k-mers from a u-mer, given a k-mer size k
    uLength = len(umer)
    i = 0
    while i <= uLength - k:
        containsDegenerate = False
        for j in range(i, i + k):
            if isDegenerate(umer[j]):
                containsDegenerate = True
                i = j + 1  # Skip to the position after the degenerated nt
                break
        if containsDegenerate:
            # Skip this k-mer if it contains a degenerate nucleotide
            continue
        kmer = umer[i:i + k]
        # Convert characters to uppercase before outputting
        kmer = kmer.upper()
        print(f"K-mer {i + 1}:\t{kmer}")
        i += 1


##########################################

# A ADAPTER // Récupérer les u-mers pour SNV
def get_SNV_umers(sequence:str, snpPos:int, kmerSize:int, snpAlt:list) -> list:
    """Pour un SNV : récupère les u-mers dans une liste, pour chaque variation possible.
    Ces u-mers sont utilisés pour générer tous les k-mers des SNP.

    Parameters:
        sequence    (str):          Séquence où fouiller
        snpPos      (int):          Position du SNP dans la séquence de référence
        kmerSize    (int):          Taille du k-mer
        snpAlt      (list(str)):    Liste contenant les variations du SNP

    Returns:
        umerList    (list):         Liste des u-mers pour chaque variation
    """
    umerList = []
    kmerPos = snpPos - kmerSize + 1
    l_kmer = sequence[snpPos - kmerSize + 1 : snpPos]
    r_kmer = sequence[snpPos+1 : snpPos + kmerSize]
    for alt in snpAlt :
        umer = l_kmer + alt + r_kmer
        # Ici pour supprimer les umers avec N
        #umerList.append((umer, kmerPos))
        # TEST - AJOUT DU ALT
        umerList.append((umer, kmerPos))
    return umerList

# A ADAPTER // Récupérer le u-mer pour DEL :
def get_DEL_umers(sequence:str, snpPos:int, kmerSize:int, snpRef:str) -> list:
    """Pour une délétion : récupère le u-mer dans une liste.
    Ce u-mer est utilisé pour générer tous les k-mers du SNP.

    Parameters:
        sequence    (SeqIO):    Séquence où rechercher
        snp_pos     (int):      Position du SNP dans la séquence de référence
        kmer_size   (int):      Taille du k-mer
        snp_ref     (str):      SNP de référence

    Returns:
        umerList    (list):     Liste des u-mers pour chaque variation    
    """
    umerList = []
    kmer_pos = int(snpPos - kmerSize + 1)
    l_kmer = sequence[snpPos - kmerSize + 1 : snpPos]
    snp = sequence[snpPos]
    r_kmer = sequence[snpPos + (len(snpRef)) : snpPos + kmerSize + len(snpRef) -1]
    umer = l_kmer + snp + r_kmer
    umerList.append((umer, kmer_pos))
    return umerList

# A ADAPTER // Récupérer les u-mers pour INS
def get_INS_umers(sequence:str, snpPos:int, kmerSize:int, snpAlt:list) -> list:
    """Pour une insertion : récupère les u-mers dans une liste, pour chaque variation possible.
    Ces u-mers max sont utilisés pour générer tous les k-mers des SNP.

    Parameters:
        sequence    (str):          Séquence où rechercher
        snpPos      (int):          Position du SNP dans la séquence de référence
        kmerSize    (int):          Taille du k-mer
        snpAlt      (list(str)):    Liste contenant les variations du SNP

    Returns:
        umerList    (list):         Liste des u-mers pour chaque variation
    """
    umerList = []
    for alt in snpAlt:
        if len(alt) < kmerSize:
            kmerPos = snpPos - kmerSize + len(alt)
            l_kmer = sequence[snpPos - kmerSize + len(alt) : snpPos]
            r_kmer = sequence[snpPos + 1 : snpPos + kmerSize - (len(alt)-1)]
            umer = l_kmer + alt + r_kmer
            umerList.append((umer, kmerPos))
        else:
            #print(f"INS TOO LONG FOR KMER SIZE : {len(alt)}")
            return umerList
    return umerList

# A ADAPTER // Récupérer les u-mers pour MNV
def get_MNV_umers(sequence:str, snpRef:str, snpPos:int, kmerSize:int, snpAlt:list) -> list:
    """Pour un MNV : récupère les u-mers dans une liste, pour chaque variation possible.
    Ces u-mers sont utilisés pour générer tous les k-mers des SNP.

    Parameters:
        sequence    (str):          Séquence où chercher
        snpRef      (str):          SNP de référence
        snpPos      (int):          Position du SNP dans la séquence de référence
        kmerSize    (int):          Taille du k-mer
        snpAlt      (list(str)):    Liste contenant les variations du SNP

    Returns:
        umerList    (list):         Liste des u-mers pour chaque variation
    """
    umerList = []
    if len(snpRef) < kmerSize:
        kmerPos = snpPos - kmerSize + len(snpRef)
        l_kmer = sequence[snpPos - kmerSize + len(snpRef): snpPos]
        r_kmer = sequence[snpPos + len(snpRef) : snpPos + kmerSize]
        for alt in snpAlt:
            umer = l_kmer + alt + r_kmer
            umerList.append((umer, kmerPos))
    else :
        #print(f"MNV TOO LONG FOR KMER SIZE : {len(snp_ref)}")
        return umerList
    return umerList

# Récupérer les u-mers pour INDEL
def get_INDEL_umers(sequence:str, snpRef:str, snpPos:int, snpAlt:list, kmerSize:int) -> list:
    """Pour un INDEL : récupère les u-mers dans une liste, pour chaque variation possible.
    Ces u-mers sont utilisés pour générer tous les k-mers des SNP.

    Parameters:
        sequence    (str):          Séquence où chercher
        snp_ref     (str):          SNP de référence
        snp_pos     (int):          Position du SNP dans la séquence de référence
        snp_alt     (list(str)):    Liste contenant les variations du SNP
        kmer_size   (int):          Taille du k-mer

    Returns:
        umerList    (list):         Liste des u-mers pour chaque variation
    """
    umerList = []
    if len(snpRef) > kmerSize :
        #print("INDEL REF IS TO LONG FOR K-MER SIZE")
        return umerList
    if len(snpRef) == 1:
        for alt in snpAlt:
            if len(alt) == len(snpRef):
                umerList += get_SNV_umers(sequence, snpPos, kmerSize, [alt])
            else :
                umerList += get_INS_umers(sequence, snpPos, kmerSize, [alt])
    else :
        for alt in snpAlt:
            if len(alt) == 1 :
                umerList += get_INS_umers(sequence, snpPos, kmerSize, [alt])
            elif len(snpRef) < len(alt) :
                umerList += get_INS_umers(sequence, snpPos, kmerSize, [alt])
            elif len(snpRef) > len(alt):
                max_kmer_list += get_INS_umers(sequence, snpPos, kmerSize, [alt])
    return max_kmer_list

# A ADAPTER // Extraction des u-mers et envoi à la fonction spécifique
def getUmerFromPos(sequence:str, pos:int, variantClass:str, kmerSize:int, snpRef:str, snpAlt:list) -> list:
    """Génère le u-mer en fonction de la classe de variant en appelant une fonction spécifique.
    Retourne une liste contenant les umers.

    Parameters:
        sequence        (str):          Séquence où chercher
        pos             (int):          Position du SNP dans la séquence de référence
        variantClass    (str):          Classe de variation pour l'appel à la fonction de génération de umer
        kmerSize        (int):          Taille du k-mer
        snpRef          (str):          SNP de référence
        snpAlt          (list(str)):    Liste contenant les variations du SNP

    Returns:
        umerList        (list):         Liste des u-mers pour chaque variation.

    """
    snpPos = pos - 1
    umerList = []

    # Exclusion des kmers contenant des points
    if snpRef == ".":
        return umerList
    for alt in snpAlt:
        if alt == ".":
            return umerList

    if variantClass == "SNV":
        return get_SNV_umers(sequence, snpPos, kmerSize, snpAlt)
    elif variantClass == "DEL":
        return get_DEL_umers(sequence, snpPos, kmerSize, snpRef)
    elif variantClass == "INDEL":
        return get_INDEL_umers(sequence, snpRef, snpPos, snpAlt, kmerSize)
    elif variantClass == "INS":
        return get_INS_umers(sequence, snpPos, kmerSize, snpAlt)
    elif variantClass == "MNV":
        return get_MNV_umers(sequence, snpRef, snpPos, kmerSize, snpAlt)
    
# Découper le u-mers pour récupérer les kmers et leurs positions
def kmerGenerator(kmerSize:int, umerToCut:str) -> list:
    """Génère les k-mers de taille voulue à partir des u-mers.
    Retourne une liste contenant tous les k-mers.

    Parameters :
        kmerSize    (int):  Taille du k-mer voulu
        umerToCut   (str):  Séquence du u-mer à découper en k-mers

    Returns:
        kmerList    (list(tuple)): Liste contenant les k-mers et leur position
    """
    #print(kmer_to_cut)
    kmerList = []
    for i in range(0, kmerSize, 1):
        # kmer[0] : séquence ; kmer[1] : position du kmer
        kmer = (umerToCut[0][i : i + kmerSize].upper(), int(umerToCut[1]) + i)
        if len(kmer[0]) == kmerSize and "N" not in kmer[0]: # pour garder des kmer de taille voulue
            kmerList.append(kmer)
    #pprint(f"kmerList : \n{kmerList}")
    return kmerList

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

###########
# LE MAIN ########################################################################################
###########
def main():

    # Gestion des arguments
    parser = argparse.ArgumentParser()

    # Création des arguments
    parser = argparse.ArgumentParser(description='Creates an index of k-mer carrying a variation, from dbSNP .vcf file and the reference sequence HG38.')
    parser.add_argument("-f", "--fasta", dest="fasta_file", help="fasta input file")
    parser.add_argument("-fai", "--faidx", dest="faidx", help="faidx index")
    parser.add_argument("-v", "--vcf", dest="dbsnp", help="dbSNP vcf file")
    parser.add_argument("-k", "--kmer_size", dest="kmer_size", default=21, type=int, help="(Optional) Select k-mer size")
    # Revenir là dessus : ça doit être le dossier à créer qui va contenir les fichiers
    parser.add_argument("-o", "--ouput", dest="output_dir", default="index", help="(Recommanded / Optional) Output folder name")

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    fastaFile = args.fasta_file
    faiPath = args.faidx
    kmerSize = args.kmer_size
    outputDirectory = args.output_dir
    vcfFile = args.dbsnp

    # Pour pas avoir à foutre des arguments à la con qui me gonflent à chaque test
    vcfFile = "/home/remycosta/phd/Anon/data/dbsnp156/GCF_000001405.40.gz"
    fastaFile = "/home/remycosta/phd/Anon/data/grch38p14/GCF_000001405.40_GRCh38.p14_genomic.fna"
    faiPath = "/home/remycosta/phd/Anon/data/grch38p14/GCF_000001405.40_GRCh38.p14_genomic.fna.fai"

    # Ouvrir les fichiers sources
    bcf_in = VariantFile(vcfFile)
    fasta = pysam.FastaFile(fastaFile, faiPath)

    # Initialize empty DataFrame with column names
    columns = ['kmer', 'rsID', 'num', 'in_genome']
    df = pd.DataFrame(columns=columns)

    # On compte les lignes pour tqdm et pour voir si ça marche
    """print("Counting lines...")
    line_count = 0
    with gzip.open(vcfFile, 'rt') as file:
        for line in file:
            line_count += 1

        print(f"Number of lines in {vcfFile}: {line_count}")"""

    # On rentre dedans
    for rec in bcf_in:
    #print(rec)
    #print(rec.alts)
    #print(rec.info['FREQ'])
    # truc à récupérer
    #chrom, snp_ref, snp_pos, rs_id, snp_alt, vc
        chrom = rec.chrom   # str
        snpRef = rec.ref    # str
        snpPos = rec.pos    # int
        rsId = rec.id       # string
        snpAlt = rec.alts   # tuple
        vc = rec.info["VC"] # string

        """print(type(chrom))
        print(type(snpRef))
        print(type(snpPos))
        print(type(rsId))
        print(type(snpAlt))
        print(type(vc))"""

        print(f"{rsId}\n\t{chrom}\n\t{snpPos}\n\t{snpRef}\n\t{snpAlt}\n\t{vc}\n")




    print("le chat")

if __name__ == '__main__':
    main()