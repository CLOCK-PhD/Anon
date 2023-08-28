#!/usr/bin/python3

"""
umer_gen.py : u-mer Generator

Programme qui va générer un fichier au format fasta contenant une séquence appelée u-mer.
Chaque u-mer est la séquence de longueur 2k-1 située autour de la position d'un SNP indiqué dans dbSNP.
Ce fichier pourra ensuite être utilisé par Jellyfish, un programme de comptage de k-mers, pour générer tous les k-mers possibles contenant un SNP.

DÉFINITIONS :

    K-MERS :
    Nucléotide de taille k.

    U-mers :
    Les k-mers dont on a besoin sont générés à partir de (2k-1)mers ;
    Ils seront nommés u-mer (k=11, u=2*11-1).

    Variants :
    Les variants (de la classe Variant), sont des k-mers identiques à d'autres,
    mais porteur d'une variation différente.

"""


import re
import os
import glob
import argparse

from itertools import product
from pprint import pprint
from Bio import SeqIO
from typing import OrderedDict
from tqdm import tqdm
from sys import getsizeof

#from variant import Variant

# Récupérer les informations contenues dans le VCF
# On pourrait en faire un objet
def getVcfLineInfo(line)-> tuple:
    """Récupère les informations contenues dans chaque ligne du fichier VCF de SNPdb.
    Retourne un tuple qui contient dans l'ordre :
        - chrom     (str):          Nom du chromosome
        - snp_ref   (str):          SNP de référence
        - snp_pos   (int):          Position du SNP dans le chromosome
        - rs_id     (str):          Identifiant du SNP
        - snp_alt   (list(str)):    Liste contenant toutes les variations connues du SNP
        - vc        (str):          Variation Class, le type de variations 

    Parameters:
        line:       Ligne du fichier vcf

    Returns:
        tuple:      chrom(str), snpRef(str), snpPos(int), rsid(str), snpAlt(list(str)), vc(str)
    """
    description = line.split("\t")
    chrom = description[0]
    chromRes = re.search("NC_0*(.*)\.", chrom)
    if chromRes:
        chrom = chromRes.group(1)
        if chrom == "24":
            chrom = "Y"
        if chrom == "23":
            chrom = "X"
    snpPos = description[1]
    rsid = description[2]
    snpRef = description[3]
    snpAlt = description[4].split(",")
    #qual = description[5]
    #filter = description[6]
    info = description[7]
    res = re.search("VC=(\w*)", info)
    if res:
        vc = res.group(1)
    else:
        vc = ""
        print(rsid)
    return chrom, snpRef, int(snpPos), rsid, snpAlt, vc

# Récupérer les u-mers pour SNV :
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
    #kmerPos = snpPos - kmerSize + 1 # Position du k-mer dans le génome
    l_kmer = sequence[snpPos - kmerSize + 1 : snpPos]
    r_kmer = sequence[snpPos+1 : snpPos + kmerSize]
    for alt in snpAlt :
        umer = l_kmer + alt + r_kmer
        # Ici pour supprimer les umers avec N
        #umerList.append((umer, kmerPos))
        # TEST - AJOUT DU ALT
        umerList.append((umer, alt))
    return umerList

# Récupérer le u-mer pour DEL :
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
    #kmer_pos = int(snpPos - kmerSize + 1)
    l_kmer = sequence[snpPos - kmerSize + 1 : snpPos]
    snp = sequence[snpPos]
    r_kmer = sequence[snpPos + (len(snpRef)) : snpPos + kmerSize + len(snpRef) -1]
    umer = l_kmer + snp + r_kmer
    alt = "DEL"
    umerList.append((umer, alt))
    return umerList

# Récupérer les u-mers pour INS
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
            #kmerPos = snpPos - kmerSize + len(alt)
            l_kmer = sequence[snpPos - kmerSize + len(alt) : snpPos]
            r_kmer = sequence[snpPos + 1 : snpPos + kmerSize - (len(alt)-1)]
            umer = l_kmer + alt + r_kmer
            umerList.append((umer, alt))
        else:
            #print(f"INS TOO LONG FOR KMER SIZE : {len(alt)}")
            return umerList
    return umerList

# Récupérer les u-mers pour MNV
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
        #kmerPos = snpPos - kmerSize + len(snpRef)
        l_kmer = sequence[snpPos - kmerSize + len(snpRef): snpPos]
        r_kmer = sequence[snpPos + len(snpRef) : snpPos + kmerSize]
        for alt in snpAlt:
            umer = l_kmer + alt + r_kmer
            umerList.append((umer, alt))
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

# Extraction des u-mers
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


    """Génère les k-mers de taille voulue au format objet, à partir des u-mers.
    Retourne une liste contenant tous les k-mers.

    Parameters:
        kmerSize    (int):  Taille voulue du k-mer
        umerToCut   (str):  u-mer à découper en k-mer
        rsid        (str):  Identifiant du SNP
        chrom       (str):  Numéro du chromosome
        snp         (str):  Valeur du SNP
        snpPos      (int):  Position du SNP

    Returns:
        kmerList    (list(Variant)):    Liste contenant les k-mers et les variants possibles.

    """
    kmerList = []
    for i in range(0, kmerSize, 1):
        # kmer[0] : séquence ; kmer[1] : position du kmer
        relPos = int(umerToCut[1]) + i
        kmer = umerToCut[0][i : i + kmerSize].upper()
        var = Variant(rsid, chrom, snp, snpPos, relPos, 0, 0)
        if len(kmer) == kmerSize and "N" not in kmer:
            # : pour garder des kmer de taille voulue
            # : pour exclure les kmers avec un "N"
            kmerList.append([kmer, var])
    #pprint(f"kmerList : \n{kmerList}")
    # Ajuster le kmerCount de chaque variant
    for k in kmerList :
        k[1].kmersCount = len(kmerList)
        #print(f"{k[0]} - {k[1].rsid} {k[1].relPos} {k[1].kmersCount} {k[1].ambiguousKmersCount}")
    #pprint(kmerList)
    return kmerList

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

    # Gestion des arguments
    parser = argparse.ArgumentParser()

    # Création des arguments
    parser = argparse.ArgumentParser(description='Creates an index of k-mer carrying a variation, from dbSNP .vcf file and the reference sequence HG38.')
    parser.add_argument("-f", "--fasta_file", dest="fasta_file", help="fasta input file")
    parser.add_argument("-i", "--dbsnp", dest="dbsnp", help="dnSNP vcf file")
    parser.add_argument("-k", "--kmer_size", dest="kmer_size", default=21, type=int, help="(Optional) Select k-mer size")
    #parser.add_argument("-o", "--ouput", dest="output_dir", default="index", help="(Recommanded / Optional) Output folder name")

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    fastaFile = args.fasta_file
    kmerSize = args.kmer_size
    #outputDirectory = args.output_dir
    vcfFile = args.dbsnp

    # Adapter le nom du fichier de sortie au chromosome analysé
    res = re.search("^.*/(.*).fasta$", fastaFile)
    if res :
        outputFileName = "VarUmer" + "_chr" + res.group(1) + ".fasta"
    try :
        f = open(outputFileName, "x")
        f.close()
    except FileExistsError:
        outputFileName = uniquify(outputFileName)
        f = open(outputFileName, "x")
        f.close()
    

    # Lecture du fichier fasta pour générer les k-mers porteurs de SNP
    seq = []
    with open(fastaFile) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq.upper())

    # Lecture du fichier vcf
    vcfLineCount = 0
    with open(vcfFile, "r") as vcf:
        for line in vcf:
            vcfLineCount += 1
    print(f"Nombre de lignes dans le fichier VCF : {vcfLineCount}")

    # Parcours du fichier vcf
    print("Creating u-mer fasta file")
    lineCount = 0
    umerNumber = 0
    pbar = tqdm(total = vcfLineCount)
    with open(outputFileName, "a") as umerFasta :
        with open(vcfFile, "r") as vcf:
                for line in vcf:
                    pbar.update(1)
                    if lineCount < vcfLineCount :
                        # 3.1 Récupération des informations
                        chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = getVcfLineInfo(line)

                        
                        # Création de la liste des u-mers
                        umerList = getUmerFromPos(seq, snp_pos, vc, kmerSize, snp_ref, snp_alt)

                        # pour vkg.py : création des k-mers à partir des u-mers
                        for umer in umerList :
                            header = ">" + rs_id + "-" + vc + "-" + snp_ref + ":" + umer[1] + "\n"
                            umer_seq = umer[0] + "\n"
                            umerFasta.write(header)
                            umerFasta.write(umer_seq)
                            umerNumber += 1

                    else:
                        break
    pbar.close()
    print(f"{umerNumber} u-mers created in {outputFileName}")



if __name__ == '__main__':
    main()