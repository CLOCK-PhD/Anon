#!/usr/bin/python3

"""
varnum_kmer.py

Un programme comptant le nombre de sites de variation possibles qu'on peut retrouver dans un k-mer de taille donnée.

Le nombre de variations par kmer est stockée dans un dictionnaire.
    Clé : Nombre de variations pour un kmer
    Valeur : Occurences (nombre de fois qu'on trouve un kmer avec ce nombre de variations)

Entrée :
    Fichier vcf de dbSNP

Sortie :
    Histogramme du nombre de variations au sein d'un k-mer

    AUTRE METHODE - liste
"""


# RÉSOLUTION POSSIBLE DU PROBLÈME : 
# Lire une première fois le fichier, compter le nombre de ligne, arrêter la boucle à len(file) - kmerSize

import re
import argparse
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

def get_vcf_line_info(line)-> tuple:
    """Récupère les informations contenues dans chaque ligne du fichier VCF de SNPdb.
    Retourne un tuple qui contient dans l'ordre :
        - chrom (str):          Nom du chromosome
        - snp_ref (str):        SNP de référence
        - snp_pos (int) :       Position du SNP dans le chromosome
        - rs_id (str):          Identifiant du SNP
        - snp_alt (list(str)):  Liste contenant toutes les variations connues du SNP
        - vc (str):             Variation Class, le type de variations 

    Parameters :
        line: Ligne du fichier vcf

    Returns:
        tuple: chrom(str), snp_ref(str), snp_pos(int), rs_id(str), snp_alt(list(str)), vc(str)
    """
    description = line.split("\t")
    chrom = description[0]
    chrom_res = re.search("NC_0*(.*)\.", chrom)
    if chrom_res:
        chrom = chrom_res.group(1)
        if chrom == "24":
            chrom = "Y"
        if chrom == "23":
            chrom = "X"
    snp_pos = description[1]
    rs_id = description[2]
    snp_ref = description[3]
    snp_alt = description[4].split(",")
    #qual = description[5]
    #filter = description[6]
    info = description[7]
    res = re.search("VC=(\w*)", info)
    if res:
        vc = res.group(1)
    else:
        vc = ""
        print(rs_id)
    return chrom, snp_ref, int(snp_pos), rs_id, snp_alt, vc


def main():

    # Gestion des arguments
    parser = argparse.ArgumentParser()

    # Création des arguments
    parser = argparse.ArgumentParser(description='Draws a bar plot showing the number of variations in a k-mer from a vcf file frome dbSNP')
    parser.add_argument("-i", "--input", dest="input", help="Path to the directory containing the vcf file")
    parser.add_argument("-k", "--kmer_size", dest="kmer_size", default=31, type=int, help="Select k-mer size")

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    dbsnp_vcf = args.input          # Dossier contenant l'index
    kmerSize = args.kmer_size          # Taille du k-mer

    #dbsnp_vcf = "../data/snp_latest/common_14.vcf"                     # Fichier .vcf dbSNP du chromosome
    #kmerSize = 31

    # Dictionnaire des kmers et leurs variations :
    varnumDict = {}

    # Récupérer le nombre de lignes
    lineNumber = 0
    with open(dbsnp_vcf, "r") as vcf :
        for line in vcf :
            lineNumber += 1

    # Ouverture du fichier vcf
    pbar = tqdm(total=lineNumber - kmerSize)
    with open(dbsnp_vcf, "r") as vcf:

        count = 0
        position = 0    # Définir la position de départ du curseur
        line = vcf.readline()
        while line and count < lineNumber - kmerSize :
            pbar.update(1)
            proximity_count = 0
            vcf.seek(position)
            chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = get_vcf_line_info(line) # Récupérer les infos de la ligne
            
            line  = vcf.readline()  # Lire la ligne suivante
            position = vcf.tell()   # Donner la position de la nouvelle ligne
            count += 1
            chrom2, snp_ref2, snp_pos2, rs_id2, snp_alt2, vc2 = get_vcf_line_info(line)

            # Comparaison et lecture des lignes suivantes :
            if line and snp_pos != snp_pos2 : # test pour le cas de la dernière ligne qui fait une erreur
                while snp_pos2 < snp_pos + kmerSize :
                    #proximity_count += 1
                    line2 = vcf.readline()
                    if line2 :
                        chrom2, snp_ref2, snp_tmp, rs_id2, snp_alt2, vc2 = get_vcf_line_info(line2)
                        if snp_tmp == snp_pos2 :
                            continue
                        else :
                            proximity_count += 1
                            snp_pos2 = snp_tmp

                # Remplissage du dictionnaire :
                try :
                    varnumDict[proximity_count] += 1    # Incrémentation si la valeur des variations existe déjà
                except KeyError :
                    varnumDict[proximity_count] = 1     # Création de la clé et initialisation de sa valeur à 1
            #print(f"\tCompteur de SNP à proximité : {proximity_count}")

    pbar.close()

    print(varnumDict)

    # Histogramme

    fig, ax = plt.subplots()
    y = varnumDict.values() # y
    x = list(varnumDict.keys()) # y
    ax.bar(x, y, log = True, ec="k", color="red")
    ax.set_xlabel(f"Nombre de variations par {kmerSize}-mer")
    ax.set_ylabel("Nombre de cas")
    plt.show()

    print("le chat")

if __name__ == '__main__':
    main()