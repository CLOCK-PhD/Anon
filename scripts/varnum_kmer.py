#!/usr/bin/python3

"""
varnum_kmer.py

Un programme comptant le nombre de sites de variation possibles qu'on peut retrouver dans un k-mer de taille donnée.
Les sites de variations sont définis sur la séquence de référence grâce aux données dbSNP, pour un chromosome.

Entrée :
    Séquence de référence du chromosome
    Fichier vcf de dbSNP

Sortie :
    Histogramme du nombre de variations au sein d'un k-mer
"""

import re
import argparse
import os
from Bio import SeqIO

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

    ref_seq = "../data/grch38p13/GCF_000001405.39/chrY.fna"   # Séquence de référence du chromosome
    dbsnp_vcf = "../data/snp_latest/24.vcf"                                   # Fichier .vcf dbSNP du chromosome
    kmerSize = 31

    # Récupérer la séquence du chromosome en mémoire - Aucune idée de pourquoi j'ai mis ça ici
    """seq = []
    with open(ref_seq) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = record.seq"""

    # Ouverture du fichier vcf
    with open(dbsnp_vcf, "r") as vcf:

        count = 0   # Pour stoper la boucle pendant les tests

        position = 0
        line = vcf.readline()
        while line and count < 20:  # Retirer le compteur à la fin des tests
            proximity_count = 0
            # Début de la boucle
            print()
            print(f"Boucle {count}")
            vcf.seek(position)
            print(f"Position du curseur : {vcf.tell()}")
            #print(line)
            # Récupérer les infos de la ligne :
            chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = get_vcf_line_info(line)
            print(f"\tChrom : {chrom}\tPosition : {snp_pos}")
            
            
            # Lire la ligne suivante :
            line  = vcf.readline()
            # Donner la position de la nouvelle ligne
            position = vcf.tell()
            #print(position)
            count += 1  # Pour les tests, à supprimer
            chrom2, snp_ref2, snp_pos2, rs_id2, snp_alt2, vc2 = get_vcf_line_info(line)
            print(f"\tNouvelle position : {snp_pos2}")


            # Test pour la comparaison et la lecture des lignes suivantes :
            while snp_pos2 < snp_pos + kmerSize :
                proximity_count += 1
                line2 = vcf.readline()
                chrom2, snp_ref2, snp_pos2, rs_id2, snp_alt2, vc2 = get_vcf_line_info(line2)
                #print(snp_pos2)

            print(f"\tCompteur de SNP à proximité : {proximity_count}")


            

    print("le chat")

if __name__ == '__main__':
    main()