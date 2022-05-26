#!/usr/bin/python3

"""
Extraire des k-mers de taille voulue (k=21 par défaut) depuis la séquence de référence du chromosome directement, 
à partir des infos vcf du fichier ref snp

Génère des fichiers .tsv contenant 100000 kmers (valeur par défaut) triés par ordre lexicographique
avec les informations suivantes :
Kmer_seq   rs_id   chromosome  snp_position    kmer_position

A cause du nombre de fichiers qui peuvent être produits, tous les 1000 fichiers .tsv,
un merge est réalisé pour les trier dans un nouveau fichier.

Tous les fichiers sont ensuite à nouveau triés et réunis dans un fichier final_merge.tsv
"""

"""
# MODIFICATIONS A APPORTER :
# Intégrer un index de k-mer en préfixe/suffixe.
# On va commencer les tests avec un préfixe de longueur 5 (soit 1024 fichiers outputs) sur le chrY

Idée : Une fois les merges partiels réalisés, on va réaliser l'index au lieu de produire le fichier unique
qui est normalement crée.

Problème à régler : dans le cas où on n'a qu'un seul fichier à la fin (par exemple si on a un nombre de fichiers
kmers moins grand que le la limite du merge, comme c'est le cas du chrY)
Dans ce cas, on devrait générer à partir des fichiers qui restent sans faire de merge.
Dans les autres cas, on fait le merge sur les derniers fichiers qui restent, puis on génère l'index à partir des fichiers
de merge.

Autre problème à gérer : l'écriture dans les fichiers. On ne pourra pas les ouvrir tous en même temps,
donc il faut pouvoir ajouter les k-mers les uns après les autres.

Idée : On peut générer les fichiers en amont avec des fonctions de recombinaisons qui vont faire une liste contenant
tous les noms possibles des suffixes, puis ajouter les kmers dans ces fichiers si leur préfixe correspond.
Une fois le merge final réalisé, on va chercher les fichiers vides et les supprimer.

Idée 2 : Créer les fichiers à la volée en fonction de leur préfixe et les remplir comme ça.

# A FAIRE : Docstring
# A FAIRE : Lister les k-mers en double et leur nombre
# ? A FAIRE : Créer un fichier contenant tous les rs_id des kmers supprimés
# ? A FAIRE : logs pour les kmers rejetés
# ? A FAIRE : lire tous les chromosomes pour lancer le programme en une fois

    - Détecter les chromosomes différents
    - Créer un dossier séparé pour chaque chromosome
    - Exécuter le programme pour chaque
    Idée :
    - Paralléliser les exécutions pour faire plusieurs chromosomes en même temps
    - Implique de parcourir le fichier vcf une première fois pour lister tous les chromosomes différents

# ? A FAIRE : paralléliser le merge pendant la génération de kmer
# IDÉE : Barre de progression
"""

import re
import sys
import argparse
import os
import heapq
from Bio import SeqIO
from typing import OrderedDict
from pprint import pprint
from os.path import isfile, join

# Récupérer les informations contenues dans le VCF
def get_vcf_line_info(line)-> tuple:
    description = line.split("\t")
    chrom = description[0]
    snp_pos = description[1]
    rs_id = description[2]
    snp_ref = description[3]
    snp_alt = description[4].split(",")
    #qual = description[5]
    #filter = description[6]
    info = description[7]
    #vc = "pas détecté"
    res = re.search("VC=(\w*)", info)
    if res:
        vc = res.group(1)
    else:
        vc = ""
        print(rs_id)
    return chrom, snp_ref, int(snp_pos), rs_id, snp_alt, vc

# Récupérer les kmer_max pour SNV :
def get_SNV_kmer_max(sequence:SeqIO, snp_pos:int, kmer_size:int, snp_alt:list) -> list:
    kmer_max_list = []
    kmer_pos = snp_pos - kmer_size + 1
    l_kmer = sequence[snp_pos - kmer_size + 1 : snp_pos]
    r_kmer = sequence[snp_pos+1 : snp_pos + kmer_size]
    for alt in snp_alt :
        kmer_max = l_kmer + alt + r_kmer
        kmer_max_list.append((kmer_max, kmer_pos))
    #print(f"Kmer pos :{snp_pos} - {kmer_size} + 1 = {kmer_pos}")
    return kmer_max_list

# Récupérer le kmer_max pour DEL :
def get_DEL_kmer_max(sequence:SeqIO, snp_pos:int, kmer_size:int, snp_ref:str) -> list:
    #print("DEL")
    kmer_max_list = []
    kmer_pos = int(snp_pos - kmer_size + 1)
    l_kmer = sequence[snp_pos - kmer_size + 1 : snp_pos]
    snp = sequence[snp_pos]
    r_kmer = sequence[snp_pos + (len(snp_ref)) : snp_pos + kmer_size + len(snp_ref) -1]
    kmer_max = l_kmer + snp + r_kmer
    kmer_max_list.append((kmer_max, kmer_pos))
    #print(f"Kmer pos :{snp_pos} - {kmer_size} + 1 = {kmer_pos}")
    #print(kmer_max_list)
    return kmer_max_list

# Récupérer les kmer_max pour INS
def get_INS_kmer_max(sequence:SeqIO, snp_pos:int, kmer_size:int, snp_alt:list) -> list:
    kmer_max_list = []
    for alt in snp_alt:
        if len(alt) < kmer_size:
            kmer_pos = snp_pos - kmer_size + len(alt)
            l_kmer = sequence[snp_pos - kmer_size + len(alt) : snp_pos]
            r_kmer = sequence[snp_pos + 1 : snp_pos + kmer_size - (len(alt)-1)]
            kmer_max = l_kmer + alt + r_kmer
            kmer_max_list.append((kmer_max, kmer_pos))
            #print(f"Kmer pos :{snp_pos} - {kmer_size} + {len(alt)} = {kmer_pos}")
        else:
            #print(f"INS TOO LONG FOR KMER SIZE : {len(alt)}")
            return kmer_max_list
    return kmer_max_list

# Récupérer les kmer_max pour MNV
def get_MNV_kmer_max(sequence:SeqIO, snp_ref:str, snp_pos:int, kmer_size:int, snp_alt:list) -> list:
    kmer_max_list = []
    if len(snp_ref) < kmer_size:
        kmer_pos = snp_pos - kmer_size + len(snp_ref)
        l_kmer = sequence[snp_pos - kmer_size + len(snp_ref): snp_pos]
        r_kmer = sequence[snp_pos + len(snp_ref) : snp_pos + kmer_size]
        for alt in snp_alt:
            kmer_max = l_kmer + alt + r_kmer
            kmer_max_list.append((kmer_max, kmer_pos))
            #print(f"Kmer pos :{snp_pos} - {kmer_size} + {len(snp_ref)} = {kmer_pos}")
    else :
        #print(f"MNV TOO LONG FOR KMER SIZE : {len(snp_ref)}")
        return kmer_max_list
    return kmer_max_list

# Récupérer les kmer_max pour INDEL
def get_INDEL_kmer_max(sequence:SeqIO, snp_ref:str, snp_pos:int, snp_alt:list, kmer_size:int) -> list:
    max_kmer_list = []
    if len(snp_ref) > kmer_size :
        #print("INDEL REF IS TO LONG FOR K-MER SIZE")
        return max_kmer_list
    if len(snp_ref) == 1:
        for alt in snp_alt:
            if len(alt) == len(snp_ref):
                max_kmer_list += get_SNV_kmer_max(sequence, snp_pos, kmer_size, [alt])
            else :
                max_kmer_list += get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])
    else :
        for alt in snp_alt:
            if len(alt) == 1 :
                max_kmer_list += get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])
            elif len(snp_ref) < len(alt) :
                max_kmer_list += get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])
            elif len(snp_ref) > len(alt):
                max_kmer_list += get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])
    return max_kmer_list

# Extraction des kmer_max
def get_kmer_from_pos(sequence:SeqIO, pos:int, variant_class:str, kmer_size:int, snp_ref:str, snp_alt:list) -> list:
    snp_pos = pos - 1
    kmer_max_list = []

    # Exclusion des kmers contenant des points
    if snp_ref == ".":
        return kmer_max_list
    for alt in snp_alt:
        if alt == ".":
            return kmer_max_list

    if variant_class == "SNV":
        return get_SNV_kmer_max(sequence, snp_pos, kmer_size, snp_alt)
    elif variant_class == "DEL":
        return get_DEL_kmer_max(sequence, snp_pos, kmer_size, snp_ref)
    elif variant_class == "INDEL":
        return get_INDEL_kmer_max(sequence, snp_ref, snp_pos, snp_alt, kmer_size)
    elif variant_class == "INS":
        return get_INS_kmer_max(sequence, snp_pos, kmer_size, snp_alt)
    elif variant_class == "MNV":
        return get_MNV_kmer_max(sequence, snp_ref, snp_pos, kmer_size, snp_alt)

# Découper le kmer_max pour récupérer les kmers et leurs positions
def kmer_generator(kmer_size:int, kmer_to_cut:SeqIO) -> SeqIO:
    #print(kmer_to_cut)
    kmer_list = []
    for i in range(0, kmer_size, 1):
        kmer = (kmer_to_cut[0][i : i + kmer_size].upper(), int(kmer_to_cut[1]) + i)
        if len(kmer[0]) == kmer_size: # pour garder des kmer de taille voulue
            kmer_list.append(kmer)
    #pprint(f"kmer_list : \n{kmer_list}")
    return kmer_list

# Merge des fichiers kmers
def merge_kmers(output_dir:str, merged_file_number, kmer_files:list)-> str:
    files = [open(filename, "r") for filename in kmer_files]
    output_merged_file_name = f"{output_dir}/{str(merged_file_number)}_merge.tsv"
    merged = heapq.merge(*files)
    prev_line = ""
    print("\tMerging files...")
    with open(output_merged_file_name, "w", encoding="utf-8") as f:
        for line in merged:
            kmer = line.split("\t")[0]
            if kmer == prev_line :
                prev_line = kmer
            else:
                f.write(line)
                prev_line = kmer
    print("\tMerging done")
    
    # Suppression des fichiers
    print("\t\tDeleting kmer files...")
    for kmer_file_name in kmer_files :
        os.remove(f"{kmer_file_name}")
    print("\t\tk-mer files deleted")

    return output_merged_file_name

def main() :
    # Gestion des arguments
    parser = argparse.ArgumentParser()

    # Création des arguments
    parser = argparse.ArgumentParser(description='Extract kmers at SNP positions from a reference VCF and fasta file')
    parser.add_argument("-i", "--input", dest="fasta_file", help="fasta input file")
    parser.add_argument("-r", "--reference", dest="ref", help="VCF SNP reference file")
    parser.add_argument("-k", "--kmer_size", dest="kmer_size", default=21, help="Select k-mer size")
    parser.add_argument("-n", dest="kmers_per_output_file", default=100000, help="Number of kmers per output file for the heap merge")
    parser.add_argument("-o", "--ouput", dest="output_dir", default="generated_kmers", help="Output folder name")
    parser.add_argument("-b", "--batch_size", dest="batch_size", default=1000, type=int, help="Limit of k-mer files to merge during execution")

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    input_file = args.fasta_file
    kmer_size = int(args.kmer_size)
    kmers_per_file = int(args.kmers_per_output_file)
    output_dir = args.output_dir
    ref = args.ref
    batch_size = args.batch_size

    # Créer le dossier de sortie
    os.makedirs(output_dir)
    # Dictionnaire contenant les kmers
    kmers = {}
    # Numéro du fichier de sortie
    file_number = 0
    # Variables pour le merge
    output_file_list = []   # Liste des fichiers de kmers à merge
    merged_file_number = 0  # Numéro du fichier mergé
    merged_file_list_for_final_merge = []

    # Récupérer la séquence du chromosome en mémoire
    seq = []
    with open(input_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = record.seq

    # Ouverture et parcours du fichier vcf de référence
    print("Generating k-mers...")
    print(f"Generating k-mer files : batch {merged_file_number}")
    with open(ref, "r") as vcf:
        for line in vcf:
            # Merge et suppression des fichiers
            if len(output_file_list) == batch_size :
                print(f"\tMerging k-mers batch {merged_file_number}")
                # Merge les fichiers de output_file_list
                merged_file_list_for_final_merge.append(merge_kmers(output_dir, merged_file_number, output_file_list))
                # Incrémentation du nom du fichier des kmers mergés
                merged_file_number += 1
                # Réinitialisation de la liste des fichiers à merge
                output_file_list = []
                print(f"Generating k-mer files : batch {merged_file_number}")

            chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = get_vcf_line_info(line)
                
            # Générer les kmers depuis les kmers_max
            kmer_max_list = get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt)
            for kmer_max in kmer_max_list:
                kmer_list = kmer_generator(kmer_size, kmer_max)

            # Place les kmers générés dans le dictionnaire
            for kmer in kmer_list :
                if len(kmers) < kmers_per_file :
                    kmers[kmer[0]] = (rs_id, chrom, snp_pos, kmer[1])
                if len(kmers) == kmers_per_file :
                    output_file_name = f"{output_dir}/{str(file_number)}_snp_k{kmer_size}.tsv"
                    kmers = OrderedDict(sorted(kmers.items()))
                    with open(output_file_name, "w", encoding="utf-8") as f:
                        for kmer, values in kmers.items() :
                            line_output = f"{kmer}\t{values[0]}\t{values[1]}\t{values[2]}\t{values[3]}\n"
                            f.write(line_output)
                    # Réinitialisation du dictionnaire des kmers
                    kmers={}
                    # Ajout du nom de fichier dans la liste des fichiers à merge
                    output_file_list.append(output_file_name)
                    file_number += 1

    # Boucle pour les tests
    """with open(ref, "r") as vcf:
        for line in vcf:
            # Boucle pour les tests :
            # ChY : 2375594
            
            # Merge et suppression des fichiers
            if len(output_file_list) == 10 :
                # Merge les fichiers de output_file_list
                merged_file_list_for_final_merge.append(merge_kmers(output_dir, merged_file_number, output_file_list))
                # Incrémentation du nom du fichier des kmers mergés
                merged_file_number += 1
                # Réinitialisation de la liste des fichiers à merge
                output_file_list = []
            if count <= 100000:
                chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = get_vcf_line_info(line)
                
                # Test pour générer les kmers depuis les kmers_max
                kmer_max_list = get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt)
                for kmer_max in kmer_max_list:
                    kmer_list = kmer_generator(kmer_size, kmer_max)

                # Place les kmers générés dans le dictionnaire :
                # Output = tsv : Séquence k-mer, ID, chromosome, position du SNP, position DU KMER
                for kmer in kmer_list :
                    if len(kmers) < kmers_per_file :
                        kmers[kmer] = (rs_id, chrom, snp_pos)
                    if len(kmers) == kmers_per_file :
                        output_file_name = f"{output_dir}/{str(file_number)}_snp_k{kmer_size}.tsv"
                        kmers = OrderedDict(sorted(kmers.items()))
                        with open(output_file_name, "w", encoding="utf-8") as f:
                            for kmer, values in kmers.items() :
                                line_output = f"{kmer[0]}\t{values[0]}\t{values[1]}\t{values[2]}\t{kmer[1]}\n"
                                f.write(line_output)
                        # Réinitialisation du dictionnaire des kmers
                        kmers={}
                        # Ajout du nom de fichier dans la liste des fichiers à merge
                        output_file_list.append(output_file_name)
                        file_number += 1
                count += 1
            else:
                break"""
    
    # Ajouter les derniers fichiers de kmers dans la liste du merge final
    merged_file_list_for_final_merge += output_file_list
    # Exporter les derniers kmers dans un fichier
    output_file_name = f"{output_dir}/{str(file_number)}_snp_k{kmer_size}.tsv"
    kmers = OrderedDict(sorted(kmers.items()))
    with open(output_file_name, "w", encoding="utf-8") as f:
        for kmer, values in kmers.items() :
            line_output = f"{kmer}\t{values[0]}\t{values[1]}\t{values[2]}\t{values[3]}\n"
            f.write(line_output)
        # Ajout du dernier fichier à la liste du merge final
        merged_file_list_for_final_merge.append(output_file_name)

    # Merge des derniers fichiers
    print("Final merge...")
    merge_kmers(output_dir, "final", merged_file_list_for_final_merge)

    print("All k-mers merged")

if __name__ == '__main__':
    main()