#!/usr/bin/python3

# Extract k-mers around SNP in rs_ch.fas files from dbSNP
# ftp.ncbi.nih.gov/snp/organisms/human_9606/rs_fasta/

# Positions des éléments dans les entêtes :
# 0 = gnl; 1=dbSNP; 5=taxid; 6=mol; 10=suspect?
# 2 = rs
# 3 = position du SNP dans la séquence
# 4 = longueur de la séquence
# 7 = class
# 8 = alleles
# 9 = build

# SNP class description :
# 1 = SNV; 2=DIV; 3=HETERIZYGOUS; 4=STR; 5=NAMED; 6=NO VARIATION; 7=MIXED; 8=MNV
# https://www.ncbi.nlm.nih.gov/projects/SNP/snp_legend.cgi?legend=snpClass

# Outp222ut = tsv : Séquence k-mer, ID, chromosome, position du SNP, position DU KMER
# Problème avec les fichiers rs_ch.fas : pas de position du snp et du kmer

# A FAIRE : Faire fonctionner avec plusieurs inputs (essayer parallélisation ?)
#       IDEE : fournir un dossier en input et traiter tous les fichiers rs_ch.fas 
# A FAIRE : LAST : Tri de tas des fichiers (autre programme ?)

import re
import sys
import argparse
import os
from Bio import SeqIO
from numpy import equal
from matplotlib.pyplot import close
from pprint import pprint

# Gestion des arguments
parser = argparse.ArgumentParser()

# Fichier d'entrée
parser.add_argument("-i", "--input", dest="rs_fasta_file", help = "rs_ch.fas input file")
parser.add_argument("-k", "--kmer_size", dest="kmer_size", default=21, help="Select k-mer size")
parser.add_argument("-n", dest="kmers_per_output_file", default=100000, help="Number of kmers per output file for the heap merge")
parser.add_argument("-o", "--ouput", dest="output_dir", default="generated_kmers", help="Output folder name")

args = parser.parse_args()
input_file = args.rs_fasta_file
kmer_size = int(args.kmer_size)
kmers_per_file = int(args.kmers_per_output_file)
output_dir = args.output_dir

snp_count = 0
kmers = []
file_number = 0
#output_dir = "generated_kmers"
os.makedirs(output_dir)

print(input_file)

# Extraire les variants possibles :
def extract_snp_var(seq_info):
    res = re.search("alleles=\"./(.*)\"$", seq_info)
    if res:
        snp_var = res.group(1)
    else:
        return 0
    
    if len(snp_var) == 1:
        x = [snp_var]
        return x
    else:
        x = snp_var.split("/")
        return x

# Sélectionner la séquence à diviser en k-mers :
def make_max_kmer(kmer_size, seq, snp_pos, seq_len):
    if snp_pos >= kmer_size : # cas ideal où len(maxkmer) = 2*kmer_size) -1
        max_kmer = seq[snp_pos - kmer_size : snp_pos + kmer_size -1]
        return max_kmer
    else :
        max_kmer = seq[0 : kmer_size + (snp_pos - 1)]
        return max_kmer

# Fonction pour remplacer le snp du kmer par ses variations connues :
def max_kmer_variations(max_kmer, snp_pos, snp_var, kmer_size):
    max_kmers_list = []
    
    if(snp_var != 0):
        for snp in snp_var :
            if len(snp) >= kmer_size:
                snp_var.remove(snp)
        if len(max_kmer) == 2*kmer_size - 1:
            for snp in snp_var:
                if len(snp) == 1:   # Cas où le SNP est de taille 1 et qu'on a un kmermax
                    max_kmer_var = max_kmer[:kmer_size - 1] + snp + max_kmer[kmer_size:]
                    max_kmers_list.append(max_kmer_var)
                else : # cas du snp long
                    max_kmer_var = max_kmer[len(snp)-1:kmer_size-1] + snp + max_kmer[kmer_size:len(max_kmer)-(len(snp)-1)]
                    max_kmers_list.append(max_kmer_var)
        else :
            for snp in snp_var:
                if len(snp) == 1 :
                    max_kmer_var = max_kmer[:snp_pos-1] + snp + max_kmer[kmer_size - snp_pos:]
                    max_kmers_list.append(max_kmer_var)
                else:
                    max_kmer_var = max_kmer[:snp_pos-1] + snp + max_kmer[kmer_size - snp_pos:len(max_kmer)-(len(snp)-1)]
                    max_kmers_list.append(max_kmer_var)
    else :
        return max_kmers_list
    return max_kmers_list

# Générer les kmers des snp
def kmer_generator(kmer_size, kmer_to_cut):
    kmer_list = []
    for i in range(0, kmer_size, 1):
        kmer = kmer_to_cut[i : i + kmer_size]
        if len(kmer) == kmer_size: # pour garder des kmer de taille voulue
            kmer_list.append(kmer)

    return kmer_list

# Extraire les kmers à partir des données rs_fasta
with open(input_file) as handle :
    for record in SeqIO.parse(handle, "fasta"):
        # Séparer les informations de l'entête :
        seq_info = record.description.split("|")
        seq_rs = seq_info[2].split(" ")[0]
        snp_pos = int(seq_info[3].split("=")[1])
        seq_len = int(seq_info[4].split("=")[1])
        
        seq_snp_var = extract_snp_var(seq_info[8])
        
        # Sélection du grand k-mer à découper en k-mer de taille voulue
        max_kmer = make_max_kmer(kmer_size, record.seq, snp_pos, seq_len)
        max_kmers_list = max_kmer_variations(str(max_kmer), snp_pos, seq_snp_var, kmer_size)

        # Génération des kmers à partir de la liste des variations de kmermax :
        record_kmer_list = []
        for var in max_kmers_list:
            kmer_list = kmer_generator(kmer_size, str(var))
            for kmer in kmer_list:
                record_kmer_list.append(kmer)

        
        for kmer in record_kmer_list :
            # Ajouter les kmers à la liste des kmers
            if len(kmers) < kmers_per_file :
                kmer_id = (kmer, seq_rs)
                kmers.append(kmer_id)
            # Exporter la liste quand on atteint un nombre de kmers dans la liste :
            if len(kmers) == kmers_per_file :
                output_file_name = f"{output_dir}/{str(file_number)}_snp_k {kmer_size}.tsv"
                kmers.sort()
                with open(output_file_name, "w", encoding="utf-8") as f:
                    for kmer in kmers :
                        line_output = f"{kmer[0]}\t{kmer[1]}\n"
                        f.write(line_output)
                kmers=[]
                file_number += 1

        # Exporter les derniers kmers dans un fichier :        
        output_file_name = f"{output_dir}/{str(file_number)}_snp_k {kmer_size}.tsv"
        kmers.sort()
        with open(output_file_name, "w", encoding="utf-8") as f:
            for kmer in kmers :
                line_output = f"{kmer[0]}\t{kmer[1]}\n"
                f.write(line_output)
        
        snp_count += 1


print(f"Nombre total de SNP dans le chromosome 1 : {snp_count}")