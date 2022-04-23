#!/usr/bin/python3

# Extract k-mers around SNP in rs_fasta files from dbSNP
# ftp.ncbi.nih.gov/snp/organisms/human_9606/rs_fasta/

# Positions des éléments dans les entêtes :
# 0 = gnl
# 1 = dbSNP
# 2 = rs
# 3 = position du SNP dans la séquence
# 4 = longueur de la séquence

from pprint import pprint
from Bio import SeqIO
import re
from matplotlib.pyplot import close

from numpy import equal

#input_file = "../data/example.fas"
input_file = "../data/rs_ch1.fas"

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
    max_kmer = seq[snp_pos - kmer_size : snp_pos + kmer_size -1]
    return max_kmer

# Générer les kmers des snp
def kmer_generator(kmer_size, max_kmer, snp_var):
    for snp in snp_var:
        a=1
        
    return 0

kmer_size = 2 # GERER LES SNP_POS < KMER_SIZE !!!!!

count = 0

origin_snp = ["N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "A", "T", "C", "G"]
# refaire sans ATCG et enregistrer les sorties pour ces nt
# structure de k-mers voulue : Séquence, ID, chromosome, position du SNP, position DU KMER

# Extraire les kmers à partir des données rs_fasta
outputFile = open("../data/weird_fucks.txt", "w")

with open(input_file) as handle :
    for record in SeqIO.parse(handle, "fasta"):
        print(record.description)
        
        # Séparer les informations de l'entête :
        seq_info = record.description.split("|")
        snp_pos = int(seq_info[3].split("=")[1])
        seq_len = int(seq_info[4].split("=")[1])
        
        seq_snp_var = extract_snp_var(seq_info[8])
        print(f"snp_pos : {snp_pos}\t seq_len : {seq_len}\t snp_var : {seq_snp_var}")
        #print(seq_info)
        #print(record.seq)
        print(f"Le SNP : {record.seq[int(snp_pos)-1]}")
        max_kmer = make_max_kmer(kmer_size, record.seq, snp_pos, seq_len)
        print(max_kmer)
        print(len(max_kmer))
        print(max_kmer[kmer_size-1])
        #print(len(record))
        count += 1

        # vérification des char utilisés pour les SNP:
        if record.seq[int(snp_pos)-1] not in origin_snp :
            outputFile.write(record.description)
            outputFile.write("\n")
            break

outputFile.close()

print(f"Nombre total de SNP dans le chromosome 1 : {count}")

# Même chose avec une liste des positions pour chercher directement 
# dans la primary assembly par chromosome
# A modifier pour que ça marche avec chaque chromosome
"""
chrom1_snp_pos = [10019, 10039]
kmer_snp = []
# parcourir les séquences
with open("grch38p7_prim_assembly.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        print(record.description)
        print(len(record))
        for e in chrom1_snp_pos:
            print(record.seq[e - kmer_size +1 : e + kmer_size])
        break
"""