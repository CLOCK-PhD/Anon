#!/usr/bin/python3

# Extract k-mers around SNP in rs_fasta files from dbSNP
# ftp.ncbi.nih.gov/snp/organisms/human_9606/rs_fasta/

# Positions des éléments dans les entêtes :
# 0 = gnl
# 1 = dbSNP
# 2 = rs
# 3 = position du SNP dans la séquence
# 4 = longueur de la séquence

# structure de k-mers voulue : Séquence, ID, chromosome, position du SNP, position DU KMER

# A FAIRE : GERER LES KMERMAX < KMERSIZE

from pprint import pprint
from Bio import SeqIO
import re
from matplotlib.pyplot import close

from numpy import equal

input_file = "../data/example.fas"
#input_file = "../data/rs_ch1.fas"

# Extraire les variants possibles :
def extract_snp_var(seq_info):
    res = re.search("alleles=\"./(.*)\"$", seq_info)
    if res:
        snp_var = res.group(1)
    else:
        # Cas où on a un résultat non gérable
        # exemple : alleles="(A)8/10"
        # GERER CE CAS EN AVAL
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
    if len(max_kmer) == 2*kmer_size - 1:
        for snp in snp_var:
            if len(snp) == 1:   # Cas où le SNP est de taille 1 et qu'on a un kmermax
                max_kmer_var = max_kmer[:kmer_size - 1] + snp + max_kmer[kmer_size:]
                max_kmers_list.append(max_kmer_var)
            else : # cas du snp long
                max_kmer_var = max_kmer[len(snp)-1:kmer_size-1] + snp + max_kmer[kmer_size:len(max_kmer)-(len(snp)-1)]
                max_kmers_list.append(max_kmer_var)
    else :
        # il est nécessaire de noter la position du snp dans ces kmers
        print("bite")
        for snp in snp_var:
            if len(snp) == 1 :
                print("kmermax réduit, snp taille 1")
                max_kmer_var = max_kmer[:snp_pos-1] + snp + max_kmer[kmer_size - snp_pos:]
                max_kmers_list.append(max_kmer_var)
            else:
                print("kmermax réduit, snp taille >1")
    return max_kmers_list

# Générer les kmers des snp  
# Fonctionne pour les kmers idéaux + kmer snp_pos<kmer_size
def kmer_generator(kmer_size, kmer_to_cut):
    kmer_list = []
    for i in range(0, kmer_size, 1):
        kmer = kmer_to_cut[i : i + kmer_size]
        if len(kmer) == kmer_size: # pour garder des kmer de taille voulue
            kmer_list.append(kmer)

    pprint(kmer_list)
    return kmer_list

kmer_size = 11
# GERER LES SNP_POS < KMER_SIZE !!!!!
# taille de la snp_pos mini : 2
# COUILLE : quand on a snp_var = 0

count = 0

"""
# Vérifications terminées
#origin_snp = ["N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]
# refaire sans ATCG et enregistrer les sorties pour ces nt
# ,"A", "T", "C", "G"
# Fichier output
#outputFile = open("../data/weird_fucks.txt", "w")
"""

# Extraire les kmers à partir des données rs_fasta
with open(input_file) as handle :
    for record in SeqIO.parse(handle, "fasta"):
        print(record.description)
        
        # Séparer les informations de l'entête :
        seq_info = record.description.split("|")
        #print(seq_info)
        #print(record.seq)
        snp_pos = int(seq_info[3].split("=")[1])
        seq_len = int(seq_info[4].split("=")[1])
        
        seq_snp_var = extract_snp_var(seq_info[8])
        
        print(f"snp_pos : {snp_pos}\t seq_len : {seq_len}\t snp_var : {seq_snp_var}")
        print(f"Le SNP : {record.seq[int(snp_pos)-1]}")
        
        # Sélection du grand k-mer à découper en k-mer de taille voulue
        max_kmer = make_max_kmer(kmer_size, record.seq, snp_pos, seq_len)
        print(max_kmer)
        print(len(max_kmer))
        print(max_kmer[kmer_size-1]) # Afficher le SNP dans le kmer_max
        max_kmers_list = max_kmer_variations(str(max_kmer), snp_pos, seq_snp_var, kmer_size)
        pprint(max_kmers_list)

        # kmer_generator(kmer_size, str(max_kmer))
        # Génération des kmers à partir de la liste des variations de kmermax :
        for var in max_kmers_list:
            kmer_generator(kmer_size, str(var))

        #print(len(record))
        
        count += 1

        """
        # vérification des char utilisés pour les SNP:
        #if record.seq[int(snp_pos)-1] not in origin_snp :
            #print(record.description)
            #print(seq_snp_var)
            #outputFile.write(record.description)
            #outputFile.write("\n")
            #break
        """
#outputFile.close()

print(f"Nombre total de SNP dans le chromosome 1 : {count}")

"""
# Même chose avec une liste des positions pour chercher directement 
# dans la primary assembly par chromosome
# A modifier pour que ça marche avec chaque chromosome
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