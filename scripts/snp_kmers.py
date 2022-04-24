#!/usr/bin/python3

# Extract k-mers around SNP in rs_fasta files from dbSNP
# ftp.ncbi.nih.gov/snp/organisms/human_9606/rs_fasta/

# Positions des éléments dans les entêtes :
# 2 = rs
# 3 = position du SNP dans la séquence
# 4 = longueur de la séquence

# Output = tsv : Séquence k-mer, ID, chromosome, position du SNP, position DU KMER
# Problème avec les fichiers rs_fasta : pas de position du snp et du kmer

# A FAIRE : Output des kmers -> 1 fichier pour les kmers de 10k seq
# A FAIRE : LAST : Tri de tas des fichiers (autre programme ?)

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
        #pprint(f"kmer max épurés : {snp_var}")
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
                    #print("kmermax réduit, snp taille 1")
                    max_kmer_var = max_kmer[:snp_pos-1] + snp + max_kmer[kmer_size - snp_pos:]
                    max_kmers_list.append(max_kmer_var)
                else:
                    #print("kmermax réduit, snp taille >1")
                    max_kmer_var = max_kmer[:snp_pos-1] + snp + max_kmer[kmer_size - snp_pos:len(max_kmer)-(len(snp)-1)]
                    max_kmers_list.append(max_kmer_var)
    else :
        #print("Machin anormal")
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

kmer_size = 21
count = 0
kmers = []

# Extraire les kmers à partir des données rs_fasta


with open(input_file) as handle :
    for record in SeqIO.parse(handle, "fasta"):
        #print(record.description)
        
        # Séparer les informations de l'entête :
        seq_info = record.description.split("|")
        seq_rs = seq_info[2].split(" ")[0]
        snp_pos = int(seq_info[3].split("=")[1])
        seq_len = int(seq_info[4].split("=")[1])
        
        seq_snp_var = extract_snp_var(seq_info[8])
        
        #print(f"snp_pos : {snp_pos}\t seq_len : {seq_len}\t snp_var : {seq_snp_var}")
        #print(f"Le SNP : {record.seq[int(snp_pos)-1]}")
        
        # Sélection du grand k-mer à découper en k-mer de taille voulue
        max_kmer = make_max_kmer(kmer_size, record.seq, snp_pos, seq_len)
        #print(max_kmer)
        #print(len(max_kmer))
        #print(max_kmer[kmer_size-1]) # Afficher le SNP dans le kmer_max
        max_kmers_list = max_kmer_variations(str(max_kmer), snp_pos, seq_snp_var, kmer_size)
        #pprint(max_kmers_list)

        # Génération des kmers à partir de la liste des variations de kmermax :
        record_kmer_list = []

        for var in max_kmers_list:
            kmer_list = kmer_generator(kmer_size, str(var))
            for kmer in kmer_list:
                record_kmer_list.append(kmer)

        for kmer in record_kmer_list :
            kmer_id = (kmer, seq_rs)
            kmers.append(kmer_id)
        count += 1

# tri lexicographique :
kmers.sort() # vraiment très dur
pprint(kmers)

print(f"Nombre total de SNP dans le chromosome 1 : {count}")