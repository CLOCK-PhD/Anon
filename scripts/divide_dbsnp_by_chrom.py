#!/usr/bin/python3

"""
Divise et découpe le fichier vcf de dbSNP (récupéré avec fetch_dbsnp_vcf.py) en plusieurs fichiers pour faciliter l'utilisation de kmer_snp_gen_index.py

Découpe le fichier vcf en fichiers :
    - header
    - chromosomes
"""

# A FAIRE : LES ARGUMENTS POUR LES FICHIERS D'ENTREE ET DE SORTIE

import re
from tqdm import tqdm

# Conversion du nom du chromosome
def convertChromName(c:str)->str:
    res = re.search("NC_00+([0-9]{1,2}).*$", c)
    if res :
        return res.group(1)
    else :
        return c

input_file = "../data/data"

headerOutputFile = "../data/snp_latest/dbsnp_header.vcf"

chrom_list = []

count = 0
total_count = 0

# Ecrire le header, parcourir tous les chromosomes présents
with open(input_file, "r") as vcf:
    print("Header : ")
    for line in vcf:       
        # Récupérer le header
        if line.startswith("#"):
            #print(line)
            count += 1
            with open(headerOutputFile, "a") as outputFile :
                outputFile.write(line)
        elif line.startswith("NT"):
            break
        else :
            total_count += 1
            count += 1
            chrom = line.split("\t")[0]
            if chrom not in chrom_list :
                print(f"\t{count-1} lignes")
                print(f"Ajout de {chrom}")
                chrom_list.append(chrom)
                count = 0
                count += 1


print(f"{total_count} lignes au total")

print(f"{len(chrom_list)} chromosomes détectés :")

chromListForFiles = []
for e in chrom_list :
    fileName = convertChromName(e)
    fileName = fileName + ".vcf"
    chromListForFiles.append(fileName)
    print(f"\t{e}")

# 2. Découper le fichier par chromosomes :

pbar = tqdm(total=total_count)

with open(input_file, "r") as vcf :
    for line in vcf :
        if line.startswith("NT"):
                break
        else : 
            info = line.split("\t")
            currentChrom = info[0]
            convertChromName(currentChrom)
            info[0] = convertChromName(currentChrom)
            if currentChrom in chrom_list:
                outputFile = "../data/snp_latest/" + info[0] + ".vcf"
                with open(outputFile, "a") as of :
                    of.write("\t".join(info))
        pbar.update(1)

pbar.close()
