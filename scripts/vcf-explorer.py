#!/usr/bin/python3

# VCF Explorer
# Pour parcourir les fichiers VCF et en sortir des infos
# Modifications au fur et à mesure des trucs débiles que j'aurai à faire avec des vcf

from pprint import pprint
import re

snp_per_chromosome_count = {}

with open("../data/snp_latest/data") as f :
#with open("test.vcf") as f :    
    current_start = ""

    for l in f :
        x = l.split("\t")
        #print(x)
        #print("len x : " + str(len(x)))
        if len(x) >= 7:
            if x[0] == current_start:
                #print('yay ' + l)
                snp_per_chromosome_count[current_start] += 1
            else:
                current_start = x[0]
                #print("changement\n" + l)
                snp_per_chromosome_count[current_start] = 1

#pprint(snp_per_chromosome_count)

total_snp_test = 0
for cle, valeur in snp_per_chromosome_count.items():
    print("Chromosome", cle, "\t", valeur)
    total_snp_test += valeur
print(f"Snp total : {total_snp_test}")

total_snp = 0
with open("snp_per_chromosome_count", "w") as f:
    for cle, valeur in snp_per_chromosome_count.items():
        line_to_write = "Chromosome " + cle + "\t" + str(valeur) + "\n"
        f.write(line_to_write)
        total_snp += valeur
    f.write("Total : " + str(total_snp))