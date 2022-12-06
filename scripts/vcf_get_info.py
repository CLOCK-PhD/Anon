#!/usr/bin/python3

"""
Programme pour récupérer toutes les informations différentes contenues dans les fichiers vcf de dbSNP
"""

input_file = "../data/snp_latest/rs_10.vcf"

with open(input_file) as f:
    for line in f:
        info = line.split("\t")[7].split(";")
        #sep_info = info.split(";")
        print(info[5:])
