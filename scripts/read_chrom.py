#!/usr/bin/python3

"""
Extraire des kmers depuis la séquence de référence du chromosome directement, 
à partir des infos vcf du fichier ref snp

Output : génère un dossier contenant des fichiers .tsv d'une taille définie (default = 100k):
##Kmer_seq   rs_id   chromosome  snp_position    kmer_position
"""

# A FAIRE : Intégrer snp_position et kmer_position
# A FAIRE : logs pour les kmers rejetés
# A FAIRE : Barre de progression
# A FAIRE : Docstring
# A FAIRE : Préciser les types des arguments des fonctions
# A FAIRE : s'occuper des kmers avec des N

import re
import sys
import argparse
import os
from Bio import SeqIO
from typing import OrderedDict
from pprint import pprint

# Récupérer les informations contenues dans le VCF
def get_vcf_line_info(line):
    description = line.split("\t")
    chrom = description[0]
    snp_pos = description[1]
    rs_id = description[2]
    snp_ref = description[3]
    snp_alt = description[4].split(",")
    #qual = description[5]
    #filter = description[6]
    info = description[7]
    vc = "pas détecté"
    res = re.search("VC=(\w*)", info)
    if res:
        vc = res.group(1)
    else:
        print(rs_id)
    return chrom, snp_ref, int(snp_pos), rs_id, snp_alt, vc

# Récupérer les kmer_max pour SNV :
def get_SNV_kmer_max(sequence, snp_pos:int, kmer_size:int, snp_alt) -> list:
    kmer_max_list = []
    # Cas classique : 1 nt changé en 1 autre nt
    #seq = sequence[snp_pos - kmer_size + 1: snp_pos + kmer_size]
    #print(f"{seq} : {len(seq)}")
    #print(variant_class)
    l_kmer = sequence[snp_pos - kmer_size + 1 : snp_pos]
    r_kmer = sequence[snp_pos+1 : snp_pos + kmer_size]
    for alt in snp_alt :
        kmer_max = l_kmer + alt + r_kmer
        kmer_max_list.append(kmer_max)
    return kmer_max_list

# Récupérer le kmer_max pour DEL :
def get_DEL_kmer_max(sequence, snp_pos:int, kmer_size:int, snp_ref):
    # Les nt après le premier sont supprimés ; il ne reste que le nt à la position indiquée
    # les régions supprimées peuvent être très longues
    # on ne modifie que la partie r_kmer
    # une seule variation possible
    kmer_max_list = []
    l_kmer = sequence[snp_pos - kmer_size + 1 : snp_pos]
    snp = sequence[snp_pos]
    r_kmer = sequence[snp_pos + (len(snp_ref)) : snp_pos + kmer_size + len(snp_ref) -1]
    """print(f"Debut : {l_kmer} : {len(l_kmer)}")
    print(f"SNP : {sequence[snp_pos]}")
    print(f"Délétion : {snp_ref} - longueur : {len(snp_ref)}")
    print(f"Fin : {r_kmer} : {len(r_kmer)}")"""
    kmer_max = l_kmer + snp + r_kmer
    #print(kmer_max)
    kmer_max_list.append(kmer_max)
    return kmer_max_list

# Récupérer les kmer_max pour INS
# A FAIRE : Ajouter un log pour les cas rejetés
def get_INS_kmer_max(sequence, snp_pos, kmer_size, snp_alt):
    # nt ajoutés après snp_pos
    # la taille peut être longue
    # la taille ne doit pas être plus grande que kmer_size
    # il peut y avoir plusieurs variations
    # Les variations peuvent être de taille différente
    # fontcionne pour les DEL, INS +, DEL -
    kmer_max_list = []
    for alt in snp_alt:
        if len(alt) < kmer_size:
            l_kmer = sequence[snp_pos - kmer_size + len(alt) : snp_pos]
            #snp = sequence[snp_pos]
            r_kmer = sequence[snp_pos + 1 : snp_pos + kmer_size - (len(alt)-1)]
            """print(f"Séquence : {sequence[snp_pos - 20 : snp_pos + 1 + 20]}")
            print(f"Debut : {l_kmer} : {len(l_kmer)}")
            print(f"SNP : {sequence[snp_pos]} - Position : {snp_pos}")
            print(f"Insertion(s) : {alt} - longueur : {len(alt)}")
            print(f"{alt} - longueur : {len(alt)}")
            print(f"Fin : {r_kmer} : {len(r_kmer)}")"""
            kmer_max = l_kmer + alt + r_kmer
            #print(f"{kmer_max} - longueur : {len(kmer_max)}")
            kmer_max_list.append(kmer_max)
        else:
            #print(f"INS TOO LONG FOR KMER SIZE : {len(alt)}")
            return kmer_max_list
    return kmer_max_list

# Récupérer les kmer_max pour MNV
# A FAIRE : Ajouter un log pour les cas rejetés
def get_MNV_kmer_max(sequence, snp_ref, snp_pos, kmer_size:int, snp_alt):
    # plusieurs NT changés par un nombre égal de NT
    # toujours de la même taille
    # il peut y avoir plusieurs variants
    # la taille ne doit pas être plus grande que kmer_size
    # la modif de taille se fait du côté l_kmer
    #print(variant_class)
    kmer_max_list = []
    if len(snp_ref) < kmer_size:
        l_kmer = sequence[snp_pos - kmer_size + len(snp_ref): snp_pos]
        r_kmer = sequence[snp_pos + len(snp_ref) : snp_pos + kmer_size]
        #snp = sequence[snp_pos]
        for alt in snp_alt:
            """print(f"Debut : {l_kmer} : {len(l_kmer)}")
            print(f"SNP : {alt}")
            print(f"Fin : {r_kmer} : {len(r_kmer)}")"""
            kmer_max = l_kmer + alt + r_kmer
            #print(kmer_max)
            kmer_max_list.append(kmer_max)
    else :
        #print(f"MNV TOO LONG FOR KMER SIZE : {len(snp_ref)}")
        return kmer_max_list
    return kmer_max_list

# Récupérer les kmer_max pour INDEL
# A FAIRE : Ajouter un log pour les cas rejetés
def get_INDEL_kmer_max(sequence, snp_ref, snp_pos, snp_alt, kmer_size:int):
    # Pire des cas, on peut tout avoir
    #print(variant_class)
    # La ref peut être plus grande que kmer_size
    # les alt peuvent être plus grands que kmer_size
    # on peut avoir des insertions et des délétions pour l'alt d'un même snp
    # il existe des délétions simples
    # il existe des insertions simples
    # on ne retrouve pas de cas de MNV, et rarement des SNV
    # Cas INS + : len(ref)>1 + len(alt)>len(ref)
    # Cas INS + : peut se prendre comme une insertion depuis le snp_pos en fait
    # Cas DEL - : len(ref)>1 + len(alt)<len(ref)
    # Cas DEL - : Délétion avec ALT à la place de snp_ref dans le kmer_max
    max_kmer_list = []
    if len(snp_ref) > kmer_size :
        #print("INDEL REF IS TO LONG FOR K-MER SIZE")
        return max_kmer_list
    if len(snp_ref) == 1:
        for alt in snp_alt:
            if len(alt) == len(snp_ref):
                #print('\tSNV')
                max_kmer_list += get_SNV_kmer_max(sequence, snp_pos, kmer_size, [alt])
            else :
                #print("\tINS")
                max_kmer_list += get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])
    else :
        for alt in snp_alt:
            if len(alt) == 1 :
                #print("\tDEL")
                #print(f"Séquence : {sequence[snp_pos - 20 : snp_pos + 1 + 20]}")
                max_kmer_list += get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])
            elif len(snp_ref) < len(alt) :
                #print("\tINS +")
                max_kmer_list += get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])
            elif len(snp_ref) > len(alt):
                #print("\tDEL -")
                """l_kmer = sequence[snp_pos - kmer_size + len(alt) : snp_pos]
                r_kmer = sequence[snp_pos + len(snp_ref) : snp_pos + len(l_kmer) + len(snp_ref)]
                print(f"Séquence : {sequence[snp_pos - 20 : snp_pos + 1 + 20]}")
                print(f"Debut : {l_kmer} : {len(l_kmer)}")
                print(f"SNP : {sequence[snp_pos]}")
                print(f"Insertion(s) : {alt} - longueur : {len(alt)}")
                print(f"{alt} - longueur : {len(alt)}")
                print(f"Fin : {r_kmer} : {len(r_kmer)}")
                print(f"kmer_max : {l_kmer}{alt}{r_kmer}")
                # essai avec INS :
                print(f"TEST : {get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])}")"""
                max_kmer_list += get_INS_kmer_max(sequence, snp_pos, kmer_size, [alt])
            """elif len(snp_ref) == len(alt):
                #print("\tMNV")
                print("bite")
                max_kmer_list += get_MNV_kmer_max(sequence, snp_pos, kmer_size, [alt])"""
    return max_kmer_list

# Extraction des kmer_max
def get_kmer_from_pos(sequence, pos:int, variant_class:str, kmer_size:int, snp_ref, snp_alt):
    snp_pos = pos - 1
    kmer_max_list = []

    # Pour essayer de virer les kmers contenant des points
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

# Découper le kmer_max
def kmer_generator(kmer_size:int, kmer_to_cut):
    kmer_list = []
    for i in range(0, kmer_size, 1):
        kmer = kmer_to_cut[i : i + kmer_size].upper()
        if len(kmer) == kmer_size: # pour garder des kmer de taille voulue
            kmer_list.append(kmer)
    return kmer_list

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

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    input_file = args.fasta_file
    kmer_size = int(args.kmer_size)
    kmers_per_file = int(args.kmers_per_output_file)
    output_dir = args.output_dir
    ref = args.ref

    # Créer le fichier de sortie
    os.makedirs(output_dir)

    # Récupérer la séquence du chromosome en mémoire
    seq = []
    with open(input_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = record.seq
            #print(record.description)

    print(f"Longueur de la séquence : {len(seq)}")

    count = 0 # Pour les testss

    # Pour vérifier les différents cas possibles dans les INDELS (partie 1) - OK
    """n_ref1_vs_alt1 = 0      # SNV       A   ->  T
    n_ref_vs_alt1 = 0       # DEL       AGG ->  A
    n_ref1_vs_alt = 0       # INS       G   ->  GT
    n_ref_vs_alt_petit = 0  # mini_del  AGG ->  AG
    n_ref_vs_alt_grand = 0  # ins+      GT  ->  GTT
    n_ref_eq_alt = 0        # MNV       ACG ->  AC
    indel_count = 0"""

    # Dictionnaire contenant les kmers
    kmers = {}

    # Numéro du fichier de sortie :
    file_number = 0

    # Ouverture et parcours du fichier vcf de référence
    with open(ref, "r") as vcf:
        for line in vcf:
            # Boucle pour les tests :
            # ChY : 2375594
            if count <= 2375594:
                chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = get_vcf_line_info(line)
                #print(f"{chrom}\t{snp_ref}\t{snp_pos}\t{rs_id}\t{snp_alt}\t{vc}")
                #print(rs_id)
                #kmer_max = get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt)
                
                # Test pour générer les kmers depuis les kmers_max
                kmer_max_list = get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt)
                for kmer_max in kmer_max_list:
                    kmer_list = kmer_generator(kmer_size, kmer_max)

                # Place les kmers générés dans le dictionnaire :
                # Output = tsv : Séquence k-mer, ID, chromosome, position du SNP, position DU KMER
                """for kmer in kmer_list :
                    kmers[kmer] = (rs_id, chrom)"""

                for kmer in kmer_list :
                    # Ajouter les kmers à la liste des kmers
                    if len(kmers) < kmers_per_file :
                        kmers[kmer] = (rs_id, chrom)
                    # Exporter la liste quand on atteint un nombre de kmers dans la liste :
                    if len(kmers) == kmers_per_file :
                        output_file_name = f"{output_dir}/{str(file_number)}_snp_k{kmer_size}.tsv"
                        kmers = OrderedDict(sorted(kmers.items()))
                        with open(output_file_name, "w", encoding="utf-8") as f:
                            for kmer, values in kmers.items() :
                                line_output = f"{kmer}\t{values[0]}\t{values[1]}\n"
                                f.write(line_output)
                        kmers={}
                        file_number += 1

                # Test général sur les kmer max - OK
                """print(f"{vc} - {rs_id}")
                print(f"Référence : {snp_ref} - Position : {snp_pos}")
                print(f"Variants : {snp_alt}")
                print(get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt))
                print("---------------------------------------------")"""
                #print(kmer_max)
                
                # TEST POUR INS - OK ?
                """if vc == "INS":
                    print("---------------------------------------------")
                    print(f"{vc} - {rs_id}")
                    print(f"Référence : {snp_ref}")
                    print(f"Variants : {snp_alt}")
                    print(get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt))"""

                # Test pour les INDELS - OK ?
                """if vc == "INDEL":
                    print("---------------------------------------------")
                    print(f"{vc} - {rs_id}")
                    print(f"Référence : {snp_ref}")
                    print(f"Variants : {snp_alt}")
                    print(get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt))"""

                # Vérifications des possiblités entre ref et alt - OK
                """if vc == "INDEL":
                    indel_count += 1
                    print(rs_id)
                    if len(snp_ref) == 1:
                        for alt in snp_alt:
                            if len(alt) == len(snp_ref):
                                n_ref1_vs_alt1 += 1
                            else :
                                n_ref1_vs_alt += 1
                    else :
                        for alt in snp_alt:
                            if len(alt) == 1 :
                                n_ref_vs_alt1 += 1
                            elif len(ref) < len(alt) :
                                n_ref_vs_alt_grand += 1
                            elif len(ref) > len(alt):
                                n_ref_vs_alt_petit += 1
                            elif len(ref) == len(alt):
                                n_ref_eq_alt += 0

                    print("---------------------------")
                    #break"""
                count += 1
            else:
                break
    
    #pprint(kmers)
    # Exporter les derniers kmers dans un fichier :        
    output_file_name = f"{output_dir}/{str(file_number)}_snp_k{kmer_size}.tsv"
    kmers = OrderedDict(sorted(kmers.items()))
    with open(output_file_name, "w", encoding="utf-8") as f:
        for kmer, values in kmers.items() :
            line_output = f"{kmer}\t{values[0]}\t{values[1]}\n"
            f.write(line_output)

    # Pour vérifier les différents cas possibles dans les indels (partie 2) - OK
    """print(f"INDELS : {indel_count}")
    print(f"SNV : {n_ref1_vs_alt1}\nDEL : {n_ref_vs_alt1}\nINS : {n_ref1_vs_alt}")
    print(f"Mini-DEL : {n_ref_vs_alt_petit}\nINS+ : {n_ref_vs_alt_grand}\nMNV : {n_ref_eq_alt}")"""

    # Voir tous les VC existantes dans le fichier de référence
    """var_dic = {}
    with open(ref, "r") as vcf:
        for line in vcf:
            a,b,c,d,e, vc = get_vcf_line_info(line)
            var_dic[vc]=1

    pprint(var_dic)"""



if __name__ == '__main__':
    main()