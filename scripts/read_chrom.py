#!/usr/bin/python3

"""
Extraire des kmers depuis la séquence de référence du chromosome directement, 
à partir des infos vcf du fichier ref snp
"""

# Output = tsv : Séquence k-mer, ID, chromosome, position du SNP, position DU KMER

import re
import sys
import argparse
import os
from Bio import SeqIO
from pprint import pprint

# OK
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


# Récupérer le kmer_max pour SNV :
def get_SNV_kmer_max(sequence, snp_pos:int, kmer_size:int, snp_alt):
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

# Récupérer le kmer_max pour INS
# A FAIRE : Ajouter un log pour les cas rejetés
def get_INS_kmer_max(sequence, snp_pos, kmer_size, snp_alt):
    # nt ajoutés après snp_pos
    # la taille peut être longue
    # la taille ne doit pas être plus grande que kmer_size
    # il peut y avoir plusieurs variations
    # Les variations peuvent être de taille différente
    kmer_max_list = []
    for alt in snp_alt:
        if len(alt) < kmer_size:
            l_kmer = sequence[snp_pos - kmer_size + len(alt) : snp_pos]
            #snp = sequence[snp_pos]
            r_kmer = sequence[snp_pos + 1 : snp_pos + kmer_size - (len(alt)-1)]
            """print(f"Debut : {l_kmer} : {len(l_kmer)}")
            print(f"SNP : {sequence[snp_pos]}")
            print(f"Insertion(s) : {alt} - longueur : {len(alt)}")
            print(f"{alt} - longueur : {len(alt)}")
            print(f"Fin : {r_kmer} : {len(r_kmer)}")"""
            kmer_max = l_kmer + alt + r_kmer
            #print(f"{kmer_max} - longueur : {len(kmer_max)}")
            kmer_max_list.append(kmer_max)
        else:
            print(f"INS TOO LONG FOR KMER SIZE : {len(alt)}")
    return kmer_max_list

# Récupérer le kmer_max pour MNV
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
        print(f"MNV TOO LONG FOR KMER SIZE : {len(snp_ref)}")
    return kmer_max_list


# EN COURS - Extraction des kmer_max
"""
SNV : DONE
DEL : DONE
INDEL : en cours sa mère.
INS : DONE
MNV : OK
Penser à donner les bonnes valeurs aux return
"""
def get_kmer_from_pos(sequence, pos:int, variant_class:str, kmer_size:int, snp_ref, snp_alt):
    snp_pos = pos - 1
    kmer_max_list = []
    if variant_class == "SNV":
        return get_SNV_kmer_max(sequence, snp_pos, kmer_size, snp_alt)
    elif variant_class == "DEL":
        return get_DEL_kmer_max(sequence, snp_pos, kmer_size, snp_ref)
    elif variant_class == "INDEL":
        # Pire des cas, on peut tout avoir
        #print(variant_class)
        # La ref peut être plus grande que kmer_size
        # les alt peuvent être plus grands que kmer_size
        # on peut avoir des insertions et des délétions pour l'alt d'un même snp
        # il existe des délétions simples
        # il existe des insertions simples
        # on ne retrouve pas de cas de MNV ou de SNV
        if len(snp_ref) < kmer_size:
            for alt in snp_alt:
                if len(alt) < kmer_size:
                    print("le champ du possible")
                else :
                    print(f"INDEL ALT TOO LONG FOR KMER SIZE : {len(alt)}")
        else :
            print(f"INDEL SNP TOO LONG FOR KMER SIZE : {len(snp_ref)}")
        return "INDEL"
    elif variant_class == "INS":
        return get_INS_kmer_max(sequence, snp_pos, kmer_size, snp_alt)
    elif variant_class == "MNV":
        return get_MNV_kmer_max(sequence, snp_ref, snp_pos, kmer_size, snp_alt)

def main() :
    # Gestion des arguments
    parser = argparse.ArgumentParser()

    # Création des arguments
    parser = argparse.ArgumentParser(description='Extract kmers at SNP positions from a reference VCF and fasta file')
    parser.add_argument("-i", "--input", dest="fasta_file", help="fasta input file")
    parser.add_argument("-r", "--reference", dest="ref", help="VCF SNP reference file")
    parser.add_argument("-k", "--kmer_size", dest="kmer_size", default=21, help="Select k-mer size")
    #parser.add_argument("-nt","--nucleotide_postition", default=1, dest="nt_pos", help="Position du nucléotide à partir duquel générer le kmer")

    # Récupération des valeurs des arguments et attribution
    args = parser.parse_args()
    input_file = args.fasta_file
    kmer_size = int(args.kmer_size)
    ref = args.ref
    #nt_pos = int(args.nt_pos) - 1

    # Récupérer la séquence du chromosome en mémoire
    seq = []
    with open(input_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = record.seq
            #print(record.description)

    #snp_pos_list = [10002, 10003, 10007, 10008]
    #print(f"Position du SNP : {seq[nt_pos]}")
    #print(f"kmer : {seq[nt_pos - kmer_size//2 : nt_pos + kmer_size//2]}")
    print(f"Longueur de la séquence : {len(seq)}")

    # Récupère le kmer à la position indiquée dans la séquence
    """for i in snp_pos_list:
        snp_pos = i -1
        print(f"Position du snp : {i}")
        print(f"SNP : {seq[snp_pos]}")
        print(f"kmer : {seq[snp_pos - kmer_size//2 : snp_pos + kmer_size//2 + 1]}")"""

    ### EN COURS DE DEV : PACOURS DU VCF
    # lire le vcf ref snp
    count = 0

    """n_ref1_vs_alt1 = 0      # SNV       A   ->  T
    n_ref_vs_alt1 = 0       # DEL       AGG ->  A
    n_ref1_vs_alt = 0       # INS       G   ->  GT
    n_ref_vs_alt_petit = 0  # mini_del  AGG ->  AG
    n_ref_vs_alt_grand = 0  # ins+      GT  ->  GTT
    n_ref_eq_alt = 0        # MNV       ACG ->  AC
    indel_count = 0"""

    with open(ref, "r") as vcf:
        #print("bla")
        for line in vcf:
            # ChY : 2375594
            if count <= 1000:
                chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = get_vcf_line_info(line)
                #print(f"{chrom}\t{snp_ref}\t{snp_pos}\t{rs_id}\t{snp_alt}\t{vc}")
                #print(rs_id)
                #kmer_max = get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt)
                print(vc)
                print(get_kmer_from_pos(seq, snp_pos, vc, kmer_size, snp_ref, snp_alt))
                #print(kmer_max)

                # Vérifications des possiblités entre ref et alt :
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

    """print(f"INDELS : {indel_count}")
    print(f"SNV : {n_ref1_vs_alt1}\nDEL : {n_ref_vs_alt1}\nINS : {n_ref1_vs_alt}")
    print(f"Mini-DEL : {n_ref_vs_alt_petit}\nINS+ : {n_ref_vs_alt_grand}\nMNV : {n_ref_eq_alt}")"""

    """var_dic = {}
    with open(ref, "r") as vcf:
        for line in vcf:
            a,b,c,d,e, vc = get_vcf_line_info(line)
            var_dic[vc]=1

    pprint(var_dic)"""



if __name__ == '__main__':
    main()