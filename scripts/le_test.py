#!/usr/bin/python3

"""
Scripts de tests en cours
"""

from itertools import product
from math import prod
from operator import index
from pprint import pprint
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
from numpy import identity
from tqdm import tqdm

# Générer un tableau de taille 4^k pour des préfixes de taille k (pref_size) - OK
def gen_prefixes(pref_size=5)->list:
    return [''.join(s) for s in product("ACGT", repeat=pref_size)]

def encode(pref:str, prefix_list:list)->int:
    try :
        return prefix_list.index(pref)
    except ValueError :
        return "nope"
    print("le chaton")

def find_kmer(suff:str, prefix_file):
    print("le chatonnage")

def main():

    # Dossier contenant les fichiers :
    # Pire cas avec l'index déjà purifié où on ne retrouve normalement rien
    #output_dir = "../data/snp_chr_name_test/purified_index"
    # Cas avec l'index non purifié
    output_dir = "../data/snp_chr_name_test"
    # Lister tous les fichiers :
    kmer_files = [f for f in listdir(output_dir) if isfile(join(output_dir, f))]
    kmer_files.sort()
    #print(kmer_files)

    # Charger la séquence en mémoire :
    seq_file = "../data/grch37p13/NC_000024.9_Homo_sapiens_chromosome_Y.fasta"
    print("Loading sequence...")
    with open(seq_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
    print(f"Sequence length : {len(seq)}")

    #la_liste = gen_prefixes()

    #f_desc = [open(filename, "r") for filename in f"{output_dir}/{kmer_files}"]

    f_desc = [open(f"{output_dir}/{filename}", "r") for filename in kmer_files]

    """pref = "TTTAG"
    id_pref = encode(pref, kmer_files)
    print(f"Nom du préfixe : {pref} - Nom du fichier {kmer_files[id_pref]}")
    suff = "AAAAAAATGTAAATCA"

    for line in f_desc[id_pref]:
        print(line)
        if line.split("\t")[0] == suff:
            print("TROUVÉ")
            break

    suff = "AAAATGTTAACTCAGG"
    for line in f_desc[id_pref]:
        print(line)
        if line.split("\t")[0] == suff:
            print("TROUVÉ")
            break"""

    """for line in f_desc[1028] :
        print(line)

    print(kmer_files[1028])"""

    # Parcourir les kmers de la séquence
    kmers_found = 0
    found_kmer_dict = {}
    ksize = 21
    prefix_size = 5
    n_kmers = len(seq) - ksize + 1
    print(f"Number of k-mers to analyse : {n_kmers}")
    print("Searching sequence k-mers in the index...")
    pbar = tqdm(total=n_kmers)
    for i in range(n_kmers):
        #print("=========== NOUVEAU KMER ===============")
        pbar.update(1)
        kmer = str(seq[i:i+ksize])
        prefix = kmer[:prefix_size]
        suffix = kmer[prefix_size:]
        
        id_pref = encode(prefix, kmer_files)

        #print(f"\t{f_desc[id_pref]}")
        #print(f"\t +++++ {id_pref} +++++")

        if id_pref != "nope" :
            #print("\tPREFIXE DANS LA LISTE : RECHERCHE DU SUFFIXE")
            f_desc[id_pref].seek(0)
            for line in f_desc[id_pref]:
                #print(line)
                if line.split("\t")[0] == suffix:
                    #print(f"\tTROUVE !!!!!!!!!!!!!!")
                    try:
                        found_kmer_dict[(prefix+line.split("\t")[0], line.split("\t")[1])].append(i)
                    except KeyError:
                        found_kmer_dict[(prefix+line.split("\t")[0], line.split("\t")[1])] = [i]
                    break

    pbar.close()
    for f in f_desc:
        f.close()

    print(len(found_kmer_dict))

if __name__ == '__main__':
    main()