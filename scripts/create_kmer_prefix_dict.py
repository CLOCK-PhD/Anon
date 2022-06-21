#!/usr/bin/python3

"""
Créer un dictionnaire qui sert d'index préfixe pour les k-mers d'une séquence :
    clés : préfixe(str)
    valeur : liste qui contient tous les suffixes associés au préfixe (str)

Étapes :
    1. Charger la séquence en mémoire
    2. Parcourir chaque k-mer de la séquance
    3. Pour chaque k-mer : découper en préfixe et suffixe
    4. Remplissage du dictionnaire :
        Si le préfixe n'existe pas, l'ajouter nouvelle clé du dictionnaire
        Ajouter le suffixe dans une liste en tant que valeur (en liste)
        Si le préfixe existe, ajouter le suffixe dans la valeur (sans doublon)
    5. Ordonner le dictionnaire
    6. Ordonner chaque liste
    7. Faire des recherches dans l'index


    A FAIRE :
        - Arrêt de la recherche : faire un stop si on dépasse l'ordre lexico à partir duquel on ne peut plus trouver le mot
            pour éviter d'aller jusqu'à la fin du fichier quand on ne le trouve pas

        Exemple :
                    str1 = "AAACG"
                    str2 = "AAAAA"
                    str1 > str2
                    True
                    str2 > str1
                    False
        
        Donc on peut essayer de mettre un break pour arrêter la recherche dès qu'on dépasse notre SUFFIXE
"""
import sys
from typing import OrderedDict
from Bio import SeqIO
from tqdm import tqdm
from pprint import pprint
from os import listdir
from os.path import isfile, join

def encode(pref:str, prefix_list:list)->int:
    try :
        return prefix_list.index(pref)
    except ValueError :
        return "nope"
    print("le chaton")

def main():

    # 0. Création des variables
    kmers_dict = {}
    ksize = 21
    prefix_size = 5

    # 1. Charger la séquence en mémoire
    seq_file = "../data/grch37p13/NC_000024.9_Homo_sapiens_chromosome_Y.fasta"
    print("Loading sequence...")
    with open(seq_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
    print(f"Sequence length : {len(seq)}")

    # 2. Parcourir les k-mers de la séquence
    # Parcourir les kmers de la séquence
    n_kmers = len(seq) - ksize + 1
    print(f"Number of k-mers to create : {n_kmers}")
    print("Creating prefix dictionnary...")
    pbar = tqdm(total=n_kmers)
    
    # 3. Pour chaque k-mer, découper en préfixe et suffixe (range(n_kmers))
    for i in range(n_kmers):
        pbar.update(1)
        kmer = str(seq[i:i+ksize])
        prefix = kmer[:prefix_size]
        suffix = kmer[prefix_size:]

        # 4. Remplissage du dictionnaire
        try :
            if suffix not in kmers_dict[prefix]:
                kmers_dict[prefix].add(suffix)
        except KeyError:
            kmers_dict[prefix] = {suffix}

    pbar.close()

    # 5. Ordonner les clés du dictionnaire :
    print("Ordering dictionnary...")
    kmers_dict = OrderedDict(sorted(kmers_dict.items()))
    print("\tDone.")
    
    # 6. Ordonner les valeurs du dictionnaire pour chaque clé :
    print("Ordering values...")
    """for value in kmers_dict.values():
        value = list(value)
        value.sort()
        #value = sorted(value)"""
    for key in kmers_dict :
        kmers_dict[key] = list(kmers_dict[key])
        kmers_dict[key].sort()
    print("\tDone.")

    # Parcourir chaque élément de la liste en valeur pour chaque clé
    kmer_count = 0
    for key in kmers_dict:
        for value in kmers_dict[key]:
            kmer_count += 1

    

    # 7. Faire les recherches dans l'index

    # Dossier contenant les fichiers :
    output_dir = "../data/snp_chr_name_test"
    # Lister tous les fichiers :
    kmer_files = [f for f in listdir(output_dir) if isfile(join(output_dir, f))]
    kmer_files.sort()

    """# Ouvrir tous les fichiers - plus nécessaire ?
    f_desc = [open(f"{output_dir}/{filename}", "r") for filename in kmer_files]

    print("Parcours des k-mers")
    for key in kmers_dict:
        id_pref = encode(key, kmer_files)
        if id_pref != "nope":
            for value in kmers_dict[key]:
                suffix = value
                for line in f_desc[id_pref]:
                    if line.split("\t")[0] == suffix:
                        print(f"Trouvé : {key}{value} à la position {len(line)}")"""

    # Test en parcourant les fichiers l'un après l'autre (vu que c'est ordonné)
    print("Recherche dans l'index...")
    not_in_index = []
    count = 0
    pbar2 = tqdm(total=len(kmers_dict))
    for key in kmers_dict:
        pbar2.update(1)
        current_file = f"{output_dir}/{key}"
        reading_position = 0
        try :
            with open(current_file, "r") as f:
                #print(f"\nOuverture du fichier {key}")
                f.seek(reading_position)
                for value in kmers_dict[key]:
                    f.seek(reading_position)
                    curr_suffix = value
                    line = f.readline()
                    if line.split("\t")[0] == curr_suffix :
                            #print(f"\tTrouvé : {key}{value} à la position {f.tell()}")
                            reading_position = f.tell() # Pour enregistrer la position et y retourner directement
                            count += 1
                            #continue
                    while line:
                        #print(f.tell())
                        line = f.readline()
                        if line.split("\t")[0] == curr_suffix :
                            #print(f"\tTrouvé : {key}{value} à la position {f.tell()}")
                            count += 1
                            reading_position = f.tell()
                            #continue
                        #print(line)
        except :
            #print(f"Le fichier {key} n'existe pas")
            not_in_index.append(key)
    pbar2.close()

    # Test en parcourant les fichiers l'un après l'autre (vu que c'est ordonné)
    """for key in kmers_dict:
        current_file = f"{output_dir}/{key}"
        try :
            with open(current_file, "r") as f:
                for value in kmers_dict[value]:
                    curr_suffix = value
                    line = f.readline()
                    if line.split("\t")[0] == curr_suffix :
                            print(f"Trouvé : {key}{value} à la position {f.tell()}")
                    while line:
                        print(f.tell())
                        line = f.readline()
                        print(line)
        except :
            print(f"Le fichier {key} n'existe pas")"""






    # Fermer tous les fichiers
    """for f in f_desc:
        f.close()"""
    
    print(f"Taille de l'index :                         {len(kmer_files)}")
    print(f"Nombre de préfixes trouvés :                {len(kmers_dict)}")
    print(f"Nombre de kmers dans le dictionnaire :      {kmer_count}")
    print(f"Nombre de kmers en commun avec l'index :    {count}")
    print("le chat")

if __name__ == '__main__':
    main()