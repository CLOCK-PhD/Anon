#!/usr/bin/python3

"""
Genome K-mer Generator

Programme de test pour générer les k-mers uniques de tout le génome
dans un dictionnaire.
"""

import glob

from sys import getsizeof
from Bio import SeqIO
from tqdm import tqdm
from pprint import pprint

def main():

    kmerSize = 31
    kmersInGenome = {}
    genomeDirectory = "/home/remycosta/phd/Anon/data/grch38p13/"

    # TEST : Lecture de tous les fichiers du dossier
    #print(glob.glob(genomeDirectory + "*.fasta"))

    fastaFiles = glob.glob(genomeDirectory + "*.fasta")
    print(len(fastaFiles))
    pprint(fastaFiles)

    for f in fastaFiles :
        #print(f)
        seq = []
        with open(f) as handle :
            print()
            print(f"Now reading file : {f}")

            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq.upper())

            # Découpage de la séquence en k-mers
            # Nombre de k-mers possibles à partir de la séquence de référence
            n_kmers = len(seq) - kmerSize + 1
            print(f"Sequence length : {len(seq)}")
            print(f"Number of k-mers : {n_kmers}")
            print()
            new_kmers = 0
            existing_kmers = 0

            pbar = tqdm(total=n_kmers)
            for i in range(n_kmers):
                pbar.update(1)
                gKmer = seq[i:i+kmerSize]
                if len(gKmer) == kmerSize and "N" not in gKmer :
                    #print(gKmer)
                    try :
                        kmersInGenome[gKmer].append(i)
                        existing_kmers += 1
                    except :
                        kmersInGenome[gKmer] = [i]
                        new_kmers += 1
            pbar.close()
            print(f"\t{new_kmers} new k-mers added to the dictionnary")
            print(f"\t{existing_kmers} k-mers already existed")
            print(f"There are now {len(kmersInGenome)} entries in the dictionnary.")
            print(f"\tDictionnary size is now : {getsizeof(kmersInGenome)} bytes.")

            pbar.close()

if __name__ == '__main__':
    print("les chats")
    main()