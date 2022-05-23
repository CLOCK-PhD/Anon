#!/usr/bin/python3

"""
Programme pour effectuer une recherche exacte des k-mers produits par le programme kmer_snp_generator.py
sur une séquence de référence.

EN COURS DE DEVELOPPEMENT

Objectifs :
1. Retrouver les k-mers dans la référence
2. Identifier les k-mers uniques
3. Récupérer les k-mers uniques dans un fichier de sortie
Objectifs supplémentaires :
- Eliminer les k-mers multi-match (définir ce qu'on en fait)
"""

import ahocorasick
import re
import sys
import argparse
import os
import heapq
import linecache
from Bio import SeqIO
from typing import OrderedDict
from pprint import pprint
from os.path import isfile, join

def main():

    # Charger la séquence de référence en mémoire
    # A modifier plus tard pour l'input
    # Test sur le chromosome 22
    seq_file = "../data/grch37p13/NC_000022.10_Homo_sapiens_chromosome_22.fasta"
    kmer_file = "../data/kmer_snp/ch22/final_merge.tsv"
    seq = ""
    with open(seq_file) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq)

    # Aho-Corasick : créer un automaton
    automaton = ahocorasick.Automaton()

    # Utiliser la classe Automaton comme un trie
    # En gros on lui file les mots qu'on veut trouver ici
    """for idx, key in enumerate("AAAAAAAAAAAAAAAAACCTC TTTTTTTTTTTTTTTTTTTGG AAAAAAAAAAAAAAAAAAAAC".split()):
        automaton.add_word(key, (idx, key))"""

    # Vérifier si les mots sont dans l'automaton
    #print("he" in automaton)
    #print("HER" in automaton)

    # Recherche dans l'automaton
    #print(automaton.get("he"))
    #print(automaton.get("she"))

    # Convertir le trie en un automaton Aho-Corasick pour permettre la recherche Aho-Corasick
    #automaton.make_automaton()

    #haystack = "AAAAAAAAAATTTTCCCAAAAATTGCCCCGTTTCAAACCCAACTAGGGACAGCAGTTAGCAT"

    # Rechercher toutes les occurences des clés dans l'input string (haystack)
    """
    Here we print the results and just check that they are correct. 
    The Automaton.iter() method return the results as two-tuples of the end index 
    where a trie key was found in the input string and the associated value for this key.
    Here we had stored as values a tuple with the original string and its trie insertion order
    """
    """for end_index, (insert_order, original_value) in automaton.iter(seq):
        start_index = end_index - len(original_value) + 1
        print((start_index, end_index, (insert_order, original_value)))
        assert seq[start_index:start_index + len(original_value)] == original_value
"""

    #print("le chat")


    # Fichier kmer
    
    # Test sur l'ouverture du fichier k-mer
    count = 0
    with open(kmer_file, "r") as f :
        for line in f:
            if count < 10000000 :
                kmer = line.split("\t")[0]
                #print(f"Ligne {count} - kmer : {kmer}")
                automaton.add_word(kmer, (count, kmer))
                count += 1
            else : break

    automaton.make_automaton()
    for end_index, (insert_order, original_value) in automaton.iter(seq):
        start_index = end_index - len(original_value) + 1
        print((start_index, end_index, (insert_order, original_value)))
        assert seq[start_index:start_index + len(original_value)] == original_value
    #automaton.remove_word(kmer)
    #print(automaton.get_stats())


    #the_file = open(kmer_file, "r")

    

    """for line in the_file:
        print(line.split("\t")[0])"""

if __name__ == "__main__" :
    main()