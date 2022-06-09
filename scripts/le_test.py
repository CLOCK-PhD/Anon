#!/usr/bin/python3

"""
Scripts de tests en cours
"""

from itertools import product
from math import prod
from pprint import pprint

# Générer un tableau de taille 4^k pour des préfixes de taille k (pref_size) - OK
def gen_prefixes(pref_size=5)->list:
    pref_list = []
    pref_list = [''.join(s) for s in product("ACGT", repeat=pref_size)]
    return pref_list

def encode():
    print("le chaton")

def main():
    print("le chat")
    la_liste = gen_prefixes()
    encode()


if __name__ == '__main__':
    main()