#!/usr/bin/python3

import os
from itertools import product
from pprint import pprint
from sys import prefix
import argparse

def genPrefixes(length:int) -> list:
    pref_list = []
    for pref in list(product("ACGT", repeat=length)):
        pref_list.append("".join(pref))
    return pref_list

def make_prefix_files(prefix_list:list, output_dir = ""):
    if output_dir == "" :
        for prefix in prefix_list :
            os.makedirs(f"index/{prefix}")
    else :
        for prefix in prefix_list :
            os.makedirs(f"{output_dir}/{prefix}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--foo', dest="foo" ,action='store_false')
    args = parser.parse_args()
    make_index = args.foo
    
    prefix_list = genPrefixes(5)

    kmer = "ACATCGACCTTAAATTTCACA"
    print(kmer)
    suffixe = kmer[5:]
    prefixe = kmer[:5]
    print(prefixe)
    print(f"{suffixe} : longueur {len(suffixe)}")
    print()
    print(kmer)
    print(f"{prefixe}{suffixe}")

    test = "ACGT\trs_id\tchr1\tsnp_pos\tkmer_pos"
    print(test)
    new_line = suffixe + "\t" + "\t".join(test.split("\t")[1:]) + "\n"
    print(new_line)
    #print("\t".join(test.split("\t")[1:]))

    if make_index:
        print("chaton")
    else :
        print("perritus")

    """output_dir = "bla"
    os.makedirs(output_dir)
    for e in prefix_list:
        open(f"{output_dir}/{e}", "a").close()"""


    