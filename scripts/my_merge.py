#!/usr/bin/python3

"""
Programme pour merge les tsv avec suppression des doublons
"""

import heapq
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Merge multiple sorted files')
    parser.add_argument('files', metavar='files', type=str, nargs='+', help='input files')

    args = parser.parse_args()

    files = [open(filename, 'r') for filename in args.files]

    """
    # Sans suppression
    merged = heapq.merge(*files)
    for line in merged:
        sys.stdout.write(line)"""

    # Suppression des doublons
    merged = heapq.merge(*files)
    prev_line = ""
    for line in merged:
        kmer = line.split("\t")[0]
        if kmer == prev_line :
            prev_line = kmer
        else:
            sys.stdout.write(line)
            prev_line = kmer


if __name__ == '__main__':
    main()