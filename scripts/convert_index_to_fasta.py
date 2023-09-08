#!/usr/bin/python3

import os
import argparse

from tqdm import tqdm

def main():

    # Gestion des arguments
    parser = argparse.ArgumentParser()

    # Création des arguments
    parser = argparse.ArgumentParser(description="Convert a vkg.py index to a unique fasta file.")
    parser.add_argument("-i", "--index_directory", dest="indexDir", help="Path to the index directory")
    #parser.add_argument("-o", "--output", dest="output_file", help="Name of the output file")

    args = parser.parse_args()
    folder_path = args.indexDir

    # Path to the folder containing the text files
    #folder_path = "/home/remycosta/phd/Anon/data/vkg_index/"

    # Get a list of all the text files in the folder
    file_list = sorted([file for file in os.listdir(folder_path)])

    # Output file where we'll write the results
    output_file = 'output.fasta'

    pbar = tqdm(total = len(file_list))
    with open(output_file, 'w') as output:
        for file_name in file_list:
            file_path = os.path.join(folder_path, file_name)
            with open(file_path, 'r') as input_file:
                pbar.update(1)
                for line in input_file:
                    elements = line.strip().split('\t')
                    if elements:
                        header = f'>{elements[1]}\n'
                        output_line = f'{file_name}{elements[0]}\n'
                        output.write(header)
                        output.write(output_line)

    pbar.close()
    print('Done! Results written to', output_file)


if __name__ == '__main__':
    main()
