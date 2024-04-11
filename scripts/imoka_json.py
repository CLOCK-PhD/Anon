#!/usr/bin/python3

import json
from pprint import pprint

def extract_kmer_info_to_simple_dict(json_file_path):
    # Load the JSON file
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    
    # Prepare a dictionary to hold the extracted information
    kmer_info_simple_dict = {}

    # Process each k-mer in the dataset
    for kmer in data['kmers']:
        # Ensure there are alignments to iterate over
        if 'alignments' in kmer and kmer['alignments'] is not None:
            for alignment in kmer['alignments']:
                # Extracting information
                kmer_sequence = kmer['kmer']
                chromosome = alignment.get('chromosome', 'N/A')
                start = alignment.get('start', 'N/A')
                end = alignment.get('end', 'N/A')
                strand = alignment.get('strand', 'N/A')
                genes = alignment.get('genes', 'N/A')
                
                # The list contains the values in the order: chromosome, start, end, strand, genes
                info_list = [chromosome, start, end, strand, genes]
                
                # Use the k-mer sequence as the key and append the information list to the dictionary
                if kmer_sequence not in kmer_info_simple_dict:
                    kmer_info_simple_dict[kmer_sequence] = info_list

    # Print a comment indicating what is found in the list associated with each k-mer
    print("# Each list contains the following information in order: [chromosome, start, end, strand, genes]")

    # Print the simplified dictionary
    #print(json.dumps(kmer_info_simple_dict, indent=4))
    pprint(kmer_info_simple_dict)

# Specify the path to your JSON file
json_file_path = '../programmes_annexes/iMOKA-master/paper_codes/results/supplementary_data_2/TCGA_OV_aggregated.json'

# Call the function with your file path
extract_kmer_info_to_simple_dict(json_file_path)



#'../programmes_annexes/iMOKA-master/paper_codes/results/supplementary_data_2/TCGA_OV_aggregated.json'