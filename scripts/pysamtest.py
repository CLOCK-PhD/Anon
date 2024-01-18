#!/usr/bin/python3

import pysam
from pysam import VariantFile

kmerSize = 21

bcf_in = VariantFile("/home/remycosta/phd/Anon/data/dbsnp156/GCF_000001405.40.gz")

fasta_path = "/home/remycosta/phd/Anon/data/grch38p14/GCF_000001405.40_GRCh38.p14_genomic.fna"
fai_path = "/home/remycosta/phd/Anon/data/grch38p14/GCF_000001405.40_GRCh38.p14_genomic.fna.fai"

# Ouvrir le fichier fasta
fasta = pysam.FastaFile(fasta_path, fai_path)

linecount = 0
for rec in bcf_in:
    #print(rec)
    #print(rec.alts)
    #print(rec.info['FREQ'])

    if rec.info["VC"] == "SNV":
        #print(rec.ref)
        #print(rec.pos)
        """if "FREQ" in rec.info.keys():
            linecount += 1"""
        
        # start et end pour choper seulement le nt
        umer = fasta.fetch(rec.chrom, rec.pos-1, rec.pos).upper()
        #print(umer.upper())

        print(f"{rec.ref} - {umer}")
        
print(linecount)

fasta.close()

"""int *len;   // Nécessaire pour faidx_fetch_seq; raison inconnue mais permet de fonctionner.
                int start = position - (kmer_size-1);
                int end = position + (kmer_size-1);
                const char *reg = chromosome_name.c_str();
                char *sequence = faidx_fetch_seq(fai, reg, start, end, len);
                if(!sequence){
                    cerr << "Error: Could not fetch the sequence" << endl;
                    return 1;
                }"""
