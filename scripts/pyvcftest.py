#!/usr/bin/python3

import vcf
import gzip

# Counting lines
def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)

print("le chat")

#vcf_unzip_path = "/home/remycosta/data/anon/snp/latest/data"
vcf_path = "/home/remycosta/phd/Anon/data/dbsnp156/GCF_000001405.40.gz"
tbi_path = "/home/remycosta/phd/Anon/data/dbsnp156/GCF_000001405.40.gz.tbi"
vcf_reader = vcf.Reader(filename=vcf_path, compressed=True)


lineCount = 0
with gzip.open(vcf_path, 'rb') as fp:
    c_generator = _count_generator(fp.read)
    # count each \n
    lineCount = sum(buffer.count(b'\n') for buffer in c_generator)
print(f"Lines: {lineCount}")


#with gzip.open(vcf_path, 'rb') as f:
#    vcf_reader = vcf.Reader(open(vcf_path, "r"))

#for record in vcf_reader:
    #print(f"Chromosome: {record.CHROM}, ALTS: {record.ALT}, INFO: {record.INFO['RS']}")
    #if "VC" in record.INFO.keys() == "SNV":
        #print(record.INFO["VC"])
    #if "FREQ" not in record.INFO.keys():
        #print(record.INFO['FREQ'])
        #print(f"No freq : {record.ID}")
#chrom = 1

#variants = vcf.reader.fetch(chrom, 10000, 20000)