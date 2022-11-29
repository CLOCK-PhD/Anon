#!/usr/bin/python3

# Fetch dbsnp VCF file

"""
Script pour télécharger la dernière version de dbSNP sur GRCH38
"""

from ftplib import FTP
import wget

# le script ici fonctionne
"""server = "ftp.ncbi.nih.gov"
path = "/snp/organisms/human_9606/VCF/"
filename = "00-All.vcf.gz"

destination = "/home/user/PhD/Anon/download"
ftp = FTP('ftp.ncbi.nih.gov')
ftp.login()
ftp.cwd(path)
ftp.retrlines('LIST')


ftp.retrbinary("RETR " + filename, open(filename, "wb").write)
ftp.quit()"""

# là c'est beaucoup plus court
#link = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz"
link = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz"
wget.download(link)