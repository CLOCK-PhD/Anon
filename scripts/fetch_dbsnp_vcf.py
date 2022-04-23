#!/usr/bin/python3

# Fetch dbsnp VCF file

import sys
import os
from ftplib import FTP
import time
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
link = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz"
wget.download(link)