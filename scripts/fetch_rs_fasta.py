#!/usr/bin/python3
from fileinput import filename
from ftplib import FTP
from pprint import pprint
import os

# Fetch rs_fasta files for all chromosomes from human_9606 in dbSNP
server = "ftp.ncbi.nih.gov"
path = "/snp/organisms/human_9606_b151_GRCh37p13/rs_fasta/"

# Création de l'objet FTP avec le nom du server
ftp = FTP('ftp.ncbi.nih.gov')

# Login
ftp.login()

# Accès au chemin du dossier dans le server
ftp.cwd(path)

# Lister les noms 
#ftp.retrlines('NLST')
ftp_files = ftp.nlst()

# Afficher les fichiers
#pprint(ftp_files)

# Télécharger les fichiers
for fileName in ftp_files:
    print(f"Dowloading : {fileName}")
    ftp.retrbinary(f"RETR {fileName}", open(fileName, "wb").write)

ftp.quit()