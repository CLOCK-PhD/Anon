#!/usr/bin/python3
from fileinput import filename
from ftplib import FTP
from pprint import pprint
import os

# Fetch rs_fasta files for all chromosomes from human_9606 in dbSNP
server = "ftp.ensembl.org"
path = "/pub/release-108/fasta/homo_sapiens/dna/"

# Création de l'objet FTP avec le nom du server
ftp = FTP(server)

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