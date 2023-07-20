#!/usr/bin/python3

# Programme pour retirer des colonnes de l'index généré par ksg
# Crée en attendant de refaire ksg.py, afin de faire des tests

import os
import pandas as pd

sourceFolder = "../data/index_full_genome"
destinationFolder = "../data/index_lite"

# Create the destination folder if it doesn't exist
if not os.path.exists(destinationFolder):
    os.makedirs(destinationFolder)

for filename in os.listdir(sourceFolder):
    filePath = os.path.join(sourceFolder, filename)
    destinationPath = os.path.join(destinationFolder, filename)

    df = pd.read_csv(filePath, delimiter="\t", header=None)
    df = df.drop(df.columns[2:6], axis=1)
    df.to_csv(destinationPath, sep="\t", index=False, header=None)