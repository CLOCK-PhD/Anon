#!/bin/bash

source_folder="../data/index_full_genome"
destination_folder="../data/index_lite_redux"

# Crée le dossier s'il n'existe pas
mkdir -p "$destination_folder"

# Parcourir chaque fichier du dossier source
for file in "$source_folder"/*; do
    if [ -f "$file" ]; then
        # Prendre le nom du fichier sans le chemin
        filename=$(basename "$file")

        # Copier les 500 000 premières lignes dans le fichier de destination
        head -n 200000 "$file" > "$destination_folder/$filename"

        echo "Copied $filename"
    fi
done