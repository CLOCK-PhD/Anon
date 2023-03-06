#!/usr/bin/python3

"""
ksg.py : K-mer SNP Generator

Version orientée objet de kmer_snp_gen_index.py

On copie plein de trucs de l'ancien fichier et on voit comment ça se passe ici.

DÉFINITIONS :

    K-MERS :
    Nucléotide de taille k.

    U-mers :
    Les k-mers sont générés à partir de (2k-1)mers ;
    Ils seront nommés u-mer (k=11, u=2*11-1).

    Variants :
    Les variants (de la classe Variant), sont des k-mers identiques à d'autres,
    mais qui proviennent d'un autre SNP.

"""

# A FAIRE : gestions des arguments
# sortir les k-mers in_gen

# objectif : le pouvoir de discrimination des k-mers

# IMPORTANT : Rédiger les objectifs avant de coder.

# IMPORTANT: Penser à faire des readme sur les données utilisées

# 

"""
Notes :
    - Les k-mers contenant un N ne sont pas gardés dans kmerGenerator()
    - La classe Kmer ne sera pas utilisée parce qu'elle complique le tri
    - On va garder la méthode avec le dictionnaire et on aura :
        clé = séquence ; valeur = [Variant]
        Si on retrouve une clé (donc la seq d'un kmer),
        on lui ajoute le "variant"
        Dans le cas contraire, on crée la clé avec son "variant".
    - Vérifier l'histoire du SNP dans le variant :
        je crois que je mets l'original et pas le variant
"""

# 

import re
import os

from variant import Variant
#from kmer import Kmer

from itertools import product
from pprint import pprint
from Bio import SeqIO
from typing import OrderedDict
from tqdm import tqdm

# Récupérer les informations contenues dans le VCF
def getVcfLineInfo(line)-> tuple:
    """Récupère les informations contenues dans chaque ligne du fichier VCF de SNPdb.
    Retourne un tuple qui contient dans l'ordre :
        - chrom     (str):          Nom du chromosome
        - snp_ref   (str):          SNP de référence
        - snp_pos   (int):          Position du SNP dans le chromosome
        - rs_id     (str):          Identifiant du SNP
        - snp_alt   (list(str)):    Liste contenant toutes les variations connues du SNP
        - vc        (str):          Variation Class, le type de variations 

    Parameters:
        line:       Ligne du fichier vcf

    Returns:
        tuple:      chrom(str), snpRef(str), snpPos(int), rsid(str), snpAlt(list(str)), vc(str)
    """
    description = line.split("\t")
    chrom = description[0]
    chromRes = re.search("NC_0*(.*)\.", chrom)
    if chromRes:
        chrom = chromRes.group(1)
        if chrom == "24":
            chrom = "Y"
        if chrom == "23":
            chrom = "X"
    snpPos = description[1]
    rsid = description[2]
    snpRef = description[3]
    snpAlt = description[4].split(",")
    #qual = description[5]
    #filter = description[6]
    info = description[7]
    res = re.search("VC=(\w*)", info)
    if res:
        vc = res.group(1)
    else:
        vc = ""
        print(rsid)
    return chrom, snpRef, int(snpPos), rsid, snpAlt, vc

# Récupérer les u-mers pour SNV :
def get_SNV_umers(sequence:str, snpPos:int, kmerSize:int, snpAlt:list) -> list:
    """Pour un SNV : récupère les u-mers dans une liste, pour chaque variation possible.
    Ces u-mers sont utilisés pour générer tous les k-mers des SNP.

    Parameters:
        sequence    (str):          Séquence où fouiller
        snpPos      (int):          Position du SNP dans la séquence de référence
        kmerSize    (int):          Taille du k-mer
        snpAlt      (list(str)):    Liste contenant les variations du SNP

    Returns:
        umerList    (list):         Liste des u-mers pour chaque variation
    """
    umerList = []
    kmerPos = snpPos - kmerSize + 1
    l_kmer = sequence[snpPos - kmerSize + 1 : snpPos]
    r_kmer = sequence[snpPos+1 : snpPos + kmerSize]
    for alt in snpAlt :
        umer = l_kmer + alt + r_kmer
        # Ici pour supprimer les umers avec N
        umerList.append((umer, kmerPos))
    return umerList

# Récupérer le u-mer pour DEL :
def get_DEL_umers(sequence:str, snpPos:int, kmerSize:int, snpRef:str) -> list:
    """Pour une délétion : récupère le u-mer dans une liste.
    Ce u-mer est utilisé pour générer tous les k-mers du SNP.

    Parameters:
        sequence    (SeqIO):    Séquence où rechercher
        snp_pos     (int):      Position du SNP dans la séquence de référence
        kmer_size   (int):      Taille du k-mer
        snp_ref     (str):      SNP de référence

    Returns:
        umerList    (list):     Liste des u-mers pour chaque variation    
    """
    umerList = []
    kmer_pos = int(snpPos - kmerSize + 1)
    l_kmer = sequence[snpPos - kmerSize + 1 : snpPos]
    snp = sequence[snpPos]
    r_kmer = sequence[snpPos + (len(snpRef)) : snpPos + kmerSize + len(snpRef) -1]
    umer = l_kmer + snp + r_kmer
    umerList.append((umer, kmer_pos))
    return umerList

# Récupérer les u-mers pour INS
def get_INS_umers(sequence:str, snpPos:int, kmerSize:int, snpAlt:list) -> list:
    """Pour une insertion : récupère les u-mers dans une liste, pour chaque variation possible.
    Ces u-mers max sont utilisés pour générer tous les k-mers des SNP.

    Parameters:
        sequence    (str):          Séquence où rechercher
        snpPos      (int):          Position du SNP dans la séquence de référence
        kmerSize    (int):          Taille du k-mer
        snpAlt      (list(str)):    Liste contenant les variations du SNP

    Returns:
        umerList    (list):         Liste des u-mers pour chaque variation
    """
    umerList = []
    for alt in snpAlt:
        if len(alt) < kmerSize:
            kmerPos = snpPos - kmerSize + len(alt)
            l_kmer = sequence[snpPos - kmerSize + len(alt) : snpPos]
            r_kmer = sequence[snpPos + 1 : snpPos + kmerSize - (len(alt)-1)]
            umer = l_kmer + alt + r_kmer
            umerList.append((umer, kmerPos))
        else:
            #print(f"INS TOO LONG FOR KMER SIZE : {len(alt)}")
            return umerList
    return umerList

# Récupérer les u-mers pour MNV
def get_MNV_umers(sequence:str, snpRef:str, snpPos:int, kmerSize:int, snpAlt:list) -> list:
    """Pour un MNV : récupère les u-mers dans une liste, pour chaque variation possible.
    Ces u-mers sont utilisés pour générer tous les k-mers des SNP.

    Parameters:
        sequence    (str):          Séquence où chercher
        snpRef      (str):          SNP de référence
        snpPos      (int):          Position du SNP dans la séquence de référence
        kmerSize    (int):          Taille du k-mer
        snpAlt      (list(str)):    Liste contenant les variations du SNP

    Returns:
        umerList    (list):         Liste des u-mers pour chaque variation
    """
    umerList = []
    if len(snpRef) < kmerSize:
        kmerPos = snpPos - kmerSize + len(snpRef)
        l_kmer = sequence[snpPos - kmerSize + len(snpRef): snpPos]
        r_kmer = sequence[snpPos + len(snpRef) : snpPos + kmerSize]
        for alt in snpAlt:
            umer = l_kmer + alt + r_kmer
            umerList.append((umer, kmerPos))
    else :
        #print(f"MNV TOO LONG FOR KMER SIZE : {len(snp_ref)}")
        return umerList
    return umerList

# Récupérer les u-mers pour INDEL
def get_INDEL_umers(sequence:str, snpRef:str, snpPos:int, snpAlt:list, kmerSize:int) -> list:
    """Pour un INDEL : récupère les u-mers dans une liste, pour chaque variation possible.
    Ces u-mers sont utilisés pour générer tous les k-mers des SNP.

    Parameters:
        sequence    (str):          Séquence où chercher
        snp_ref     (str):          SNP de référence
        snp_pos     (int):          Position du SNP dans la séquence de référence
        snp_alt     (list(str)):    Liste contenant les variations du SNP
        kmer_size   (int):          Taille du k-mer

    Returns:
        umerList    (list):         Liste des u-mers pour chaque variation
    """
    umerList = []
    if len(snpRef) > kmerSize :
        #print("INDEL REF IS TO LONG FOR K-MER SIZE")
        return umerList
    if len(snpRef) == 1:
        for alt in snpAlt:
            if len(alt) == len(snpRef):
                umerList += get_SNV_umers(sequence, snpPos, kmerSize, [alt])
            else :
                umerList += get_INS_umers(sequence, snpPos, kmerSize, [alt])
    else :
        for alt in snpAlt:
            if len(alt) == 1 :
                umerList += get_INS_umers(sequence, snpPos, kmerSize, [alt])
            elif len(snpRef) < len(alt) :
                umerList += get_INS_umers(sequence, snpPos, kmerSize, [alt])
            elif len(snp_ref) > len(alt):
                max_kmer_list += get_INS_umers(sequence, snpPos, kmerSize, [alt])
    return max_kmer_list

# Extraction des u-mers
def getUmerFromPos(sequence:str, pos:int, variantClass:str, kmerSize:int, snpRef:str, snpAlt:list) -> list:
    """Génère le u-mer en fonction de la classe de variant en appelant une fonction spécifique.
    Retourne une liste contenant les umers.

    Parameters:
        sequence        (str):          Séquence où chercher
        pos             (int):          Position du SNP dans la séquence de référence
        variantClass    (str):          Classe de variation pour l'appel à la fonction de génération de umer
        kmerSize        (int):          Taille du k-mer
        snpRef          (str):          SNP de référence
        snpAlt          (list(str)):    Liste contenant les variations du SNP

    Returns:
        umerList        (list):         Liste des u-mers pour chaque variation.

    """
    snpPos = pos - 1
    umerList = []

    # Exclusion des kmers contenant des points
    if snpRef == ".":
        return umerList
    for alt in snpAlt:
        if alt == ".":
            return umerList

    if variantClass == "SNV":
        return get_SNV_umers(sequence, snpPos, kmerSize, snpAlt)
    elif variantClass == "DEL":
        return get_DEL_umers(sequence, snpPos, kmerSize, snpRef)
    elif variantClass == "INDEL":
        return get_INDEL_umers(sequence, snpRef, snpPos, snpAlt, kmerSize)
    elif variantClass == "INS":
        return get_INS_umers(sequence, snpPos, kmerSize, snpAlt)
    elif variantClass == "MNV":
        return get_MNV_umers(sequence, snpRef, snpPos, kmerSize, snpAlt)

# Découper le u-mers pour récupérer les kmers et leurs positions
def kmerGenerator(kmerSize:int, umerToCut:str) -> list:
    """Génère les k-mers de taille voulue à partir des u-mers.
    Retourne une liste contenant tous les k-mers.

    Parameters :
        kmerSize    (int):  Taille du k-mer voulu
        umerToCut   (str):  Séquence du u-mer à découper en k-mers

    Returns:
        kmerList    (list(tuple)): Liste contenant les k-mers et leur position
    """
    #print(kmer_to_cut)
    kmerList = []
    for i in range(0, kmerSize, 1):
        # kmer[0] : séquence ; kmer[1] : position du kmer
        kmer = (umerToCut[0][i : i + kmerSize].upper(), int(umerToCut[1]) + i)
        if len(kmer[0]) == kmerSize and "N" not in kmer[0]: # pour garder des kmer de taille voulue
            kmerList.append(kmer)
    #pprint(f"kmerList : \n{kmerList}")
    return kmerList

def kmerObjGen(kmerSize:int, umerToCut:str, rsid:str, chrom:str, snp:str, snpPos:int):
    """Génère les k-mers de taille voulue au format objet, à partir des u-mers.
    Retourne une liste contenant tous les k-mers.

    Parameters:
        kmerSize    (int):  Taille voulue du k-mer
        umerToCut   (str):  u-mer à découper en k-mer
        rsid        (str):  Identifiant du SNP
        chrom       (str):  Numéro du chromosome
        snp         (str):  Valeur du SNP
        snpPos      (int):  Position du SNP

    Returns:
        kmerList    (list(Variant)):    Liste contenant les k-mers et les variants possibles.

    """
    kmerList = []
    for i in range(0, kmerSize, 1):
        # kmer[0] : séquence ; kmer[1] : position du kmer
        relPos = int(umerToCut[1]) + i
        kmer = umerToCut[0][i : i + kmerSize].upper()
        var = Variant(rsid, chrom, snp, snpPos, relPos, 0, 0)
        if len(kmer) == kmerSize and "N" not in kmer:
            # : pour garder des kmer de taille voulue
            # : pour exclure les kmers avec un "N"
            kmerList.append([kmer, var])
    #pprint(f"kmerList : \n{kmerList}")
    # Ajuster le kmerCount de chaque variant
    for k in kmerList :
        k[1].kmersCount = len(kmerList)
        #print(f"{k[0]} - {k[1].rsid} {k[1].relPos} {k[1].kmersCount} {k[1].ambiguousKmersCount}")
    #pprint(kmerList)
    return kmerList

# Création de l'index des suffixes des k-mers
def createIndexFromDict(d:dict, prefixSize:int, outputDir:str):
    """Génère l'index des k-mers n'ayant qu'un seul variant et qui ne sont pas retrouvés dans le génome.

    Parameters :
        d           (dict): Dictionnaire de k-mers
        prefixSize  (int):  Taille du préfixe pour la création de l'index
        ouputDir    (str):  Nom du fichier de sortie
    """

    print("Creating index...")

    # Créer les suffixes possibles:
    prefixes = product("ACGT", repeat=prefixSize)
    prefList = []
    for pref in prefixes :
        prefList.append(''.join(pref))

    # On va ensuite aller chercher tous les kmers qui ont ce suffixe,
    # les pop pour les mettre dans un nouveau dico de kmers des suffixes,
    # trier le dico, puis l'écrire dans un fichier au nom du préfixe.
    pbar = tqdm(total=len(prefList))
    for pref in prefList:
        pbar.update(1)
        index = {}
        popList = []
        for k in d:
            if k[:5] == pref:
                popList.append(k)
        for pop in popList:
            index[pop[5:]] = d.pop(pop)
        with open(f"{outputDir}/{pref}", "a") as f:
            for suff, var in index.items() :
                line = f"{suff}\t{var[0].variantProperties}\n"
                f.write(line)
    pbar.close()
    print("\tDone.")

# Générer un nom de dossier unique
def uniquify(path:str) -> str:
    """
    Génère un nom de dossier unique pour FileExistsError

    Parameters :
        path    (str):  Nom du chemin du dossier à créer

    Returns :
        path    (str):  Nouveau nom du chemin du dossier
    """
    filename, extension = os.path.splitext(path)
    counter = 1

    while os.path.exists(path):
        path = filename + "_" + str(counter) + extension
        counter += 1

    return path

def main():

    fastaFile = "../data/grch38p13/Y.fasta"
    vcfFile = "../data/snp_latest/Y_common_snv.vcf"
    kmerSize = 31
    kmers_per_file = 100000

    # Création des variables
    kmers = {}                              # Dictionnaire contenant les kmers
    dupKmers = {}                           # Dico des kmers multiples
    kmersObj = []                           # Test - Liste des objets Kmer
    file_number = 0                         # Numéro du fichier de sortie
    output_file_list = []                   # Liste des fichiers de kmers à merge
    merged_file_number = 0                  # Numéro du fichier mergé
    merged_file_list_for_final_merge = []   # Liste des fichiers pour le merge final
    prefixSize = 5                          # A modifier plus tard par une fonction qui trouvera la taille idéale du préfixe.
    totalKmers = 0
    outputDirectory = "../data/ksg_test/kmer_snp_index"

    try :
        os.makedirs(outputDirectory)
    except FileExistsError:
        outputDirectory = uniquify(outputDirectory)
        os.makedirs(outputDirectory)
    

    # Lecture du fichier fasta
    seq = []
    with open(fastaFile) as handle :
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq.upper())

    # Lecture du fichier vcf
    vcfLineCount = 0
    with open(vcfFile, "r") as vcf:
        for line in vcf:
            vcfLineCount += 1
    print(f"Nombre de lignes dans le fichier VCF : {vcfLineCount}")

    # Parcours du fichier vcf
    print(f"Creating {kmerSize}-mers dict...")
    lineCount = 0
    pbar = tqdm(total = vcfLineCount)
    with open(vcfFile, "r") as vcf:
            for line in vcf:
                if lineCount < vcfLineCount :
                    # 3.1 Récupération des informations
                    chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = getVcfLineInfo(line)

                    # 3.2 Génération des kmers à partir des umers
                    umerList = getUmerFromPos(seq, snp_pos, vc, kmerSize, snp_ref, snp_alt)
                    for umer in umerList:
                        kmerList = kmerObjGen(kmerSize, umer, rs_id, chrom, snp_ref, snp_pos)
                    lineCount += 1
                    pbar.update(1)

                    totalKmers += len(kmerList)
                    # Pas forcément nécessaire de trier à cette étape
                    #kmerList.sort(key=lambda x:x[0])
                    for kmer in kmerList:
                        # kmer[0] : séquence ; kmer[1] : variants
                        # si la clé existe, on ajoute le variant
                        try :
                            kmers[kmer[0]].append(kmer[1])
                        # sinon, on ajoute la clé avec son variant
                        except KeyError:
                            kmers[kmer[0]] = [kmer[1]]                

                else:
                    break
    pbar.close()
    print(f"\t{len(kmers)} keys added to the dictionnary")
    """print(len(kmersObj))
    kmersObj.sort(key=lambda x: x.sequence)"""

    # Marquage des kmers génomiques
    kmersInDict = len(kmers)
    inGenomeCount = 0
    n_kmers = len(seq) - kmerSize + 1
    print(f"Number of k-mers to analyse : {n_kmers}")
    print(f"Looking for genmics {kmerSize}-mers in the dictionnary...")
    pbar2 = tqdm(total=n_kmers)
    for i in range(n_kmers):
        pbar2.update(1)
        gKmer = seq[i:i+kmerSize]
        if gKmer in kmers :
            in_genome = True
            if kmers[gKmer][-1] != True :
                kmers[gKmer].append(in_genome)
                inGenomeCount += 1
    pbar2.close()


    # PURGE
    # Détecter les k-mers répétés
    # (On ne peut pas modifier un dictionnaire pendant sa lecture)
    dupCount = 0
    dupList = []
    #genKmersCount = 0 # nombre de kmers du génome qu'on retrouve
    for k, v in kmers.items():
        if len(v) > 1:
            dupCount += 1
            dupList.append(str(k))

    print(f"Purge du dictionnaire")
    pbar4 = tqdm(total=len(dupList))
    for k in dupList :
        if k in kmers:
            pbar4.update(1)
            dupKmers[k] = kmers.pop(k)
    pbar4.close()

    print("Sorting dict...")
    kmers = OrderedDict(sorted(kmers.items()))
    print("\tDone.")

    # Afficher le dictionnaires des kmers triés:
    """for k, v in kmers.items():
        print(f"{k}\t{v[0].variantProperties}")"""

    print(f"Total de k-mers générés : {totalKmers}")
    print(f"Nombre de kmers différents dans le dictionnaire : {kmersInDict}")
    print(f"Nombre de k-mers répétés : {dupCount}/{kmersInDict}")
    print(f"Nombres de k-mers dans le génome : {inGenomeCount}")
    print(f"Nombre de kmers conservés : {len(kmers)}")

    # Création de l'index
    createIndexFromDict(kmers, prefixSize, outputDirectory)

if __name__ == '__main__':
    main()