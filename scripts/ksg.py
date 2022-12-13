#!/usr/bin/python3

"""
ksg.py : K-mer SNP Generator

Version orientée objet de kmer_snp_gen_index.py

On copie plein de trucs de l'ancien fichier et on voit comment ça se passe ici.

Les k-mers sont générés à partir de (2k-1)mers ; ils seront nommés u-mer (k=11, u=2*11-1).
"""

"""
Notes :
    - Les k-mers contenant un N ne sont pas gardés dans kmerGenerator()
"""
import re

from variant import Variant
from kmer import Kmer

from pprint import pprint
from Bio import SeqIO
from typing import OrderedDict

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

# Extraction des kmer_max
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

# Découper le kmer_max pour récupérer les kmers et leurs positions
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
    kmerList = []
    for i in range(0, kmerSize, 1):
        # kmer[0] : séquence ; kmer[1] : position du kmer
        relPos = int(umerToCut[1]) + i
        kmer = Kmer(umerToCut[0][i : i + kmerSize].upper())
        var = Variant(rsid, chrom, snp, snpPos, relPos, 0, 0)
        if len(kmer.sequence) == kmerSize and "N" not in kmer.sequence: # pour garder des kmer de taille voulue
            kmer.addVariant(var)
            kmerList.append(kmer)
    #pprint(f"kmerList : \n{kmerList}")
    return kmerList

def main():

    fastaFile = "../data/grch38p13/Y.fasta"
    vcfFile = "../data/snp_latest/Y_common_snv.vcf"
    kmerSize = 31
    kmers_per_file = 100000

    # Création des variables
    kmers = {}                              # Dictionnaire contenant les kmers
    kmersObj = []                           # Test - Liste des objets Kmer
    file_number = 0                         # Numéro du fichier de sortie
    output_file_list = []                   # Liste des fichiers de kmers à merge
    merged_file_number = 0                  # Numéro du fichier mergé
    merged_file_list_for_final_merge = []   # Liste des fichiers pour le merge final
    prefix_size = 5                         # A modifier plus tard par une fonction qui trouvera la taille idéale du préfixe.
    
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
    lineCount = 0
    with open(vcfFile, "r") as vcf:
            for line in vcf:
                if lineCount < 5 :
                    # 3.1 Récupération des informations
                    chrom, snp_ref, snp_pos, rs_id, snp_alt, vc = getVcfLineInfo(line)

                    # 3.2 Génération des kmers à partir des umers
                    umerList = getUmerFromPos(seq, snp_pos, vc, kmerSize, snp_ref, snp_alt)
                    for umer in umerList:
                        kmerList = kmerObjGen(kmerSize, umer, rs_id, chrom, snp_ref, snp_pos)
                    lineCount += 1

                    kmerList.sort(key=lambda x:x.sequence)

                    for k in kmerList:
                        if k.sequence in kmersObj 0
                            print("bla")
                        else :
                            kmersObj.append(k)
                

                else:
                    break
    
    print(len(kmersObj))
    kmersObj.sort(key=lambda x: x.sequence)

    for k in kmersObj:
        print(k.sequence)


if __name__ == '__main__':
    main()
    """print("le chat")
    # test pour voir si l'import des classes fonctionne
    var1 = Variant("rs1", "X", "T", 1, 1, 0, 0)
    print(var1.rsid)
    k1 = Kmer("ATCG")
    k1.addVariant(var1)
    k1.show_kmer()"""