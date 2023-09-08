#!/usr/bin/python3

"""
Programme d'analyse des fichiers vcf de dbSNP.

Analyse la source de chaque fréquence enregistrée, et de corriger l'allélisme pour les SNP disposant de variants qui n'ont pas de fréqunce.
Analyse les SNPs qui sont à proximité des autres pour une taille de k-mer donnée.
Fournit une liste des SNP qui sont à la même position.

Dresse un histogramme des distributions pour chaque fichier, puis pour tous les fichiers.
"""

# A FAIRE
"""
* Gestion des arguments :
    - Dossier des vcf
    - taille des k-mers

* OK : Décompter les SNP di/tri/tétra alléliques
    - OK : dans le fichier original
    - OK : après correction en supprimant les alt qui ont une fréquence nulle

* Extraction des fréquences
    - OK : Récupérer les fréquences et leurs sources
    - Définir comment sélectionner la meilleure fréquence
    - Sélection

* Vérifications des positions des SNP :
    - OK : Positions identiques
        - Faire un dictionnaire qui contient le nom du chromosome et le nombre de snp a la mm position
    - OK : Décompte des SNP à proximité (avant-après) (->liste avec toutes les positions)

* NOUVEAU PROBLÈME :
    Certains SNPs ont les même rsid.
"""

# TACHES OPTIONNELLES :
#   - OK : tqdm : importer le lecteur rapide de fichier et faire la barre de progression.


import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from pprint import pprint
from tqdm import tqdm
from math import ceil

class Snp :

    def __init__(self, chr:str, ref:str, pos:int, rsid:str, alt:list, varClass:str, freq:list) -> None:
        """
        Générateur de l'objet dbSNP.
        Contient les informations contenues dans chaque ligne du fichier vcf.

        Parameters:
        - chrom     (str):          Nom du chromosome
        - snpRef    (str):          SNP de référence
        - pos    (int):          Position du SNP dans le chromosome
        - rsid      (str):          Identifiant du SNP
        - alt    (list(str)):    Liste contenant toutes les variations connues du SNP
        - vc        (str):          Variation Class, le type de variations
        - freq      (list(tuples))  Fréquences : source et fréquences observées
        """

        self._chr = chr
        self._ref = ref
        self._pos = pos
        self._rsid = rsid
        self._alt = alt
        self._varClass = varClass
        self._freq = freq
    
    @property
    def chr(self):
        return self._chr
    @chr.setter
    def chr(self, c:str):
        self._chr = c

    @property
    def ref(self):
        return self._ref
    @ref.setter
    def ref(self, r:str):
        self._ref = r

    @property
    def pos(self):
        return self._pos
    @pos.setter
    def pos(self, p:int):
        self._pos = p

    @property
    def rsid(self):
        return self._rsid
    @rsid.setter
    def rsid(self, id:str):
        self._rsid = id

    @property
    def alt(self):
        return self._alt
    @alt.setter
    def alt(self, a:list):
        self._alt = a

    @property
    def varClass(self):
        return self._varClass
    @varClass.setter
    def varClass(self, vc:str):
        self._varClass = vc

    @property
    def freq(self):
        return self._freq
    @freq.setter
    def freq(self, f:list):
        self._freq = f

class FreqInfo:

    def __init__(self, source:str, freqRef:float, freqAlt1:float, freqAlt2=0.0, freqAlt3=0.0) -> None:
        self._source = source
        self._freqRef = freqRef
        self._freqAlt1 = freqAlt1
        self._freqAlt2 = freqAlt2
        self._freqAlt3 = freqAlt3
        
    @property
    def source(self):
        return self._source
    @source.setter
    def source(self, s:str):
        self._source = s

    @property
    def freqRef(self):
        return self._freqRef
    @freqRef.setter
    def freqRef(self, f):
        self._freqRef = f

    @property
    def freqAlt1(self):
        return self._freqAlt1
    @freqAlt1.setter
    def freqAlt1(self, f):
        self._freqAlt1 = f

    @property
    def freqAlt2(self):
        return self._freqAlt2
    @freqAlt2.setter
    def freqAlt2(self, f):
        self._freqAlt2 = f

    @property
    def freqAlt3(self):
        return self._freqAlt3
    @freqAlt3.setter
    def freqAlt3(self, f):
        self._freqAlt3 = f

    @property
    def positive_frequencies(self):
        pf = 0
        if self._freqAlt1 > 0.0:
            pf += 1
        if self._freqAlt2 > 0.0:
            pf += 1
        if self._freqAlt3 > 0.0:
            pf += 1
        return pf

# Counting lines
def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)

# Correcting Frenquencies
def floatFreq(freq:str)->int:
    if freq == ".":
        f = 0
        return float(f)
    else:
        return float(freq)

# Extract info from vcf file
def getVcfLineInfo(line)-> Snp:
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
        Snp:        chrom(str), snpRef(str), snpPos(int), rsid(str), snpAlt(list(str)), vc(str)
    """
    # Split the vcf file line
    description = line.split("\t")

    # Get chromosome
    chrom = description[0]
    chromRes = re.search("NC_0*(.*)\.", chrom)
    if chromRes:
        chrom = chromRes.group(1)
        if chrom == "24":
            chrom = "Y"
        if chrom == "23":
            chrom = "X"
    # Get SNP position
    snpPos = description[1]

    # Get RSID
    rsid = description[2]

    # Get SNP reference
    snpRef = description[3]

    # Get SNP Alt(s)
    snpAlt = description[4].split(",")

    # Qual and filters are empty in dbSNP vcf file
    #qual = description[5]
    #filter = description[6]

    # Extract INFO
    info = description[7]
    # Get Variant Class and FREQ info
    res = re.search("VC=(\w*).*FREQ=(.*?);", info)
    if res:
        vc = res.group(1)
        frequencies_info = res.group(2).split("|")
    else:
        vc = ""
        print(rsid)

    # Get Frequencies
    #print()
    #print(rsid)
    freqList = []
    for f in frequencies_info :
        res = re.search("(.*):(.*)", f)
        if res :
            sourceName = res.group(1)
            #print(sourceName)
            freqs = res.group(2).split(",")
            #print(freqs)
            #print(len(freqs))
            if len(freqs) == 2:
                freq = FreqInfo(sourceName, floatFreq(freqs[0]), freqAlt1=floatFreq(freqs[1]))
                freqList.append(freq)
            if len(freqs) == 3:
                freq = FreqInfo(sourceName, floatFreq(freqs[0]), floatFreq(freqs[1]), freqAlt2=floatFreq(freqs[2]))
                freqList.append(freq)
            if len(freqs) == 4:
                freq = FreqInfo(sourceName, floatFreq(freqs[0]), floatFreq(freqs[1]), floatFreq(freqs[2]), freqAlt3=floatFreq(freqs[3]))
                freqList.append(freq)
        #print(len(freqList))
        #print(freqList)

    snp = Snp(chrom, snpRef, int(snpPos), rsid, snpAlt, vc, freqList)
    return snp

# Correct SNP allelism by deleting the ALT that doesn't have a frequency
def correctSnpAllelism(snp:Snp)->Snp:
    return snp

def snpVicinityCount(snpPosList:list, kmerSize:int)->dict:
    snpCounts = {}
    # Iterate through SNP positions
    print("Checking SNPs' proximity")
    pbarVicinity = tqdm(total=len(snpPosList))
    for i in range(len(snpPosList)):
        pbarVicinity.update(1)
        count = 0
        # Counting SNPs in next positions
        currentPos = i
        while True:
            currentPos += 1
            try :
                if (snpPosList[currentPos] < snpPosList[i] + kmerSize):
                    count += 1
                else :
                    break
            except IndexError:
                break
        # Counting SNPs in previous positions
        currentPos = i
        while True:
            currentPos -= 1
            try:
                if (snpPosList[currentPos] > snpPosList[i] - kmerSize) and currentPos != -1:
                    count += 1
                else:
                    break
            except IndexError:
                break

        # Update the dictionary
        if count in snpCounts:
            snpCounts[count] += 1
        else:
            snpCounts[count] = 1
    
    pbarVicinity.close()
    return snpCounts

# Create the graph for the SNP proximity Count
def snpVicinityGraph(snp_counts:dict, graphFileName:str):
    plt.figure(figsize=[19.2, 10.8])
    counts = list(snp_counts.keys())
    frequency = list(snp_counts.values())
    total_counts = sum(frequency)
    # Calculate the percentage for each count
    percentages = [(count / total_counts) * 100 for count in frequency]
    # Create colormap
    cmap = get_cmap('plasma')
    # Normalize the data to map it to the colormap
    norm = Normalize(vmin=min(counts), vmax=max(counts))
    bar_width = 0.8 # Set the width of the bars (adjust as needed)
    bars = plt.bar(counts, frequency, color=cmap(norm(counts)), edgecolor='black', linewidth=1, width=bar_width)
    # Set y-axis to a logarithmic scale
    plt.yscale('log')
    # Add percentages as text on top of each bar
    for bar, percentage in zip(bars, percentages):
        font_size = min(12, bar.get_height() * 0.1)  # Adjust 0.3 as needed for font size scaling
        plt.text(bar.get_x() + bar.get_width() / 2, 
                 bar.get_height(), 
                 f'{percentage:.2f}%',
                 ha='center', 
                 va='bottom')
    # Customize the plot
    plt.xlabel("Number of Nearby SNPs")
    plt.ylabel("Frequency (log scale)")
    plt.title("Distribution of nearby SNPs")
    plt.grid(axis='y', linestyle="--", alpha=0.5)
    #plt.show()
    plt.savefig(f"{graphFileName}_SNP_proximity_.png", bbox_inches="tight", dpi=200)
    plt.close()

def main():

    inputDir = "/home/remycosta/phd/Anon/data/snp_latest/"
    vcfFiles = sorted(glob.glob(inputDir + "*.vcf"))
    print(f"Opening directory : {len(vcfFiles)} files found to analyse")
    # Test : 1 file
    #vcfFiles = ["/home/remycosta/phd/Anon/data/snp_latest/snv_common_chr21.vcf"]
    # Test : 2 files
    #vcfFiles = ["/home/remycosta/phd/Anon/data/snp_latest/snv_common_chr21.vcf", "/home/remycosta/phd/Anon/data/snp_latest/snv_common_chr22.vcf"]
    # Pour les tests :
    kmerSize = 21

    # Global records:
    # Allelic type count : original
    global_di_snps = 0                      # Total number of di-allelic SNPs in dbSNP
    global_tri_snps = 0                     # Total number of tri-allelic SNPs in dbSNP
    global_tetra_snps = 0                   # Total number of tetra-allelic SNPs in dbSNP
    # Allelic type count : corrected
    global_corrected_di_snps = 0            # Total number of di-allelic SNPs, after correction
    global_corrected_tri_snps = 0           # Total number of tri-allelic SNPs, after correction
    global_corrected_tetra_snps = 0         # Total number of tetra-allelic SNPs, after correction
    global_corrected_no_freq_alt_snp = 0    # Total number of SNPs with no frequencies recorded, after correction
    global_empty_freq = []                  # List of SNPs with no frequency recorded
    # Proximity
    global_proximity_distribution = {}      # Dictionnary with the number of SNPs being at proximity of the k-mer as keys, and the number of time it was recorded as values
    global_freq_sources = {}                # Dictionary with all the sources as keys and the number of time they appear in dbSNP as values
    sameSnpPosRecord = []                   # List of SNPs' rsids being at the same position

    for inputFile in vcfFiles :
        print("\n---------------------------------------------------------")
        print(f"Reading {inputFile}")
        graphFileName = re.search("_([A-Za-z0-9]*).vcf$", inputFile).group(1)
        # Allelic Type Count
        di_snps = 0                     # Number of di-allelic SNPs
        tri_snps = 0                    # Number of tri-allelic SNPs
        tetra_snps = 0                  # Number of tetra-allelic SNPs
        # Décompte SNP corrigé
        corrected_no_freq_alt_snp = 0   # Number of SNPs with no frequencies recorded, after correction
        corrected_di_snps = 0           # Number of di-allelic SNPs, after correction
        corrected_tri_snps = 0          # Number of tri-allelic SNPs, after correction
        corrected_tetra_snps = 0        # Number of tetra-allelic SNPs, after correction
        # Weird SNPs
        empty_freq = []                 # List of SNPs' rsids with no allelic frequencies
        # Proximité des SNPS
        snpPositions = []               # List of all SNPs positions for proximity analysis
        # Frequences : 
        freq_sources = {}               # Key=sources, Values = number of times the source has been counted

        # Counting number of lines for progression bar
        lineCount = 0
        with open(inputFile, 'rb') as fp:
            c_generator = _count_generator(fp.raw.read)
            # count each \n
            lineCount = sum(buffer.count(b'\n') for buffer in c_generator)

        # Opening file
        pbar = tqdm(total=lineCount)
        with open(inputFile, "r") as f:
            last_snp = Snp("0", "N", 0, "rs0", [], "", [])
            samePosSnpCount = 0
            # Processing file line
            for line in f:
                pbar.update(1)
                # Create the SNP object
                snp = getVcfLineInfo(line)
                #print(snp.rsid, len(snp.alt), snp.alt)
                # Add the snp position to the list snp position list for proximity analysis
                snpPositions.append(snp.pos)
                #print(f"{snp.pos} - {last_snp.pos}")

                # Allelic Count
                if len(snp.alt) == 1:
                    di_snps +=1
                if len(snp.alt) == 2:
                    tri_snps += 1
                if len(snp.alt) == 3:
                    tetra_snps += 1

                # Correction allélisme
                allele_correction = []
                # Check frequencies
                for f in snp.freq:
                    allele_correction.append(f.positive_frequencies)
                    #print("\t", f.source, "\t",  f.freqRef, "\t", f.freqAlt1, "\t", f.freqAlt2,"\t", f.freqAlt3)
                    #print(f"\tAllelism : {len(snp.alt)}\tTrue allelism : {f.positive_frequencies}")
                    try :
                        freq_sources[f.source] += 1
                    except KeyError :
                        freq_sources[f.source] = 1
                
                # Correcting allelism
                corrected_allelism = 0
                for e in allele_correction :
                    corrected_allelism += e
                try :
                    corrected_allelism = ceil(corrected_allelism / len(allele_correction))
                except ZeroDivisionError :
                    empty_freq.append(snp.rsid)
                #print(f"{len(snp.alt)} - {corrected_allelism}")
                #print(f"{len(allele_correction)}")
                # Record corrected allelism
                if corrected_allelism == 0:
                    corrected_no_freq_alt_snp += 1
                elif corrected_allelism == 1:
                    corrected_di_snps += 1
                elif corrected_allelism == 2:
                    corrected_tri_snps += 1
                elif corrected_allelism == 3:
                    corrected_tetra_snps += 1

                # Check if SNP is at the same position than the last
                if snp.pos == last_snp.pos:
                    #print(f"{snp.rsid}:{snp.pos}\t{last_snp.rsid}:{last_snp.pos}")
                    samePosSnpCount += 1
                    sameSnpPosRecord.append(f"{graphFileName}\t{snp.rsid}\t{last_snp.rsid}")
                # Keep last SNP in memory
                last_snp = snp              
        pbar.close()

        # Compter les sources des fréquences
        #pprint(freq_sources)

        # SNP proximity distribution
        snpCounts = snpVicinityCount(snpPositions, kmerSize)
        snpVicinityGraph(snpCounts, graphFileName)

        # Affichage analyses
        print(f"\nAnlysis results for {graphFileName}:")
        print(f"SNP allelism")
        print(f"di\ttri\ttetra\ttotal")
        print(f"{di_snps}\t{tri_snps}\t{tetra_snps}\t{di_snps+tri_snps+tetra_snps}")
        print(f"{corrected_di_snps}\t{corrected_tri_snps}\t{corrected_tetra_snps}\t{corrected_tetra_snps+corrected_di_snps+corrected_tri_snps+corrected_no_freq_alt_snp}")
        print(f"Empty frequencies:")
        pprint(empty_freq)
        print(f"Number of SNPs at the same position: {samePosSnpCount}")

        ###########################
        # ADD DATA TO GLOBAL DATA #
        ###########################
        # SOURCES
        for source, freq in freq_sources.items():
            if source in global_freq_sources:
                global_freq_sources[source] += freq
            else:
                global_freq_sources[source] = freq
        # ALLELIC TYPE
        # Original
        global_di_snps += di_snps
        global_tri_snps += tri_snps
        global_tetra_snps += tetra_snps
        # Corrected
        global_corrected_di_snps += corrected_di_snps
        global_corrected_tri_snps += corrected_tri_snps
        global_corrected_tetra_snps += corrected_tetra_snps
        global_corrected_no_freq_alt_snp += corrected_no_freq_alt_snp
        # Weird SNPS :
        for e in empty_freq:
            global_empty_freq.append(e)
        # PROXIMITÉ
        for key, value in snpCounts.items():
            if key in global_proximity_distribution:
                global_proximity_distribution[key] += value
            else:
                global_proximity_distribution[key] = value

        # SOURCE DISTRIBUTION HISTOGRAM
        categories = list(freq_sources.keys())
        values = list(freq_sources.values())
        plt.figure(figsize=[19.2, 10.8])
        plt.bar(categories, values, color='red', edgecolor='black', linewidth=1.2)
        plt.xlabel('Source project')
        plt.ylabel('Number of SNPs recorded')
        plt.title('SNP distribution by project')
        plt.xticks(rotation=75)
        # Calculate and add percentage values as labels
        for i in range(len(categories)):
            percentage = (values[i] / lineCount) * 100
            plt.text(categories[i], values[i], f'{percentage:.2f}%', ha='center', va='bottom')
        #plt.show()
        plt.grid(axis='y', linestyle="--", alpha=0.5)
        plt.savefig(f"source_distribution_{graphFileName}.png", bbox_inches="tight", dpi=200)
        plt.close()

        # ALLELIC TYPE HISTOGRAM
        original_data = [di_snps, tri_snps, tetra_snps]
        corrected_data = [corrected_di_snps, corrected_tri_snps, corrected_tetra_snps]
        x_labels = ["diallelic", "triallelic", "tetra-allelic"]
        # Create an array of x-values for the bars
        x = np.arange(len(x_labels))
        # Set the width of the bars
        bar_width = 1
        # Create the figure and two subplots
        plt.figure(figsize=[19.2, 10.8])
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
        # Define colors and edgecolors for the bars
        colors = ['purple', 'orange', 'yellow']
        edgecolors = 'black'
        # Plot the "Original" histogram in the first subplot
        ax1.bar(x, original_data, bar_width, label='Original', color=colors, edgecolor=edgecolors)
        # Set x-axis labels and customize the plot
        ax1.set_xticks(x)
        ax1.set_xticklabels(x_labels)
        #ax1.set_xlabel('Allelic Types')
        ax1.set_title('dbSNP original data')
        ax1.yaxis.grid(True, alpha=0.5, linestyle="--")
        # Plot the "Corrected" histogram in the second subplot
        ax2.bar(x, corrected_data, bar_width, label='Corrected', color=colors, edgecolor=edgecolors)
        # Set x-axis labels and customize the plot
        ax2.set_xticks(x)
        ax2.set_xticklabels(x_labels)
        #ax2.set_xlabel('Allelic Types')
        ax2.set_title('Curated version')
        ax2.yaxis.grid(True, alpha=0.5, linestyle="--")
        plt.tight_layout()
        # Show the plot
        #plt.show()
        plt.savefig(f"Allelic_types_correction_{graphFileName}.png", bbox_inches="tight", dpi=200)
        plt.close()

    ###############################
    # AFFICHAGE ANALYSES GLOBALES #
    ###############################
    global_snp_number = global_di_snps + global_tri_snps + global_tetra_snps
    print("\n---------------------")
    print("GLOBAL ANALYSIS")
    print("---------------------")
    print(f"SNP allelism")
    print(f"di\ttri\ttetra\ttotal")
    print(f"{global_di_snps}\t{global_tri_snps}\t{global_tetra_snps}\t{global_di_snps+global_tri_snps+global_tetra_snps}")
    print(f"{global_corrected_di_snps}\t{global_corrected_tri_snps}\t{global_corrected_tetra_snps}\t{global_corrected_tetra_snps+global_corrected_di_snps+global_corrected_tri_snps+global_corrected_no_freq_alt_snp}")
    print(f"Empty frequencies:")
    pprint(global_empty_freq)
    print(f"Total number of SNPs at the same position: {len(sameSnpPosRecord)}")
    with open("same_snp_record.tsv", "a") as f:
        for e in sameSnpPosRecord:
            f.write(f"{e}\n")

    #####################
    # GLOBAL HISTOGRAMS #
    #####################
    # GLOBAL SOURCE DISTRIBUTION HISTOGRAM
    categories = list(global_freq_sources.keys())
    values = list(global_freq_sources.values())
    plt.figure(figsize=[19.2, 10.8])
    plt.bar(categories, values, color='red', edgecolor='black', linewidth=1.2)
    plt.xlabel('Source project')
    plt.ylabel('Number of SNPs recorded')
    plt.title('SNP distribution by project')
    plt.xticks(rotation=75)
    # Calculate and add percentage values as labels
    for i in range(len(categories)):
        percentage = (values[i] / global_snp_number) * 100
        plt.text(categories[i], values[i], f'{percentage:.2f}%', ha='center', va='bottom')
    #plt.show()
    plt.grid(axis='y', linestyle="--", alpha=0.5)
    plt.savefig(f"global_source_distribution.png", bbox_inches="tight", dpi=200)
    plt.close()

    # ALLELIC TYPE DISTRIBUTION HISTOGRAM - LOG SCALE
    original_data = [global_di_snps, global_tri_snps, global_tetra_snps]
    corrected_data = [global_corrected_di_snps, global_corrected_tri_snps, global_corrected_tetra_snps]
    x_labels = ["diallelic", "triallelic", "tetra-allelic"]
    # Create an array of x-values for the bars
    x = np.arange(len(x_labels))
    # Set the width of the bars
    bar_width = 1
    # Create the figure and two subplots
    plt.figure(figsize=[19.2, 10.8])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    # Define colors and edgecolors for the bars
    colors = ['purple', 'orange', 'yellow']
    edgecolors = 'black'
    # Plot the "Original" histogram in the first subplot
    ax1.bar(x, original_data, bar_width, label='Original', color=colors, edgecolor=edgecolors)
    # Set x-axis labels and customize the plot
    ax1.set_xticks(x)
    ax1.set_xticklabels(x_labels)
    #ax1.set_xlabel('Allelic Types')
    ax1.set_title('dbSNP original data')
    ax1.yaxis.grid(True, alpha=0.5, linestyle="--")
    ax1.set_yscale('log')
    # Plot the "Corrected" histogram in the second subplot
    ax2.bar(x, corrected_data, bar_width, label='Corrected', color=colors, edgecolor=edgecolors)
    # Set x-axis labels and customize the plot
    ax2.set_xticks(x)
    ax2.set_xticklabels(x_labels)
    #ax2.set_xlabel('Allelic Types')
    ax2.set_title('Curated version')
    ax2.yaxis.grid(True, alpha=0.5, linestyle="--")
    ax2.set_yscale('log')
    plt.tight_layout()
    # Show the plot
    #plt.show()
    plt.savefig(f"Global_allelic_types_correction_log.png", bbox_inches="tight", dpi=200)
    plt.close()

    # ALLELIC TYPE DISTRIBUTION HISTOGRAM
    original_data = [global_di_snps, global_tri_snps, global_tetra_snps]
    corrected_data = [global_corrected_di_snps, global_corrected_tri_snps, global_corrected_tetra_snps]
    x_labels = ["diallelic", "triallelic", "tetra-allelic"]
    # Create an array of x-values for the bars
    x = np.arange(len(x_labels))
    # Set the width of the bars
    bar_width = 1
    # Create the figure and two subplots
    plt.figure(figsize=[19.2, 10.8])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    # Define colors and edgecolors for the bars
    colors = ['purple', 'orange', 'yellow']
    edgecolors = 'black'
    # Plot the "Original" histogram in the first subplot
    ax1.bar(x, original_data, bar_width, label='Original', color=colors, edgecolor=edgecolors)
    # Set x-axis labels and customize the plot
    ax1.set_xticks(x)
    ax1.set_xticklabels(x_labels)
    #ax1.set_xlabel('Allelic Types')
    ax1.set_title('dbSNP original data')
    ax1.yaxis.grid(True, alpha=0.5, linestyle="--")
    # Plot the "Corrected" histogram in the second subplot
    ax2.bar(x, corrected_data, bar_width, label='Corrected', color=colors, edgecolor=edgecolors)
    # Set x-axis labels and customize the plot
    ax2.set_xticks(x)
    ax2.set_xticklabels(x_labels)
    #ax2.set_xlabel('Allelic Types')
    ax2.set_title('Curated version')
    ax2.yaxis.grid(True, alpha=0.5, linestyle="--")
    plt.tight_layout()
    # Show the plot
    #plt.show()
    plt.savefig(f"Global_allelic_types_correction.png", bbox_inches="tight", dpi=200)
    plt.close()

    # GLOBAL PROXIMITY DISTRIBUTION
    snpVicinityGraph(global_proximity_distribution, "global")


if __name__ == '__main__':
    main()