# Récapitulatif des résultats et opérations effectuées

## Contexte

Le programme vkg.py permet de créer l'index des k-mers uniques porteurs des variations connues contenues dans dbSNP.
Ce programme génère tous les k-mers possibles qui portent une variation enregistrée dans dbSNP.
Une fois les k-mers générés, le programme trie les k-mers pour ne garder que les k-mers uniques.
Le programme supprime ensuite tous les k-mers générés qui peuvent apparaître dans le génome de référence.
On obtient ainsi un index de k-mers uniques n'apparaissant pas dans le génome de référence et qui peuvent servir de marqueur d'identifcation.

On souhaite repérer ces k-mers uniques dans les données brutes de séquençage afin de déceler des biomarqueurs susceptibles de permettre la réidentification d'un individu.

Cette recherche se fait en parallèle du programme principal en cours de développement et qui devrait assumer à lui seul toutes les étapes de recherches dans les données brutes de séquençage.



## Création de l'index :

### Données
- Variations : dbSNP (snp_latest), SNV uniquement, COMMON uniquement
- Génome de référence : grch38p13 (autosomes et chromosomes sexuels, primary assembly)
- Taille des k-mers : 21

#### Fichiers locaux des données :
- SNP : `/home/remycosta/phd/Anon/data/snp_latest`

#### Nombre de SNV COMMON :
`wc -l *.vcf`

    889628 10_common_snv.vcf
    847948 11_common_snv.vcf
    814264 12_common_snv.vcf
    579026 13_common_snv.vcf
    537986 14_common_snv.vcf
    522375 15_common_snv.vcf
    602163 16_common_snv.vcf
    535673 17_common_snv.vcf
    479952 18_common_snv.vcf
    454888 19_common_snv.vcf
   1353715 1_common_snv.vcf
    451455 20_common_snv.vcf
    305364 21_common_snv.vcf
    309350 22_common_snv.vcf
   1382556 2_common_snv.vcf
   1135944 3_common_snv.vcf
   1161593 4_common_snv.vcf
   1015784 5_common_snv.vcf
   1034358 6_common_snv.vcf
   1061607 7_common_snv.vcf
    883347 8_common_snv.vcf
    854904 9_common_snv.vcf
    699526 X_common_snv.vcf
    242395 Y_common_snv.vcf
  18155801 total

On a un total de 18 155 801 SNP.

### Résultats
#### Fichiers locaux :
- Index :       `/home/remycosta/phd/Anon/data/vkg_index`
- Index Fasta : `/home/remycosta/phd/Anon/data/test_jellyfish/vkg_index_21mers.fasta`

On dénombre 131 052 736 21-mers uniques porteurs d'une variation générés à partir des 18 155 801 de SNV contenus dans la sélection des SNP de dbSNP pour le génome humain, en excluant les k-mers qui peuvent être présents naturellement dans le génome humain.



## Vérifier que l'index ne contient bien que des k-mers uniques

Le programme vkg.py développé assure normalement de n'obtenir que des k-mers uniques qui n'apparaissent pas dans le génome de référence afin de conférer la qualité de biomarqueur potentiel aux k-mers.

On souhaite néanmoins verifier qu'ils sont bien uniques, et l'on va utiliser le programme jellyfish (https://github.com/gmarcais/Jellyfish) pour l'étape de vérification suplémentaire.

### Méthode
- Utilisation de jellyfish pour le comptage des k-mers du génome.
- Conversion de l'index au format fasta pour l'utiliser avec la commande `query` de jellyfish (`convert_index_to_fasta.py`)

La création d'une matrice de comptage des 21 mers du génome complet est réalisée avec jellyfish.
(Les chromosomes indiqués comme _primary assembly_ dans le fihcier fasta de grch38)

#### Commande pour réaliser le comptage :
`jellyfish count -m 21 -s 100M -t 10 full_genome.fasta`

#### Commande pour réaliser la requête :
`time jellyfish query ../grch38_mer_counts.jf -s ../vkg_index_21mers.fasta > kmer_count.txt`

    real	2m24,908s
    user	2m3,546s
    sys 	0m6,075s

#### Commande pour vérifier qu'aucun k-mer n'a un comptage supérieur à 0 :
`grep -v "0" kmer_count.txt`


### Résultats de la vérification :
#### Fichiers locaux :
- Comptage :    `/home/remycosta/phd/Anon/data/test_jellyfish/grch38_mer_counts.jf`
- Requête :     `/home/remycosta/phd/Anon/data/test_jellyfish/verification_kmers_uniques_index_jellyfish/kmer_count`

La requête des k-mers contenus dans l'index au format fasta sur le comptage de k-mers du génome grch38 a été réalisée.
Aucun k-mer n'a un comptage supérieur à 0, ce qui confirme que l'index ne contient bien que des k-mers uniques.
Pas de modification à apporter au programme, il génère correctement les k-mers uniques.



## Retrouver les k-mers de l'index dans un fichier de données de séquençage brut

On souhaite maintenant compter le nombre de fois que les k-mers porteurs de variations susceptibles d'être identifiantes apparaissent dans les données brutes de séquençage issues des génomes contenus dans le projet 1000genomes.

Pour cela, on va récupérer les données brutes de séquençage d'un génome dans la base de données du projet 1000genomes, au format fastq.
On va réaliser un comptage des k-mers contenus dans le génome avec jellyfish.
On va ensuite effectuer une requête des k-mers de l'index au format fasta sur le comptage issu de jellyfish.
On analysera enfin le résultat de la requête.

### Données :
- Génome : HG00501 (ERR020236_1.fastq)
- Index : vkg_index_21mers.fasta

#### Détail du génome HG00501 :
- Sex:              female
- Populations:      Han Chinese South, East Asian Ancestry (https://www.internationalgenome.org/data-portal/population/CHS)
- Biosample ID:     SAME123362 (https://www.ebi.ac.uk/biosamples/samples/SAME123362)
- Cell line source: HG00501 at Coriell (https://catalog.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=HG00501&Product=CC)

#### Related samples :
- Child :           HG00502 (https://www.internationalgenome.org/data-portal/sample/HG00502)
- Sibling:          HG00512 (https://www.internationalgenome.org/data-portal/sample/HG00512)
- Second Order :    HG00514 (https://www.internationalgenome.org/data-portal/sample/HG00514)
- Sibling :         HG00524 (https://www.internationalgenome.org/data-portal/sample/HG00524)
- Second Order :    HG00526 (https://www.internationalgenome.org/data-portal/sample/HG00526)

#### Remarques :
On dispose donc de génomes de la même famille pour des analyses ultérieures

L'execution de jellyfish avec l'index converti en fichier fasta sur le génome ERR20236 a pris 151 minutes.
Le fichier de sortie est appelé 'the_count'.

### Méthode
- Comptage des k-mers du génome HG00501
- Requête des k-mers de l'index au format fasta sur le fichier de comptage.

#### Commande pour réaliser le comptage :
`jellyfish count -m 21 -s 100M -t 10 ERR020236_1.fastq`

Remarque : Pas sûr pour les options -s et -t. Refaire pour vérification et confirmation des résultats.

#### Commande pour réaliser la requête :
`time jellyfish query mer_counts.jf -s ../test_jellyfish/vkg_index_21mers.fasta > the_count`
    real	151m57,195s
    user	2m46,899s
    sys     16m29,619s


### Résultats
#### Fichiers locaux :
- Comptage :  `/home/remycosta/phd/Anon/data/raw_seq_data/mer_counts.jf`
- Requête : `/home/remycosta/phd/Anon/data/raw_seq_data/the_count`

On dénombre 131 052 578 k-mers analysés ; soit 158 k-mers disparus.

## Analyse 

On souhaite désormais observer quels k-mers ont été retrouvés dans les données de séquençage brutes, soit les résultats de la requête effectuées dans l'étape précédente.

Pour cela, on a développé un programme python : analyse_jellyfish_query.py
Ce programme permet de séparer les k-mers en fonction du nombre de fois où ils ont été détectés, et de les séparer dans des fichiers différents, qui sera utile pour retrouver à quel rsid sont associés les k-mers (information contenue dans l'index original), dans la prochaine étape.
Il réalise également un histogramme permettant de visualiser les résultats.

### Données :
- Requête : `/home/remycosta/phd/Anon/data/raw_seq_data/the_count`

### Méthode :
- Répartition des k-mers dans des fichiers différents selon leur comptage.

### Résultats :
#### Fichiers locaux :
- `/home/remycosta/phd/Anon/data/raw_seq_data/kmers_repartition`

#### Décompte réalisé avec le programme analyse_jellyfish_query.py :

Number of k-mers not detected: 111363116 (84.97590638774004%)
Number of k-mers counted more than 1 time: 19689462 (15.024093612259959%)
Number of k-mers counted [1, 10] times : 19270734 (14.704582156331178%)
Number of k-mers counted [11, 20] times: 268394 (0.20479871826710652%)
Number of k-mers counted [21, 30] times: 57401 (0.04379997774633628%)
Number of k-mers counted [31, 40] times: 25080 (0.019137357221618335%)
Number of k-mers counted [41, 50] times: 14212 (0.010844502425583723%)
Number of k-mers counted [51, 60] times: 9203 (0.007022372348905643%)
Number of k-mers counted [61, 70] times: 6738 (0.005141447885138132%)
Number of k-mers counted [71, 80] times: 4879 (0.0037229332489743163%)
Number of k-mers counted [81, 90] times: 3610 (0.002754619600081427%)
Number of k-mers counted [91, 100] times: 3004 (0.00229220977247773%)
Number of k-mers counted [101, 500] times: 20938 (0.015976793680472275%)
Number of k-mers counted [501, 1000] times: 2606 (0.001988514869200055%)
Number of k-mers counted [1001, 5000] times: 2058 (0.0015703620878026528%)
Number of k-mers counted [5001, 10000] times: 501 (0.00038228931291988777%)
Number of k-mers counted more than 10000 times: 104 (7.935746216301064e-05%)

Environ 85% des k-mers ne sont pas détectés.
14.70% des k-mers apparaissent entre 1 et 10 fois.
On constate que certains k-mers sont énormément présents, pouvant apparaitre plusieurs centaines voire milliers de fois.

Il convient donc de recherche à quels SNP sont associés ces k-mers, et essayer de comprendre d'où peut provenir une telle présence.



## Interprétation des résultats d'analyse et identification des SNP

En recherchant les k-mers issus du comptage de jellyfish dans l'index, certains ne sont pas retrouvés : pourquoi ? Séquence inverse ?

En y regardant de plus près, on peut observer ceci :

`jellyfish query mer_counts.jf "CATTCCATTCCGGATGATTCC"`
CATTCCATTCCGGATGATTCC 19978

Jellyfish compte le k-mer CATTCCATTCCGGATGATTCC 19 978 fois.

Si on vérifie directement le nombre d'occurence de ce k-mer dans le fichier brut de séquençage :
`grep "CATTCCATTCCGGATGATTCC" ERR020236_1.fastq | wc -l`
9865
On peut voir qu'on ne retrouve le k-mer que 9865 fois.

En faisant la vérifiction pour son complément GGAATCATCCGGAATGGAATG :
`grep "GGAATCATCCGGAATGGAATG" ERR020236_1.fastq | wc -l`
10099

Bien que restant dans le même ordre de grandeur, on ne retrouve que 19964, soit 14 k-mers disparus.

En revanche, en regardant dans le fichier mer_count.jf :
`jellyfish query mer_counts.jf CATTCCATTCCGGATGATTCC`
CATTCCATTCCGGATGATTCC 19978

`jellyfish query mer_counts.jf GGAATCATCCGGAATGGAATG`
ATTCCATTCCGGATGATTCC 19978

On voit bien qu'il s'agit du même k-mer.

Ceci étant posé, il faudrait s'assurer de ne pas avoir de comptage de k-mers canoniques.
Il y a un comptage des k-mers canoniques.



# dbSNP 156 et GRCH38p14

18 343 772 SNPs répondent aux même choix de sélection, cependant, cela comprend le génome mitochondrial.
