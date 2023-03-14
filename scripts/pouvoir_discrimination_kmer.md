# Pouvoir de discrimination des k-mers

## Introduction

- démocratisation des HTS
- explosion de leur utilisation en santé
- quantités énormes de données disponibles
- nécessité de partager les données
- partage extrêmement réglementé et contrôlé
- risque majeur : ré-identification du patient
- atteindre l'anonymat ou le pseudonymat permet le partage
- Possible de partager les données en cas de dé-identification
- La présence de marqueurs biologiques permet la ré-identification d'une personne
- Il est nécessaire de pouvoir estimer quelles données brutes seront à risques d'être partagées
- On propose une approche basées sur les k-mers pour estimer ce risque de partage des données brutes
- Ces k-mers seront porteurs de marqueurs et permettront d'identifier les *reads* suceptibles de présenter un risque lors du partage des données.

## Objectif

On souhaite utiliser les k-mers comme des outils d'évaluation de risque pour la protection des données paratagées sur les données brutes de séquençage humain.

Les k-mers peuvent être porteurs, ou non, de marqueurs biologiques, notamment les SNP, qui permettent de réidentifier une personne.

Ici, nous déveleppons une approche basée sur les k-mers porteurs de SNP, et donc porteurs d'un biomarqueur, qui une fois détectés sur les *reads*, nous permettraient d'estimer le risque associé au partage d'un tel *read* issu des données brutes de séquençage.

Si un *read* est considéré comme étant "à risque", c'est-à-dire présentant des k-mers porteurs de marqueurs, il sera donc possible de supprimer celui-ci des données brutes afin de préserver l'anonymat du patient.

## Méthode

Afin de créer les k-mers porteurs de marqueurs d'identification, il va être nécessaire d'aller chercher dans les banques de données disponibles ces biomarqueurs.

Une fois sélectionnés, la création des k-mers pourra être faite et la construction d'une structure de donnée sera alors nécessaire pour pouvoir par la suite procéder à l'estimation du risque, en sélectionnait les k-mers ayant une pertinence pour à la ré-identification, c'est-à-dire ceux apparaissant de manière unique (ou peu de fois) dans le génome.


### Sélection des marqueurs

#### Biomarqueurs et réidentification

- Plusieurs marqueurs permettent la ré-identification
- STR, SNP, ET
- On va s'intéresser aux SNP qui sont les plus grands représentants du polymorphisme

Les SNP :
- plus de 90% des variations du génome
- reprendre les descriptions des SNP d'une présentation

#### Sélection des SNP

- dbSNP
- SNV : toutes les entrées de dbSNP ne sont pas des SNP en eux-même mais des marqueurs de polymorphismes de nature différentes. On a décidé de se concentrer sur les SNV en premier lieu.
- COMMON => SNP identifiés grâce au projet 1000genomes
- Pourquoi ? Parce qu'ils contiennent une fréquence observée dans la population, contrairement à beaucoup d'autres SNP inclus dans dbSNP, et qu'ils sont nombreux à provenir de cette source. Elle est donc plus fiable que d'autres. D'autres sources de SNP pourront être ajoutées par la suite.

### Construction des k-mers

On cherche à créer k-mers porteurs de SNP précédemment cités.
On va créer tous les k-mers possibles pouvant porter les variations possibles des SNP.
Pour chaque SNP, on aura donc tous les k-mers qui porteront chaque variation du SNP à chaque position dans le k-mer.
Les k-mers porteurs de SNP apparaissant plusieurs fois seront considérés comme des "variants".
Les k-mers variants seront à différencier des k-mers porteurs de SNP uniques, qui sont les porteurs de risque le plus élevé pour la ré-identfication.

- Sélection du SNP
- recherche du u-mer portant le SNP au milieu de sa séquence
- découpe du u-mer en k-mers
- génère tous les k-mers à partir du u-mer
- les k-mers sont placés dans un dictionnaire python : clé = séquence, valeur = informations relatives au k-mer porteur de SNP

### Fabrication de l'index des SNP

Afin de faciliter le tri et la recherche des k-mers porteurs de SNP crées, ceux-ci seront placés dans un index.

L'index principal contiendra uniquement les k-mers porteurs de SNP qui n'ont pas d'autres variants qu'eux-même, et qu'on ne retrouve pas dans les k-mers constitutifs du génome.

Cet index est un index des préfixes, donc chaque fichier sera les combinaisons possibles des 5 premiers nucléotides du k-mer. Chaque fichier contiendra ensuite les suffixes, ordonnées lexicographiquement, et contentant les informations relatives au SNP à partir duquel ils ont été créés.

Les autres k-mers seront ceux considérés comme provenant de différentes sources, c'est à dire les k-mers non uniques qui sont apparus plus d'une fois, que ce soit lors de la création des k-mers à partir des SNP ou alors comme existant déjà en tant que k-mer constitutif du génome.

On appelera de tels k-mers des k-mers ambigus.Ces k-mers auront une pertinence moindre pour estimer le risque, mais il est important de ne pas tous les écarter.

Les k-mers constitutifs du génome ne portant pas de variation de SNP n'ont pas de pertinence et ne seront pas conservés.

### Évaluation du risque de ré-identification

Une fois nos données de ré-identification créés, il faudra aller chercher dans les *reads* de nos données brutes s'ils contiennent les k-mers porteurs de SNP créés.

En fonction du nombre de k-mers porteurs de SNP trouvés, et de leurs caractéristiques (qu'ils soient uniques, ou en fonction du nombre de variants qu'ils présentent), on pourra définir les risques qu'un *read* soit suceptible de permettre une ré-identification.

Pour ce faire, on va utiliser les données brutes de séquençage disponibles du projet 1000genomes.