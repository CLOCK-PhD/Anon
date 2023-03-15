# Pouvoir de discrimination des k-mers

ébauche
- ajouter plus de détails et les sources
- refaire en latex quand ce sera plus complet

## Objectif

On souhaite utiliser les k-mers comme des outils d'évaluation de risque du partage des données brutes de séquençages humains.

Les k-mers peuvent être porteurs, ou non, de marqueurs biologiques, notamment les SNP, qui permettent de réidentifier une personne.

Ici, nous développons une approche basée sur les k-mers porteurs de SNP, et donc porteurs d'un biomarqueur, qui une fois détectés sur les *reads* issu des données brutes de séquençage de patient, nous permettraient d'estimer le risque associé au partage d'un tel *read*.

Si un *read* est considéré comme étant "à risque", c'est-à-dire présentant des k-mers porteurs de marqueurs, il sera donc possible de supprimer celui-ci des données brutes afin de contribuer à l'anonymat du patient.

## Introduction

La démocratisation des technologies de séquençage haut débit (HTS) dans la recherche et la santé a permis ces dernières années une explosion des données génomiques et transcriptomiques. 
Si elles ont permis de grandes avancées dans la rechecher et le diagnostique, ces données restent extrêmement sensibles. 
En effet, les données génétiques contiennent des informations qui permettent l'identification des individus, ainsi que des traits ou maladies qui leurs sont associés (ou qui pourront être détéctés à l'avenir avec les avancées des connaissances scientifiques).
Les problématiques associées au partage de ces données sont donc au centre de l'attention des réglementations afin de protéger la vie privée et les droits des patients.

Le partage des données est un élément crucial dans la recherche scientifique, que ce soit pour des raisons de reproductibilité ou encore pour réaliser de nouvelles découvertes permettant de faire avancer les méthodes existantes de diagnostique ou la compréhension des maladies.
Il est donc indispensable de pouvoir partager ces données en respectant les droits des patients et assurer la sécurité de leurs informations génétiques.

Si pour l'instant, l'approche privilégiée consiste à restreindre les accès aux données à des personnes autorisées, ces procédures longues et laborieuses limitent drastiquement l'accès et le partage, les données restant limitées aux seuls laboratoires ou hopitaux disposant des autorisations nécessaires pour pouvoir les stocker et les utiliser.

En Europe, le RGPD a imposé de nouvelles mesures concernant les données génétiques.
Celles-ci sont définies comme "les données à caractère personnel relatives aux caractéristiques génétéiques héréditaires ou acquises d'une personne physique qui donnent des information uniques sur la physiologie ou l'état de santé de cette personne physique, et qui résultent, notamment, d'une analyse d'un échantillon biologique de la personne en question".
Le règlement permet leur partage sous certaines conditions dans un but de recherche scientifique. La pseudonymisation des données (Article 4, paragraphe 5), qui consiste à faire en sorte que les données à caractères personnels "ne puissent pas être attribuées à une personne précise sans avoir recours à des informations complémentaires", est dans le champ du RGPD.

En revanche, les données anonymes "à savoir les informations ne concernant pas une personne physique identifiée ou identifiable", ou les données anonymisées "de telle manière que la personne concernée ne soit pas ou plus identifiable", ne sont pas concernées par le règlement, et de telles informations peuvent être utilisées "à des fins statistiques ou de recherche" (26).
Les données dont auront été extraites toutes les informations pouvant permettre l'identification d'une personne pourront être utilisées.

Or, il existe de nombreux marqueurs biologiques permettant d'identifier une personne, ou à défaut de pouvoir l'identifier de manière unique, de connaître son origine biogéographique, ses liens de parentés et sa généalogie, ses traits phénotypiques, la présence de maladies génétiques, ou encore obtenir des informations sur comment il peut réagir à un traitement (pharmacogénétique).

Dans le but de pouvoir partager les données génétiques issues du séquençage humain, il est donc indispensable de pouvoir définir si elles sont identifiantes, ou alors dans quel cas elles peuvent considérées comme pseudonymisées ou anonymisées.

De nombreuses études ont par le passé démontré qu'un certain nombre de marqueurs génétiques étaient nécessaires pour ré-identifier une personne, le type et la quantité étant dépendante du contexte.
Ainsi, il est nécessaire de pouvoir estimer quelles données brutes seront à risque d'être partagées, et comment elles présentent ces risques.

Ceci aura pour but de mieux identifier les risques et ainsi de permettre une dé-identification des données génétiques de manière efficace et irreversible, tout en conservant l'information nécessaire à leur traitement.

Nous proposons ici une approche basée sur les k-mers pour estimer ce risque de partage des données brutes issues du séquançage humain

[plus de détail sur les k-mers et où ils sont utilisés]

Les k-mers que nous utiliseront dans notre approche seront porteurs de marqueurs et permettront d'identifier les *reads* suceptibles de présenter un risque lors du partage des données.

## Sélection des marqueurs

Afin de créer les k-mers porteurs de marqueurs d'identification, il va être nécessaire d'aller chercher dans les banques de données disponibles les biomarqueurs.

Une fois sélectionnés, la création des k-mers pourra être faite et la construction d'une structure de donnée sera alors nécessaire pour pouvoir par la suite procéder à l'estimation du risque, en sélectionnait les k-mers ayant une pertinence pour à la ré-identification, c'est-à-dire ceux apparaissant de manière unique (ou peu de fois) dans le génome.

### Biomarqueurs et réidentification

[Intro rapide sur les différents marqueurs qu'on pourrait utiliser]
- Plusieurs marqueurs permettent la ré-identification
- STR, SNP, ET
- On va s'intéresser aux SNP qui sont les plus grands représentants du polymorphisme

Les SNP sont des variation (polymorphisme) d'un nucléotide spécifique à un endroit précis du génome, et constitutent la première source de polymorphisme en contribuant à 90% des variations observées dans le génome. 
Ils sont présents partout dans le génome, et peuvent être associés à des maladies et des phénotypes, ce qui en fait des marqueurs de choix pour la ré-identification.

De précédents études ont démontré qu'un ensemble de 45 SNP choisis permettraient avec une erreur de 10^-15 d'identifier la majeure partie de la population mondiale. Trois cents SNO sélectionnés aléatoirement fourniraient quant à eux suffisamment d'information pour identifier de manière unique une personne. Ainsi, le risque de ré-identification augmente avec le nombre de SNP.

### Sélection des SNP

La *Single Nucleotide Polymorphism database* (dbSNP) est la plus grande archive publique contenant les polymorphismes connus.
Cette collection comprend les substitutions de nucléotides simples (SNP), mais également des délétions ou insertions de petite taille (*Deletion Insertion Polymorphisms*, DIP), des insertions d'éléments transposables, ainsi que des variations de répétition de séquence microsatellites.

Les SNP répertoriés dans dbSNP proviennent de nombreuses sources différentes et n'ont pour certains pas d'information relatives à leur fréquence dans la population qui pourrait nous permettre de définir leur pouvoir de discrimination afin de réidentifier un individu.

Il est donc nécessaire d'effectuer un tri sur toutes les données disponibles pour se concentrer sur les SNP qui présentent un risque pour l'anonymat.

Pour cela, nous avons décidé dans un premier temps de nous concentrer sur les seuls SNV (Single Nucleotide Variations) et qui ont été identifiés grâce au projet 1000 génomes (qu'on peut distinguer par la présence du terme "COMMON" dans le fichier vcf).
Ces derniers contiennent une fréquence observée dans la population, et consitute une source fiable. D'autres sources de SNP pourront être ajoutées par la suite.


## Construction des k-mers

Afin de créer les k-mers porteurs de SNP, il va être nécessaire de créer tous les k-mers possibles pouvant porter toutes les variations possibles de chaque SNP
Pour chaque SNP, on aura donc tous les k-mers qui porteront chaque variation du SNP à chaque position dans le k-mer.
Les k-mers porteurs de SNP apparaissant plusieurs fois seront considérés comme des "variants".
Les k-mers variants seront à différencier des k-mers porteurs de SNP uniques, qui sont les porteurs de risque le plus élevé pour la ré-identfication.

Création du k-mer :
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

## Évaluation du risque de ré-identification

Une fois nos données de ré-identification créés, il faudra aller chercher dans les *reads* de nos données brutes s'ils contiennent les k-mers porteurs de SNP créés.

En fonction du nombre de k-mers porteurs de SNP trouvés, et de leurs caractéristiques (qu'ils soient uniques, ou en fonction du nombre de variants qu'ils présentent), on pourra définir les risques qu'un *read* soit suceptible de permettre une ré-identification.

Pour ce faire, on va utiliser les données brutes de séquençage disponibles du projet 1000genomes.