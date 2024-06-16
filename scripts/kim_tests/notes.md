# Historique des tests mis en place pour la lecture du fichier VCF et du génome de référence pour la génération de k-mers porteurs de variants.


## Création des fichiers de tests

L'objectif était de réaliser des fichiers de petites tailles appropriés pour tester la lecture d'un fichier au format VCF qui renvoie sur un génome de référence, ainsi que le fichier FASTA correspondant à ce génome de référence.

Ce génome fictif est consituté de quatre chromosomes de 250 bases chacun appelés A, B, C et D.

## Ajouts

Énumérer tous les cas qu'on a créé.
Cas à la con : gérer les événements en position 1 (INDEL)

### Fichier VCF des variants

Nous avons réalisé un fichier VCF qui simule le contenu de dbSNP. Le fichier de test vcf `snp_test.vcf` contient 10 variants pour chaque chromosomes, les variants étant placés aléatoirement sur le chromosome. Il reflète la diversité du contenu de dbSNP et respecte les spécifications du format. On a donc les chromosomes qui sont ordonnés, et chaque variant est identifié par identifiant unique. Le champs INFO est formaté selon celui de dbSNP, et le header est le même.

Nous avons donc un total de 40 variants avec leurs informations associées (chr, pos, id, ref, alt, qual, filter, info)

### Fichier FASTA

Les séquences du fichier FASTA ont été réalisée grâce à un petit script Python qui lit le fichier VCF, enregistre la position des nucléotides dans la séquence de référence, et génère une séquence aléatoire qui intègre ces variants à la position correspondante.

Il contient donc quatre chromosomes de 250 bases, dont les positions indiquées par le fichier VCF correspondent aux nucléotides du champs "REF".
Si le VCF indique qu'on a un A en 10ème position, on aura un A en dixième position sur le chromosome ; les autres nucléotides sont attribués aléatoirement.
Si on utilise le programme à nouveau, on aura une autre séquence, mais les nucléotides indiqués par le VCF seront toujours les mêmes.

## Objectif

- Utilisation de vcfpp 
- Récupération des champs d'intérêt pour KIM
- Test pour la lecture du fichier FASTA pour la génération des k-mers

Nous souhaitions observer si l'utilisation de htslib pour la lecture du FASTA pouvait prendre en compte le fait que les chromosomes n'était pas ordonnés comme dans le fichier VCF.
Notre fichier VCF contient les chromosomes ordonnés, et chaque variant est ordonné en fonction de sa position dans le génome de référence.
Normalement, le fichier FASTA a les séquences qui sont également ordonnées. Ici nous avons modifié cet ordre dans le cadre du test.

Ordre du VCF : A, B, C, D ; variants ordonnés par odre croissant en fonction de leur position.
Ordre du FASTA : C, A, B, D

## Programme de test

Nom du programme : `verif.cpp`

### Commande de Compilation

```sh
g++ -std=c++17 -o verif verif.cpp -lhts
```

### Détails

Nous souhaitons utiliser vcfpp plutôt que de htslib/vcf.h pour plus de clarté dans le code ; cette librairie offre une surcouche plus compréhensible et documentée que htslib.
Pour être utilisée, VCFPP nécessite que le vcfpp.h soit présent, et htslib installée, donc la compilation requiert toujours de faire le lien en utilisant `-lhts`

#### Récupéartion des données

VCFPP propose des fonctions plus simples permettant d'éviter la lourdeur du code de htslib.

1. **Initialisation du fichier VCF:**

VCFPP:
```cpp
    // Ouvrir le fichier VCF
    vcfpp::BcfReader vcf(vcf_file.c_str());
    vcfpp::BcfRecord record(vcf.header);
```

HTSLIB:
```cpp
bcf_hdr_t *vcf_header = bcf_hdr_read(vcf_file);
bcf1_t *vcf_record = bcf_init();
```

2. **Parcours du fichier VCF :**

VCFPP:
```cpp
while (vcf.getNextVariant(record)) {
    ...
}
```

HTSLIB:
```cpp
while (bcf_read(vcf_file, vcf_header, vcf_record) >= 0) {
    ...
}
```

3. **Récupération des champs**

Ici, on est dans la boucle `while`.

`vcfpp::BcfRecord`

VCFPP

```cpp
std::string chrom = record.CHROM();     // CHROM
int pos = record.POS();                 // POS
std::string ref_allele = record.REF();  // REF
std::string alt_allele = record.ALT();  // ALT
std::string info = record.allINFO();    // INFO
```

À propos de record.ALT() : c'est une string, donc si on a plusieurs ALT, on aura une string ou chacun est séparé par une virgule. Il faut donc faire un split ou un truc dans le genre pour accéder au ALT voulu et extraire les infos.

Pour record.allINFO(), il existe une fonction spécifique pour retrouver la valeur d'un champs, mais cela nécessite de connaître les champs qui sont dans INFO. Cependant, VCFPP ne propose pas de récupérer les sous-champs de INFO qui sont contenus dans le header. Pour ça, il faut passer par htslib.

C'est pourquoi nous avons implémenté une fonction `inline` à vcfpp.h qui permet de récupérer ces champs.

```cpp
    /** 
     * @brief Extract all INFO field from header
     * 
     * This function iterates through the header records and extracts the IDs
     * of all INFO fields. It returns these IDs as a vector of strings.
     * 
     * @return std::vector<std::string> A vector containing the IDs of all INFO fields
    */
    inline std::vector<std::string> getINFO_IDs() {
        std::vector<std::string> ids;
        ids.reserve(hdr->nhrec);
        for (int i = 0; i < hdr->nhrec; ++i) {
            bcf_hrec_t* hrec = hdr->hrec[i];
            if (hrec->type == BCF_HL_INFO) {
                for (int j = 0; j < hrec->nkeys; ++j) {
                    if (std::string(hrec->keys[j]) == "ID") {
                        ids.push_back(hrec->vals[j]);
                    }
                }
            }
        }
        return ids;
    }
```

HTSLIB

Commande obligatoire pour parcourir le VCF :

```cpp
bcf_unpack(vcf_record, BCF_UN_ALL);
```

Récupération des champs :


```cpp

```


4. **Récupération des champs INFO**


VCFPP
```cpp
std::string info_vc;
record.getINFO("VC", info_vc)
```
On récupère bien le champs, mais seulement s'il existe.

Si le champs n'existe pas, on a l'erreur suivante :
```sh
terminate called after throwing an instance of 'std::runtime_error'
  what():  RS has to be of string type

Abandon (core dumped)
```
C'est parce que RS est un int et pas un string, problème régé.

Pour pseudogène par exemple :

```sh
terminate called after throwing an instance of 'std::runtime_error'
  what():  there is no PSEUDOGENEINFO tag in INFO of this variant.

Abandon (core dumped)
```

Cette erreur est retrouvée pour tous les autres champs s'ils sont absents.
Il est donc nécessaire de lever une exception pour ne pas faire planter le programme.

```cpp
std::string info_vc;
try{
    record.getINFO("VC", info_vc);
} catch (const std::runtime_error& e) {
    std::cerr << "Warning :" << e.what();
}
```


## Modifications apportées à vcfpp.h

### Résumé
1. Extraction des champs INFO du header
2. Gestion des sous-champs de type Flag dans INFO
3. Extraction des allèles dans un vecteur de strings plutôt qu'un string (rawALT), et changement de la fonction "ALT()" pour générer la string à partir du vecteur.

### Extraction des champs INFO du header

```cpp
    /** 
     * @brief Extract all INFO field from header
     * 
     * This function iterates through the header records and extracts the IDs
     * of all INFO fields. It returns these IDs as a vector of strings.
     * 
     * @return std::vector<std::string> A vector containing the IDs of all INFO fields
    */
    inline std::vector<std::string> getINFO_IDs() {
        std::vector<std::string> ids;
        ids.reserve(hdr->nhrec);
        for (int i = 0; i < hdr->nhrec; ++i) {
            bcf_hrec_t* hrec = hdr->hrec[i];
            if (hrec->type == BCF_HL_INFO) {
                for (int j = 0; j < hrec->nkeys; ++j) {
                    if (std::string(hrec->keys[j]) == "ID") {
                        ids.push_back(hrec->vals[j]);
                    }
                }
            }
        }
        return ids;
    }
```

### Gestion des sous-champs de type Flag dans le champs INFO

```cpp
  /**
     * @brief get tag value in INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input is std::string
     * @return bool
     * */
    template<typename T>
    isString<T> getINFO(std::string tag, T & v)
    {
        info = bcf_get_info(header->hdr, line.get(), tag.c_str());
        if(!info) throw std::runtime_error("there is no " + tag + " tag in INFO of this variant.\n");
        // if(info->type == BCF_BT_CHAR) (origin)
        if((info->type == BCF_BT_CHAR) || (info->type == BCF_HT_FLAG)) // Handle Flags
            v = std::string(reinterpret_cast<char *>(info->vptr), info->vptr_len);
        else
            throw std::runtime_error(tag + " has to be of string type\n");
    }
```

Cet ajout permet maintenant d'extraire un sous-champ de INFO qui est de type Flag, et qui n'a généralement pas de valeur associée (même si c'est possible).

### Extraction des allèles dans un vecteur de strings 

vcfpp.h permet d'obtenir tous les allèles alternatifs avec la fonction `record.ALT()`. Cette fonction retourne une string qui donne tous les ALT, séparés par une virgule. Bien que pratique et avec une syntaxe claire, ce format de donnée n'est pas très pratique puisqu'il demande de séparer la chaîne pour en récupérer les ALT.

On a donc ajouté une nouvelle fonction `rawALT()` qui extrait les ALT dans un vecteur de string, ce qui est beaucoup plus facile pour y accéder. Ensuite, on a modifié la fonction `ALT()` pour qu'elle génère la string des allèles alternatifs à partir du vecteur rawALT. On a conservé le nom original de cette fonction.

```cpp
/** @brief return raw ALT alleles as vector of strings */
    inline std::vector<std::string> rawALT() const
    {
        std::vector<std::string> v;
        v.reserve(line->n_allele);
        for(int i = 1; i < line->n_allele; i++)
        {
            v.push_back(line->d.allele[i]);
        }
        return v;
    }

    /** @brief return raw ALT alleles as string */
    inline std::string ALT() const
    {
        std::vector<std::string> v = rawALT();
        std::string s;
        for(size_t i = 0; i < v.size(); i++)
        {
            s += v[i] + ",";
        }
        if(s.length() > 1) s.pop_back();
        return s;
    }
```

Les deux fonctions sont maintenant disponibles.



## Problèmes rencontrés avec htslib

### Extraction des valeurs contenues dans INFO

Pour le cas du sous-champs "VC" de dbSNP, j'ai rencontré le problème suivant sur mon fichier de test:

```sh
VC=INDEL/Project1:0.88,0.08,0.04|Project2:0.85,0.10,0.05
VC=INDEL%Project1:0.91,0.09|Project2:0.89,0.11GENEF:100287106
VC=SNV(Project2:0.94,0.06,0|Project4:0.95,0.5,.
VC=SNVProject2:0.90,0.05,0.05
VC=INS8Project1:0.80,0.20|Project2:0.75,0.25|Project4:0.78,0.22PseudogeneA:100287102
```

Au lieu de simplement extraire la valeur de VC, le programme récupérait le champ suivant en ajoutant un caractère à la place du ";". On a même un cas où deux champs sont récupérés.
On retrouve : "/", "%", "(", "8", "22", ou rien.
Cette erreur n'apparait que si le champs VC n'est pas suivi d'un champ de type flag. Cela peut être dû à mon programme ou je supprimais des caractères invisibles. Toujours est-il que c'est un peu pourri si les caractères invisibles empêche le bon fonctionnement du programme, leur suppression assui.



Pour les "copiers/collers":
```cpp

```

```sh

```

