/*
OK : TROUVER COMMENT DETECTER "COMMON" -> marche pour tout type de champ sans valeur dans INFO avec la fonction checkInfoField
OK : DETECTER SI REF EST DEGEN
OK : SEPARER LES DIFFERENTS PROJETS AVEC LEURS FREQUENCES RESPECTIVES
OK : SELECTIONNER UN SEUL PROJET
EN COURS : EXTRAIRE LE UMER - ATTENTION AUX EXTREMITES DE LA SEQUENCE !!
    OK : SNV
    A FAIRE : DEL
    A FAIRE : INS
    A FAIRE : MNV
    A FAIRE : INDEL
EN COURS:  : GENERER LES ALT-UMER
    OK : SNV + Fréquences nulles
    A FAIRE : DEL
    A FAIRE : INS
    A FAIRE : MNV
    A FAIRE : INDEL
OK : GENERER LES V-KMERS
*/

// ATTENTION : J'ai modifié un peu le VCF pour les tests, et les séquences peuvent ne plus correspondre. Il faudrait les générer à nouveau pour être sur.
// En vrai ça concerne seulement rs5 (A:200) où j'ai changé REF en N.

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include "vcfpp.h"
#include <htslib/faidx.h>

const int K = 21;  // Taille du k-mer
const int UMER_LENGTH = 2 * K - 1;  // Taille du u-mer

// Fonction pour diviser une chaîne en utilisant un délimiteur
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Fonction pour récupérer les fréquences des projets et les stocker dans une map
std::map<std::string, std::vector<double>> get_all_freqs(const std::string& info_freq) {
    std::map<std::string, std::vector<double>> frequencies_map;
    std::vector<std::string> projets = split(info_freq, '|');
    for (const std::string& projet : projets) {
        std::vector<std::string> parts = split(projet, ':');
        if (parts.size() != 2) continue; // Vérifier qu'on a bien un projet et des fréquences

        std::string project_name = parts[0];
        std::vector<std::string> freqs_str = split(parts[1], ',');
        std::vector<double> freqs;

        for (std::string& freq : freqs_str) {
            if (freq == ".") {
                freq = "0";
            }
            freqs.push_back(std::stod(freq));
        }

        frequencies_map[project_name] = freqs;
    }
    return frequencies_map;
}

// Fonction pour obtenir les fréquences d'un projet donné
std::vector<double> get_freq_from_project_name(const std::map<std::string, std::vector<double>>& frequencies_map, const std::string& project_name) {
    auto it = frequencies_map.find(project_name);
    if (it != frequencies_map.end()) {
        return it->second;
    } else {
        std::cerr << "Warning:" << project_name << " doesn't exist for this variant" << std::endl;
        return {};
    }
}

// Function to check if a base is degenerated (not A, C, G, or T)
bool isDegenerate(char base) {
    return (base == 'R' || base == 'Y' || base == 'S' || base == 'W' ||
            base == 'K' || base == 'M' || base == 'B' || base == 'D' ||
            base == 'H' || base == 'V' || base == 'N');
}

// Fonction pour vérifier si un champ sans valeur est présent dans INFO (VCFPP ne propose pas cette option)
bool checkInfoField(const std::string& info, const std::string& field) {
    std::vector<std::string> parts = split(info, ';');
    for (const std::string& part : parts) {
        if (part == field) {
            return true;
        }
    }
    return false;
}


// Fonction pour extraire le u-mer autour de la position POS
std::string extract_snv_umer(const faidx_t* fai, const std::string& chrom, int pos, int& variant_pos_in_umer) {
    int start = pos - K;  // Calcul du début de la séquence u-mer
    int end = pos + K - 2;  // Calcul de la fin de la séquence u-mer

    if (start < 0) {  // Gestion des indices négatifs
        start = 0;
    }

    // Calcul de la longueur réelle de la séquence à extraire
    int len;
    char* seq = faidx_fetch_seq(fai, chrom.c_str(), start, end, &len);
    if (seq == nullptr) {
        std::cerr << "Erreur lors de la récupération de la séquence pour " << chrom << ":" << pos << std::endl;
        return "";
    }

    // Calculer la position du variant dans le u-mer
    variant_pos_in_umer = pos - start; // Position exacte sans ajustement

    std::string umer(seq, len);
    free(seq);
    return umer;
}

//Fonction pour générer les k-mers
std::vector<std::string> generate_kmers(const std::string& umer, int k) {
    std::vector<std::string> kmers;
    for (size_t i = 0; i <= umer.size() - k; ++i) {
        kmers.push_back(umer.substr(i, k));
    }
    return kmers;
}

// Fonction pour générer les alt-umers des SNV, avec correction de l'allélisme (implique d'avoir sélectionné un projet)
std::vector<std::string> generate_alt_umers(const std::string& umer, const std::vector<std::string>& alt_alleles, int variant_pos_in_umer, std::vector<double> specfic_project_freqs) {
    std::vector<std::string> alt_umers;
    if (specfic_project_freqs.size()>0){
        std::string left = umer.substr(0, variant_pos_in_umer);
        std::string right = umer.substr(variant_pos_in_umer + 1);

        for(size_t i = 0; i < alt_alleles.size(); i++){
            if(specfic_project_freqs[i+1] == 0){
                std::cerr << "ALT ("<< alt_alleles[i] <<") freq is " << specfic_project_freqs[i+1] << ", skipping k-mer generation" << std::endl;
            } else {
                std::cout << alt_alleles[i] << '\t' << specfic_project_freqs[i+1] << std::endl;
                std::string alt = alt_alleles[i];
                std::string alt_umer = left + alt + right;
                alt_umers.push_back(alt_umer);
            }
        }
    } else {
        std::cerr << "Warning: project name not valid for this variant" << std::endl;
    }
    return alt_umers;
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <VCF file> <FASTA file>" << std::endl;
        return 1;
    }

    std::string vcf_file = argv[1];
    std::string fasta_file = argv[2];

    // Ouvrir le fichier VCF
    vcfpp::BcfReader vcf(vcf_file.c_str());
    vcfpp::BcfRecord record(vcf.header);

    // Obtenir les sous-champs INFO de l'en-tête - C'est gadget mais c'était pour voir
    // C'est pas possible de faire ça avec VCFPP, c'est dommage.
    // Du coup j'ai ajouté un truc dans vcfpp.h
    /*bcf_hdr_t* header = vcf.header.get_hdr();
    std::cout << "INFO available INFO fields:" << std::endl;
    for (int i = 0; i < header->nhrec; ++i) {
        bcf_hrec_t* hrec = header->hrec[i];
        if (hrec->type == BCF_HL_INFO) {
            for (int j = 0; j < hrec->nkeys; ++j) {
                if (strcmp(hrec->keys[j], "ID") == 0) {
                    std::cout << hrec->vals[j] << " ";
                }
            }
        }
    }*/

    // Récupérer les champs 
    std::vector<std::string> infos_ids = vcf.header.getINFO_IDs();
    std::cout << "Available INFO fields:" << std::endl;
    for (size_t i = 0; i < infos_ids.size(); ++i) {
        std::cout << infos_ids[i] << " ";
    }
    std::cout << std::endl;

    // Ouvrir le fichier FASTA
    faidx_t* fai = fai_load(fasta_file.c_str());
    if (!fai) {
        std::cerr << "Erreur lors de l'ouverture du fichier FASTA: " << fasta_file << std::endl;
        return 1;
    }

    // Parcourir chaque variant dans le fichier VCF
    while (vcf.getNextVariant(record)) {
        std::cout << "-----------" << std::endl;
        ///////////////////////////////////////////
        // Récupérer les informations des champs //
        ///////////////////////////////////////////

        std::string chrom = record.CHROM();         // CHROM
        int pos = record.POS();                     // POS
        std::string id = record.ID();               // ID
        std::string ref_allele = record.REF();      // REF
        std::string alt_alleles_str = record.ALT(); // ALT
        std::string info = record.allINFO();        // INFO (tout le champ INFO)

        // Vérifier si REF contient des nucléotides dégénérés
        bool degenerate_found = false;
        for (char base : ref_allele) {
            if (isDegenerate(base)) {
                degenerate_found = true;
                std::cerr << "Erreur : Le nucléotide " << base << " dans REF est dégénéré pour " << id << " " << chrom << ":" << pos << std::endl;
                break;
            }
        }
        if (degenerate_found) {
            continue; // Passer au variant suivant
        }

        // Fréquences: Map, clé = nom du projet, valeurs = fréquences.
        std::map<std::string, std::vector<double>> frequencies_map;

        // Obtenir le nucléotide de référence à la position spécifiée dans le fichier FASTA
        int len;
        char* seq = faidx_fetch_seq(fai, chrom.c_str(), pos - 1, pos - 1 + ref_allele.size() - 1, &len);
        if (seq == nullptr) {
            std::cerr << "Erreur lors de la récupération de la séquence pour " << chrom << ":" << pos << std::endl;
            fai_destroy(fai);
            return 1;
        }
        
        // Comparer l'allèle de référence du VCF avec le nucléotide de la séquence FASTA
        std::string fasta_allele(seq, len);
        std::string result = (ref_allele == fasta_allele) ? "Check" : "Not found";

        // Afficher les résultats
        std::cout <<"CHROM:" << chrom << "\t" 
                    << "POS:" << pos << "\t" 
                    << "ID:" << id << "\t" 
                    << "REF:" << ref_allele << "\t" 
                    << "ALT:" << alt_alleles_str << "\t" 
                    << "In FASTA: " << result << std::endl;
        //std::cout << info << std::endl;
        
        // Affichage ALT après extraction de la chaine
        // std::vector<std::string> alt_alleles = split(alt_alleles_str, ',');
        std::vector<std::string> alt_alleles = record.rawALT(); // Nouvelle fonction
        std::cout << "ALTS: ";
        for (const std::string& alt : alt_alleles) {
            std::cout << alt << "\t";
        }
        std::cout << std::endl;
        
        // Récupération des champs contenus dans INFO :
        std::cout << "\nINFO field values:" << std::endl;
        std::string info_vc;
        try{
            record.getINFO("VC", info_vc);
        } catch (const std::runtime_error& e) {
            std::cerr << "Warning: " << e.what();
        }
        
        std::string info_freq;
        try{
            record.getINFO("FREQ", info_freq);
        } catch (const std::runtime_error& e) {
            std::cerr << "Warning: " << e.what();
        }
        
        if (!info_vc.empty()){std::cout << "VC:\t" << info_vc << std::endl;}
        if (!info_freq.empty()){std::cout << "FREQ:\t" << info_freq << std::endl;}

        // Vérifier si le champ COMMON est présent dans INFO
        bool has_common = checkInfoField(record.allINFO(), "COMMON");
        if (has_common) {
            std::cout << "COMMON" << std::endl;
        } else {
            std::cout << "Warning: COMMON field is not present in INFO" << std::endl;
        }

        int info_rs;
        try{
            record.getINFO("RS", info_rs);
            std::cout << "RS:\t" << info_rs << std::endl;
        } catch (const std::runtime_error& e) {
            std::cerr << "Warning: " << e.what();
        }

        // Récupération du flag COMMON dans INFO, nouvelle version (modif vcfpp.h)
        std::string info_common;
        try{
            record.getINFO("COMMON", info_common);
            std::cout << "COMMON:\t'" << info_common << "'" << std::endl;
        } catch (const std::runtime_error& e) {
            std::cerr << "Warning: " << e.what();
        }

        // Séparer les projets et leurs valeurs
        if (!info_freq.empty()) {
            std::cout << "FREQ by project:" << std::endl;
            frequencies_map = get_all_freqs(info_freq);
            for (const std::pair<std::string, std::vector<double>>& entry : frequencies_map) {
                std::cout << "- " << entry.first << ": ";
                for (const double& freq : entry.second) {
                    std::cout << freq << " ";
                }
                std::cout << std::endl;
            }
        }

        // Afficher les fréquences d'un projet spécifique
        std::string nom_du_projet_specifique = "Project2";
        std::vector<double> selected_project_freqs;
        std::cout << "\nSelected project for k-mers generation: " << nom_du_projet_specifique << std::endl;
        if (!info_freq.empty()){
            selected_project_freqs = get_freq_from_project_name(frequencies_map, nom_du_projet_specifique);
            for (const double& freq : selected_project_freqs) {
                std::cout << freq << " ";
            }
            std::cout << std::endl;
        }

        // A FAIRE : EXTRAIRE LE UMER - ATTENTION AUX EXTREMITES DE LA SEQUENCE !!
        std::cout << "\nGenerating k-mers..." << std::endl;
        // Extraire le u-mer
        if (info_vc == "SNV") {
            int variant_pos_in_umer;
            std::string umer = extract_snv_umer(fai, chrom, pos, variant_pos_in_umer);
            if (umer.empty()) {
                continue; // Passer au variant si l'extraction a échoué
            }
            std::cout << "u-mer pour " << id << " (" << chrom << ":" << pos << ") : " << umer << std::endl;

            // Générer les ALT-UMER
            std::vector<std::string> alt_umers = generate_alt_umers(umer, alt_alleles, variant_pos_in_umer - 1, selected_project_freqs);
            for (size_t i = 0; i < alt_umers.size(); ++i) {
                std::cout << "ALT-UMER pour " << alt_alleles[i] << " : " << alt_umers[i] << std::endl;

                // Générer les k-mers à partir de chaque alt_umer
                std::vector<std::string> kmers = generate_kmers(alt_umers[i], K);
                std::cout << "k-mers pour ALT-UMER " << alt_alleles[i] << ":" << std::endl;
                for (const std::string& kmer : kmers) {
                    std::cout << kmer << std::endl;
                }
            }
        } else {
            if(!info_vc.empty()){
                std::cout << "We don't know how to handle " << info_vc << "... yet." << std::endl;
            } else {
                std::cout << "Can't generate k-mers. Check prior errors." << std::endl;
            }
        }


        free(seq);
    }

    fai_destroy(fai);
    return 0;
}
