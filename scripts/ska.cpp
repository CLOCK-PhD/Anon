/*
SKA : SNP K-mer Analyzer

Programme test pour détecter les SNP présents au sein d'un k-mer.

(fait à partir de vkg.cpp où j'avais déjà fait plein de manips pour explorer dbSNP)

Entrée :
Fichier (json ?) de k-mers contenant :
- La séquence du k-mer
- Le chromosome sur lequel se situe le k-mer
- La position de départ "start" sur le génome de référence
- La position de fin "end" sur le génome de référence
- Le brin (voir comment on gère ce cas plus tard)
*/

/* A FAIRE :
    0. Lire le fichier de k-mers, récupérer les données nécessaires.
    1. Vérifier l'ouverture et la lecture du fichier vcf dbSNP
    2. Récupérer tous les SNPs entre les positions start et end
    3. Récupérer la séquence de référence aux positions start et end pour le chromosome approprié
    4. Comparer la séquence du k-mer et le k-mer extrait de la référence
    5. Identifier les SNPs présents ou non
    6. Récupérer les fréquences des SNP détectés
    7. Calculer le score
*/

/* APRÈS :
    - Extraire les infos associées au SNP sur PheGenI (Phenotype-Genotype Integrator)
    https://www.ncbi.nlm.nih.gov/gap/phegeni
*/

/* PLAN :
    - faire 1
    - avant de faire avec le json, faire un test avec des  variables directement dans le code pour les k-mers (2 ou 3)
    - 

*/

/* NOTES :

    1. Le programme de base est fait pour lire dbSNP puis aller chercher la séquence de référence.
    On devrait pouvoir procéder autrement : d'abord récupérer la séquence de référence et comparer les deux.
    ATTENTION : penser à s'assurer que les positions correspondent.
    En conservant les positions, on devrait pouvoir identifier quelles sont les variations présentes dans le  k-mer.

    La première étape devrait donc être de récupérer la séquence et de faire une comparaison entre les deux.

    Ensuite, on devrait pouvoir faire des requêtes dans dbSNP une fois qu'on aura identifié les variations présentes.
    Si on a des SNP, ça devrait être simple, mais il faudra prendre en compte les autres types de variation possibles.

    On va commencer par ça.

*/


#include <vector>
#include <cstdlib>
#include <iostream>
#include <htslib/hts.h>
#include <htslib/vcf_sweep.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/faidx.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <cctype> // for toupper

using namespace std;


    // Noms des chromosomes GRCh38p14
    /*
    chr1 = NC_000001.11
    chr2 = NC_000002.12
    chr3 = NC_000003.13
    chr4 = NC_000004.12
    chr5 = NC_000005.10
    chr6 = NC_000006.12
    chr7 = NC_000007.14
    chr8 = NC_000008.11
    chr9 = NC_000009.12
    chr10 = NC_000010.11
    chr11 = NC_000011.10
    chr12 = NC_000012.12
    chr13 = NC_000013.11
    chr14 = NC_000014.9
    chr15 = NC_000015.10
    chr16 = NC_000016.10
    chr17 = NC_000017.11
    chr17 = NC_000018.10
    chr18 = NC_000019.10
    chr20 = NC_000020.11
    chr21 = NC_000021.9
    chr22 = NC_000022.11
    chr23 = NC_000023.11
    chr24 = NC_000024.10
    mt = NC_012920.1
    */ 

///////////////////////////////
// FONCTIONS DEVELOPPEES ICI //
///////////////////////////////
void printDifferences(const std::string& str1, const std::string& str2) {
    // Check if the strings have the same length
    if (str1.length() != str2.length()) {
        std::cerr << "Error: Strings do not have the same length." << std::endl;
        return; // Exit the function if lengths differ
    }

    // Iterate through each character since the lengths are the same
    for (size_t i = 0; i < str1.length(); ++i) {
        // If characters at the current position are different
        if (str1[i] != str2[i]) {
            std::cout << "Difference at position " << i << ": "
                      << "'" << str1[i] << "' (str1) vs "
                      << "'" << str2[i] << "' (str2)" << std::endl;
        }
    }
}

std::string toUpperCase(const std::string& input) {
    std::string result = input; // Make a copy of the input
    for (char& c : result) { // Use reference to modify the string in-place
        c = std::toupper(c); // Convert each character to uppercase
    }
    return result;
}

bool isPrintable(char c) {
    // Check if the character is a printable ASCII character
    return (c >= 32 && c <= 126);
}

std::string removeNonPrintableChars(const std::string& input) {
    // Remove hidden characters from a string
    std::string result;
    for (char c : input) {
        if (isPrintable(c)) {
            result += c;
        }
    }
    return result;
}

std::string displayHiddenChars(const std::string& input) {
    // Show hidden characters in a string
    std::string result;
    for (char c : input) {
        if (c >= 32 && c <= 126) {
            result += c; // Printable characters are added as is
        } else {
            // Non-printable characters are displayed as hexadecimal values
            result += "\\x" + std::to_string(static_cast<unsigned int>(c) & 0xFF);
        }
    }
    return result;
}

//  FONCTIONS ----------------------------------------------------------------------
bool isDegenerate(char base) {
    // Define a function to check if a base is degenerated
    return (base == 'R' || base == 'Y' || base == 'S' || base == 'W' || base == 'K' ||
            base == 'M' || base == 'B' || base == 'D' || base == 'H' || base == 'V' || base == 'N');
}

// Pas utile dans ce programme
/*
void kmer_generator(const string &umer, int k){
    // Generate k-mers from a u-mer, given a k-mer size k
    int uLength = umer.length();
    int i = 0; // compteur pour la boucle
    int kmer_count = 0; // compter de k-mers générés
    while(i <= uLength - k){
        bool containsDegenerate = false;
        for (int j = i; j < i + k; j++) {
            if (isDegenerate(umer[j])) {
                containsDegenerate = true;
                i = j+1; // Skip to the position after the degenerated nt
                break;
            }
        }
        if (containsDegenerate) {
            // Skip this k-mer if it contains a degenerate nucleotide
            continue;
        }
        std::string kmer = umer.substr(i, k);
        kmer_count ++;
        // Convert characters to uppercase before outputting
        for (char& c : kmer) {
            c = std::toupper(c);
        }
        std::cout << "K-mer " << i + 1 << ":\t" << kmer << std::endl;
        i++;        
    }
    cout << "K-mers generated: " << kmer_count << endl;

}*/

std::vector<std::string> split(const std::string &s, char delimiter) {
    // Initialize a vector of strings to store tokens
    std::vector<std::string> tokens;

    // Convert the string to a stream for parsing
    std::istringstream stream(s);

    // Initialize a string to store each token
    std::string token;

    // Loop through the stream and extract tokens
    while (std::getline(stream, token, delimiter)) {
        // Add each token to the vector
        tokens.push_back(token);
    }

    // Return the vector of tokens
    return tokens;
}

bool check_source(char* my_source, const std::string &input){

    char delimiter = '|'; // Split each source and its frequencies from the others

    std::vector<std::string> tokens = split(input, delimiter);
    std::vector<float> dbgap_freq;

    for (size_t i=0; i < tokens.size(); i++){
        std::vector<std::string> subtokens = split(tokens[i], ':');
        std::string sourcename = subtokens[0];
        if (sourcename == my_source){
            return true;
        }
    }
    return false;
}


std::vector<float> get_dbgap_freq(const std::string &input){
    // Return a vector containing the frequencies of dbGaP_PopFreq

    char delimiter = '|'; // Split each source and its frequencies from the others

    std::vector<std::string> tokens = split(input, delimiter);
    std::vector<float> dbgap_freq;

    for (size_t i=0; i < tokens.size(); i++){
        std::vector<std::string> subtokens = split(tokens[i], ':');
        std::string sourcename = subtokens[0];
        std::string allele_frequencies = subtokens[1];
        if (sourcename != "dbGaP_PopFreq"){
            continue;
        }
        else{
            std::vector<std::string> frequencies = split(allele_frequencies, ',');
            for(size_t j=0; j < frequencies.size(); j++){
                if (removeNonPrintableChars(frequencies[j]) == "."){
                    dbgap_freq.push_back(0.0);
                }
                else{
                    dbgap_freq.push_back(std::stof(frequencies[j]));
                }
            }
        }
    }

    return dbgap_freq;
}

// FIN FONCTIONS -------------------------------------------------------------------</s>

int main() {

    ///////////////////
    // OPENING FILES /////////////////////////////////////////////////////////////////////
    ///////////////////

    // const char* source_name = "dbGaP_PopFreq";

    // OPENING FASTA FILE - OK
    const char* fasta_file_path = "/home/remycosta/phd/Anon/data/grch38p14/GCF_000001405.40_GRCh38.p14_genomic.fna";
    faidx_t* fai = fai_load(fasta_file_path);
    if(!fai){
        cerr << "Error opening fasta file: " << fasta_file_path << endl;
        return 1;
    }

    // OPENING VCF FILE - OK
    const char* vcf_file_path = "/home/remycosta/phd/Anon/data/dbsnp156/GCF_000001405.40.gz";
    htsFile *vcf_file = bcf_open(vcf_file_path, "r");
    if (!vcf_file) {
        std::cerr << "Error opening VCF file: " << vcf_file_path << std::endl;
        return 1;
    }
    ///////////////////////////
    //INITIALIZING VARIABLES ///////////////////////////////////////////////////////////////
    ///////////////////////////

    // Test grossier avec un k-mer tout moche
    //TGTGCTCAGCAAAGACAGACTTGAAGGTCCG': ['chr22', 39317497, 39317528, '+', [15]]
    /*string chrom_name = "NC_000022.11";
    int start_pos = 39317497;
    int end_pos = 39317528;
    string strand = "+";
    string kmer = "TGTGCTCAGCAAAGACAGACTTGAAGGTCCG";
    */

    //'TTAGTGCCAACAATAAGTGGCCATTCGTAAA': ['chr2', 147927985, 147928016, '+', [35]],
    // PAS DE DIFF ICI
    /*string chrom_name = "NC_000002.12";
    int start_pos = 147927985;
    int end_pos = 147928016;
    string strand = "+";
    string kmer = "TTAGTGCCAACAATAAGTGGCCATTCGTAAA";
    */

    //'AAGTGGCGTGAACTCGGCTCACTGCAACCTC': ['chr1', 41085894, 41085925, '+', [0]],
    string chrom_name = "NC_000001.11";
    int start_pos = 41085894;
    int end_pos = 41085925;
    string strand = "+";
    string kmer = "AAGTGGCGTGAACTCGGCTCACTGCAACCTC";


    cout << kmer << " " << start_pos << " " << end_pos << endl;


    //////////////////////////////////////
    // Initialize VCF header and record //
    //////////////////////////////////////

    bcf_hdr_t *vcf_header = bcf_hdr_read(vcf_file);
    bcf1_t *vcf_record = bcf_init();

    ////////////////////////////
    // GET REFERENCE SEQUENCE //
    ////////////////////////////

    int len;   // Mandatory to use faidx_fetch_seq()
    //int start = position - (kmer_size-1);
    //int end = position + (kmer_size-1);
    int start = start_pos;
    int end = end_pos -1 ;
    //const char *reg = chromosome_name.c_str();
    const char *reg = chrom_name.c_str();
    //char *ref_sequence = faidx_fetch_seq(fai, reg, start, end, &len);
    char *s = faidx_fetch_seq(fai, reg, start, end, &len);
    if(!s){
        cerr << "WARNING: Could not fetch the sequence" << endl;
    }
    string ref_kmer = s;
    ref_kmer = toUpperCase(ref_kmer);
    free(s);
    switch (len) {
        case 1:
            cerr << "WARNING: unable to retrieve the sequence corresponding to chromosome "
                 << reg << " between " << start << " and " << end << " positions" << endl;
            break;
        case 2:
            cerr << "WARNING: invalid chromosome name (sequence not found)" << endl;
            break;
        default:
            len = 0;
    }
    cout << ref_kmer << endl;

    printDifferences(ref_kmer, kmer);


    // OK - Vérification pour la position du SNP :
    int start_verif = start_pos; // c'est bien la valeur indiquée par la fonction avec le test bête
    int end_verif = start_verif;
    char *s2 = faidx_fetch_seq(fai, reg, start_verif, end_verif, &len);
    string verif_base = s2;
    free(s2);
    cout << "Verifying base at position "<<start_verif<<": "<<verif_base<<endl;
    int target_pos = start_verif ;

    // TEST POUR LIRE LE VCF ET RETROUVER LE SNP
    cout << "On essaie de retrouver le SNP" << endl;
    while (bcf_read(vcf_file, vcf_header, vcf_record) == 0) {
        // AJOUT DE CODE POUR LES INFOS SUPPLEMENTAIRES
        bcf_unpack(vcf_record, BCF_UN_ALL);
        // Get ID
        char *rsid = vcf_record->d.id;
        // Get REF (vcf_record->d.allele[0] is REF other is ALT)     
        char* ref = vcf_record->d.allele[0];
        // Get ALT
        vector<string> alts(vcf_record->n_allele -1);
        for (int i = 1; i < vcf_record->n_allele; ++i){
            alts[i-1] = vcf_record->d.allele[i];
        }
        // VC (Variant Class)
        bcf_info_t *info_vc = bcf_get_info(vcf_header, vcf_record, "VC");






        int32_t pos = vcf_record->pos +1 ; // Positions are 0-based in htslib
        const char *chrom = bcf_hdr_id2name(vcf_header, vcf_record->rid);



        // Check if the current record matches the target chromosome and position
        if (chrom_name == chrom && target_pos == pos) {
            std::cout << "Match found at " << chrom << ":" << pos << std::endl;
            // Add code here to process or print the record details
            cout << chrom << "\t" << rsid  << "\t" << pos << "\t" << ref << "\t";
            for (size_t i = 0; i < alts.size(); i++){
            if(i == alts.size()-1){
                cout << alts[i] << endl;
            } else {
                cout << alts[i] << ", ";
            }
        }
            break; // Exit loop after finding the match
        }
    }



    // Read each record in the VCF file
    /*while (bcf_read(vcf_file, vcf_header, vcf_record) >= 0) {
        // UNPACK - MANDATORY
        //bcf_unpack(vcf_record, BCF_UN_STR);   // Pour unpack 
        //bcf_unpack(vcf_record, BCF_UN_INFO);  // Pour unpack jusqu'à INFO
        bcf_unpack(vcf_record, BCF_UN_ALL); // Unpack everything

        /////////////////////////////////////////////////////////////////////////////////////////
        // Accessing VCF fields - 1:CHROM, 2:POS, 3:ID, 4:REF, 5;ALT, 6:QUAL, 7:FILTER, 8:INFO //
        /////////////////////////////////////////////////////////////////////////////////////////
        
        // Get chromosome (no hidden char)
        string chromosome_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
        // EXCLUSION : NW et NT
        // On vérifie si le chromosome commencent par NC
        if (chromosome_name.compare(0, 6, "NC_000") != 0) {
            //cerr << "Not Primary assembly sequence : " << chromosome_name << endl;
            continue;
        }
        
        // Get position
        int position = vcf_record->pos; // 1-based position in VCF
        
        // Get ID
        char *rsid = vcf_record->d.id;
        
        // Get REF (vcf_record->d.allele[0] is REF other is ALT)     
        char* ref = vcf_record->d.allele[0];
        // EXCLUSION : Nucléotide dégénéré
        bool ref_checkpoint = true;
        for (const char* ptr = ref; *ptr !='\0'; ++ptr){
            if (isDegenerate(*ptr)) {
                cerr << rsid << ":\tDegenerated nucleotide found at pos "<< position+1 <<" of "<< chromosome_name << " (" << ref << ")" <<endl;
                ref_checkpoint = false;
                break;
            }
        }
        if(!ref_checkpoint){
            continue;
        }
        
        // Get ALT
        vector<string> alts(vcf_record->n_allele -1);
        for (int i = 1; i < vcf_record->n_allele; ++i){
            alts[i-1] = vcf_record->d.allele[i];
        }

        ///////////////////////////////
        // Accessing INFO_FIELD_NAME //
        ///////////////////////////////

        // COMMON
        bcf_info_t *info_common = bcf_get_info(vcf_header, vcf_record, "COMMON");
        // EXCLUSION - COMMON
        if(!info_common){
            continue;
        }

        // VC (Variant Class)
        bcf_info_t *info_vc = bcf_get_info(vcf_header, vcf_record, "VC");
        // Get Variant Class (VC) - Warning: Contains hiddent characters
        string var_class = (char *)(info_vc->vptr);
        var_class = removeNonPrintableChars(var_class);
        // EXCLUSION - VARIANT CLASS
        if (var_class != "SNV"){
            continue;
        }

        // FREQ - Warning: Contains hidden characters
        bcf_info_t *info_freqs = bcf_get_info(vcf_header, vcf_record, "FREQ");
        // Get frequency project source (FREQ) - OK
        const char *freqs;
        if (info_freqs && info_freqs->type == BCF_BT_CHAR){
            freqs = (char*)(info_freqs->vptr);
        }
        // CREATE DBGAP FREQUENCIES VECTOR
        vector<float> dbgap_freqs = get_dbgap_freq(freqs);
        if (dbgap_freqs.size() == 0){
            //cerr << rsid << ": no dbGaP_PopFreq" << endl;
            continue;
        }
        
        /////////////////////
        // AFFICHAGE INFOS ///////////////////////////////////////////////////////////////////////
        /////////////////////
        cout << "-------------------------------------" << endl;
        // Print the data
        cout << chromosome_name << "\t" << rsid  << "\t" << position << "\t" << ref << "\t";
        for (size_t i = 0; i < alts.size(); i++){
            if(i == alts.size()-1){
                cout << alts[i] << endl;
            } else {
                cout << alts[i] << ", ";
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////

        /////////////////////////////
        // RETRIEVE FASTA SEQUENCE // PROBLEME AVEC LE NOMBRE DE STRING
        /////////////////////////////

        // Dealing with SNVs
        int len;   // Mandatory to use faidx_fetch_seq()
        int start = position - (kmer_size-1);
        int end = position + (kmer_size-1);
        const char *reg = chromosome_name.c_str();
        char *s = faidx_fetch_seq(fai, reg, start, end, &len);
        if(!s){
            cerr << "WARNING: Could not fetch the sequence" << endl;
            continue;
        }
        string umer = s;
        free(s);
        switch (len) {
            case 1:
                cerr << "WARNING: unable to retrieve the sequence corresponding to chromosome "
                     << reg << " between " << start << " and " << end << " positions" << endl;
                break;
            case 2:
                cerr << "WARNING: invalid chromosome name (sequence not found)" << endl;
                break;
            default:
                len = 0;
        }
        if (len) {
            continue;
        }

        // EXCLUSION - REF != NT in the sequence
        if (ref[0] != toupper(umer[kmer_size-1])){
            cerr << "WARNING : " << chromosome_name << " " << rsid << endl;
            cerr << "\t" << ref[0] << "\t" << umer[kmer_size-1] << endl;
            //cerr << "\t" << ref << "\t" << sequence << endl;
            continue;
        }

        /////////////////////////////
        // GENERATING K-MERS : WIP ///////////////////////////////////////////////////////////////
        /////////////////////////////

        // Si on est là, c'est qu'on a passé tous les tests
        selected_snps_count++;

        // TEST GENERATION ALT_UMERS:
        if (umer.length() < kmer_size){
            cerr << "The length of the SNP is less than the size of the k-mers (" << kmer_size << "). Skipping this variant..." << endl;
            continue;
        }

        // SUPPRESSION DES ALT A FREQUENCE NULLE
        // Lorsque l'on veut comparer deux entiers non signés A et B à 1 près (typiquement si A == B - 1),
        // il faut préférer comparer A + 1 == B car 0 - 1 n'est pas possible pour un entier non signé et
        // vaut le plus grand entier possible.
        if (alts.size() + 1 != dbgap_freqs.size()){
            continue;
        } else {
            for (size_t i=0; i < alts.size(); i++){
                if ((dbgap_freqs[i+1]) == 0){
                    continue;
                } else {
                    string left = umer.substr(0,kmer_size-1);   //OK
                    string right = umer.substr(kmer_size);      //OK
                    //cout << left << "X" << right << endl;
                    //cout << left << alts[i] << right << endl;
                    string alt_umer = left + alts[i] + right;
                    //cout << alt_umer << endl;
                    //cout << "K-mers for " << ref << "->" << alts[i] << ":" << endl;
                    //kmer_generator(alt_umer, kmer_size);
                }
            }
        }

    }*/

    ////////////////////////
    // Close and clean up //
    ////////////////////////

    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf_file);
    fai_destroy(fai);

    //cout << "SNPs selected:\t" << selected_snps_count << endl;

    return 0;
}