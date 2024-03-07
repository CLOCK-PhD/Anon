/*
VKG : Variant K-mer Generator

Programme de test pour générer un index de k-mer, version c++.
(Version Python : ksv.py)

Crée tous k-mers porteurs de variants possibles à partir du fichier vcf de dbSNP,
et du génome de référence HG38.

TEST PC FIXE

*/

/* A FAIRE :
Divers ----------------------
* OK - Comptage
Exclusion -------------------
* OK - Exclure sur le nom du chromosome dès le début, vu que c'est la première info qu'on a
    => récupérer seulement les NC (primary assembly)
* PRIO - EN COURS - Réorganiser les priorités pour les exclusions
* OK - Exclure les SNP qui ont REF=N
* OK - Exclure les umers de longueur < 2k-1
* OK - Exclure les cas où REF != REF dans la séquence d'origine
* OK - Voir si on peut exclure les séquences qui ne sont pas NC000...
Urgent ----------------------
* OK - Générer les umers pour chaque var
* Faire des vrais k-mers, pas juste un print
Développer ------------------
* Fonctions pour gérer tous les cas des VC (voir vkg.py)
    OK - SNV (c'est le cas le plus simple donc bon...)
    MNV
    INS
    DEL
    INDEL
* Log messages d'erreurs
    Pour ne plus avoir à afficher des trucs dans la console, et garder une trace de ce qui s'est passé
Gros truc --------------------
* Marquage des k-mers présents dans le génome de référence
* Exclusion des k-mers identiques (ou alors marquage)
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

using namespace std;

// TEST REMOVE HIDDEN CHAR - OK ---------------------------------------------------
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
// FIN TEST REMOVE HIDDEN CHAR - OK ------------------------------------------------

//  FONCTIONS ----------------------------------------------------------------------
bool isDegenerate(char base) {
    // Define a function to check if a base is degenerated
    return (base == 'R' || base == 'Y' || base == 'S' || base == 'W' || base == 'K' ||
            base == 'M' || base == 'B' || base == 'D' || base == 'H' || base == 'V' || base == 'N');
}

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

}

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

    // Ne fonctionne pas à cause du problème des strings

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

    size_t kmer_size = 21;

    // STATS
    long selected_snps_count = 0;
    //int kmers_count;

    // Initialize VCF header and record
    bcf_hdr_t *vcf_header = bcf_hdr_read(vcf_file);
    bcf1_t *vcf_record = bcf_init();

    // Read each record in the VCF file
    while (bcf_read(vcf_file, vcf_header, vcf_record) >= 0) {
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
        // EXCLUSION
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
        /* cout << "-------------------------------------" << endl;
        // Print the data
        cout << chromosome_name << "\t" << rsid  << "\t" << position << "\t" << ref << "\t";
        for (int i = 0; i < alts.size(); i++){
            if(i == alts.size()-1){
                cout << alts[i] << endl;
            } else {
                cout << alts[i] << ", ";
            }
        }*/
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
                    kmer_generator(alt_umer, kmer_size);
                }
            }
        }


        // TEST KMER_GENERATOR() - OK
        //cout << "TEST KMER_GENERATOR()" << endl;
        //cout << umer << endl;
        //kmer_generator(umer, kmer_size);

    }

    // Close and clean up
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf_file);
    fai_destroy(fai);

    cout << "SNPs selected:\t" << selected_snps_count << endl;

    return 0;
}