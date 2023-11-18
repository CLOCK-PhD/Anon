/*
VKG : Variant K-mer Generator

Programme de test pour générer un index de k-mer, version c++.
(Version Python : ksv.py)

Crée tous k-mers porteurs de variants possibles à partir du fichier vcf de dbSNP,
et du génome de référence HG38.

*/

/* A FAIRE :
Exclusion 
* OK - Exclure les SNP qui ont REF=N
* Exclure les umers de longueur < 2k-1
* OK - Exclure les cas où REF != REF dans la séquence d'origine
* Voir si on peut exclure les séquences qui ne sont pas NC000...
Urgent
* Générer les umers pour chaque var
* Faire des vrais k-mers, pas juste un print
Développer :
* Fonctions pour gérer tous les cas des VC (voir vkg.py)
    OK - SNV (c'est le cas le plus simple donc bon...)
    MNV
    INS
    DEL
    INDEL
* Log messages d'erreurs
    Pour ne plus avoir à afficher des trucs dans la console, et garder une trace de ce qui s'est passé
Gros truc :
* Comparer aux k-mers du génome de référence
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

// TEST REMOVE HIDDEN CHAR - OK
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
// FIN TEST REMOVE HIDDEN CHAR - OK

bool isDegenerate(char base) {
    // Define a function to check if a base is degenerate
    return (base == 'R' || base == 'Y' || base == 'S' || base == 'W' || base == 'K' ||
            base == 'M' || base == 'B' || base == 'D' || base == 'H' || base == 'V' || base == 'N');
}

void kmer_generator(string umer, int k){
    // Generate k-mers from a u-mer, given a k-mer size k
    int uLength = umer.length();
    int i = 0;
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
        // Convert characters to uppercase before outputting
        for (char& c : kmer) {
            c = std::toupper(c);
        }
        std::cout << "K-mer " << i + 1 << ":\t" << kmer << std::endl;
        i++;        
    }

}

int main() {

    ///////////////////
    // OPENING FILES //
    ///////////////////

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

    // VARIABLES
    int kmer_size = 21;
    // STATS
    int selected_snps_count;
    //int kmers_count;

    // Initialize VCF header and record
    bcf_hdr_t *vcf_header = bcf_hdr_read(vcf_file);
    bcf1_t *vcf_record = bcf_init();

    // Read each record in the VCF file
    while (bcf_read(vcf_file, vcf_header, vcf_record) >= 0) {
        // UNPACK - MANDATORY
        //bcf_unpack(vcf_record, BCF_UN_STR);   // Pour unpack 
        //bcf_unpack(vcf_record, BCF_UN_INFO); // Pour unpack jusqu'à INFO
        bcf_unpack(vcf_record, BCF_UN_ALL); // Unpack everything

        // Accessing INFO_FIELD_NAME
        // COMMON - OK
        bcf_info_t *info_common = bcf_get_info(vcf_header, vcf_record, "COMMON");        
        // VC (Variant Class) - OK
        bcf_info_t *info_vc = bcf_get_info(vcf_header, vcf_record, "VC");
        // FREQ - OK
        bcf_info_t *info_freqs = bcf_get_info(vcf_header, vcf_record, "FREQ");

        // Accessing vcf fields
        //1:CHROM, 2:POS, 3:ID, 4:REF, 5;ALT, 6:QUAL, 7:FILTER, 8:INFO
        // Get chromosome - OK (no hidden char)
        string chromosome_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
        // Get position - OK
        int position = vcf_record->pos; // 1-based position in VCF
        // Get ID - OK
        char *rsid = vcf_record->d.id;
        // Get REF (vcf_record->d.allele[0] is REF other is ALT) - OK
        char* ref = vcf_record->d.allele[0];
        // Get ALT - OK
        vector<string> alts(vcf_record->n_allele -1);
        for (int i = 1; i < vcf_record->n_allele; ++i){
            alts[i-1] = vcf_record->d.allele[i];
        }
        // TEST EXCLUSION : NT DEGEN
        // On ne prend que le premier char mais normalement c'est le seul cas où on trouve N
        if (isDegenerate(ref[0])){
            cout << "DEGEN : " << ref << " for " << rsid << " in " << chromosome_name << endl;
            continue;
        }
        // Faire le test où on parcourt chaque char de REF, et on stop dès qu'un N est rencontré
        
        // Print the data
        //std::cout << "Chromosome: " << chromosome_name << ", Position: " << position << std::endl;
        //std::cout << "ID : " << rsid << std::endl;
        //std::cout << "REF " << ref << std::endl;
        // Print ALT - OK
        //std::cout << vcf_record->n_allele << std::endl; // nombre de ref + alleles
        /*std::cout << "ALT ";
        for (int i = 1; i < vcf_record->n_allele; ++i){
            std::cout << vcf_record->d.allele[i] << " ";
        }
        std::cout << std::endl;
        for (int i = 0; i < alts.size(); i++){
            cout << alts[i] << " ";
        }
        cout << endl;*/
        
        // check if is snp - OK
        //std::cout << "Is SNP : " << bcf_is_snp(vcf_record) << std::endl;

        // GET INFO - WIP (1/2)
        //std::cout << "INFOS" << std::endl;

        // Get Variant Class (VC) - Warning : Contains hiddent characters - OK
        string var_class = (char *)(info_vc->vptr);
        var_class = removeNonPrintableChars(var_class);
        /*if (info_vc && info_vc->type == BCF_BT_CHAR){
            string vc_string = (char*)(info_vc->vptr);
            //cout << "VARIANT CLASS :\t" <<vc_string << endl;
        }*/

        // EXCLUSION : NOT SNV - OK
        if (var_class != "SNV"){
            continue;
        }
        
        // Get frequency project source (FREQ) - OK
        const char *freqs;
        if (info_freqs && info_freqs->type == BCF_BT_CHAR){
            freqs = (char*)(info_freqs->vptr);
        }

        // IS COMMON - OK
        /*if (info_common){
            cout << "COMMON : True" << endl;
        } else {
            cout << "COMMON : False" << endl;
        }*/
        // EXCLUSION - NOT COMMON - OK
        if(!info_common){
            continue;
        }

        /////////////////////////////
        // RETRIEVE FASTA SEQUENCE // OK : PROBLEME AVEC LE NOMBRE DE STRING
        /////////////////////////////

        // Dealing with SNVs
        int *len;   // Nécessaire pour faidx_fetch_seq; raison inconnue mais permet de fonctionner.
        int start = position - (kmer_size-1);
        int end = position + (kmer_size-1);
        const char *reg = chromosome_name.c_str();
        char *sequence = faidx_fetch_seq(fai, reg, start, end, len);
        if(!sequence){
            cerr << "Error: Could not fetch the sequence" << endl;
            return 1;
        }

        // TEST EXCLUSION : NW et NT - OK, par conservation des NC uniquement.
        // On vérifie si les séquences commencent par NC
        if (chromosome_name.compare(0, 2, "NC") != 0) {
            //cerr << "Not Primary assembly sequence : " << chromosome_name << endl;
            free(sequence);
            continue;
        }
        // il faudrait penser à s'occuper du cas de la séquence MT

        // EXCLUSION - REF != NT in the sequence
        if (ref[0] != toupper(sequence[kmer_size-1])){
            cerr << "WARNING : " << chromosome_name << " " << rsid << endl;
            cerr << "\t" << ref[0] << "\t" << sequence[kmer_size-1] << endl;
            cerr << "\t" << ref << "\t" << sequence << endl;
            free(sequence);
            continue;
        }
        //cout << "u-mer :\t" << endl;
        //cout << sequence << endl;


        /////////////////////////////
        // GENERATING K-MERS : WIP //
        /////////////////////////////

        // Si on est là, c'est qu'on a passé tous les tests
        selected_snps_count++;

        /////////////////////
        // AFFICHAGE INFOS //
        /////////////////////

        cout << "-------------------------------------" << endl;
        // Print the data
        cout << chromosome_name << "\t" << rsid  << "\t" << position << "\t" << ref << "\t";
        for (int i = 0; i < alts.size(); i++){
            if(i == alts.size()-1){
                cout << alts[i] << endl;
            } else {
                cout << alts[i] << ", ";
            }
        }

        // Print FREQ
        cout << freqs << endl;
        // FIN AFFICHAGE INFOS

        // GESTION FREQ
        // PB : CORE DUMPED SI LA DERNIERE FREQUENCE SE TERMINE PAR "."
        // Récupérer les différentes fréquences, sélectionner une référence, ajuster les alt

        /*// Create a stringstream from the input string
        istringstream iss(freqs);
        string token;
        // Vector to store extracted source names and values
        std::vector<std::pair<std::string, std::vector<float>>> sources;
        // Tokenize the input string by '|' character
        while (std::getline(iss, token, '|')) {
            // Tokenize each part by ':' character
            std::istringstream partStream(token);
            std::string sourceName;
            std::vector<float> values;
            std::getline(partStream, sourceName, ':');

            // Tokenize the values part by ','
            while (std::getline(partStream, token, ',')) {
                // Replace "." with "0" if it's not a valid number
                bool isValidNumber = true;
                for (char c : token) {
                    if (!isdigit(c) && c != '.' && c != '-')
                        isValidNumber = false;
                }
                if (token == "." || (token.empty() && isValidNumber))
                    token = "0";

                // Convert the string to a float
                float value = std::stof(token);
                values.push_back(value);
            // Convert the string to a float
            //values.push_back(std::stof(token));
            }
            // Add the source name and values to the vector
            sources.push_back(std::make_pair(sourceName, values));
        }

        // AFFICHAGE TEST FREQ
        for (const auto& source : sources) {
            std::cout << "Source Name: " << source.first << std::endl;
            std::cout << "Values: ";
            for (const auto& value : source.second) {
                std::cout << value << " ";
            }
        std::cout << std::endl;
        }*/

        // TEST GENERATION ALT_UMERS:
        //cout << "TEST ALT UMERS" << endl;
        string umer = sequence;
        //cout << umer << endl;
        for (int i=0; i < alts.size(); i++){
            string left = umer.substr(0,kmer_size-1);   //OK
            string right = umer.substr(kmer_size);      //OK
            //cout << left << "X" << right << endl;
            //cout << left << alts[i] << right << endl;
            string alt_umer = left + alts[i] + right;
            //cout << alt_umer << endl;g++ vkg.cpp -std=c++17 -Wall -Wextra -o truc -lhts            //cout << "K-mers for " << ref << "->" << alts[i] << ":" << endl;
            kmer_generator(alt_umer, kmer_size);
        }

        // TEST KMER_GENERATOR() - OK
        //cout << "TEST KMER_GENERATOR()" << endl;
        //cout << umer << endl;
        //kmer_generator(umer, kmer_size);

        // Free memory space
        free(sequence);

        /* STRINGS SACRIFIABLES
        Pour une raison que j'ignore, le nombre de strings dans le programme entraine un core dumped.
        Cela peut arriver si on ajoute un string ou un char* (moins fréquent), ou parfois qu'on en enlève.
        Certaines conditions nécessitent plusieurs sacrifices.
        Une "solution" pour contourner le problème est :
        - mettre une partie du code en commentaire (lecture de fasta, en général)
        - compiler
        - ajouter des strings, compiler
        - Décommenter
        - Compiler
        - Quand on rencontre un problème : on sacrifie une string du stock
        - Quand on arrive au bout : on recommence.
        */
        string sacrifice1;    
        string sacrifice2;    
        string sacrifice3;    
        string sacrifice4;   
        string sacrifice5;   
        string sacrifice6;
        string sacrifice7;
        string sacrifice8;
        string sacrifice9;
        string sacrifice10;
        string sacrifice11;
        string sacrifice12;
        string sacrifice13; // FREQ affichage : sacrifice 1, mais affiche la première ligne avant le core dumped
        string sacrifice14; // Là j'ai dû la rajouter quand j'ai supprimé une autre string;
        string sacrifice15;
        string sacrifice16;
        string sacrifice17;
        string sacrifice18;
        //string sacrifice19;
        //string sacrifice20;
        //string sacrifice21;
        //string sacrifice22;
    }

    // Close and clean up
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf_file);
    fai_destroy(fai);

    cout << "SNPs selected:\t" << selected_snps_count << endl;

    return 0;
}