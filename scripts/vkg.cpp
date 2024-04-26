/*
OVKG : Object Variant K-mer Generator

Programme pour extraire les SNV COMMON du vcf de dbSNP et générer les k-mers porteurs de variant correspondants.

*/

/* A FAIRE :
Divers ----------------------
* OK - Comptage
Exclusion -------------------
* OK - Exclure sur le nom du chromosome dès le début, vu que c'est la première info qu'on a
    => récupérer seulement les NC (primary assembly)
* OK - Réorganiser les priorités pour les exclusions
* OK - Exclure les SNP qui ont REF=N
* OK - Exclure les umers de longueur < 2k-1
* OK - Exclure les cas où REF != REF dans la séquence d'origine
* OK - Voir si on peut exclure les séquences qui ne sont pas NC000...
Urgent ----------------------
* OK - Générer les umers pour chaque var
* OK - Faire des vrais k-mers, pas juste un print
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
* Exclusion des k-mers identiques (ou alors marquage) - géré par kim ?
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

class Kmer{
    public:
        string sequence;
        string chromosome;
        string rsid;
        int kmer_start;
        int kmer_end;
        int snp_position_genome;
        int snp_position_kmer;
        char ref_allele;
        std::string alt_allele;
        float reference_allele_frequency;
        float alternative_allele_frequency;

        // CONSTRUCTOR
        Kmer(string seq, string chrom, string id, int start, int end, int snp_gen_pos, int snp_kmer_pos, char ref, std::string alt, float ref_freq, float alt_freq)
            :sequence(seq), 
            chromosome(chrom), 
            rsid(id), 
            kmer_start(start), 
            kmer_end(end),
            snp_position_genome(snp_gen_pos), 
            snp_position_kmer(snp_kmer_pos),
            ref_allele(ref),
            alt_allele(alt),
            reference_allele_frequency(ref_freq),
            alternative_allele_frequency(alt_freq)
        {}

        // Display k-mer information
        void display(){
            cout << "K-mer sequence: " << sequence << "\n"
                << chromosome << "\t"
                << kmer_start << "\t" 
                << kmer_end << "\n"
                << rsid << "\t"
                << "REF: " << ref_allele << " (" << reference_allele_frequency << ")" << "\t"
                << "ALT: " << alt_allele << " (" << alternative_allele_frequency <<")" <<"\n"
                << "SNP Position in Genome:\t" << snp_position_genome << "\n"
                << "SNP Position in K-mer:\t" << snp_position_kmer << "\n"
                << std::endl;
        };
};

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

//  FUNCTIONS ----------------------------------------------------------------------
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

std::vector<Kmer> generate_kmers(const std::string &umer, int k, const std::string &chromosome, const std::string &rsid, int snp_position_genome, char ref_allele, std::string alt_allele, float ref_frequency, float alt_frequency){
    int uLength = umer.length();
    int i = 0;
    std::vector<Kmer> kmer_objects;
    int snp_position_umer = k; // SNP position in the u-mer, center position

    while (i <= uLength - k){
        bool containsDegenerate = false;
        for(int j = i; j < i + k; j++){
            if(isDegenerate(umer[j])){
                containsDegenerate = true;
                i= j + 1; // Skip to the position after the degenerate nucleotide
                break;
            }
        }

        if (containsDegenerate) {
            // Skip this k-mer if it contains a degenerate nucleotide
            continue;
        } else {
            // Creating k-mer and convert characters to uppercase
            std::string kmer = umer.substr(i, k);
            std::transform(kmer.begin(), kmer.end(), kmer.begin(),
                            [](unsigned char c){ return std::toupper(c); });
            // Creating the last values for the k-mer object
            int kmer_start = snp_position_genome - snp_position_umer + i;
            int kmer_end = kmer_start + k;
            int snp_position_kmer = snp_position_umer - i;  // Corrected to account for SNP's relative position
            // Creating k-mer object
            Kmer kmer_obj(kmer, chromosome, rsid, kmer_start, kmer_end, snp_position_genome, snp_position_kmer, ref_allele, alt_allele, ref_frequency, alt_frequency);
            kmer_objects.push_back(kmer_obj);
        }
        i++; // Move to the next base for the next k-mer
    }

    // Print each Kmer object in the vector
    for (Kmer& kmer : kmer_objects) {
        kmer.display();
    }

    return kmer_objects;
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


int main() {

    ///////////////////
    // OPENING FILES //
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
    //INITIALIZING VARIABLES //
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
        // Recap info : chromosome_name (str), position (int), rsid (char*), ref (char*),
        // alts (vector, string),var_class (str), freqs (vector, float), umer (str), alt_umer (str)
        
        // Get chromosome (no hidden char)
        string chromosome_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
        // EXCLUSION : NW et NT
        // Check if chromosome begins with NC ; attention : ne prend pas l'ADNmt avec 'NC_000'
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
        // EXCLUSION : Degenerated nt
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
        // Get Variant Class (VC) - Warning: Contains hidden characters
        string var_class = (char *)(info_vc->vptr);
        var_class = removeNonPrintableChars(var_class);
        // EXCLUSION - VARIANT CLASS
        if (var_class != "SNV"){
            continue;
        }

        // FREQ - Warning: Contains hidden characters
        bcf_info_t *info_freqs = bcf_get_info(vcf_header, vcf_record, "FREQ");
        // Get frequency project source (FREQ)
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
        // AFFICHAGE INFOS //
        /////////////////////
        /*cout << "-------------------------------------" << endl;
        // Print the data
        cout << chromosome_name << "\t" << rsid  << "\t" << position << "\t" << ref << "\t";
        for (size_t i = 0; i < alts.size(); i++){
            if(i == alts.size()-1){
                cout << alts[i] << endl;
            } else {
                cout << alts[i] << ", ";
            }
        }*/

        /////////////////////////////
        // RETRIEVE FASTA SEQUENCE //
        /////////////////////////////

        // Dealing with SNVs
        // Get sequence from fasta file
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
                    cout << alts[i] << " has a null frequency. Skipping..." << endl;
                    continue;
                } else {
                    string left = umer.substr(0,kmer_size-1);
                    string right = umer.substr(kmer_size);
                    string alt_umer = left + alts[i] + right;
                    //kmer_generator(alt_umer, kmer_size);
                    generate_kmers(alt_umer, kmer_size, chromosome_name, rsid, position, ref[i], alts[i], dbgap_freqs[0], dbgap_freqs[i+1]);
                }
            }
        }

    }

    // Close and clean up
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf_file);
    fai_destroy(fai);

    cout << "SNPs selected:\t" << selected_snps_count << endl;

    return 0;
}