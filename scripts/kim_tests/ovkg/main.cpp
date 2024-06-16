/* Version pour les tests avec mes fichiers.
Je vire le traitement des infos pour le chromosome pour pas faire d'erreur

Compilation : g++ -std=c++17 -Wall -Wextra main.cpp Kmer.cpp -o gatin -lhts
*/


#include "Kmer.h"
#include <vector>
#include <cstdlib>
#include <fstream>
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


// Fonction pour convertir le chromosome en un entier (au lieu de NC_000...)
int extractChromosomeNumber(const std::string &chromosome) {
    // Trouver la position du premier underscore '_'
    size_t firstUnderscore = chromosome.find('_');
    if (firstUnderscore == std::string::npos) {
        std::cerr << "Erreur : Premier underscore non trouvé dans " << chromosome << std::endl;
        return -1; // Si pas de '_', retourner une valeur indicative d'erreur
    }

    // Trouver la position du premier point '.'
    size_t firstDot = chromosome.find('.', firstUnderscore);
    if (firstDot == std::string::npos) {
        std::cerr << "Erreur : Premier point non trouvé dans " << chromosome << std::endl;
        return -1; // Si pas de '.', retourner une valeur indicative d'erreur
    }

    // Extraire le numéro du chromosome entre le premier underscore et le premier point
    std::string number = chromosome.substr(firstUnderscore + 1, firstDot - firstUnderscore - 1);

    // Convertir en entier pour supprimer les zéros en tête
    try {
        int chromosomeNumber = std::stoi(number);
        return chromosomeNumber;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Erreur : Conversion en entier échouée pour " << number << " dans " << chromosome << std::endl;
        return -1; // Si la conversion échoue, retourner une valeur indicative d'erreur
    }
}

int main() {
    ///////////////////
    // OPENING FILES //
    ///////////////////

    // OPENING FASTA FILE - OK
    //const char* fasta_file_path = "/home/remycosta/phd/Anon/data/grch38p14/GCF_000001405.40_GRCh38.p14_genomic.fna";
    const char* fasta_file_path = "../ref_vcf_test.fasta";
    faidx_t* fai = fai_load(fasta_file_path);
    if(!fai){
        cerr << "Error opening fasta file: " << fasta_file_path << endl;
        return 1;
    }

    // OPENING VCF FILE - OK
    //const char* vcf_file_path = "/home/remycosta/phd/Anon/data/dbsnp156/GCF_000001405.40.gz";
    const char* vcf_file_path = "../snp_test.vcf";
    htsFile *vcf_file = bcf_open(vcf_file_path, "r");
    if (!vcf_file) {
        std::cerr << "Error opening VCF file: " << vcf_file_path << std::endl;
        return 1;
    }

    // Test - OUTPUT TSV
    // Path to the output TSV file
    std::string outputFilePath = "kmers.tsv";
    // Open an output file stream
    std::ofstream outputFile(outputFilePath);

    // Check if file is open
    if (!outputFile.is_open()) {
        std::cerr << "Failed to open the file for writing: " << outputFilePath << std::endl;
        return 1;
    }

    // Write the header to the TSV file
    //outputFile << "#KMER\tID\tCHR\tSTART\tEND\tSNP_POS\tIN_KMER\tREF\tFREQ\tALT\tMAF\n";
    // Reduced TSV
    outputFile << "#KMER\tID\tALT\tMAF\tCHROM:POS\n";



    ///////////////////////////
    //INITIALIZING VARIABLES //
    ///////////////////////////

    size_t kmer_size = 5;

    // STATS
    size_t selected_snps_count = 0;
    size_t kmers_count = 0;

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
        
        std::cout << "-------------" << endl;
        // Get chromosome (no hidden char)
        string chromosome_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
        std::cout << "CHROM:" <<chromosome_name << std::endl; // test
        // EXCLUSION : NW et NT
        // Check if chromosome begins with NC ; attention : ne prend pas l'ADNmt avec 'NC_000'
        /*if (chromosome_name.compare(0, 6, "NC_000") != 0) {
            //cerr << "Not Primary assembly sequence : " << chromosome_name << endl;
            continue;
        }*/
        
        // Get position
        size_t position = vcf_record->pos; // 1-based position in VCF
        std::cout << "POS:" << position << std::endl; // test
        
        // Get ID
        char *rsid = vcf_record->d.id;
        // Convert rsID from char* to size_t by removing "rs":
        string rsid_str(rsid);
        if (rsid_str.substr(0, 2) == "rs") { // Check if rsID begins with rs
        rsid_str = rsid_str.substr(2); // Remove "rs" prefix
        }
        size_t rsid_num = stoul(rsid_str);
        std::cout << "RSID:" << rsid << endl; // test
        
        // Get REF (vcf_record->d.allele[0] is REF other is ALT)     
        char* ref = vcf_record->d.allele[0];
        std::cout << "REF:" << ref << endl;
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
        std::cout << "ALTS:" << endl; // test
        vector<string> alts(vcf_record->n_allele -1);
        for (int i = 1; i < vcf_record->n_allele; ++i){
            alts[i-1] = vcf_record->d.allele[i];
            std::cout << alts[i-1] << '\t'; //test
        }
        std::cout << endl; // test

        ///////////////////////////////
        // Accessing INFO_FIELD_NAME //
        ///////////////////////////////

        // COMMON
        bcf_info_t *info_common = bcf_get_info(vcf_header, vcf_record, "COMMON");
        // EXCLUSION - COMMON
        if(!info_common){
            std::cout << "COMMON non présent" << endl;
            //continue;
        }

        // VC (Variant Class)
        bcf_info_t *info_vc = bcf_get_info(vcf_header, vcf_record, "VC");
        // Get Variant Class (VC) - Warning: Contains hidden characters
        string var_class = (char *)(info_vc->vptr);
        var_class = removeNonPrintableChars(var_class);
        std::cout << "VC=" << var_class << std::endl;
        // EXCLUSION - VARIANT CLASS
        if (var_class != "SNV"){
            cerr << "VC is not SNV" << endl;
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
                // check the frequency, change value to exclude frequencies 
                //"== 0" for ALT with null frequencies
                // "< 0.1" will exclude all ALT with a frequency lower than 0.1)
                if ((dbgap_freqs[i+1]) < 0.1){
                    //cout << alts[i] << " has a null frequency. Skipping..." << endl;
                    continue;
                } else {
                    string left = umer.substr(0,kmer_size-1);
                    string right = umer.substr(kmer_size);
                    string alt_umer = left + alts[i] + right;
                    //generate_kmers(alt_umer, kmer_size, chromosome_name, rsid, position, ref[i], alts[i], dbgap_freqs[0], dbgap_freqs[i+1]);
                    // Create a vector containing all the k-mers from the u-mer
                    std::vector<Kmer> kmers = generate_kmers(alt_umer, kmer_size, chromosome_name, rsid_num, position, ref[i], alts[i], dbgap_freqs[0], dbgap_freqs[i+1]);

                    // Iterate over each Kmer object and output the data
                    // Version lourde
                    /*for (const auto& kmer : kmers) {
                        outputFile << kmer.sequence << '\t'
                        << kmer.rsid << '\t'
                        << kmer.chromosome << '\t'
                        << kmer.kmer_start << '\t'
                        << kmer.kmer_end << '\t'
                        << kmer.snp_position_genome << '\t'
                        << kmer.snp_position_kmer << '\t'
                        << kmer.ref_allele << '\t'
                        << kmer.reference_allele_frequency << '\t'
                        << kmer.alt_allele << '\t'
                        << kmer.alternative_allele_frequency << '\n';
                    }*/
                    // Version allégée
                    kmers_count += kmers.size();

                    for (const auto& kmer : kmers) {
                        outputFile << kmer.sequence << '\t'
                        << kmer.rsid << '\t'
                        << kmer.alt_allele << '\t'
                        //<< kmer.alternative_allele_frequency << '\n';
                        << kmer.alternative_allele_frequency << '\t'
                        //<< extractChromosomeNumber(kmer.chromosome) << ':' << kmer.snp_position_genome << '\n';
                        << kmer.chromosome << ':' << kmer.snp_position_genome << '\n';
                    }
                }
            }
        }

    }

    // Close and clean up
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf_file);
    fai_destroy(fai);

    std::cout << "SNPs selected:\t" << selected_snps_count << endl;
    std::cout << "K-mers generated:\t" << kmers_count << endl;

    outputFile.close();

    return 0;
}