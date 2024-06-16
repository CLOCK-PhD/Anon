#include "Kmer.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <htslib/hts.h>
#include <htslib/vcf_sweep.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/faidx.h>

Kmer::Kmer(std::string seq, std::string chrom, size_t id, size_t start, size_t end, size_t snp_gen_pos, size_t snp_kmer_pos, char ref, std::string alt, float ref_freq, float alt_freq)
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

void Kmer::display(){
        std::cout << "K-mer sequence: " << sequence << "\n"
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

bool isDegenerate(char base) {
    // Define a function to check if a base is degenerated
    return (base == 'R' || base == 'Y' || base == 'S' || base == 'W' || base == 'K' ||
            base == 'M' || base == 'B' || base == 'D' || base == 'H' || base == 'V' || base == 'N');
}

std::vector<Kmer> generate_kmers(const std::string &umer, size_t &k, const std::string &chromosome, const size_t &rsid, size_t &snp_position_genome, char &ref_allele, std::string &alt_allele, float &ref_frequency, float &alt_frequency){
    size_t uLength = umer.length();
    size_t i = 0;
    std::vector<Kmer> kmer_objects;
    size_t snp_position_umer = k; // SNP position in the u-mer, center position

    while (i <= uLength - k){
        bool containsDegenerate = false;
        for(size_t j = i; j < i + k; j++){
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
            size_t kmer_start = snp_position_genome - snp_position_umer + i;
            size_t kmer_end = kmer_start + k;
            size_t snp_position_kmer = snp_position_umer - i;
            // Creating k-mer object
            Kmer kmer_obj(kmer, chromosome, rsid, kmer_start, kmer_end, snp_position_genome, snp_position_kmer, ref_allele, alt_allele, ref_frequency, alt_frequency);
            kmer_objects.push_back(kmer_obj);
        }
        i++; // Move to the next base for the next k-mer
    }

    // Print each Kmer object in the vector
    /*for (Kmer& kmer : kmer_objects) {
        kmer.display();
    }*/

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