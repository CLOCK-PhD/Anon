// Kmer.h

#ifndef KMER_H
#define KMER_H

#include <string>
#include <vector>

/**
 * @brief Represents a k-mer with associated SNP information.
 * 
 * This class stores detailed information about a k-mer, including its sequence,
 * associated chromosome, SNP details, and allele frequencies.
 */
class Kmer {
public:
    std::string sequence; ///< K-mer nucleotide sequence
    std::string chromosome; ///< Chromosome where the k-mer is located
    size_t rsid; ///< Reference SNP ID from dbSNP, without "rs" prefix
    size_t kmer_start; ///< Start position of the k-mer in the genome
    size_t kmer_end; ///< End position of the k-mer in the genome
    size_t snp_position_genome; ///< Position of the SNP within the genome
    size_t snp_position_kmer; ///< Position of the SNP within the k-mer sequence
    char ref_allele; ///< Reference allele at the SNP position
    std::string alt_allele; ///< Alternative allele at the SNP position
    float reference_allele_frequency; ///< Frequency of the reference allele in the population
    float alternative_allele_frequency; ///< Frequency of the alternative allele in the population

    /**
     * @brief Construct a new Kmer object
     * 
     * @param seq K-mer sequence
     * @param chrom Chromosome
     * @param id SNP ID
     * @param start Start position of the k-mer in the genome
     * @param end End position of the k-mer in the genome
     * @param snp_gen_pos Position of the SNP in the genome
     * @param snp_kmer_pos Position of the SNP in the k-mer
     * @param ref Reference allele
     * @param alt Alternative allele
     * @param ref_freq Frequency of the reference allele
     * @param alt_freq Frequency of the alternative allele
     */
    Kmer(std::string seq, std::string chrom, size_t id, size_t start, size_t end, size_t snp_gen_pos, size_t snp_kmer_pos, char ref, std::string alt, float ref_freq, float alt_freq);

    /**
     * @brief Display the k-mer information to standard output.
     */
    void display();
};

/**
 * @brief Check if a character is a printable ASCII character.
 * 
 * @param c Character to check
 * @return true If the character is printable
 * @return false If the character is not printable
 */
bool isPrintable(char c);

/**
 * @brief Remove non-printable ASCII characters from a string.
 * 
 * @param input String to be processed
 * @return std::string Processed string with non-printable characters removed
 */
std::string removeNonPrintableChars(const std::string& input);

/**
 * @brief Display non-printable characters in a string as hexadecimal values.
 * 
 * @param input String to be processed
 * @return std::string String with non-printable characters shown as hexadecimal
 */
std::string displayHiddenChars(const std::string& input);

/**
 * @brief Determine if a DNA base character is a degenerate nucleotide.
 * 
 * @param base DNA base to check
 * @return true If the base is degenerate
 * @return false If the base is not degenerate
 */
bool isDegenerate(char base);

/**
 * @brief Generate a vector of Kmer objects from a given u-mer sequence.
 * 
 * @param umer U-mer sequence
 * @param k Length of k-mers to generate
 * @param chromosome Chromosome of the u-mer
 * @param rsid SNP ID associated with the u-mer
 * @param snp_position_genome Genomic position of the SNP
 * @param ref_allele Reference allele at the SNP position
 * @param alt_allele Alternative allele at the SNP position
 * @param ref_frequency Frequency of the reference allele
 * @param alt_frequency Frequency of the alternative allele
 * @return std::vector<Kmer> Vector of generated k-mers
 */
std::vector<Kmer> generate_kmers(const std::string &umer, size_t &k, const std::string &chromosome, const size_t &rsid, size_t &snp_position_genome, char &ref_allele, std::string &alt_allele, float &ref_frequency, float &alt_frequency);

/**
 * @brief Split a string into a vector of strings based on a specified delimiter.
 * 
 * @param s String to split
 * @param delimiter Character used as the delimiter for splitting
 * @return std::vector<std::string> Vector of split strings
 */
std::vector<std::string> split(const std::string &s, char delimiter);

/**
 * @brief Check if a given source string is found within an input string.
 * 
 * @param my_source Source string to find
 * @param input Input string to check within
 * @return true If the source string is found
 * @return false If the source string is not found
 */
bool check_source(char* my_source, const std::string &input);

/**
 * @brief Extract and return allele frequencies from a formatted string.
 * 
 * @param input String containing allele frequency data
 * @return std::vector<float> Vector of allele frequencies
 */
std::vector<float> get_dbgap_freq(const std::string &input);

#endif
