/*
VKG : Variant K-mer Generator

Programme de test pour générer un index de k-mer, version c++.
(Version Python : ksv.py)

Crée tous k-mers porteurs de variants possibles à partir du fichier vcf de dbSNP,
et du génome de référence HG38.

*/

#include <vector>
#include <cstdlib>
#include <iostream>
#include <htslib/vcf_sweep.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

std::vector<std::string> extractAltRef(char *als) {
    std::vector<std::string> res;
    std::string str = "";
    int i = 0;
    int j = 0;
    while(j != 2) {
        if(als[i] == '\0') {
            res.push_back(str);
            str.clear();
            j++;
        } else {
            str += als[i];
        }
        i++;
    }
    return res;
}

int main() {
    // Test sur le chromosome Y
    const char* vcf_file_path = "../data/test.vcf";

    // Open VCF file for reading
    htsFile *vcf_file = bcf_open(vcf_file_path, "r");
    if (!vcf_file) {
        std::cerr << "Error opening VCF file: " << vcf_file_path << std::endl;
        return 1;
    }

    // Initialize VCF header and record
    bcf_hdr_t *vcf_header = bcf_hdr_read(vcf_file);
    bcf1_t *vcf_record = bcf_init();

    // Read each record in the VCF file
    bcf_unpack(vcf_record, BCF_UN_STR);
    while (bcf_read(vcf_file, vcf_header, vcf_record) == 0) {
        // Access VCF fields using bcf_get functions
        // For example, to get the chromosome (CHROM) and position (POS):
        std::string chromosome = bcf_hdr_id2name(vcf_header, vcf_record->rid); // Chromosome - OK
        int position = vcf_record->pos + 1; // 1-based position in VCF - Position - OK
        //std::vector<std::string> altRef = extractAltRef(vcf_record->d.als);
        char* refalt = vcf_record->d.als;

        // Process other VCF fields as needed
        // ...
        // Print or process the data
        std::cout << "Chromosome: " << chromosome << ", Position: " << position << std::endl;
        //std::cout << "ID : " << rs_id << std::endl;
        //std::cout << "REF " << altRef[0] << std::endl;
        //std::cout << "ALT " << altRef[1] << std::endl;
        std::cout << refalt << std::endl;
    }

    // Close and clean up
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf_file);

    return 0;
}
