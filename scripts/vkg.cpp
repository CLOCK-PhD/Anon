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
#include <htslib/hts.h>
#include <htslib/vcf_sweep.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/faidx.h>

using namespace std;

int main() {
    // Test vcf complet
    const char* vcf_file_path = "/home/remycosta/phd/Anon/data/dbsnp156/GCF_000001405.40.gz";
    // test chromosome Y
    //const char* vcf_file_path = "/home/remycosta/phd/Anon/data/snp_latest/divers/test.vcf.gz";

    // TEST FASTA
    int kmer_size = 21;
    const char* fasta_file_path = "/home/remycosta/phd/Anon/data/grch38p14/GCF_000001405.40_GRCh38.p14_genomic.fna";
    //cout << "before load" << endl;
    //system("ls -al /home/remycosta/phd/Anon/data/grch38p14/");
    faidx_t* fai = fai_load(fasta_file_path);
    if(!fai){
        cerr << "Error opening fasta file: " << fasta_file_path << endl;
        return 1;
    }
    //cout << "after load" << endl;
    //system("ls -al /home/remycosta/phd/Anon/data/grch38p14/");

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
    while (bcf_read(vcf_file, vcf_header, vcf_record) == 0) {
        // UNPACK - MANDATORY
        //bcf_unpack(vcf_record, BCF_UN_STR);
        //bcf_unpack(vcf_record, BCF_UN_INFO);
        bcf_unpack(vcf_record, BCF_UN_ALL);

        // Accessing INFO_FIELD_NAME
        // COMMON
        bcf_info_t *info_common = bcf_get_info(vcf_header, vcf_record, "COMMON");
        
        // VC (Variant Class)
        bcf_info_t *info_vc = bcf_get_info(vcf_header, vcf_record, "VC");

        // Accessing vcf fields
        // Get chromosome
        string chromosome_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
        
        // Get position
        int position = vcf_record->pos; // 1-based position in VCF
        
        // Get REF (vcf_record->d.allele[0] is REF other is ALT)
        string ref = vcf_record->d.allele[0];
        
        // Get ID
        char *rsid = vcf_record->d.id;

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
        std::cout << std::endl;*/

        // check if is snp - OK
        //std::cout << "Is SNP : " << bcf_is_snp(vcf_record) << std::endl;

        // TEST RECUPERATION INFO - OK SI AVANT FASTA - EN FAIT C'EST AUTRE CHOSE MAIS CEST OKER MAINTENANT
        //std::cout << "INFOS" << std::endl;
        // print VC CLASS - OK
        //string vc;
        //string vc_string = (char*)(info_vc->vptr); // ICI ÇA MARCHE PAS
        //cout << "VARIANT CLASS: " << info_vc->vptr << endl;
        const char* var_class = (char *)(info_vc->vptr); // ICI ÇA MARCHE
        //cout << var_class << endl;
        /*if (info_vc && info_vc->type == BCF_BT_CHAR){
            //cout << "bonjour" << endl;
            string vc_string = (char*)(info_vc->vptr);
            //cout << "VARIANT CLASS :\t" <<vc_string << endl;
        }*/

        // FREQUENCIES - OK
        /*bcf_info_t *info_freqs = bcf_get_info(vcf_header, vcf_record, "FREQ");
        if (info_freqs && info_freqs->type == BCF_BT_CHAR){
            cout << "FREQUENCIES : " << endl;
            const char *freqs = (char*)(info_freqs->vptr);
            cout << freqs << endl;
        }*/

        // IS COMMON - OK
        /*if (info_common){
            cout << "COMMON : True" << endl;
        } else {
            cout << "COMMON : False" << endl;
        }*/

        // TEST FASTA - OK ! Mais bloque VC si on le met avant alors faut laisser ça à la fin.
        //cout << "TEST PRINT SEQUENCE" << endl;
        int *len;
        // TEST POS
        int start = position;
        int end = position;
        //int start = position-(kmer_size-1);   // for TEST POS
        //int end = position+(kmer_size+1);     // for TEST POS
        const char *reg = chromosome_name.c_str();
        char *sequence = faidx_fetch_seq(fai, reg, start, end, len);
        if(!sequence){
            cerr << "Error: Could not fetch the sequence" << endl;
            return 1;
        }

        // TEST POS
        //cout << "VERIF : " << ref << " " << sequence << endl;   // for TEST POS
        //cout << sequence << endl;
        string seq = sequence;
        if (ref[0] != toupper(seq[0])){
            cerr << "WARNING" << endl;
            cerr << "\t" << ref[0] << "\t" << seq[0] << endl;
            cerr << "\t" << ref << "\t" << seq << endl;
            cerr << chromosome_name << " " << rsid << endl;
        }


        // Free memory space
        free(sequence);
    }

    // Close and clean up
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf_file);
    fai_destroy(fai);

    return 0;
}
