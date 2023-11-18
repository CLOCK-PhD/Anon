#include <vector>
#include <cstdlib>
#include <iostream>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <iostream>

void read_vcf(char *fname)
{
    // open vcf file
    htsFile *fp = hts_open(fname, "rb");
    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec = bcf_init();
    // save for each vcf record
    while (bcf_read(fp, hdr, rec) >= 0)
    {
        // unpack for read REF,ALT,INFO,etc
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);
        // read CHROM
        std::cout << bcf_hdr_id2name(hdr, rec->rid) << std::endl;
        // read POS
        std::cout << rec->pos << std::endl;
        // read ID(rsid)
        std::cout << rec->d.id << std::endl;
        // rec->d.allele[0] is REF other is ALT
        std::cout << rec->n_allele << std::endl;
        for (int i = 0; i < rec->n_allele; ++i)
            std::cout << rec->d.allele[i] << std::endl;
        // check if is snp
        std::cout << "ttt " << bcf_is_snp(rec) << std::endl;
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ((ret = hts_close(fp)))
    {
        std::cerr << "hts_close(" << fname << "): non-zero status " << ret << std::endl;
        exit(ret);
    }
}

int main(){
    char* fname = "../../data/dbsnp156/GCF_000001405.40.gz";

    // open vcf file
    htsFile *fp = hts_open(fname, "rb");
    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec = bcf_init();
    // save for each vcf record
    while (bcf_read(fp, hdr, rec) >= 0)
    {
        // unpack for read REF,ALT,INFO,etc
        //bcf_unpack(rec, BCF_UN_STR);
        //bcf_unpack(rec, BCF_UN_INFO);
        bcf_unpack(rec, BCF_UN_ALL);
        // read CHROM
        std::cout << bcf_hdr_id2name(hdr, rec->rid) << std::endl;
        // read POS
        std::cout << rec->pos << std::endl;
        // read ID(rsid)
        std::cout << rec->d.id << std::endl;
        // rec->d.allele[0] is REF other is ALT
        std::cout << rec->n_allele << std::endl;
        for (int i = 0; i < rec->n_allele; ++i)
            std::cout << rec->d.allele[i] << std::endl;
        // check if is snp
        std::cout << "ttt " << bcf_is_snp(rec) << std::endl;
    }
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ((ret = hts_close(fp)))
    {
        std::cerr << "hts_close(" << fname << "): non-zero status " << ret << std::endl;
        exit(ret);
    }

    return 0;
}