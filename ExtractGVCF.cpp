#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <cstring>
#include "ExtractGVCF.hpp"

using namespace std;

// Define a function to extract the metrics from GVCF file.
int extract_vqsr_metrics(const char* filename, BGZF* output) {
    htsFile *vcf = bcf_open(filename, "r"); // Open the VCF file.
    if (!vcf) {
        cerr << "Error opening VCF file: " << filename << endl;
        return 1;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(vcf); // Read the header of the VCF file.
    bcf1_t *rec = bcf_init(); // Initialize a BCF record.

    while (bcf_read(vcf, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        vector<string> alts_vec; // Get ALTs
        for (int i=1; i < rec -> n_allele; ++i) {
            alts_vec.push_back(rec->d.allele[i]);
        }

        if (alts_vec[0] != "<NON_REF>") {
            const char* chrom = bcf_hdr_id2name(hdr, rec->rid); // Get CHROM
            int pos = rec->pos + 1; // Get POS
            string ref = rec->d.allele[0]; // Get REF

            float qual = rec->qual; // Get QUAL

            // Get RAW_MQandDP
            int32_t* raw_mq_dp = nullptr;
            int n_raw_mq_dp = 0;
            int32_t raw_mq_val;
            int32_t raw_dp_val;
            int ret = bcf_get_info_int32(hdr, rec, "RAW_MQandDP", &raw_mq_dp, &n_raw_mq_dp);
            if (ret > 0 && n_raw_mq_dp >= 2) {
                raw_mq_val = raw_mq_dp[0];
                raw_dp_val = raw_mq_dp[1];
            }
            else {
                raw_mq_val = -9999; // Default value if not present
                raw_dp_val = -9999; // Default value if not present
            }

            // Get MQRankSum
            float* mq_ranksum = nullptr;
            int n_mq_ranksum = 0;
            float mq_ranksum_val;
            ret = bcf_get_info_float(hdr, rec, "MQRankSum", &mq_ranksum, &n_mq_ranksum);
            if (ret > 0 && n_mq_ranksum > 0) {
                mq_ranksum_val = mq_ranksum[0];
            }
            else {
                mq_ranksum_val = -9999.0; // Default value if not present
            }

            // Get ReadPosRankSum
            float* read_pos_ranksum = nullptr;
            int n_read_pos_ranksum = 0;
            float read_pos_ranksum_val;
            ret = bcf_get_info_float(hdr, rec, "ReadPosRankSum", &read_pos_ranksum, &n_read_pos_ranksum);
            if (ret > 0 && n_read_pos_ranksum > 0) {
                read_pos_ranksum_val = read_pos_ranksum[0];
            }
            else {
                read_pos_ranksum_val = -9999.0; // Default value if not present
            }

            // Get SB
            int32_t* sb = nullptr;
            int sb_count = 0;
            vector<int32_t> sb_vec;
            ret = bcf_get_format_int32(hdr, rec, "SB", &sb, &sb_count);
            if (ret > 0 && sb) {
                for (int i=0; i<4; ++i) {
                    sb_vec.push_back(sb[i]);
                }
            }

            // Output
            string line;
            line = string(chrom) + "\t" + to_string(pos) + "\t" + string(ref) + "\t";

            for (size_t i = 0; i < alts_vec.size(); ++i) {
                if (i) line = line + ",";
                line = line + string(alts_vec[i]);
            }
            
            line = line + "\t" + to_string(qual) + "\t" + to_string(raw_dp_val) + "\t" + to_string(raw_mq_val) + "\t" + to_string(mq_ranksum_val) + "\t" + to_string(read_pos_ranksum_val) + "\t";

            for (size_t i = 0; i < sb_vec.size(); ++i) {
                if (i != 0) {
                    line = line + "," + to_string(sb_vec[i]);
                }
                else{
                    line = line + to_string(sb_vec[i]);
                }
            }

            line = line + "\n";
            
            if (bgzf_write(output, line.c_str(), line.size()) < 0) {
                std::cerr << "Failed to write " << line << " to BGZF file.\n";
                return 1;
            }

            free(mq_ranksum);
            free(read_pos_ranksum);
            free(raw_mq_dp);
            free(sb);
        }
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);

    return 0;
}