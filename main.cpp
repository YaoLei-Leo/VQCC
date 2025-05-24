#include <iostream>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include "ExtractGVCF.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <vcf_file> <out_file>" << endl;
        return 1;
    }

    const char* vcf_path = argv[1];
    const char* out_path = argv[2];
    BGZF* output = bgzf_open(out_path, "w");

    std::string line = "CHROM\tPOS\tREF\tALT\tQUAL\tRAW_DP\tRAW_MQ\tMQRankSum\tReadPosRankSum\tSB\n";
    
    if(bgzf_write(output, line.c_str(), line.size()) < 0) {
        std::cerr << "Failed to write header to BGZF file.\n";
        bgzf_close(output);
        return 1;
    }

    extract_vqsr_metrics(vcf_path, output);

    bgzf_close(output);

    return 0;
}