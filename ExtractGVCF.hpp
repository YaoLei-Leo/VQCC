#ifndef EXTRACT_GVCF
#define EXTRACT_GVCF

int extract_vqsr_metrics(const char* filename, BGZF* output);

#endif