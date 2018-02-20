#ifndef _REF_IO   /* Include guard */
#define _REF_IO

void print_selective(const char *format, ...);

void printTime();

void chomp(char* c);

void cut_at_first_space(char* c);

void readToResult(resultsVector *rv, char * seq, unsigned int pos, unsigned int length);

void revCompReadToResult(resultsVector *rv, char * seq, unsigned int pos, unsigned int length);

void readToResult_quality_aware(resultsVector *rv, char* seq, char* q, unsigned int pos, unsigned int length);

void qualityToResult(resultsVector *rv, char * seq, unsigned int pos, unsigned int length);

void revQualityToResult(resultsVector *rv, char * seq, unsigned int pos, unsigned int length);

void postProcessResults(setting arg, resultsVector *rv);

void printCSV(FILE* file, resultsVector rv);

void printCSV_quality_aware_bases(FILE *file, resultsVector rv);

void printHtml(FILE *file, setting s, resultsVector rv);

void printGnuplotDat(FILE* file, resultsVector rv);

void printGnuplotDat_quality_aware_bases(FILE* file, resultsVector rv);

void writeConsensusReferenceFile(FILE* file, resultsVector rv, char*);

void plotQualityProfile(setting arg, int length, char* tempFileName);

void plotCoverage(setting arg, char* tempFileName);

void plotBaseFrequencies(setting arg, char* tempFileName);

void plotMajorBase(setting arg, char* tempFileName);

void printAvgCoverage(setting arg, resultsVector rv);

void printSamLine(int position, int position_offset, char* referenceName, char* name, char* seq, char* qual, int lenght,  unsigned int isComplemented, unsigned int mapping_quality, Placement* placement);

void printSamHeader(char* nameReference, unsigned int referenceLength);

void printSamHeader_additional_reference(char* nameReference, unsigned int referenceLength);

char* get_reference_name_by_position(globalVariables g, int pos);

int get_reference_offset_by_position(globalVariables g, int pos);
#endif
