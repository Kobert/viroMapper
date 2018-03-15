#ifndef _SIM_QC   /* Include guard */
#define _SIM_QC

#include <stdio.h>



// #define _singleUnsigned
#define _multiUnsigned
#define _numMulti 1


//setup for A==0001, C==0010, G==0100, T==1000
// #define basesPerByte 2
// #define bitsPerBase 4

//setup for A==00, C==01, G==10, T==11, DOES NOT ACCOUNT FOR DEGENERATE BASES!
#define basesPerByte 4
#define bitsPerBase 2

//TODO This is experimental. See if it works or remove...
#define _limitBasesPerWindow 8
#define _perfectList


#ifdef _perfectList
#define _hashTableSize 65537

#else
//#define _hashTableSize 10
//#define _hashTableSize 1327
//#define _hashTableSize 11777
//#define _hashTableSize   1002017 
//#define _hashTableSize 1252159
// #define _hashTableSize 100003
// #define _hashTableSize 10030037
#define _hashTableSize 100300037

#endif

#define _hashValue 37
//#define _hashValue 1000861



//#ifdef _32BIT
FILE* fopen64(const char *filename, const char *type);
//#endif

int good_hit;
int bad_hit;
int multi_hit;

typedef struct
{
unsigned int verbose;

//contains the name of the first reference file (only for historic reasons)
char* referenceFile;

//contains all reference files
char** multi_referenceFiles;
int num_multi_references;

char* readsFile;

char* dataInFile;
char* outFilePrefix;
int qFloor;

unsigned int printSAM;

unsigned int doCsvFile;
unsigned int doGnuplotFile;
unsigned int doJSFile;
unsigned int concealResults;

unsigned int executeReferenceOnly;

unsigned int writeConsensus;
	
unsigned int numPlots;
//double pValue;

unsigned int doTrim;
unsigned int doExtendedTesting;

unsigned int storeReads;

unsigned int mapOnly;
unsigned int primes;

unsigned int windows;
unsigned int minFracs;

// Variable denotes usage of old and depreciated output layouts
unsigned int LEGACY;
}setting;


#ifdef _singleUnsigned
typedef struct
{
unsigned int value;
unsigned int repeats;

unsigned int *hitList;
unsigned int allocatedHitList;
unsigned int inHitList;

int next;

}hashEntry;
#endif
#ifdef _multiUnsigned
typedef struct
{
unsigned int *value;
unsigned int repeats;

unsigned int *hitList;
unsigned int allocatedHitList;
unsigned int inHitList;

int next;

}hashEntry;
#endif



typedef struct
{
char * seq;
unsigned int pos;

unsigned int length;

char matchesReference;
char isSubstring;

}Read;


typedef struct
{
    int number_of_placements;
    
    unsigned int has_indel;
    unsigned int has_insertion;
    unsigned int has_deletion;
    
    int start_leftmost;
    int start_clipping;
    
    int start_indel;
    int indel_length;
    
    int end_clipping;
    int end_rightmost;
    
    int num_left_match;
    int num_right_match;
    
    
}Placement;


typedef struct
{
char * seq;

unsigned int hits;
unsigned int * isEndFor;
unsigned int * isStartFor;
unsigned int * weight;

//unsigned int length;
  
}Haplotype;

typedef struct
{
  char * referenceSequence;
  unsigned int referenceSequenceLength; 
  
  char ** reference_names;
  unsigned int * references_individual_lengths;
  unsigned int * references_individual_start;
  int number_of_references;
  
  int * break_points; //Points where one reference ends and another begins
  

 //hash table of POINTERS! for entries.
 hashEntry** hashTable;
 //Table holding the actual entrier.
 hashEntry* entryTable;
 unsigned int itemsInTable;

 
 Read * read;
 unsigned int numReads;
 unsigned int allocatedNumReads;
 
 Haplotype * altSequences;
 unsigned int numAltSequences;
 unsigned int allocatedNumAltSequences;
 
 char *** brokenSequences;
 unsigned int * numBrokenSequences;
 unsigned int * allocatedNumBrokenSequences;
 
  char *** brokenSequencesShifted;
 unsigned int * numBrokenSequencesShifted;
 unsigned int * allocatedNumBrokenSequencesShifted;
 
 long long unsigned trimmed;
 long long unsigned kept;
 
 double avgRatioTrimmed;
 unsigned int avgTrimmedOf;
  
 unsigned int number_of_too_short_reads;
 
}globalVariables;


typedef struct
{
  unsigned int A;
  unsigned int C;
  unsigned int G;
  unsigned int T;
  unsigned int N;
  
  double qA;
  double qC;
  double qG;
  double qT;
  double qN;
  
  double ML_A;
  double ML_C;
  double ML_G;
  double ML_T;
  double ML_N;
  
  unsigned int coverage;
  unsigned int qFloorCoverage;
  
  unsigned int * indivError;

  unsigned int indivErrorLength;
  
  char majorBase;
  char secondBase_char;
  unsigned int majorBaseCount;
  unsigned int secondBaseCount;
  double frequencyMajorBase;
  double frequencySecondBase;
  
  double frequencyMajorBase_quality_corrected;
  double frequencySecondBase_quality_corrected;
  
  
  unsigned int meanPhError;
  unsigned int medianPhError;
  unsigned int lowerPhQuartile;
  unsigned int upperPhQuartile;
  unsigned int minPhError;
  unsigned int maxPhError;
  
  double meanError;
  double medianError;
  double minError;
  double maxError;
  
  unsigned int numIndels;
  
  

}result;

typedef struct
{
  result * results;
   unsigned int assignedLength;
   
   double averageCoverage;
   double averageQFloorCoverage;
   
   unsigned int indels;
   
   
   unsigned int numReferenceMatches;

   unsigned int hit;
   unsigned int miss;
   
}resultsVector;

//int populateHashTableWithKey(setting * arg, unsigned int hashValue, unsigned int hashTableSize);

typedef struct
{
  unsigned int * e;
  unsigned int num;  
  unsigned int allocatedNum;
  
  int maxLength;
}MultiList;





#endif
