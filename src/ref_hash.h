#ifndef _HASH   /* Include guard */
#define _HASH



void freeTable(hashEntry* entryList, unsigned int length);

//void addToHash(unsigned int *slidingWindow, unsigned int sequencePosition, hashEntry **hashTable, hashEntry *entryTable, unsigned int *itemsInTable, unsigned int *hit, unsigned int *chained, unsigned int hashValue, unsigned int hashTableSize);
int testReference(globalVariables* g, setting s);

int populateHashTable(setting * arg, globalVariables* globalVar, char *seq, hashEntry ***hashTable, hashEntry **entryTable, unsigned int *itemsInTable);

int populateHashTableWithKey(setting * arg, globalVariables* globalVar, char *seq,  hashEntry ***hashTable, hashEntry **entryTable, unsigned int *itemsInTable, unsigned int hashValue, unsigned int hashTableSize);

int getNextNonEmptyLine(FILE * file, char ** read, size_t * bytes);

int getNextRead(FILE * file, char ** read, unsigned int * bytes, unsigned int * firstPosition);

int getNextQuality(FILE * file, char ** read, unsigned int * bytes, unsigned int * firstPosition);

void hashMapReads(setting arg, globalVariables *globalVar, resultsVector *rv);
#endif
