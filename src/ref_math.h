#ifndef _REF_MATH
#define _REF_MATH

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

long unsigned binom(unsigned int n,unsigned int k);

double Q2P(unsigned int q);

double cQ2P(char c);

unsigned int P2Q(double p);

double probRandom(setting s, unsigned int k, unsigned int refLength, unsigned int length);

unsigned int mapping_quality(setting s, unsigned int k, unsigned int refLength, unsigned int length);
#endif
