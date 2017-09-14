#ifndef _REF_MATH
#define _REF_MATH


long unsigned binom(unsigned int n,unsigned int k);

double Q2P(unsigned int q);

unsigned int P2Q(double p);

double probRandom(setting s, unsigned int k, unsigned int refLength, unsigned int length);

unsigned int mapping_quality(setting s, unsigned int k, unsigned int refLength, unsigned int length);
#endif
