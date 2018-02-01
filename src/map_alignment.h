
#ifndef _MAP_ALIGN   /* Include guard */
#define _MAP_ALIGN

#define DIAG         1
#define UPOPEN       2
#define LEFTOPEN     4
#define UP           8
#define LEFT        16
#define UPEXT       32
#define LEFTEXT     64

typedef struct cstack_s
{
  char c;
  struct cstack_s * next;
} cstack_t;

int alignSingleDeletion(unsigned int length, char* ref1, char* ref2, char* read);

int alignSingleInsertion(unsigned int length, char* ref1, char* ref2, char* read, int offset);

int alignAffine(unsigned int length_a, const char * a, unsigned int length_b, const char * b);
#endif