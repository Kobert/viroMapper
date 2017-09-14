//Functions to test combinations of primenumbers for efficient hashing techniques


#include <stdio.h>
#include <stdlib.h>


#include "referenceAssembly.h"
#include "ref_arg.h"
#include "ref_hash.h"
#include "ref_primes.h"


unsigned int primesHuge(unsigned int pos)
{
 unsigned int prime[] = {
1000003,	1000033,	1000037,	1000039,	1000081,	1000099,	1000117,	1000121,	1000133,	1000151,
1000159,	1000171,	1000183,	1000187,	1000193,	1000199,	1000211,	1000213,	1000231,	1000249,
1000253,	1000273,	1000289,	1000291,	1000303,	1000313,	1000333,	1000357,	1000367,	1000381,
1000393,	1000397,	1000403,	1000409,	1000423,	1000427,	1000429,	1000453,	1000457,	1000507,
1000537,	1000541,	1000547,	1000577,	1000579,	1000589,	1000609,	1000619,	1000621,	1000639,
1000651,	1000667,	1000669,	1000679,	1000691,	1000697,	1000721,	1000723,	1000763,	1000777,
1000793,	1000829,	1000847,	1000849,	1000859,	1000861,	1000889,	1000907,	1000919,	1000921,
1000931,	1000969,	1000973,	1000981,	1000999,	1001003,	1001017,	1001023,	1001027,	1001041,
1001069,	1001081,	1001087,	1001089,	1001093,	1001107,	1001123,	1001153,	1001159,	1001173,
1001177,	1001191,	1001197,	1001219,	1001237,	1001267,	1001279,	1001291,	1001303,	1001311,
1001321,	1001323,	1001327,	1001347,	1001353,	1001369,	1001381,	1001387,	1001389,	1001401,
1001411,	1001431,	1001447,	1001459,	1001467,	1001491,	1001501,	1001527,	1001531,	1001549,
1001551,	1001563,	1001569,	1001587,	1001593,	1001621,	1001629,	1001639,	1001659,	1001669,
1001683,	1001687,	1001713,	1001723,	1001743,	1001783,	1001797,	1001801,	1001807,	1001809,
1001821,	1001831,	1001839,	1001911,	1001933,	1001941,	1001947,	1001953,	1001977,	1001981,
1001983,	1001989,	1002017,	1002049,	1002061,	1002073,	1002077,	1002083,	1002091,	1002101,
1002109,	1002121,	1002143,	1002149,	1002151,	1002173,	1002191,	1002227,	1002241,	1002247,
1002257,	1002259,	1002263,	1002289,	1002299,	1002341,	1002343,	1002347,	1002349,	1002359,
1002361,	1002377,	1002403,	1002427,	1002433,	1002451,	1002457,	1002467,	1002481,	1002487,
1002493,	1002503,	1002511,	1002517,	1002523,	1002527,	1002553,	1002569,	1002577,	1002583};
  unsigned int length=sizeof(prime)/sizeof(unsigned int);
  return prime[pos % length];
}

unsigned int primesHigh(unsigned int pos)
{
 unsigned int prime[] = {40459,  40471,  40483,  40487,  40493,  40499,  40507,  40519,  40529,  40531, \
  40543,  40559,  40577,  40583,  40591,  40597,  40609,  40627,  40637,  40639, \
  40693,  40697,  40699,  40709,  40739,  40751,  40759,  40763,  40771,  40787, \
  40801,  40813,  40819,  40823,  40829,  40841,  40847,  40849,  40853,  40867, \
  40879,  40883,  40897,  40903,  40927,  40933,  40939,  40949,  40961,  40973, \
  40993,  41011,  41017,  41023,  41039,  41047,  41051,  41057,  41077,  41081, \
  41113,  41117,  41131,  41141,  41143,  41149,  41161,  41177,  41179,  41183, \
  41189,  41201,  41203,  41213,  41221,  41227,  41231,  41233,  41243,  41257, \
  41263,  41269,  41281,  41299,  41333,  41341,  41351,  41357,  41381,  41387, \
  41389,  41399,  41411,  41413,  41443,  41453,  41467,  41479,  41491,  41507, \
  41513,  41519,  41521,  41539,  41543,  41549,  41579,  41593,  41597,  41603, \
  41609,  41611,  41617,  41621,  41627,  41641,  41647,  41651,  41659,  41669, \
  41681,  41687,  41719,  41729,  41737,  41759,  41761,  41771,  41777,  41801, \
  41809,  41813,  41843,  41849,  41851,  41863,  41879,  41887,  41893,  41897, \
  41903,  41911,  41927,  41941,  41947,  41953,  41957,  41959,  41969,  41981, \
  41983,  41999,  42013,  42017,  42019,  42023,  42043,  42061,  42071,  42073, \
  42083,  42089,  42101,  42131,  42139,  42157,  42169,  42179,  42181,  42187, \
  42193,  42197,  42209,  42221,  42223,  42227,  42239,  42257,  42281,  42283 };
  unsigned int length=sizeof(prime)/sizeof(unsigned int);
  return prime[pos % length];
}
/*
unsigned int primesHigh(unsigned int pos)
{
 unsigned int prime[] = { 
  98017,  98041,  98047,  98057,  98081,  98101,  98123,  98129,  98143,  98179, 
  98207,  98213,  98221,  98227,  98251,  98257,  98269,  98297,  98299,  98317, 
  98321,  98323,  98327,  98347,  98369,  98377,  98387,  98389,  98407,  98411, 
  98419,  98429,  98443,  98453,  98459,  98467,  98473,  98479,  98491,  98507, 
  98519,  98533,  98543,  98561,  98563,  98573,  98597,  98621,  98627,  98639, 
  98641,  98663,  98669,  98689,  98711,  98713,  98717,  98729,  98731,  98737,
  98773,  98779,  98801,  98807,  98809,  98837,  98849,  98867,  98869,  98873, 
  98887,  98893,  98897,  98899,  98909,  98911,  98927,  98929,  98939,  98947, 
  98953,  98963,  98981,  98993,  98999,  99013,  99017,  99023,  99041,  99053, 
  99079,  99083,  99089,  99103,  99109,  99119,  99131,  99133,  99137,  99139, 
  99149,  99173,  99181,  99191,  99223,  99233,  99241,  99251,  99257,  99259, 
  99277,  99289,  99317,  99347,  99349,  99367,  99371,  99377,  99391,  99397, 
  99401,  99409,  99431,  99439,  99469,  99487,  99497,  99523,  99527,  99529, 
  99551,  99559,  99563,  99571,  99577,  99581,  99607,  99611,  99623,  99643, 
  99661,  99667,  99679,  99689,  99707,  99709,  99713,  99719,  99721,  99733, 
  99761,  99767,  99787,  99793,  99809,  99817,  99823,  99829,  99833,  99839, 
  99859,  99871,  99877,  99881,  99901,  99907,  99923,  99929,  99961,  99971, 
  99989,  99991, 100003, 100019, 100043, 100049, 100057, 100069, 100103, 100109 };
  unsigned int length=sizeof(prime)/sizeof(unsigned int);
  return prime[pos % length];
}
*/


unsigned int primesLow(unsigned int pos)
{
 unsigned int prime[] = { 31,     37,     41,     43,     47,     53,     59,     61,     67,     71, 
     73,     79,     83,     89,     97,    101,    103,    107,    109,    113, 
    127,    131,    137,    139,    149,    151,    157,    163,    167,    173, 
    179,    181,    191,    193,    197,    199,    211,    223,    227,    229, 
    233,    239,    241,    251,    257,    263,    269,    271,    277,    281, 
    283,    293,    307,    311,    313,    317,    331,    337,    347,    349, 
    353,    359,    367,    373,    379,    383,    389,    397,    401,    409, 
    419,    421,    431,    433,    439,    443,    449,    457,    461,    463, 
    467,    479,    487,    491,    499,    503,    509,    521,    523,    541, 
    547,    557,    563,    569,    571,    577,    587,    593,    599,    601, 
    607,    613,    617,    619,    631,    641,    643,    647,    653,    659, 
    661,    673,    677,    683,    691,    701,    709,    719,    727,    733, 
    739,    743,    751,    757,    761,    769,    773,    787,    797,    809, 
    811,    821,    823,    827,    829,    839,    853,    857,    859,    863, 
    877,    881,    883,    887,    907,    911,    919,    929,    937,    941, 
    947,    953,    967,    971,    977,    983,    991,    997,   1009,   1013, 
   1019,   1021,   1031,   1033,   1039,   1049,   1051,   1061,   1063,   1069, 
   1087,   1091,   1093,   1097,   1103,   1109,   1117,   1123,   1129,   1151, 
   1153,   1163,   1171,   1181,   1187,   1193,   1201,   1213,   1217,   1223, 
   1229};
unsigned int length=sizeof(prime)/sizeof(unsigned int);
  return prime[pos % length];
}






void testPrimes(setting * arg, globalVariables *globalVar)
{
  setting tempArg = *arg;
  unsigned int i,j;
  
  unsigned int silent = 0;
  tempArg.verbose = silent;
  
  unsigned int hashTableSize 	= _hashTableSize;
  unsigned int hashValue	= _hashValue;
  unsigned int minValue 	= _hashValue;
  unsigned int minTableSize 	= _hashTableSize;
  
  unsigned int new;
 
  //hash table of POINTERS! for entries.
  hashEntry **hashTable; //must still be allocated in function. Where to free...?
  hashTable = (hashEntry**)calloc(1, sizeof(hashEntry*)); 
  
  //Table holding the actual entrier.
  hashEntry *entryTable; 
  entryTable = (hashEntry*)calloc(1, sizeof(hashEntry)); 
  
  unsigned int itemsInTable = 0;
    
  unsigned int min 	= populateHashTableWithKey(&tempArg, globalVar, globalVar->referenceSequence, &hashTable, &entryTable, &itemsInTable, hashValue, hashTableSize);
 
  
  for(j=0;j<200;j++)
  {
   hashTableSize = primesHigh(j);
  for(i=0; i<200; i++)
  {
    
              
	      hashValue = primesHuge(i);
	      new	= populateHashTableWithKey(&tempArg, globalVar, globalVar->referenceSequence, &hashTable, &entryTable,  &itemsInTable, hashValue, hashTableSize);
	      
	      if(new<min)
	      {
		min		= new;
		minValue	= hashValue;
		minTableSize	= hashTableSize;
	      printf("min %u, Value %u, size %u \n",min, minValue, minTableSize);
  	      }
  	   //   printf("new %u\n", new);

  }
    for(i=0; i<180; i++)
  {
    
	      hashValue = primesHigh(i);
	      new	= populateHashTableWithKey(&tempArg, globalVar, globalVar->referenceSequence, &hashTable, &entryTable,  &itemsInTable, hashValue, hashTableSize);
	      
	      if(new<min)
	      {
		min		= new;
		minValue	= hashValue;
		minTableSize	= hashTableSize;
	      printf("min %u, Value %u, size %u \n",min, minValue, minTableSize);
  	      }
  	   //   printf("new %u\n", new);
  	   
  
  }
      for(i=0; i<181; i++)
  {
    
	      hashValue = primesLow(i);
	      new	= populateHashTableWithKey(&tempArg, globalVar, globalVar->referenceSequence, &hashTable, &entryTable,  &itemsInTable, hashValue, hashTableSize);
	      
	      if(new<min)
	      {
		min		= new;
		minValue	= hashValue;
		minTableSize	= hashTableSize;
	      printf("min %u, Value %u, size %u \n",min, minValue, minTableSize);
  	      }
  	   //   printf("new %u\n", new);
  	   
  
  }
  	     printf("%u\n",j);
  }
  
  freeTable(entryTable, itemsInTable);
  free(entryTable);
  free(hashTable);
  
  printf("min %u, Value %u, size %u \n",min, minValue, minTableSize);

}
