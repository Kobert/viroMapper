//All functions that manipulate or access the hashtable or entries go here.
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "referenceAssembly.h"
#include "ref_slidingWindow.h"
#include "ref_hash.h"
#include "ref_io.h"
#include "ref_math.h"
#include "ref_arg.h"

//------------UTILITIES-FOR-ENTRIES-IN-HASHTABLE-------------------------------------------------------------------------------------------






#ifdef _singleUnsigned
void initEntry(hashEntry* entry, unsigned int * val)
{
entry->value=*val;
entry->repeats=1;
entry->next=-1;
}
#endif
#ifdef _multiUnsigned
void initEntry(hashEntry* entry, unsigned int * val)
{

unsigned int i;

entry->inHitList = 0;
entry->allocatedHitList = 10; 
// entry->allocatedHitList = 120; 
entry->hitList = (unsigned int *)calloc(entry->allocatedHitList, sizeof(unsigned int));
entry->value = (unsigned int *)calloc(_numMulti, sizeof(unsigned int));

for(i=0;i<_numMulti;i++)
{
entry->value[i]=val[i];
}

entry->repeats=1;
entry->next=-1;

}
#endif
  
#ifdef _singleUnsigned
void freeTable(hashEntry* entryList, unsigned int length)	//Not a clean solution, function should simply not be called if singleUnsigned is set...
{

}  
#endif
#ifdef _multiUnsigned
void freeTable(hashEntry* entryList, unsigned int length)
{
  
unsigned int i;
for(i=0; i<length; i++)
{
  hashEntry entry = entryList[i];
  
  free(entry.hitList);
  free(entry.value);
}
  
}
#endif  
  
#ifdef _singleUnsigned 
void printEntry(hashEntry e)
{
printWindow(&(e.value));
printf("Repeated %u times\n",e.repeats);

}
#endif
#ifdef _multiUnsigned 
void printEntry(hashEntry e)
{
//printWindow(&(e.value));
printWindow(e.value);
printf("Repeated %u times\n",e.repeats);

}
#endif

#ifdef _singleUnsigned
unsigned int sEQe(unsigned int* s,hashEntry e)	//Function to check whether a slidingwindow s matches an entry from the hashtable.
{						//Does not do much right now, but is anticipated to grow more complex.

if(*s==e.value)
return 1;

return 0;
}
#endif
#ifdef _multiUnsigned
unsigned int sEQe(unsigned int* s,hashEntry e)	//Function to check whether a slidingwindow s matches an entry from the hashtable.
{						//Does not do much right now, but is anticipated to grow more complex.

unsigned int i;

for(i=0; i<_numMulti; i++)
{
  if(s[i] != e.value[i])
    return 0;
}
return 1;

  
}
#endif
//-----------------------------------------------------------------------------------------------------------------------------





//-----------HASHING-FUNCTIONS-------------------------------------------------------------------------------------------------
#ifdef _singleUnsigned
unsigned int hashUI(unsigned int *value, unsigned int hashValue, unsigned int hashTableSize)	//Simple hashing function to hash unsigned integer values to some more managable range.
												//Must be changed to something more complex later...
{
  
unsigned long result = (((*value)*hashValue) % hashTableSize); //Use long unsigned int here to avoid problem with integer limits
return (unsigned int)result;
}
#endif
#ifdef _multiUnsigned
unsigned int hashUI(unsigned int *value, unsigned int hashValue, unsigned int hashTableSize)	//Simple hashing function to hash unsigned integer values to some more managable range.
{						//Must be changed to something more complex later...

unsigned int i;
unsigned long result=0;			//Use long unsigned int here to avoid problem with integer limits
for(i=0;i<_numMulti;i++)
{
result += value[i]*hashValue; 
}

return (unsigned int)(result % hashTableSize);

}
#endif



//If the hash position is not assigned yet, simply add the current element to it with default starting values.
void addNewToHashtable(unsigned int pos, unsigned int *slidingWindow, unsigned int sequencePosition, hashEntry **hashTable,  hashEntry *entryTable, unsigned int *itemsInTable)
{
 		hashEntry e;
		initEntry(&e,slidingWindow);
		
		e.hitList[0] = sequencePosition;
		e.inHitList++;
		entryTable[*itemsInTable]=e;
		hashTable[pos]=&entryTable[*itemsInTable];
		(*itemsInTable)++; 

}


//Check if there are chained elements. If none match, add a new entry to the list.
void chainHash(unsigned int *slidingWindow, unsigned int sequencePosition, hashEntry *entryTable, unsigned int *itemsInTable, hashEntry *e, unsigned int *hit, unsigned int *chained)
{
 		unsigned int found = 0;
		hashEntry f =*e;
		hashEntry *h= e;

			while(f.next >= 0)
			{

			h = &entryTable[f.next];
			f = entryTable[f.next];	
				if(sEQe(slidingWindow,f))	
				{
				h->repeats++;
				h->hitList[h->inHitList] = sequencePosition;
				h->inHitList++;
                                
                                        if(h->inHitList >= h->allocatedHitList)
                                        {
                                        h->hitList = (unsigned int *)realloc(h->hitList, (h->allocatedHitList * 2) * sizeof(unsigned int));
                                        h->allocatedHitList = h->allocatedHitList * 2;
        
                                        // 	assert(h->inHitList < h->allocatedHitList); //TODO This assertion must be replaced by a check for reallocating the memory!!!!!
                                        }
        
				
				(*hit)++;
				found = 1;
//				printEntry(*h);
				}	

			}


			//Case no matching pattern was found:
			if(found == 0)
			{
			hashEntry g;
			initEntry(&g,slidingWindow);
			g.hitList[0] = sequencePosition;
			g.inHitList++;//TODO can be set to 1?
			assert(g.inHitList < g.allocatedHitList); 
                        //TODO really h not g? TODO Can be removed without fault since set to 1. Otherwise this assertion must be replaced by a check for reallocating the memory!!!!! 
                
			
			entryTable[*itemsInTable]=g;
			h->next	 = *itemsInTable;		
			(*itemsInTable)++;
			(*chained)++;
			} 
}


//Hashing procedure. 
//-Get hash position
//-check if position is empty
//  -add to position
//-check for equality with first element
//  -add to frequency
//-check for equality with chained elements
//  -add to frequency
//-No match
//  -add to list and establish chain
void addToHash(unsigned int *slidingWindow, unsigned int sequencePosition, hashEntry **hashTable, hashEntry *entryTable, 
	       unsigned int *itemsInTable, unsigned int *hit, unsigned int *chained, unsigned int hashValue, unsigned int hashTableSize)
{
  
unsigned int pos = hashUI(slidingWindow, hashValue, hashTableSize);
	if(hashTable[pos]==NULL)
	{
		  
	  addNewToHashtable(pos, slidingWindow, sequencePosition, hashTable,  entryTable, itemsInTable);
		  
	}else{
		  
	hashEntry *e=hashTable[pos];


	if(sEQe(slidingWindow,*e))	//If we find a matching entry, simply count one more repeat
        {
	
//	printEntry(*e);
        e->repeats++;
        e->hitList[e->inHitList] = sequencePosition;
	e->inHitList++;
        
        if(e->inHitList >= e->allocatedHitList)
        {
        e->hitList = (unsigned int *)realloc(e->hitList, (e->allocatedHitList * 2) * sizeof(unsigned int));
        e->allocatedHitList = e->allocatedHitList * 2;
        
// 	assert(e->inHitList < e->allocatedHitList); //TODO This assertion must be replaced by a check for reallocating the memory!!!!!
        }
        
        
        
        (*hit)++;
	}else{
		  
	  chainHash(slidingWindow, sequencePosition, entryTable, itemsInTable, e, hit, chained);

	}
        //NOTE this seems to have been a gigantic mistake!
//         chainHash(slidingWindow, sequencePosition, entryTable, itemsInTable, e, hit, chained);
	  
	}
	
}
//---------------------------------------------------------------------------------------------------------



//----Setup-Hashtable---------------------------------------------------------------------------------------


int populateHashTable(setting * arg, globalVariables* globalVar, char *seq, hashEntry ***hashTable, hashEntry **entryTable, unsigned int *itemsInTable)
{
  return populateHashTableWithKey(arg, globalVar, seq, hashTable, entryTable, itemsInTable, _hashValue, _hashTableSize);
}

int populateHashTableWithKey(setting * arg, globalVariables* globalVar, char *seq,  hashEntry ***hashTable, hashEntry **entryTable, unsigned int *itemsInTable, unsigned int hashValue, unsigned int hashTableSize)
{
  //Test that no mutually exclusive variables have been set
  int _assertion = 0;
  #ifdef _singleUnsigned
  _assertion++;
  #endif
  #ifdef _multiUnsigned
  _assertion++;
  #endif
  assert(_assertion==1);
  assert(basesPerByte * bitsPerBase == 8);
   #ifdef _singleUnsigned
  unsigned int _numMulti=1;
   #endif
  assert(_numMulti);
  
  unsigned int verbose = arg->verbose;



int i;


unsigned int *slidingWindow;                        	//This slidingWindow will store 'sizeOfWindow' consecutive 
							//bases from the sequence.
#ifdef _singleUnsigned
  	      slidingWindow = (unsigned int *)calloc(1,sizeof(unsigned int));
	      *slidingWindow = 0;							
#endif							
#ifdef _multiUnsigned
  	      slidingWindow = (unsigned int *)calloc(_numMulti,sizeof(unsigned int));
	      for(i=0;i<_numMulti;i++)
	      slidingWindow[i] = 0;								      
#endif
unsigned int sizeOfWindow  = basesPerWindow();	 	//Each base is encoded in bitsPerBase bits, thus basesPerByte the number of
							// (8bit) bytes is the size the window.

if(verbose > 0)
print_selective("\n\t Number of bases in the window: %u\n",sizeOfWindow );


*hashTable = realloc(*hashTable, hashTableSize * sizeof(hashEntry*));//DO NOT FORGET TO FREE THIS WHERE APPROPRIATE!!! where to put *s ?

for(i=0;i<hashTableSize;i++)
(*hashTable)[i]=NULL;

//Table holding the actual entries.
//Allocate enough memory to be sure to store all information.
//Size could be reduced again later (but for some reason that didn't work well...)
*entryTable = realloc(*entryTable, globalVar->referenceSequenceLength * sizeof(hashEntry)); 
*itemsInTable=0;

unsigned int 	hit	= 0,
		chained = 0;




initWindow(slidingWindow, seq);


  if(verbose > 1)
  printf("%s\n",seq);
  
	//for(i= sizeOfWindow-1 ; i<strlen(seq)-1 ; i++)//Go through the whole sequence. -1 because we do not want to read the newline or EOF character
	for(i= sizeOfWindow-1 ; i<strlen(seq) ; i++)//Go through the whole sequence. not -1 because they are trimmed
	{
	pushToWindow(slidingWindow, seq[i]);
	    if(verbose > 1)
	    printWindow(slidingWindow);

        addToHash(slidingWindow, (i - (sizeOfWindow-1)) ,  *hashTable, *entryTable, itemsInTable, &hit, &chained, hashValue, hashTableSize);

        if(i % 10000 == 0)
            print_selective(" | completion: %f", i/(double)globalVar->referenceSequenceLength);
	}




unsigned int max=0;
if(verbose > 0)
{
for(i=0;i<*itemsInTable;i++)
{
if((*entryTable)[i].repeats>max)
max=(*entryTable)[i].repeats;
}


print_selective("\t max %u\n",max);

print_selective("\t items in table: %u, hit: %u, chained %u\n", *itemsInTable,hit, chained);
}

//*entryTable = realloc(*entryTable,(*itemsInTable)* sizeof(hashEntry));//TODO find out why this didn"t work... 

free(slidingWindow);



return chained;
}


//-----------------------------------------------------------------------------------------------------------------------


//--Look-up-hashtable-entries--------------------------------------------------------------------------------------------
int getMapPostition(setting s, resultsVector* rv, unsigned int *hitList, unsigned int *whichHit, unsigned int howManyPos, unsigned int* max_fragments)
{
//       print_selective("1 gMP\n");
  unsigned int max = 0;
  int whichMax = -1;
  unsigned int i;
  
  unsigned int count = 0;
  
  unsigned int 	foundFirst  = 0,
		foundSecond = 0;
  unsigned int  multi_max   = 0;
 //unsigned int   foundThird = 0;		
//   print_selective("Pos: %d ", howManyPos);
  for(i = 0; i<howManyPos; i++)
  {
   

//       print_selective("1 for\n");
//       
//       print_selective("i: %d:\n ", i);
//       print_selective("wH d: %d:\n ", whichHit[i]);
//       print_selective("wH u: %u:\n ", whichHit[i]);
//       print_selective("hL: %d:\n ", hitList[whichHit[i]]);
//       
//       
//       
//       print_selective("%d: %d,\n ", whichHit[i], hitList[whichHit[i]]);
       if(hitList[whichHit[i]] >= s.minFracs)
    {
        count++;
      //if(foundSecond)
	//foundThird = 1;
      
     if(!foundFirst)
     {
       foundFirst=1;
     }else{
       foundSecond=1; 
       result* r = &(rv->results[whichHit[i]]);
       r->numIndels++;
     }
      
    }
    
//       print_selective("2 for\n");
    if(hitList[whichHit[i]] > max)
    {
     max =  hitList[whichHit[i]];
     whichMax = whichHit[i];
    }
    
//       print_selective("3 for\n");
  }
  
  
//       print_selective("2 gMP\n");
//   print_selective("\n");
//   if(count>2)
//   {
//         for(i = 0; i<howManyPos; i++)
//   {
//     if(hitList[whichHit[i]] >= s.minFracs)
//     {
// //         print_selective("%u (%u) ",whichHit[i], hitList[whichHit[i]]);
//     }
//   }
// //   print_selective("count: %u\n",count);
//   }
  
  
  *max_fragments = max;
 // if(foundThird)
 //   printf("------------------------------------------------------------------------------------------------------------\n");
  if(foundSecond)
  {
  //    for(i = 0; i<howManyPos; i++)
 // {
 //   printf("%u\t ", whichHit[i]);
 // }
 //   printf("\n");
        for(i = 0; i<howManyPos; i++)
  {
 //   printf("%u\t ", hitList[whichHit[i]]);
      if(hitList[whichHit[i]] == max)
    {
        multi_max++;
        if(multi_max > 1) //If the placement cannot be uniquely defined discard the read
        {
//             printf("countthis\n");
            multi_hit++;
            return -1;
        }
    }
  }
 // printf("\n");
 //   printf("\n");
 //   rv->indels++;
 //   return -1; //WARNING This means that reads with indels are simply discarded!!!!!!!!!!!!!
  }
  
//       print_selective("3 gMP\n");
 // if(foundThird)
 //   printf("------------------------------------------------------------------------------------------------------------\n");
  
  
	if(max >= s.minFracs) //TODO find actual criteria! 2 is somewhat good already. 	length ref: 10700, length frag 100, num frag 2, p value: 0.000162873
	{	      //							length ref: 10700, length frag 400, num frag 2, p value: 0.0029322845 
		      //			3 is better 			length ref: 10700, length frag 400, num frag 3, p value: 0.0000007169 
		      //             max > 3 for refLength 10700, sequencelength <= 500 gives 10^-9 (4*10^-10)
		      //		Results regarding coverage do not change much compared to max >2. However, less Zika is assembled for Dengue reference -> ++
	if(whichMax < 10724)
        good_hit++;
        else
        bad_hit++;
            
        return whichMax;
	}
  
  return -1;
}

int checkMismatches(setting s, unsigned int *hitList, unsigned int *whichHit, unsigned int howManyPos, char*seq, unsigned int maxMiss)
{
  unsigned int max = 0;
  unsigned int i;
  	
  int target;
  
  target = floor((double)strlen(seq)/ basesPerWindow() ) - maxMiss;
  
  if(target<0)
    return 0;
  
  
  for(i = 0; i<howManyPos; i++)
  {

    if(hitList[whichHit[i]] > max)
    {
     max =  hitList[whichHit[i]];
    }
  }
  

  
	if(max >= s.minFracs && (max >= target) ) 
	{
	  return 1;
	}
  
  return 0;
}

//For a fragment, check if it matched a hash entry. If it does, note down the possible positions that result from this.
int checkHash(hashEntry ** hashTable, hashEntry * entryTable, unsigned int * window,  
	      unsigned int *hitTable, unsigned int **whichHitP, unsigned int * num_allocated_whichHit, unsigned int *howManyPos, unsigned int fragmentNumber)
{
    unsigned int *whichHit = *whichHitP;
//     printf("in cH\n");
  unsigned int i;
   unsigned int pos = hashUI(window, _hashValue, _hashTableSize);
   int sequencePosition;
    
   if(hashTable[pos] != NULL)
   {
    
     	hashEntry *e=hashTable[pos];

	if(sEQe(window,*e))	
	{

	  assert(e->inHitList > 0);
	  for(i = 0; i<e->inHitList; i++)
	  {
	    sequencePosition = e->hitList[i] - fragmentNumber * basesPerWindow(); 
	    if(sequencePosition < 0)//TODO actually check overlap here!
	    { 
	      continue;
	    }
	    
// 	    int reallocated = 0;
	    
	    if(hitTable[sequencePosition ] == 0)
	    {
                while(*howManyPos >= *num_allocated_whichHit)  //A simple if() might also suffice instead of while(). But this should be safer still.
                {
                    *num_allocated_whichHit = *num_allocated_whichHit *2;
//                     printf("pointer before: %p\n", whichHit);
                    whichHit = (unsigned int *)realloc(whichHit, *num_allocated_whichHit * sizeof(unsigned int));
                   assert(whichHit && "whichHit could not be reallocated");
//                  unsigned int *_temp_whichHit = (unsigned int *)realloc(whichHit, *num_allocated_whichHit * sizeof(unsigned int));
//                      assert(_temp_whichHit && "whichHit could not be reallocated");
//                free(whichHit);
//                whichHit = _temp_whichHit;
//                       printf("pointer after: %p\n", whichHit);
//                     reallocated = 1;
//                     printf("REALLOCATED: %d, sequencePosition: %d, howManyPos: %d\n", *num_allocated_whichHit, sequencePosition, *howManyPos);
                
                    *whichHitP = whichHit;
                }
                
	     whichHit[*howManyPos] = sequencePosition;
	  
	    (*howManyPos)++; 
            
	    }
	    hitTable[sequencePosition ]++;

//             if(reallocated)
//                 printf("After reallocation element 0: %d\n", whichHit[0]);
	  }
// 	   printf("After forloop element 0: %d\n", whichHit[0]);
	return 1;
	}else{

		hashEntry f =*e;
		

			while(f.next >= 0)
			{

			
			f = entryTable[f.next];	
				if(sEQe(window,f))	
				{
				  assert(f.inHitList > 0);
				  for(i = 0; i<f.inHitList; i++)
				  {
				    sequencePosition = f.hitList[i] - fragmentNumber * basesPerWindow(); 
				  	    if(sequencePosition < 0)
					    { 
					    continue;
					    }
					    
					    if(hitTable[sequencePosition] == 0)
					    {
                                                
                                                while(*howManyPos >= *num_allocated_whichHit)  //A simple if() might also suffice instead of while(). But this should be safer still.
                                                {
                                                *num_allocated_whichHit = *num_allocated_whichHit *2;
                                                whichHit = (unsigned int *)realloc(whichHit, *num_allocated_whichHit * sizeof(unsigned int));
                                                assert(whichHit && "whichHit could not be reallocated (else clause)");
                                                *whichHitP = whichHit;    
                                                }
                                                
					    whichHit[*howManyPos] = sequencePosition;
	  
					    (*howManyPos)++;	      
					    }
				  hitTable[sequencePosition]++;
				  
				  }
				return 1;				
				}	

			}



	}
   }
   return 0;  
}

// From http://stackoverflow.com/questions/122616/how-do-i-trim-leading-trailing-whitespace-in-a-standard-way
char *trim(char *str)
{
    size_t len = 0;
    char *frontp = str;
    char *endp = NULL;

    if( str == NULL ) { return NULL; }
    if( str[0] == '\0' ) { return str; }

    len = strlen(str);
    endp = str + len;

    /* Move the front and back pointers to address the first non-whitespace
     * characters from each end.
     */
    while( isspace(*frontp) ) { ++frontp; }
    if( endp != frontp )
    {
        while( isspace(*(--endp)) && endp != frontp ) {}
    }

    if( str + len - 1 != endp )
            *(endp + 1) = '\0';
    else if( frontp != str &&  endp == frontp )
    {           
      *str = '\0';
    }
    /* Shift the string so that it starts at str so that if it's dynamically
     * allocated, we can still free it on the returned pointer.  Note the reuse
     * of endp to mean the front of the string buffer now.
     */
    endp = str;
    if( frontp != str )
    {
            while( *frontp ) { *endp++ = *frontp++; }
            *endp = '\0';
    }


    return str;
}


int getNextNonEmptyLine(FILE * file, char ** read, size_t * bytes)
{

//     print_selective("in getNextNonEnptyLine\n");
//     fflush(stderr);
    
 unsigned int i;  

  unsigned int foundBeginning = 0;
  
  
  int seqLength =  getline(read, bytes,file);
    
//   print_selective("seqLength: %d, read %s\n", seqLength, *read);

        while(!foundBeginning && seqLength > 0)
	{
//         print_selective("while getNextNonEnptyLine\n");
	i = 0;
	  while((*read)[i] == ' ' || (*read)[i] == '\t')
	  {
	  i++;
	  }
	  
	  if((*read)[i] == '\n')
	  {
	  seqLength =   getline(read, bytes,file);
	  }else{
    
	    if((*read)[i] == EOF || seqLength <= 0 )
	    return 0;
    
	  foundBeginning = 1;
	  }
	  
// 	print_selective("i %d read \"%s\", seqLength=%d\n", i, *read, seqLength);
	}
 if(seqLength<1)
     return 0;

    //Simply delete leading and trailing whitespaces
    trim(*read);

    if(*read[0] == 0)//If the next non empty string is empty, none exists
	{
	//*read[1] = 0;
	//*read[0] = ' ';
	return 0;
	}
  
  return 1;
}

int getNthLine(int n, FILE * file, char ** read, size_t * bytes)
{
//     print_selective("in getNthLine\n");
//     fflush(stderr);
    
  assert(n>0);
   unsigned int i;
   
   for(i = 0; i<n; i++)
   {
     if(!getNextNonEmptyLine(file, read, bytes))
       return 0;
   }
   
   return 1;
}

/*
int getNextRead(FILE * file, char ** read, unsigned int * bytes, unsigned int * firstPosition)
{
 
 unsigned int i;  

  unsigned int foundBeginning = 0;
  
  int seqLength =  getline(read, bytes,file);
  
  
  unsigned int isName = 0;
  unsigned int isPlus = 0;
  
  unsigned int running = 1;
  
  while(!foundBeginning)
  {
    i = 0;
  while((*read)[i] == ' ' || (*read)[i] == '\t')
  {
    i++;
  }
  if((*read)[i] == '\n')
  {
    seqLength =   getline(read, bytes,file);
  }else{
    
  if((*read)[i] == EOF || seqLength <= 0)
    return 0;
    
    
  if((*read)[i] == '@' || (*read)[i] == '>')
  {
    isName = 1;
  }
  if((*read)[i] == '+')
  {
    isPlus = 1;
  }
  foundBeginning = 1;
  }
  }
 
  
  if(isName)
  {
   return getNextRead(file, read, bytes, &i);
  }
    if(isPlus)
  {
    getNextRead(file, read, bytes, &i);
    running = getNextRead(file, read, bytes, &i);
    
  }


  trim(*read); //TODO either rewrite with trimmed string in mind, or get firstposition to work
 if(*read[0] == 0)//An empty string crashes the program somewhere downstream, so simply avoid this....
 {
   *read[1] = 0;
   *read[0] = ' ';
   return 0;
 }
  *firstPosition = i;

  return running;
  
}

int getNextQuality(FILE * file, char ** read, unsigned int * bytes, unsigned int * firstPosition)
{
 
 unsigned int i;  

  unsigned int foundBeginning = 0;
  
  int seqLength =  getline(read, bytes,file);
  
  
//  unsigned int isName = 0;
  unsigned int isPlus = 0;
//  unsigned int isQ    = 0;
  
  unsigned int running = 1;
  
  while(!foundBeginning)
  {
    i = 0;
  while((*read)[i] == ' ' || (*read)[i] == '\t')
  {
    i++;
  }
  if((*read)[i] == '\n')
  {
    seqLength =   getline(read, bytes,file);
  }else{
    
  if((*read)[i] == EOF || seqLength <= 0)
    return 0;
    
    

  if((*read)[i] == '+')
  {
    if((*read)[i + 1] == ' ' || (*read)[i + 1] == '\n' || (*read)[i + 1] == '\t')
    isPlus = 1;

  }
  foundBeginning = 1;
  }
  }
  
    if(isPlus)
  {
   return getNextQuality(file, read, bytes, &i);
  }


  trim(*read); //TODO either rewrite with trimmed string in mind, or get firstposition to work
 if(*read[0] == 0)//An empty string crashes the program somewhere downstream, so simply avoid this....
 {
   *read[1] = 0;
   *read[0] = ' ';
 }
 
  *firstPosition = i;
  return running;
  
}
*/

int placeFragments(setting s, globalVariables* g, resultsVector* rv, hashEntry ** hashTable, hashEntry * entryTable, char * seq, unsigned int firstPosition, 
			    unsigned int * hit, unsigned int *miss, unsigned int * hitPerRound, unsigned int *max_fragments)
{
    
// print_selective("in PF\n");
   unsigned int i, pos;
   
   unsigned int howManyPos = 0,  
		localHitPerRound= 0,
		length 		= 0,
		lhit 		= 0,
		lmiss 		= 0,
		fragments 	= 0;
		
    unsigned int *window;
  
    unsigned int num = 1;
    #ifdef _multiUnsigned
    num = _numMulti;
    #endif
  
    unsigned int num_allocated_whichHit = 250;
    unsigned int * whichHit = (unsigned int *)calloc(num_allocated_whichHit, sizeof(unsigned int));
    
    window = (unsigned int*)malloc(num * sizeof(unsigned int));
    
		
   unsigned int * hitList = (unsigned int *)calloc(strlen(g->referenceSequence),sizeof(unsigned int)); 
   
   
  length = strlen(seq) - firstPosition;// Right now no -1 since strings are trimmed beforehand!
  assert(length > 0);
  
  
  fragments  = floor( (double)length/(double)basesPerWindow() );
  if(fragments <= 0)
  {
    fprintf(stderr, "fragments %u, length %u\n", fragments, length);
  assert(fragments > 0); //TODO if this actually happens, the read is to short and should be ignored/discarded
  }
  
// printf("element 0 before: %d\n", whichHit[0]);
  for(i=0; i<fragments ;i++)
  {
    pos = i * basesPerWindow() + firstPosition;
    makeWindow(window, seq, pos);
    
  if(checkHash(hashTable, entryTable, window, hitList, &whichHit, &num_allocated_whichHit, &howManyPos, i))
  {
      localHitPerRound++;
   lhit++;
    }else{
     lmiss++; 
    }
//    if(i==0)
//        printf("element 0 during: %d\n", whichHit[0]);

  } 
  
// printf("element 0 after: %d\n", whichHit[0]);

  *miss 	+= lmiss;
  *hit 		+= lhit;
  *hitPerRound	 = localHitPerRound;
  
//  print_selective("1 PF\n");
  free(window);
  
assert(whichHit[0] <= strlen(g->referenceSequence) && "whichHit has too large an entry");

  int result = getMapPostition(s, rv, hitList, whichHit, howManyPos, max_fragments); 
     
  //if( !checkMismatches(s, hitList, whichHit, howManyPos, seq, 3) )
  //{
  //  result = -1;
  //}
//     print_selective("2 PF\n");
  free(hitList); //TODO this is the problem!
  
  free(whichHit);

  return result;
}






//void hashMapReadsWithKey(setting arg,  hashEntry **hashTable, hashEntry *entryTable, unsigned int itemsInTable, unsigned int hashValue, unsigned int hashTableSize);
void hashMapReadsWithKey(setting arg, globalVariables *globalVar, resultsVector *rv, unsigned int hashValue, unsigned int hashTableSize);
void cleanHashMapReadsWithKey(setting arg, globalVariables *globalVar, resultsVector *rv, unsigned int hashValue, unsigned int hashTableSize);

void hashMapReads(setting arg, globalVariables *globalVar, resultsVector *rv)
{
   cleanHashMapReadsWithKey( arg, globalVar, rv, _hashValue, _hashTableSize);
}
/*
void hashMapReadsWithKey(setting arg, globalVariables *globalVar, resultsVector *rv, unsigned int hashValue, unsigned int hashTableSize)
{
  
    FILE *readsFile;
  readsFile = fopen(arg.readsFile,"r");
  if(!readsFile){
  readsFile = fopen64(arg.readsFile,"r");//For 32 bit systems this may solve the issue
  }
  if(!readsFile){
    printf("\n[ERROR:] Could not open file %s\n", arg.readsFile);
      printf("         Please make sure the path is correct and you have read permission.\n");
      assert(readsFile);
  }
  
  
  unsigned int firstPosition = 0;

  unsigned int bytesToRead = 10;
  char * seq;
  seq = (char*)malloc(bytesToRead );
  

  int reversePosition = -1;
	
  unsigned int hit  = 0;
  unsigned int miss = 0;
  unsigned int mapped  = 0;
  unsigned int unmapped = 0;
  
  unsigned int suspects = 0;
  unsigned int hitPerRound = 0;
  unsigned int one = 0, two = 0, three = 0, four = 0;
  unsigned int numReads = 0;
  
  	//while(getNextRead(readsFile, &seq, &bytesToRead, &firstPosition))
  
  unsigned int readDiff = 2;//During the first iteration of the while loop only read the name and then the read
	while(getNthLine(readDiff, readsFile, &seq, &bytesToRead))
	{
	readDiff = 4;//Reads are expected to fall 4 non-empty lines apart. (in between are "+", quality profile, and name of next read @)
	numReads ++;

//	printf("read: %s\n", seq);
//	printf("length: %u\n", strlen(seq));	
	int whichPos = placeFragments(arg, rv, globalVar->hashTable, globalVar->entryTable, seq, firstPosition, &hit, &miss, &hitPerRound);
	int firstAlt;
	 	char * complement;
		complement = strdup(seq);
		reverseComplementSequence(complement, seq);
                firstAlt = placeFragments(arg, rv, globalVar->complementHashTable, globalVar->complementEntryTable, complement, firstPosition, &hit, &miss, &hitPerRound);	 
		if(firstAlt >= 0)
		firstAlt = globalVar->referenceSequenceLength - (firstAlt + strlen(seq));
		//if(whichPos < 0 && firstAlt >= 0)
		//{
		//printf("Yes! %u\n", whichPos);
		//}
	
			if(whichPos != firstAlt && (whichPos >= 0 && firstAlt >= 0))
		{
		  suspects++;
		  if(0 && abs(whichPos - firstAlt) > 1 )
		  {
		  printf("\nPos %d, firstAltPos %d \n", whichPos, firstAlt);
		printf("%s\n", seq);
		if(firstAlt >= 0)
		printf("%.*s\n", strlen(seq), globalVar->referenceSequence + firstAlt);
		if(whichPos >= 0)
		  printf("%.*s\n", strlen(seq), globalVar->referenceSequence + whichPos);
		  }
		//		reversePosition = placeFragments(arg, globalVar->hashTable, globalVar->entryTable, complement, firstPosition, &hit, &miss, &hitPerRound);
		}
		
		
		if( whichPos < 0 || whichPos > (globalVar->referenceSequenceLength - strlen(seq)) )//TODO -1 can also simply mean that the querry starts before the reference, may have to check this...
		{    
		  
		//char * complement;
		//complement = strdup(seq);
		//reverseComplementSequence(complement, seq);
		
		whichPos = placeFragments(arg, rv, globalVar->complementHashTable, globalVar->complementEntryTable, seq, firstPosition, &hit, &miss, &hitPerRound);  
		reversePosition = -1;
		
		if(whichPos > 0)
		reversePosition = globalVar->referenceSequenceLength - (whichPos + strlen(seq));
		
	//	int altPos = placeFragments(arg, rv, globalVar->hashTable, globalVar->entryTable, complement, firstPosition, &hit, &miss, &hitPerRound);
	//	if(altPos != reversePosition && (reversePosition >=0 && altPos >= 0 ))
	//	{
	//	  suspects++;
	//	  if(0 && abs(reversePosition - altPos) > 1 )
	//	  {
	//	  printf("\nRevPos %d, altPos %d \n", reversePosition, altPos);
	//	printf("%s\n", complement);
	//	if(altPos >= 0)
	//	printf("%.*s\n", strlen(seq), globalVar->referenceSequence + altPos);
	//	if(reversePosition >= 0)
	//	  printf("%.*s\n", strlen(seq), globalVar->referenceSequence + reversePosition);
	//	  }
		//		reversePosition = placeFragments(arg, globalVar->hashTable, globalVar->entryTable, complement, firstPosition, &hit, &miss, &hitPerRound);
	//	}
		if( reversePosition < 0 || reversePosition > (globalVar->referenceSequenceLength - strlen(seq)) )//TODO carefully check this condition
		{    
		unmapped++;
		}else{
		  
//		  		printf("seq;\n%s\n",seq);
//		printf("complement;\n%s\n",complement);
//				printf("reference:\n%.*s\n", strlen(seq), globalVar->referenceSequence + reversePosition);
				
		//revCompReadToResult(rv, seq, reversePosition, strlen(seq));
		  		  if(arg.doTrim)
		  {
		  unsigned int an,bn;
		  trimRegions(arg, globalVar, complement, reversePosition, &an, &bn);
		  readToResult(rv, &(complement[an]), reversePosition+an, bn -an); 
		  }else{
		readToResult(rv, complement, reversePosition, strlen(seq)); 
		  }
		  
		readDiff = 2; // Set to two because 2 lines will be read now!
		getNthLine(2, readsFile, &seq, &bytesToRead);
		revQualityToResult(rv, seq, reversePosition, strlen(seq));
		
		mapped++; 
		//printf("numreads %u , whichPos %d, revPos %d\n", numReads, whichPos, reversePosition);
		//printf("[ReVCompl] fraction of hits within sequence%f\n", evalComplementStatisticExpensive(&(globalVar.referenceSequence[reversePosition]), seq, strlen(seq) -1));
		 		

		}
		

		free(complement);
		
		}else{
		//  printf("numreads %u , whichPos %d, revPos %d, ref length %u, seq length %u\n", numReads, whichPos, reversePosition, globalVar.referenceSequenceLength, strlen(seq));
		//  printf("[Regular] fraction of hits within sequence%f\n", evalStatisticExpensive(&(globalVar.referenceSequence[whichPos]), seq, strlen(seq) -1));
		  if(arg.doTrim)
		  {
		  unsigned int an,bn;
		  trimRegions(arg, globalVar, seq, whichPos, &an, &bn);
		  readToResult(rv, &(seq[an]), whichPos+an, bn -an);
		  
		  }else{
		  
		  readToResult(rv, seq, whichPos, strlen(seq));
		  }
		readDiff = 2; // Set to two because 2 lines will be read now!
		getNthLine(2, readsFile, &seq, &bytesToRead);
		qualityToResult(rv, seq, whichPos, strlen(seq));
		
		mapped++; 
		}
  
		if(hitPerRound > 3)
		four++;
		if(hitPerRound > 2)
		three++;
		if(hitPerRound > 1)
		two++;
		if(hitPerRound > 0)
		one++;
  
	}

	printf("Number of suspects %u of %u (%f)\n", rv->indels, numReads, (double)rv->indels/numReads);
 
  printf("fraction of mapped %f (%u / %u)\n",mapped/(double)(mapped+unmapped), mapped, (mapped+unmapped));
  printf("fraction of hits %f (%u / %u)\n",hit/(double)(hit+miss), hit, (hit+miss));
  printf("fraction of one hits %f , two %f, three %f, four %f\n",one/(double)(numReads), two/(double)(numReads), three/(double)(numReads), four/(double)(numReads));
 

  free(seq);
  fclose(readsFile);
}

*/





void cleanHashMapReadsWithKey(setting arg, globalVariables *globalVar, resultsVector *rv, unsigned int hashValue, unsigned int hashTableSize)
{
  
  //Open File and check that it worked
      FILE *readsFile;
  readsFile = fopen(arg.readsFile,"r");
  if(!readsFile){
  readsFile = fopen64(arg.readsFile,"r");//For 32 bit systems this may solve the issue
  }
  if(!readsFile){
    fprintf(stderr, "\n[ERROR:] Could not open file: \"%s\"\n", arg.readsFile);
      fprintf(stderr, "         Please make sure the path is correct and you have read permission.\n");
      assert(readsFile);
  }
  
  FILE *samFile;
  samFile = fopen("_temp_placeholder_sam.sam","w");
  printSamHeader(samFile, "Reference_Name", globalVar->referenceSequenceLength);
  
  //Allocate space for short reads and complement thereof
 // unsigned int bytesToRead = 10;
  size_t bytesToRead = 10;
  char * seq;
  seq = (char*)malloc(bytesToRead );

  size_t bytesToRead_name = 10;
  char * seq_name;
  seq_name = (char*)malloc(bytesToRead_name );
  
//unsigned int bytesToReadQ = 10;
  size_t bytesToReadQ = 10;
  char * seqQ;
  seqQ = (char*)malloc(bytesToReadQ );
  
  unsigned int 	numReads = 0,
		numMatches = 0,
		numForwardMatches = 0,
		numReverseMatches = 0,
		numShiftedMatches = 0;
  
  
  unsigned int 	hit = 0,
		miss = 0,
		hitPerRound = 0;
  
  // wichPosition stores the position of the reference genome, where the left most base of the short read is aligned.
  int whichPos = -1;
  
  // Variables to keep track of whether a matching pattern has been found, and what configuration was used. 
  unsigned int 	matchForward = 0,
		matchReverse = 0,
		matchShifted = 0;
  unsigned int match = 0;
  
  unsigned int max_fragments = 0, mappingQuality = 255;
  
//     unsigned int readDiff = 2;//During the first iteration of the while loop only read the name and then the read
  unsigned int readDiff = 1;//During the first iteration of the while loop only read the name and then the read
	// This loop goes over all short read's dna data
	while(getNthLine(readDiff, readsFile, &seq_name, &bytesToRead_name))
	{
//              printTime();
//              print_selective(" in while: Name %s\n", seq_name);
//              fflush(stderr);
	  getNthLine(readDiff, readsFile, &seq, &bytesToRead);
	  char * complement;
	  complement = strdup(seq);
	  
	  char * result;
	  
	  unsigned int shift = 0;
	  //complement = (char*)calloc(strlen(seq) + 1 , sizeof(char));  
	  //complement[strlen(seq)] = '\0';
	  
// 	readDiff = 4;	// Reads are expected to fall 4 non-empty lines apart. 
			// (in between are "+", quality profile, and name of next read starting with "@")
	readDiff = 3;	// Reads are expected to fall 4 non-empty lines apart, but we want the name next. 
	numReads ++;
	
	matchForward = 0;
	matchReverse = 0;
	matchShifted = 0;
	match = 0;
        max_fragments = 0;
        mappingQuality = 255;
	
	// Find the left most position where the read aligns to the reference.
	// The first try is performed on the regular sequence with no offset.
	whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, seq, 0, &hit, &miss, &hitPerRound, &max_fragments);
	
	  if( whichPos >= 0 && whichPos <= (globalVar->referenceSequenceLength - strlen(seq)) )
	  {
	    result = seq;
	      matchForward = 1;
	      match = 1;
	  }
		
	  // Next, test whether the reverse complemented sequence matched the reference.
	  if(!match)
	  {
	   reverseComplementSequence(complement, seq);
    
	   if(strlen(complement) != strlen(seq))
	   {
	   fprintf(stderr, "\nstrlen(complement) = %lu, strlen(seq) =%lu\n", (long unsigned)strlen(complement), strlen(seq));
	   assert(strlen(complement) == strlen(seq));
	   }
	   
	   whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, complement, 0, &hit, &miss, &hitPerRound, &max_fragments);
	  
	
	  if( whichPos >= 0 && whichPos <= (globalVar->referenceSequenceLength - strlen(seq)) )
	  {
	    result = complement;
	      matchReverse = 1;
	      match = 1;
	  }		
	  }

	  // Now, test whether shifting it such that the fragments are right aligned gives a valid matching.
	  // If it was already right aligned during the first iteration, use a centered method instead.
	  if(!match  && arg.doExtendedTesting)
	  {

	    shift = strlen(seq) % basesPerWindow();
	    if(shift == 0)
	    {
	    //  shift = floor(basesPerWindow()/2.0);
	    }
	    
	   whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, seq, shift, &hit, &miss, &hitPerRound, &max_fragments);
	   whichPos = whichPos - shift;
	   
	  if( whichPos >= 0 && whichPos <= (globalVar->referenceSequenceLength - strlen(seq)) )
	  {
	    result = seq;	    
	      matchForward = 1;
	      match = 1;
	  }		
	  }	

	  // Lastly, do the shifting for the reverse complemented sequence
	  if(!match && arg.doExtendedTesting)
	  {

	    
	   whichPos = placeFragments(arg, globalVar, rv, globalVar->hashTable, globalVar->entryTable, complement, shift, &hit, &miss, &hitPerRound, &max_fragments);
	   whichPos = whichPos - shift;
	   
	  if( whichPos >= 0 && whichPos <= (globalVar->referenceSequenceLength - strlen(seq)) )
	  {
	    result = complement;	    
	      matchReverse = 1;
	      match = 1;
	  }		
	  }	
	  
	  
	  // Finally, when all relevant configurations have been testes, or one good configuration found, evaluate the alignment
	  if(match)
	  {
              
            mappingQuality = mapping_quality(arg, max_fragments, globalVar->referenceSequenceLength, strlen(seq));
//	    printf("Position: %d\n", whichPos);
	    numMatches++;
	    
	    if(matchForward)
	    {
	      numForwardMatches++;
	    }
	    if(matchReverse)
	    {
	      numReverseMatches++; 
	    }
	    if(matchShifted)
	    {
	      numShiftedMatches++;
	    }
	    
	    
	    if(arg.doTrim)
		  {
		  unsigned int an,bn;
		  trimRegions(arg, globalVar, result, whichPos, &an, &bn);
		  readToResult(rv, &(result[an]), whichPos+an, bn -an); 
                    
                  
//             printf("%lu %d %.2f %.2f\n", strlen(seq), (bn-an), evalStatisticExpensive(&(globalVar->referenceSequence[whichPos]), result, strlen(seq)), evalStatisticExpensive(&(globalVar->referenceSequence[whichPos+ an]), &(result[an]), bn -an));
	//    printf("%.*s\n", bn -an, globalVar->referenceSequence+(whichPos+an));
//         printf("%.*s\n", (int)strlen(seq), globalVar->referenceSequence+(whichPos));
//            printf("%.*s\n\n", bn-an, result+an);
// 	     printf("%s",result);
  
//              int diff_NOTNEEDED = 151 - strlen(seq);
//              int i_NOTNEEDED;
//              for(i_NOTNEEDED = 0; i_NOTNEEDED<diff_NOTNEEDED ;i_NOTNEEDED++)
//              {
//                  printf("N");
//              }
// 	     printf("\n");   
//                 printDifferences(&(globalVar->referenceSequence[whichPos]), result, strlen(seq));
            
             readDiff = 2;
	          getNthLine(readDiff, readsFile, &seqQ, &bytesToReadQ);
		    if(matchReverse)
		    {
		    reverseSequence(seqQ, seqQ);  
		    }
	          qualityToResult(rv, &(seqQ[an]), whichPos+an, bn -an);
// 		  printf("%s\n\n", seqQ);
                  if(bn > an)
                    printSamLine( (whichPos+an), "Reference_Name", &seq_name[1], &(result[an]), &(seqQ[an]), (bn -an),  matchReverse, mappingQuality, NULL); 
//                   sumDifferentQs(&(globalVar->referenceSequence[whichPos]), result, seqQ,  strlen(seq));
                  
		    if( arg.storeReads )
		    {
		//    storeSequence(globalVar, rv, &(result[an]), whichPos+an, bn -an);
                     storeBrokenSequence(globalVar, rv, &(result[an]), whichPos+an, bn -an);
		     storeBrokenSequenceShifted(globalVar, rv, &(result[an]), whichPos+an, bn -an);
		    } 

		  
		  }else{
		readToResult(rv, result, whichPos, strlen(seq)); 
	    
	        readDiff = 2;
		getNthLine(readDiff, readsFile, &seqQ, &bytesToReadQ);
		    if(matchReverse)
		    {
		    reverseSequence(seqQ, seqQ);  
		    }
		  qualityToResult(rv, seqQ, whichPos, strlen(seq)); 
		  

                    printSamLine( whichPos, "Reference_Name", &seq_name[1], result, seqQ, strlen(result),  matchReverse, mappingQuality, NULL); 

		    if( arg.storeReads )
		    {
		  //  storeSequence(globalVar, rv, result, whichPos, strlen(seq));
                     storeBrokenSequence(globalVar, rv, result, whichPos, strlen(seq));
                     storeBrokenSequenceShifted(globalVar, rv, result, whichPos, strlen(seq));
                     
                    } 
		  
		  }
	  }else{ 
//	  printf("Unmatched Position: %d\n", whichPos);
            readDiff = 2;
            getNthLine(readDiff, readsFile, &seqQ, &bytesToReadQ);
	    printSamLine( -1, "Reference_Name", &seq_name[1], seq, seqQ, strlen(seq),  0, mappingQuality, NULL); 

	    
	  }
	  readDiff = 1; // We should be at quality and want to advance to the name.
	  free(complement);
          
//           print_selective(" end while\n");
//           fflush(stderr);
	}
  
//        print_selective(" out of while\n");
//        fflush(stderr);
       
  print_selective("\n\tNumber of reads: %u, number of matching reads: %u (%.2f%%)\n", numReads, numMatches, 100.0*(double)numMatches/numReads);
  print_selective("\t(Forward %.2f%%, reverse %.2f%% (Shifted: %.2f%%))\n", 100.0*(double)numForwardMatches/numReads,
								 100.0*(double)numReverseMatches/numReads, 
								 100.0*(double)numShiftedMatches/numReads);
  
  //Free allocated data structures
  fclose(samFile);
  
  fclose(readsFile);
  
  free(seq_name);
  free(seqQ);
  free(seq);
  
}







int testReferenceBySize(globalVariables* g, setting s, unsigned int seqSize)
{
    
int isGoodReference = 1;//One denotes a good reference sequence

int i, whichPos;
//unsigned int seqSize;

unsigned int hit=0, miss=0, hitPerRound=0;
unsigned int max_fragments = 0;
//TODO do dynamic
//seqSize = 250;

char * seq;

 resultsVector testRv;
 
 initResults(&testRv, strlen(g->referenceSequence));
 

seq = (char *)calloc(seqSize + 1, sizeof(char));

    for(i = 0 ; i< g->referenceSequenceLength - seqSize ; i++ )
    {
    memcpy( seq, &(g->referenceSequence[i]), seqSize );
    
    whichPos = placeFragments(s, g, &testRv, g->hashTable, g->entryTable, seq, 0, &hit, &miss, &hitPerRound, &max_fragments);
    if(whichPos != i && whichPos >= 0)
    {
        isGoodReference = 0;

  //   printf("Not hit %d != %d\n", i, whichPos);   
             break;
    }
    }
    
    
     
         for(i = 0 ; i< g->referenceSequenceLength - seqSize ; i++ )
    {
               result r = testRv.results[i];
    if(r.numIndels>0)
    {
        isGoodReference = 0;

    //   printf("Different insertionpoints found (%d: %u).\n", i, r.numIndels); 
               break;
    }
        
        
    }
    
    
 
     
 

     
freeResults(&testRv);

free(seq); 
 //   printf("returning %u\n", isGoodReference);
    return isGoodReference;
       
}


void printGoodNumber(setting s,unsigned int num)
{
    printf("Test concluded\n");
    printf("\nThe given reference sequence :\n%s\n", s.referenceFile);        
    printf("may be a candidate for processing with this tool.\n");    
            
    printf("Good behavior can however only be expected for short reads\n");
    printf("with up to %u characters each!\n", num);
    printf("\nLonger reads may cause suboptimal results.\n\n");
    
}


void printBadReference(setting s)
 {
        printf("Test concluded\n");
        printf("\nThe given reference sequence :\n%s\n", s.referenceFile);  
        printf("is likely a BAD candidate for this tool!\n");
        printf("Please consider using an alternative software, or proceed with caution.\n\n");  
 }
 
void printPerfectReference(setting s)
{
        printf("Test concluded\n");
        printf("\nThe given reference sequence :\n%s\n", s.referenceFile);        
       printf("is likely a GOOD candidate for processing with this tool.\n\n");    
}

int testReference(globalVariables* g, setting s)
{
    int i, isGoodReference = 0, lastGood=-1;
    unsigned int seqSize = 50;
    
    for(i = 50; i < g->referenceSequenceLength; i = i *2)
    {
        seqSize = i;
        
    if(testReferenceBySize(g, s, seqSize))
    {
     isGoodReference = 1;        
     lastGood = seqSize;
    }else{
        if(isGoodReference)
        {
         printGoodNumber(s, lastGood);
         return 1;
        }else{
         printBadReference(s);
         return 0;
        }
        
    }
    
    }
    
    printPerfectReference(s);
    
    return 1;
}









//----STRINGS-------------------------------------------------------------------------------

//From http://www.cse.yorku.ca/~oz/hash.html

//unsigned int stringHash(char *str, unsigned int hashTableSize)
//{
  
//        unsigned long hash = 5381;
//        char c;

//        while (c = *str++)
//            hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

//        return (unsigned int)(hash % hashTableSize);
//    }
//}






