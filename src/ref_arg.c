#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "referenceAssembly.h"
#include "ref_arg.h"
#include "ref_hash.h"
#include "ref_slidingWindow.h"
#include "ref_io.h" 
#include "map_arg.h"
 
 
 void help()
 {
  printf("\n"); 
  printf("Function arguments are:\n"); 
  printf("\n"); 
  printf(" -r for specifying a file containing the reference sequence (fasta format only)\n");    
  printf(" -s for specifying a file with shortreads (fastq format only)\n");    
  printf("\n");  
  printf("\n");    
  printf(" -o specify a prefix for the output names (can include a path).\n");  
  printf("    If no Prefix is provided, only generic temporary files will be written.\n");
  printf("\n");    
  printf(" -c Print results in form of a *.csv file.\n");
  printf("    A path should be provided with -o for this option.\n");    
  printf(" -g Plot results in Gnuplot.\n");
  printf("    A path can be provided with -o for this option.\n");
  printf("    If no path is provided, no *.png files will be written.\n");
  printf(" -y Write fasta file containing the consensus sequence of mapped reads.\n");
  printf("    A path can be provided with -o for this option.\n");
  printf(" -m Map reads. I.e. do not assemble haplotypes.\n");  
  printf("\n");  
  printf("\n");    
  printf(" -q Specify a lower bound for the quality scores (integer valued Phred score). \n");  
  printf("    Any base with a lower score will not count towards the coverage. \n");  
  printf("    This option only affects the coverage, all other plots remain the same! \n");     printf(" -t Trim any suffix and prefix that does not contain 8 consecutive matches.\n");    
  printf(" -h print this help messgae\n");    
  printf("\n\n");
  printf("Developers Options:\n");    
  printf(" -d for more verbose debug output\n"); 
  printf(" -f Specify how many fractions must match at least. (default 4, must be 2 or bigger)\n"); 
  printf(" -z for testing primenumber configurations (no reads file from -s is needed)\n");   
  printf("\n");
 }
 
void my_free(void* f)
{
free(f);
f = NULL;    
}


void freeSetting(setting *s)
{
    int i;
//       print_selective("1B\n");
  if(s->referenceFile)
  {    
    my_free(s->referenceFile);
  }
  
//       print_selective("1BB\n");
  for(i = 0; i < s->num_multi_references; i++)
  {
      free(s->multi_referenceFiles[i]);
  }
  
//       print_selective("2B\n");
  if(s->num_multi_references > 0)
  {
   free(s->multi_referenceFiles);   
  }
  
//       print_selective("3B\n");
  if(s->readsFile)
  {
    free(s->readsFile);
  }
  
//       print_selective("4B\n");
    if(s->outFilePrefix)
  {
    free(s->outFilePrefix);
  }
  
//       print_selective("5B\n");
}
 
 setting initArgs(int argc, char **argv)
 {
     
     good_hit= 0;
     bad_hit = 0;
     multi_hit = 0;
     
  setting s;
  
  s.verbose = 1;
  
  int index;
  int c;
  
  s.mapOnly = 0;
  s.primes  = 0;
  
  s.storeReads = 1;

#ifdef _multiUnsigned
  s.windows = _numMulti;
#else
  s.windows = 1;
#endif
  
  s.referenceFile = NULL;
  
  s.readsFile = NULL;

  s.outFilePrefix = NULL;
  
  s.doCsvFile = 0;
  s.doGnuplotFile = 0;
  
  s.numPlots = 40;
  
  s.minFracs = 4;
  //s.pValue = -1.0;
  
  s.doTrim = 0;
  s.doExtendedTesting = 0;
  
  s.writeConsensus = 0;
	
  opterr = 0;
  
  s.qFloor = -1;


  while ((c = getopt (argc, argv, "cdef:ghmr:s:tpw:o:q:yz")) != -1)
    switch (c)
      {
      case 'h':
	help();
	exit(1);
	break;
      case 'c':
	s.doCsvFile = 1;
	break;	
      case 'e':
	s.doExtendedTesting = 1;
	break;		
      case 'f':
	assert(atoi(optarg) > 1);
	s.minFracs = atoi(optarg);
	break;		
      case 'g':
	s.doGnuplotFile = 1;
	break;
      case 'm':
	s.mapOnly = 1;
	s.storeReads = 0;
	break;
      case 'z'://Check behaviour of different prime numbers
	s.primes = 1;
        break;      
//     case 'p':
//	s.pValue = atoi(optarg);//printf("%s\n", optarg);
//        break;
      case 'r':
	s.referenceFile = strdup(optarg);//printf("%s\n", optarg);
        break;
      case 'd':
        s.verbose = 2;
        break;
      case 's':
        s.readsFile = strdup(optarg);
        break;
      case 't':
        s.doTrim = 1;
        break;
      case 'o':
        s.outFilePrefix = strdup(optarg);
        break;	
      case 'q':
        s.qFloor = atoi(optarg);
        break;	
      case 'w':
        s.windows = atoi(optarg);
	fprintf(stderr, "[ERROR:] Feature '-w' not live yet! Aborting analysis.");
	assert(0);
        break;
      case 'y':
        s.writeConsensus = 1;
	break;  
      case '?':
        if (optopt == 's' || optopt == 'r' || optopt == 'w' )
	{
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	  assert(0);
	}
	else if (isprint (optopt))
	{
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	
	  assert(0);
	}
	else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return s;
      default:
        abort ();
      }

      
	for (index = optind; index < argc; index++)
	print_selective ("Non-option argument %s\n", argv[index]);

	if(!s.referenceFile)
	{
	fprintf(stderr, "Please specify reference file with -r command\n");
	assert(s.referenceFile);
	}
	
	if(!s.readsFile && !s.primes)
	{   
	 fprintf(stderr, "Please specify reads file with -s command\n");
	assert(s.readsFile);  
	}
  

 
  if(s.doGnuplotFile)
  {
	FILE * gnuplotPipe = popen("gnuplot -persistent", "w");
	if(!pclose(gnuplotPipe))
	{

	}else{
	map_help();
	fprintf(stderr, "\n[ERROR:] Could not execute gnuplot.\n");
	fprintf(stderr, "         Please make sure gnuplot is installed,\n");
	fprintf(stderr, "         (preferably gnuplot-x11) or specify another\n");
	fprintf(stderr, "         output format.\n\n");
	assert(0);
	}
  }
 
    if(s.doGnuplotFile && !s.outFilePrefix)
  {
   print_selective("\n[Note:] Gnuplot output is demanded without specifying an output path (via -o).\n ");   
   print_selective("        Plots will be opened in a seperate window after completion of the program.\n ");
   print_selective("        No *.png files will be written.\n ");
  }
  
  
    if(s.doCsvFile && !s.outFilePrefix)
  {
   print_selective("\n[Warning:] csv output is set without specifying an output path (via -o).\n ");   
   print_selective("           A generict temporary file will be written.\n");
  }
  
    if(!s.doGnuplotFile && !s.doCsvFile)//TODO add any new output formats here
  {
   print_selective("\n[Note:] No output format specified.\n ");   
   print_selective("        Only basic statistics will be printed to the screen.\n ");
   print_selective("        To see a list of available options run the program with the -h flag.\n ");
  }  
  
  return s;
 }
 
 

  //-------------Read-Multiple-Reference-Files---------------------------------------------------------------
 //Will attempt to read arg->multi_referenceFiles and concatinate its contents to globalVar->referenceSequence
 void read_multiple_ReferenceFiles(setting *arg, globalVariables *globalVar, char* fileName)
 {
   

     FILE *referenceFile;
     referenceFile = fopen(fileName,"r");
	
        //Check that file opened succesfully
	if(!referenceFile)
	{
	fprintf(stderr, "\n[ERROR:] Could not open file %s\n", fileName);
	  fprintf(stderr, "         Please make sure the path is correct and you have read permission.\n");
	assert(referenceFile);
	}
      
  unsigned int verbose = arg->verbose;    
  
  unsigned int i;
      
  int ic;
  char c;
  
  unsigned int foundName = 0;
  unsigned int determinedLF = 0;
  unsigned int determinedCR = 0;  
  
  int name_length = 0;
  char * name;
  size_t bytesToRead_name = 10;
  name = (char*)malloc(bytesToRead_name);
  
  int seqLength;
//unsigned int bytesToRead = 10;
  size_t bytesToRead = 10;
  char * seq;
  seq = (char*)malloc(bytesToRead);
  
  seqLength = getline(&seq, &bytesToRead,referenceFile);
  assert(seqLength >= 0);
  
  unsigned int sumLengths = 0;
  seqLength = getline(&seq, &bytesToRead,referenceFile);
   
	while( seqLength > 0)
	{
	sumLengths += seqLength;
	seqLength = getline(&seq, &bytesToRead,referenceFile);
	}
   
   rewind(referenceFile);
   
     
        //The getline() command used below does not support linebreaks that are denoted only by \r
        // (instead of \n or \r\n). So make sure we have a readable format.
	while(!determinedLF)
	{
	assert(ic = fgetc(referenceFile));
	c = (char)ic;
		if(ic == EOF)
		{
		fprintf(stderr, "[ERROR:] Could not find line feed character \\n in %s\n", arg->referenceFile);
			if(determinedCR)
			{
			fprintf(stderr, "         However, CR \\r was found. At present this software relies on \\n though.\n");
			fprintf(stderr, "         You may try converting the format with tools such as mac2unix of the dos2unix package.\n");
			}
		assert(ic != EOF);
		}

		if(c == '\n' )
		determinedLF = 1;
		if(c == '\r')
		determinedCR = 1;  
	}
	
    rewind(referenceFile);
    
    
   //Reallocate sufficient space for concatinating reference sequence.
   // Space allocated here will in almost all cases be strictly more than needed.
   int original_length=0;
   if(globalVar->referenceSequence)
   {
   original_length=strlen(globalVar->referenceSequence);
   globalVar->referenceSequence = (char*)realloc(globalVar->referenceSequence, (original_length+sumLengths+1) * sizeof(char));
   }else{
    globalVar->referenceSequence = (char*)malloc((sumLengths+1) * sizeof(char));
   }
   
        //Check to see whether we find a sequence name starting with either > or @
        while(!foundName)
	{
	assert(ic = fgetc(referenceFile));
	c = (char)ic;
		if(ic == EOF)
		{
		fprintf(stderr, "[ERROR:] Could not find beginning of sequence in %s\n", arg->referenceFile);
		assert(ic != EOF);
		}

		 if(c == '>' || c == '@')
		foundName=1;
	}
  

  name_length = getline(&name, &bytesToRead_name,referenceFile);
  assert(name_length);
  chomp(name);
  cut_at_first_space(name);
  
//   print_selective("Name:             %s (%d)\n\n", name, name_length);
print_selective(" Reference name:  %s\n\n", name);
//    print_selective("Name after chomp: %s (%d)\n", name, name_length);
//   exit(0);
 
  //TODO THIS HERE IS THE PROBLEM!!!!!!!!!
	globalVar->reference_names[globalVar->number_of_references] = strdup(name);
//         globalVar->reference_names[globalVar->number_of_references] = (char*)malloc((strlen(name)+1) * sizeof(char));
//         memcpy(globalVar->reference_names[globalVar->number_of_references], name, strlen(name));
//         globalVar->reference_names[globalVar->number_of_references][strlen(name)] = '\0';
	globalVar->number_of_references++;	
	
	
  seqLength = getline(&seq, &bytesToRead,referenceFile);
  sumLengths = 0;
  while(seqLength>=0)
    {

        if(verbose > 1)
        {
        printf("%s",seq);
        }
        
	for(i=0;i<strlen(seq);i++)//Go through the whole sequence and simply write down the appropriate bases.
	{
	  if( mapUnsafe(seq[i]) < 0)
	  continue;
            {
            globalVar->referenceSequence[original_length + sumLengths] = seq[i];
	
            sumLengths++;
            }
	}
	  fflush(stdout);
	 seqLength = getline(&seq, &bytesToRead,referenceFile);
    }
   
   if(verbose > 1)
   printf("\nSumLengths = %u\n", sumLengths);
   
   globalVar->referenceSequence = (char*)realloc(globalVar->referenceSequence, (original_length+sumLengths+1)*sizeof(char));
   globalVar->referenceSequence[original_length+sumLengths] = '\0';
   globalVar->referenceSequenceLength = original_length+sumLengths;
   
   
   fclose(referenceFile);
  
   free(seq);

   free(name);
 }
 
 //-------------Read-Reference-File---------------------------------------------------------------
 //Will attempt to read arg->referenceFile and write its contents to globalVar->referenceSequence
 void readReferenceFile(setting *arg, globalVariables *globalVar)
 {

     FILE *referenceFile;
     referenceFile = fopen(arg->referenceFile,"r");
	
        //Check that file opened succesfully
	if(!referenceFile)
	{
	fprintf(stderr, "\n[ERROR:] Could not open file %s\n", arg->referenceFile);
	  fprintf(stderr, "         Please make sure the path is correct and you have read permission.\n");
	assert(referenceFile);
	}
      
  unsigned int verbose = arg->verbose;    
  
  unsigned int i;
      
  int ic;
  char c;
  
  unsigned int foundName = 0;
  unsigned int determinedLF = 0;
  unsigned int determinedCR = 0;  

  int name_length = 0;
  char * name;
  size_t bytesToRead_name = 10;
  name = (char*)malloc(bytesToRead_name);
  
  int seqLength;
//unsigned int bytesToRead = 10;
  size_t bytesToRead = 10;
  char * seq;
  seq = (char*)malloc(bytesToRead);
  
  seqLength = getline(&seq, &bytesToRead,referenceFile);
  assert(seqLength >= 0);
  
  unsigned int sumLengths = 0;
  seqLength = getline(&seq, &bytesToRead,referenceFile);
   
	while( seqLength > 0)
	{
	sumLengths += seqLength;
	seqLength = getline(&seq, &bytesToRead,referenceFile);
	}
   
   rewind(referenceFile);
   
          //The getline() command used later does not support linebreaks that are denoted only by \r
        // (instead of \n or \r\n). So make sure we have a readable format.
	while(!determinedLF)
	{
	assert(ic = fgetc(referenceFile));
	c = (char)ic;
		if(ic == EOF)
		{
		fprintf(stderr, "[ERROR:] Could not find line feed character \\n in %s\n", arg->referenceFile);
			if(determinedCR)
			{
			fprintf(stderr, "         However, CR \\r was found. At present this software relies on \\n though.\n");
			fprintf(stderr, "         You may try converting the format with tools such as mac2unix of the dos2unix package.\n");
			}
		assert(ic != EOF);
		}

		if(c == '\n' )
		determinedLF = 1;
		if(c == '\r')
		determinedCR = 1;  
	}
        
rewind(referenceFile);
   
   //Allocate initial space for reference sequence.
   // Space allocated here will in almost all cases be strictly more than needed.
   globalVar->referenceSequence = (char*)malloc((sumLengths+1) * sizeof(char));
   
        //Check to see whether we find a sequence name starting with either > or @
        while(!foundName)
	{
	assert(ic = fgetc(referenceFile));
	c = (char)ic;
		if(ic == EOF)
		{
		fprintf(stderr, "[ERROR:] Could not find beginning of sequence in %s\n", arg->referenceFile);
		assert(ic != EOF);
		}

		 if(c == '>' || c == '@')
		foundName=1;
	}
  
  name_length = getline(&name, &bytesToRead_name,referenceFile);
  assert(name_length);
  chomp(name);
  cut_at_first_space(name);
  print_selective(" Reference name:  %s\n", name);
// print_selective("Name:             %s (%d)\n", name, name_length);

//   print_selective("Name after chomp: %s (%d)\n", name, name_length);
//   exit(0);
  
	globalVar->reference_names[globalVar->number_of_references] = strdup(name);
	globalVar->number_of_references++;
        
  seqLength = getline(&seq, &bytesToRead,referenceFile);
  sumLengths = 0;
  while(seqLength>=0)
{

  if(verbose > 1)
  {
  printf("%s",seq);
  }
	for(i=0;i<strlen(seq);i++)//Go through the whole sequence and simply write down the appropriate bases.
	{
	  if( mapUnsafe(seq[i]) < 0)
	  continue;
	  {
	globalVar->referenceSequence[sumLengths] = seq[i];
	
	sumLengths++;
	  }
	}
	  fflush(stdout);
	 seqLength = getline(&seq, &bytesToRead,referenceFile);
}
   
//    print_selective("1\n");
   
   if(verbose > 1)
   printf("\nSumLengths = %u\n", sumLengths);
   
   globalVar->referenceSequence = (char*)realloc(globalVar->referenceSequence, (sumLengths+1)*sizeof(char));
   globalVar->referenceSequence[sumLengths] = '\0';
   globalVar->referenceSequenceLength = sumLengths;
//       print_selective("1\n");
   globalVar->break_points[globalVar->number_of_references] = sumLengths;//place the first break between reference sequences at the length of the first sequence.
   
//    print_selective("2\n");
   
   fclose(referenceFile);
  
   free(seq);

   free(name);
   
 }
    //------------------------------------------------------------------------------------
    

    
 void freeGlobalVariables(setting *arg, globalVariables *globalVar)
 {
  unsigned int i,j;
  
//   print_selective("1\n");
    free(globalVar->break_points);
//    print_selective("3\n");
    free(globalVar->referenceSequence); 
//     print_selective("4\n");
    free(globalVar->hashTable);
//     print_selective("5\n");
    freeTable(globalVar->entryTable, globalVar->itemsInTable);
//       print_selective("6\n");
    free(globalVar->entryTable);
//       print_selective("7a\n");
   
      for(i = 0; i < globalVar->number_of_references; i++)
          free(globalVar->reference_names[i]);

    free(globalVar->reference_names);

//   print_selective("7A\n");

//       return;
    
 // /* Only needed if sequences are really to be stored.
    if(arg->storeReads)
    {
  for(i=0; i<globalVar->numReads; i++)
  {
	Read *r = &(globalVar->read[i]);
	if(!r->matchesReference)
	{
	free(r->seq);
	}

  }
//*/
  free(globalVar->read);

//   print_selective("3\n");
    for(i=0; i<globalVar->numAltSequences; i++)
  {
      free(globalVar->altSequences[i].seq);
  }
  free(globalVar->altSequences);
  
  
//   print_selective("4\n");
    
    for(j =0 ; j < globalVar->referenceSequenceLength/50 ; j++)
    {
    char ** list;
    list = globalVar->brokenSequences[j];
    
        for(i=0; i<globalVar->numBrokenSequences[j]; i++)
        {
            free(list[i]);
        }
        
       free(list);
        
    }
 
//   print_selective("5\n");
    free(globalVar->numBrokenSequences);   
    free(globalVar->allocatedNumBrokenSequences);
    free(globalVar->brokenSequences);
   
//   print_selective("6\n");
    for(j =0 ; j < globalVar->referenceSequenceLength/50 ; j++)
    {
    char ** list;
    list = globalVar->brokenSequencesShifted[j];
    
        for(i=0; i<globalVar->numBrokenSequencesShifted[j]; i++)
        {
            free(list[i]);
        }
        
        free(list);
        
    }
    
//   print_selective("7\n");
   free(globalVar->numBrokenSequencesShifted);   
   free(globalVar->allocatedNumBrokenSequencesShifted);
   free(globalVar->brokenSequencesShifted);
   
    }//close if()arg->storeReads) 
//   print_selective("8\n");
 }
 
 
    
void initGlobalVariables(setting *arg, globalVariables *globalVar)
{
 int i;
 
 globalVar->number_of_references = 0;//So far no reference has been added
 
//  globalVar->reference_names = (char**)calloc(globalVar->number_of_references, sizeof(char*));
 globalVar->reference_names = (char**)calloc(arg->num_multi_references, sizeof(char*));
 
 globalVar->break_points = (int*)calloc(arg->num_multi_references + 1, sizeof(int));//+1 because we have an initial 0
 globalVar->break_points[0] = 0; //This is redundant and only placed here for explicity
 
 globalVar->referenceSequence = NULL;
 
 print_selective("Reading File [%d]: %s\n", 0, arg->multi_referenceFiles[0]);
 readReferenceFile(arg, globalVar);      
 
    for(i = 1; i < arg->num_multi_references; i++)
    {
        print_selective("Reading File [%d]: %s\n", i, arg->multi_referenceFiles[i]);
//         print_selective("(Next File: %s)\n", arg->multi_referenceFiles[i+1]);
//          print_selective("(problem File: %s)\n", arg->multi_referenceFiles[45]);
     read_multiple_ReferenceFiles(arg, globalVar, arg->multi_referenceFiles[i]);
    }
    
 assert(globalVar->number_of_references == arg->num_multi_references);
 
  //allocate initial hash table of POINTERS! for entries.
  globalVar->hashTable = (hashEntry**)calloc(1, sizeof(hashEntry*));
  
  globalVar->entryTable = (hashEntry*)calloc(globalVar->referenceSequenceLength, sizeof(hashEntry)); 
  
  globalVar->itemsInTable = 0;
  
  //allocate initial hash table of POINTERS! for entries of complement sequence.
  //globalVar->complementHashTable = (hashEntry**)calloc(1, sizeof(hashEntry*));
  
  //globalVar->complementEntryTable = (hashEntry*)calloc(globalVar->referenceSequenceLength, sizeof(hashEntry)); 
  
  //globalVar->itemsInComplementTable = 0;
  
//   time_t timer;
//   char buffer[10];
//   struct tm* tm_info;

//   time(&timer);
//   tm_info = localtime(&timer);

//   strftime(buffer, 10, "%H:%M:%S", tm_info);
    
//   print_selective("\n[%s] Setting up hashtable for reference sequence...\n", buffer);
  print_selective("\n");
  printTime();
  print_selective(" Setting up hashtable for reference sequence...\n");
  
  populateHashTable(arg, globalVar, globalVar->referenceSequence, &(globalVar->hashTable), &(globalVar->entryTable), &(globalVar->itemsInTable));
 //   printf("\nItems in table %d\n",globalVar->itemsInTable);
  
 //Actually set up reverse complemented sequence
  //char * reverseComplement;
  //reverseComplement = (char*)calloc(strlen(globalVar->referenceSequence) +1, sizeof(char));
  
  //reverseComplementSequence(reverseComplement, globalVar->referenceSequence);

  //reverseComplement[strlen(globalVar->referenceSequence)]='\0';//setting terminal character
  
  //printf("\nSetting up hashtable for reverse complemented reference sequence...");
  //populateHashTable(arg, globalVar, reverseComplement, &(globalVar->complementHashTable), &(globalVar->complementEntryTable), &(globalVar->itemsInComplementTable));
  
  //free(reverseComplement);
 

  // Set up variables to store individual reads
if(arg->storeReads){
    globalVar->numReads = 0;
	globalVar->allocatedNumReads = 1024;
	globalVar->read = (Read*)calloc(globalVar->allocatedNumReads ,sizeof(Read));
	
        	// Set up variables to store individual alternative sequences (proto haplotypes)
	globalVar->numAltSequences = 0;
	globalVar->allocatedNumAltSequences = 4;
	globalVar->altSequences = (Haplotype*)calloc(globalVar->allocatedNumAltSequences ,sizeof(Haplotype));
        
        //TODO 50 must not be hardcoded
        globalVar->numBrokenSequences = (unsigned int *)calloc(globalVar->referenceSequenceLength/50 + 1, sizeof(unsigned int) );
        globalVar->allocatedNumBrokenSequences = (unsigned int *)calloc(globalVar->referenceSequenceLength/50 + 1, sizeof(unsigned int) );
        globalVar->brokenSequences = (char***)calloc(globalVar->referenceSequenceLength/50 + 1, sizeof(char**) );
        
        globalVar->numBrokenSequencesShifted = (unsigned int *)calloc(globalVar->referenceSequenceLength/50 + 1, sizeof(unsigned int) );
        globalVar->allocatedNumBrokenSequencesShifted = (unsigned int *)calloc(globalVar->referenceSequenceLength/50 + 1, sizeof(unsigned int) );
        globalVar->brokenSequencesShifted = (char***)calloc(globalVar->referenceSequenceLength/50 + 1, sizeof(char**) );
//	for(i=0 ; i < globalVar->allocatedNumReads; i++)
 //       {
   //         globalVar->read[i].matchesReference = 0;
     //   }
}
   
 globalVar->trimmed = 0;
 globalVar->kept = 0;
 
 globalVar->avgRatioTrimmed = 0.0;
 globalVar->avgTrimmedOf = 0;
  
//   time(&timer);
//   tm_info = localtime(&timer);

//   strftime(buffer, 10, "%H:%M:%S", tm_info);
//   print_selective("\n[%s] Setup completed!\n\n", buffer);
 print_selective("\n");
 printTime();
  print_selective(" Setup completed!\n\n");
}



void freeResults(resultsVector *rv)
{
    unsigned int i;
    
    result* r;
//    printf("rv->assignedLength  in free %u\n",rv->assignedLength);
    for(i = 0; i<rv->assignedLength; i++)
    {
      r = &(rv->results[i]);
      free(r->indivError);
      //free( (&(rv->results[i]))->indivError);
    }
  
  free(rv->results);
}

void reallocResults(resultsVector *rv, unsigned int length)
{
  unsigned int i;
  
  rv->results = realloc(rv->results, length * sizeof(result));
	     
  result* r;
		for(i = rv->assignedLength; i < length; i++)
		{
		r = &(rv->results[i]);
		*r = (const result){ 0 };
		r->indivErrorLength = 41;
		r->indivError = (unsigned int*)calloc(r->indivErrorLength, sizeof(unsigned int));
		}
    
	     rv->assignedLength = length;
  
//  printf("rv->assignedLength  in realloc %u\n",rv->assignedLength);
  
}

void initResults(resultsVector *rv, unsigned int length)
{
  unsigned int i;
  
    rv->results =(result*)calloc(length, sizeof(result));
    rv->assignedLength = length;
    
    result *r;
    for(i = 0; i<length; i++)
    {
      r = &(rv->results[i]);
      r->indivErrorLength = 41;
      r->indivError = (unsigned int*)calloc(r->indivErrorLength, sizeof(unsigned int));
      
      r->numIndels = 0;
    }
    
    rv->indels = 0;
    rv->numReferenceMatches = 0;
    
   rv->hit  = 0;
   rv->miss = 0;
}





//-----------------------------------------------------------------------------------------------------

void cleanForExit(setting * s, globalVariables * g,  resultsVector * rv)
{
    
//   print_selective("freeing results\n");
  freeResults(rv);
  
  
//   print_selective("freeing global variables\n");
  freeGlobalVariables(s, g);
  
  
//   print_selective("freeing settings\n");
  freeSetting(s);
  
//   print_selective("end\n");
}
