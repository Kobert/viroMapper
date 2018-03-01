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
 
#include  "map_arg.h"
 
 void map_help()
 {
  printf("\n"); 
  printf("This is viroMapper v. 0.0.6\n\n"); 
  printf("Function arguments are:\n"); 
  printf("\n"); 
  printf(" -r <string>       for specifying a file containing the reference sequence (fasta format only)\n");    
  printf(" -s <string>       for specifying a file with shortreads (fastq format only)\n");    
  printf("\n");  
  printf("\n");   
  printf("Optional arguments:\n\n"); 
  printf(" -o <Path/string>  specify a prefix for the output names (can include a path).\n");  
  printf("                   If no Prefix is provided, only generic temporary files will be written.\n");
  printf("\n");    
  printf("                   Output formats:\n");    
  
  printf(" -C                Write fasta file containing the consensus sequence.\n");   
  printf("                   A path should be provided with -o for this option.\n");    
  printf(" -G                Plot results in Gnuplot.\n");
  printf("                   A path can be provided with -o for this option.\n");
  printf("                   If no path is provided, no *.png files will be written.\n");  
  printf(" -S                Print alignment information in SAM format to stdout.\n");
  printf(" -T                Only test reference sequence for suitability and exit.\n");   
  printf("                   (No fastq file is needed with this otion.)\n");   
  printf(" -V                Print results in form of a *.csv file.\n");
  printf("                   A path should be provided with -o for this option.\n");    
  printf("\n");  
  printf("\n");    
  printf(" -t                Trim any suffix and prefix that does not contain 8 consecutive matches.\n"); 
  printf("                   (This option is suggested for most analyses.)\n");
  printf("\n");    
  printf(" -h                Print this help messgae\n");    
  printf("\n\n");
  printf("Developers Options:\n");    
  printf(" -d        for more verbose debug output\n"); 
  printf(" -f <INT>  Specify how many fractions must match at least. (default 4, must be 2 or larger)\n");  
  printf(" -e        For extended testing (also consider right aligned sequences).\n");  
  printf(" -L        Old output format for csv header\n");
  printf("\n");
 }
 


 setting map_initArgs(int argc, char **argv)
 {
  setting s;
  
  s.verbose = 1;
  
  int index;
  int c;
  
  
  s.mapOnly = 1;
  s.primes  = 0;
  
  s.printSAM = 0;
   
  s.storeReads = 0;

#ifdef _multiUnsigned
  s.windows = _numMulti;
#else
  s.windows = 1;
#endif
  
  s.referenceFile = NULL;
  
  s.multi_referenceFiles = NULL;
  s.num_multi_references = 0;
  
  s.readsFile = NULL;

  s.outFilePrefix = NULL;
  
  s.doCsvFile = 0;
  s.doGnuplotFile = 0;
  s.doJSFile = 0;
  s.writeConsensus = 0;
  
  s.numPlots = 40;
  
  s.minFracs = 4;
  //s.pValue = -1.0;
  
  s.doTrim = 0;
  s.doExtendedTesting = 0;
  
  s.executeReferenceOnly = 0;
  
  opterr = 0;
  
  s.qFloor = -1;
  
//   Variable that denotes the usage of old/discontinued parameters/output formats/layouts
  s.LEGACY = 0;


  while ((c = getopt (argc, argv, "acCdef:gGhjJLmr:s:StTpVw:o:O:q:xyz")) != -1)
  {          
    switch (c)
      {
      case 'h':
	map_help();
	exit(1);
	break;
      case 'a':
          print_selective("\n[NOTE:] You are using the discontinued -%c option name. \n", c);
          print_selective("        To get rid of this warning, please use -%c for the same result\n\n", 'S');          
      case 'S':
        s.printSAM = 1;  
        break;
      case 'c':
          print_selective("\n[NOTE:] You are using the discontinued -%c option name. \n", c);
          print_selective("        To get rid of this warning, please use -%c for the same result\n\n", 'V');        
      case 'V':
	s.doCsvFile = 1;
	break;	
      case 'y':
          print_selective("\n[NOTE:] You are using the discontinued -%c option name. \n", c);
          print_selective("        To get rid of this warning, please use -%c for the same result\n\n", 'C');        
      case 'C':
          s.writeConsensus = 1;
        break;
      case 'e':
	s.doExtendedTesting = 1;
	break;		
      case 'f':
	assert(atoi(optarg) > 1);
	s.minFracs = atoi(optarg);
	break;		
      case 'g':
        print_selective("\n[NOTE:] You are using the discontinued -%c option name. \n", c);
        print_selective("        To get rid of this warning, please use -%c for the same result\n\n", 'G');        
      case 'G':
	s.doGnuplotFile = 1;
	break;
      case 'j':
        print_selective("\n[NOTE:] You are using the discontinued -%c option name. \n", c);
        print_selective("        To get rid of this warning, please use -%c for the same result\n\n", 'J');          
      case 'J':
        s.doJSFile = 1;
        break;
      case 'L':
        s.LEGACY = 1;
      case 'm':
	s.mapOnly = 1;
	s.storeReads = 0;
	break;
      case 'z'://Check behaviour of different prime numbers
	fprintf(stderr, "[ERROR:] Feature '-z' not live yet! Aborting analysis.");
	assert(0);
        s.primes = 1;
        break;      
//     case 'p':
//	s.pValue = atoi(optarg);//printf("%s\n", optarg);
//        break;
      case 'r':
        if(!s.referenceFile) //This is mostly needed for historic reasons. Everything should now work with multi_referenceFiles
	s.referenceFile = strdup(optarg);
	
	s.multi_referenceFiles = (char**)realloc(s.multi_referenceFiles,sizeof(char*)*(s.num_multi_references+1) );
        
        s.multi_referenceFiles[s.num_multi_references] = strdup(optarg);
	s.num_multi_references++;
        
            while(argv[optind] && (argv[optind])[0]!='-')
            {
            s.multi_referenceFiles = (char**)realloc(s.multi_referenceFiles,sizeof(char*)*(s.num_multi_references+1) );
            s.multi_referenceFiles[s.num_multi_references] = strdup(argv[optind]);
	
            s.num_multi_references++;
       
//             fprintf(stderr, "%s ",argv[optind]);
            optind++;
            }
            
//             for(int i = 0; i<s.num_multi_references; i++)
//             {
//                 fprintf(stderr, "%s \n", s.multi_referenceFiles[i]);
//             }
            
        break;
      case 'd':
        s.verbose = 2;
        break;
      case 's':
        s.readsFile = strdup(optarg);
        break;
      case 'x':
          print_selective("\n[NOTE:] You are using the discontinued -%c option name. \n", c);
          print_selective("        To get rid of this warning, please use -%c for the same result\n\n", 'T');                  
      case 'T':
        s.executeReferenceOnly = 1;
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
      case '?':
        if (optopt == 's' || optopt == 'r' || optopt == 'w' )
	{
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	  assert(0);
	}
	else if (isprint (optopt))
	{
          fprintf (stderr, "[ERROR:] Unknown option `-%c'.\n", optopt);
	
	  assert(0);
	}
	else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return s;
      default:
        map_help();
        abort ();
      }
//       fprintf(stderr,"TEST: optind %d: %s\n", optind, argv[optind]);
 }
      
	for (index = optind; index < argc; index++)
	print_selective ("Non-option argument %s\n", argv[index]);

	if(!s.referenceFile)
	{
        map_help();
	fprintf(stderr, "[ERROR:] Please specify a reference file with the \'-r\' command.\n\n");
	assert(s.referenceFile);
	}
	
	if(!s.readsFile && !s.primes && !s.executeReferenceOnly)
	{
         map_help();   
	 fprintf(stderr, "[ERROR:] Please specify a reads file with the -s command,\n");
         fprintf(stderr, "         or use an option that does not rely on it (e.g \'-x\')\n\n");
         
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
   print_selective("           A generict temporary file will be written.\n\n");
  }
  
  if(s.doJSFile && !s.outFilePrefix)
  {
   print_selective("\n[Warning:] JS output is set without specifying an output path (via -o).\n ");   
   print_selective("           A generict temporary file will be written.\n\n");
  }
  
      if(s.writeConsensus && !s.outFilePrefix)
  {
   print_selective("\n[Warning:] fasta consensus output is requested without specifying an output path (via -o).\n ");   
   print_selective("           A generict temporary file will be written.\n");
  }
  
    if(!s.doGnuplotFile && !s.doCsvFile && !s.executeReferenceOnly && !s.writeConsensus && !s.printSAM && !s.doJSFile)//TODO add any new output formats here
  {
   print_selective("\n[Note:] No output format specified.\n ");   
   print_selective("        Only basic statistics will be printed to the screen.\n ");
   print_selective("        To see a list of available options run the program with the -h flag.\n\n");
  }  
  
  if(s.num_multi_references > 1 && s.executeReferenceOnly)
  {
   print_selective("\n[Note:] Multiple references given for \'-x\'!\n");
   
   print_selective("        Only basic numbers will be calculated.\n");
   print_selective("        If you wabt to test for suitability via \"-x\", please \n");
   print_selective("        only submit a single reference sequence.\n\n");
 
  }
  
  return s;
 }
 
 
