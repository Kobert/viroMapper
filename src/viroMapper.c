
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>


#include "referenceAssembly.h" 
#include "ref_arg.h"			// Functions to interpret user input and initialize datastructures go here.
#include "ref_primes.h"			// Functions to test combinations of primenumbers for efficient hashing techniques
#include "ref_slidingWindow.h"		// Contains all functions that Manipulate or access the slidingWindow, or actual bases themselves.
#include "ref_hash.h"			// All functions that manipulate or access the tashtable or entries go here.
#include "ref_io.h"			// Functions that hanlde input, output, as well as result formatting (should) go here.
#include "ref_math.h"			// Anything more complicated than simple multiplication goes here

#include "map_arg.h"



int main(int argc, char **argv)
{
  
//  fprintf(stderr,"FIX chainHash!!!\n");
//    exit(0);
  
 setting arg = map_initArgs(argc, argv);
 
 globalVariables globalVar;
  
 initGlobalVariables(&arg, &globalVar);
 
 resultsVector rv;
 
 initResults(&rv, strlen(globalVar.referenceSequence));
 
  fflush(stdout);
  
  if(arg.executeReferenceOnly)
  {
     testReference(&globalVar, arg);
          
     print_selective("[NOTE:]    \"-x\" set. Exiting now...\n"); 
      
     cleanForExit(&arg, &globalVar,  &rv);
  
  return 0;   
  }
  
//   Placement placement;
//   //insertion of 1
// //   char* t_c = "ACGCGGAGUAGCCTCGUGUGAGCCCC";
// //      align_to_reference_2_positions(arg, &globalVar, t_c, 99, 100);
// //      align_to_reference_2_positions(arg, &globalVar, t_c, 100, 99);
//     
// //     Insertion of 3
//        char* t_c = "ACGCGGAGUAGCCTTTCGUGUGAGCCCC";
// //        align_to_reference_2_positions(arg, &globalVar, t_c, 100, 97);
//     placement = align_to_reference_2_positions(arg, &globalVar, t_c, 97, 100);
//     
//     
//     //Deletion of 1
// //   char* t_c = "ACGCGGAGUAGCCGUGUGAGCCCC";
// //      align_to_reference_2_positions(arg, &globalVar, t_c, 100, 101);
// //      align_to_reference_2_positions(arg, &globalVar, t_c, 101, 100);
//      
// //      Deletion of 3 not working
// //   char* t_c = "ACGCGGAGUAGGUGUGAGCCCC";
// //      align_to_reference_2_positions(arg, &globalVar, t_c, 103, 100);
// //     align_to_reference_2_positions(arg, &globalVar, t_c, 100, 103);
//     
// //   Placement p = align_to_reference_2_positions(arg, &globalVar, t_c, 99, 100);
// //   for(int temrun = 0; temrun<11; temrun++)
// //   {
// //   unsigned int temwind = temrun, temlength = 250;
// //   printf("rl: %u, l: %u, f: %u\n", globalVar.referenceSequenceLength, temlength, temwind);
// //   printf("%.20f\n", probRandom(arg, temwind, globalVar.referenceSequenceLength, temlength) );
// //   printf("%u\n", mapping_quality(arg, temwind, globalVar.referenceSequenceLength, temlength) );
// //  
// //   }
//     
// 
//          print_selective("Position: %d\n", placement.start_leftmost);
//         print_selective("CIGAR: %dM", placement.num_left_match);
//         if(placement.has_insertion == 1)
//             print_selective("%dI", placement.indel_length);
//         else
//         print_selective("%dD", placement.indel_length);
//             
//           print_selective("%dM\n", placement.num_right_match);
//     
//    exit(1);
// //   


  hashMapReads(arg, &globalVar, &rv);


  print_selective("\nGood hits: %d, multi hits %d bad hits %d (bad + multi %.11f)\n", good_hit, multi_hit, bad_hit, (bad_hit+multi_hit)/(double)(good_hit+bad_hit+multi_hit));
  
  postProcessResults(arg, &rv);

  printAvgCoverage(arg, rv);
  
 
  if(arg.doGnuplotFile)
  {
  char *gnuplotFileName;
if(arg.outFilePrefix)
{
gnuplotFileName = (char*)calloc(strlen(arg.outFilePrefix)+16, sizeof(char));
sprintf(gnuplotFileName, "%s.dataPoints.dat", arg.outFilePrefix);
}else{
 gnuplotFileName = (char*)calloc(26, sizeof(char));

 sprintf(gnuplotFileName, "_TEMPORARY.dataPoints.dat");
}

FILE * gnuplotFile = fopen(gnuplotFileName,"w");

printGnuplotDat(gnuplotFile, rv);

fclose(gnuplotFile);

print_selective("\nDatafile for Gnuplot written to: \"%s\"\n",gnuplotFileName);
if(!arg.outFilePrefix)
{
 print_selective("File \"%s\" will be overwritten on the next call to this program!\n", gnuplotFileName); 
}

plotQualityProfile(arg, rv.assignedLength, gnuplotFileName);
print_selective("Qualityprofile plotted...\n");

plotCoverage(arg, gnuplotFileName);
print_selective("Coverageprofile plotted...\n");

plotMajorBase(arg, gnuplotFileName);
print_selective("Majorbases plotted...\n");

free(gnuplotFileName);
}
 

 
 
 
 if(arg.doCsvFile)
{
    char *csvFileName;
  if(arg.outFilePrefix)
{
  csvFileName = (char*)calloc(strlen(arg.outFilePrefix)+5, sizeof(char));
sprintf(csvFileName, "%s.csv", arg.outFilePrefix);
}else{
  csvFileName = (char*)calloc(10, sizeof(char));
sprintf(csvFileName, "_temp.csv");
}
FILE * csvFile = fopen(csvFileName,"w");

printCSV(csvFile, rv);
fclose(csvFile);
print_selective("\ncsv output written to:           \"%s\"\n", csvFileName);
free(csvFileName);
}



char *htmlFileName;
    if(arg.outFilePrefix)
    {
    htmlFileName = (char*)calloc(strlen(arg.outFilePrefix)+6, sizeof(char));
    sprintf(htmlFileName, "%s.html", arg.outFilePrefix);
    }else{
    htmlFileName = (char*)calloc(11, sizeof(char));
    sprintf(htmlFileName, "_temp.html");
    }
FILE * htmlFile = fopen(htmlFileName,"w");

printHtml(htmlFile, arg, rv);
fclose(htmlFile);
print_selective("\nHtml output written to:           \"%s\"\n", htmlFileName);
free(htmlFileName);





if(arg.writeConsensus)
{
    char *consensFileName;
  if(arg.outFilePrefix)
{
  consensFileName = (char*)calloc(strlen(arg.outFilePrefix)+7, sizeof(char));
sprintf(consensFileName, "%s.fasta", arg.outFilePrefix);
}else{
  consensFileName = (char*)calloc(12, sizeof(char));
sprintf(consensFileName, "_temp.fasta");
}
FILE * consensFile = fopen(consensFileName,"w");

writeConsensusReferenceFile(consensFile, rv, "consensus");
fclose(consensFile);
print_selective("\nConsensus sequence written to:       \"%s\"\n", consensFileName);
free(consensFileName);
}

if(arg.doTrim)
{
print_selective("\n\t Avg trimmed ratio:   %.4f (Average per read)\n",globalVar.avgRatioTrimmed);

print_selective("\t Total trimmed ratio: %.4f (Fraction of sites)\n\n",(double)globalVar.trimmed/(globalVar.trimmed + globalVar.kept));
} 
  

  printTime();
  print_selective(" Freeing data-structures...\n");
  cleanForExit(&arg, &globalVar,  &rv);
  
  printTime();
  print_selective(" Exiting.\n");
  
  return 0;

}



