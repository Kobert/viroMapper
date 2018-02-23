#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h>
#include <time.h>

#include "referenceAssembly.h"
#include "ref_math.h"


//---------UTILITY-FUNCTIONS---------------------------------------------------------------------------------------------

//The only reason for this function is to quickly change the output behaviour code-wide
void print_selective(const char *format, ...)
{
    va_list args;


    va_start(args,format);
    vfprintf(stderr, format, args);
    va_end(args);   

    
}

void printTime()
{
    time_t timer;
    char buffer[10];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(buffer, 10, "%H:%M:%S", tm_info);
       
    print_selective("[%s]", buffer);     
}

void chomp(char* c)
{
    int i;
    for(i = strlen(c)-1; i >= 0; i--)
    {
        if(c[i] == '\n' || c[i] == '\n')
        {
            c[i] = '\0';
        }else{
         break;   
        }
    }
}

void cut_at_first_space(char* c)
{
    
int i;
int found_non_space = 0;

    for(i=0; i<strlen(c); i++)
    {
        if(c[i] != ' ')
        found_non_space = 1;
        
        if(found_non_space && c[i] == ' ')
        {
        c[i] = '\0';
        return;   
        }
    }
    
}

//return of this function must be freed!
// This function simply chacks the string from the right and
// stops as soon as the first "/" is encountered.
// Everything left of this is assumed to be the file path.
char* remove_path_by_strdup(char* c)
{
    
int i;
char* result;

    for(i=strlen(c); i>0; i--)
    {
        if(c[i] == '/')
        {
        i++;
        break;
        }
    }

    assert(i<strlen(c) && "String seems to end wuth \"/\" (or other problem present).");
result = strdup(&c[i]);
    
return result;
}

//---------------------------------------------------------
void readToResult(resultsVector *rv, char * seq, unsigned int pos, unsigned int length)
{
  
    unsigned int i;
  
    for(i = 0; i < length; i++)
    {
    result* r = &(rv->results[pos + i]);
   
    char c = seq[i];
    

       switch (c)
      {
      case 'A':
      case 'a':
	   r->coverage ++;
	r->A++;
	break;
      case 'C':
      case 'c':
	   r->coverage ++;
	r->C++;
	break;
      case 'G':
      case 'g':
	   r->coverage ++;
	r->G++;
	break;
      case 'T':
      case 't':
	   r->coverage ++;
	r->T++;
	break;
      case 'N':
      case 'X':
      case '-':
      case 'W':
      case 'w':
      case 'S':
      case 's':
      case 'M':
      case 'm':
      case 'K':
      case 'k':
      case 'R':
      case 'r':
      case 'Y':
      case 'y':
      case 'B':
      case 'b':
      case 'd':
      case 'D':
      case 'H':
      case 'h':
      case 'V':
      case 'v':
	   r->coverage ++;
	r->N++;
	break;
      default:
	fprintf(stderr, "seg ='%s'",seq);
	assert(0);//A sequence should not contain anything but DNA characters... 
		  // TODO if this is a limiting factor, simply skip the entry but keep i small
		  // Also, make sure the same is true for qualities...
        break;
      }
  }
  
  
  
}

void revCompReadToResult(resultsVector *rv, char * seq, unsigned int pos, unsigned int length)
{
  
  unsigned int i;
  
  for(i = 0; i < length; i++)
  {
   result* r = &(rv->results[pos + i]);
   

    
    char c = seq[length -1 -i];
       switch (c)
      {
      case 'T':
      case 't':
	   r->coverage ++;
	r->A++;
	break;
      case 'G':
      case 'g':
	   r->coverage ++;
	r->C++;
	break;
      case 'C':
      case 'c':
	   r->coverage ++;
	r->G++;
	break;
      case 'A':
      case 'a':
	   r->coverage ++;
	r->T++;
	break;
      case 'N':
      case 'X':
      case '-':
      case 'W':
      case 'w':
      case 'S':
      case 's':
      case 'M':
      case 'm':
      case 'K':
      case 'k':
      case 'R':
      case 'r':
      case 'Y':
      case 'y':
      case 'B':
      case 'b':
      case 'd':
      case 'D':
      case 'H':
      case 'h':
      case 'V':
      case 'v':
        r->coverage ++;
	r->N++;
	break;
      default:
	fprintf(stderr, "seg ='%s'",seq);
	assert(0);//A sequence should not contain anything but DNA characters... 
		  // TODO if this is a limiting factor, simply skip the entry but keep i small
		  // Also, make sure the same is true for qualities...
        break;
      }
  }
  
  
  
}


void readToResult_quality_aware(resultsVector *rv, char* seq, char* q, unsigned int pos, unsigned int length)
{
  
    unsigned int i;
//     unsigned int value;
      double match, miss;
    
    for(i = 0; i < length; i++)
    {
        
//       if((int)q[i] < 33)
//         {
//         fprintf(stderr, "%s\n",q);
//         fprintf(stderr, "q[%u]='%c'\n",i, q[i]);
//         assert((int)q[i] >= 33);
//         }
//       value = (int)q[i] - 33;
      
 //    If coverage is distributed over all competing sites
//    match = (1-cQ2P(q[i]));
//    miss = 1/3.0*(1-match);
      
//    If coverage is distributed over called base and N   
      match = (1-cQ2P(q[i]));
//    This is done to retain mathematical equivalence between both mothods.
//       I.e. if the portion assigned to N is redistributed to all 4 bases equally, both methods are the same.
      match = match - 1/3.0*(1-match);
      if(match < 0.0){
       match = 0.0;   
      }
      miss = (1-match);
      
      
      
      
    result* r = &(rv->results[pos + i]);
   
    char c = seq[i];
    

       switch (c)
      {
      case 'A':
      case 'a':
	   r->coverage ++;
	r->A++;
        r->qA = r->qA + match;
//         r->qC = r->qC + miss;
//         r->qG = r->qG + miss;
//         r->qT = r->qT + miss;
	r->qN = r->qN + miss;
	break;
      case 'C':
      case 'c':
	   r->coverage ++;
	r->C++;
// 	r->qA = r->qA + miss;
        r->qC = r->qC + match;
//         r->qG = r->qG + miss;
//         r->qT = r->qT + miss;
	r->qN = r->qN + miss;
	break;
      case 'G':
      case 'g':
	   r->coverage ++;
	r->G++;
// 	r->qA = r->qA + miss;
//         r->qC = r->qC + miss;
        r->qG = r->qG + match;
//         r->qT = r->qT + miss;
	r->qN = r->qN + miss;
	break;
      case 'T':
      case 't':
	   r->coverage ++;
	r->T++;
// 	r->qA = r->qA + miss;
//         r->qC = r->qC + miss;
//         r->qG = r->qG + miss;
        r->qT = r->qT + match;
	r->qN = r->qN + miss;
	break;
      case 'N':
      case 'X':
      case '-':
      case 'W':
      case 'w':
      case 'S':
      case 's':
      case 'M':
      case 'm':
      case 'K':
      case 'k':
      case 'R':
      case 'r':
      case 'Y':
      case 'y':
      case 'B':
      case 'b':
      case 'd':
      case 'D':
      case 'H':
      case 'h':
      case 'V':
      case 'v':
	   r->coverage ++;
	r->N++;
        r->qN = r->qN + 1.0;
	break;
      default:
	fprintf(stderr, "seg ='%s'",seq);
	assert(0);//A sequence should not contain anything but DNA characters... 
		  // TODO if this is a limiting factor, simply skip the entry but keep i small
		  // Also, make sure the same is true for qualities...
        break;
      }
  }
  
  
  
}

void revCompReadToResult_quality_aware(resultsVector *rv, char * seq, char* q, unsigned int pos, unsigned int length)
{
  
  unsigned int i;
//   unsigned int value;
  double match, miss;
  
  for(i = 0; i < length; i++)
  {
  
//       if((int)q[i] < 33)
//         {
//         fprintf(stderr, "%s\n",q);
//         fprintf(stderr, "q[%u]='%c'\n",i, q[i]);
//         assert((int)q[i] >= 33);
//         }
//       value = (int)q[i] - 33;
      
      
//    If coverage is distributed over all competing sites
//    match = (1-cQ2P(q[i]));
//    miss = 1/3.0*(1-match);
      
//    If coverage is distributed over called base and N   
      match = (1-cQ2P(q[i]));
//    This is done to retain mathematical equivalence between both mothods.
//       I.e. if the portion assigned to N is redistributed to all 4 bases equally, both methods are the same.
      match = match - 1/3.0*(1-match);
      if(match < 0.0){
       match = 0.0;   
      }
      miss = (1-match);
      
      
  
   result* r = &(rv->results[pos + i]);
   

    
    char c = seq[length -1 -i];
       switch (c)
      {
      case 'T':
      case 't':
	   r->coverage ++;
	r->A++;
        r->qA = r->qA + match;
        r->qN = r->qN + miss;
//         r->qC = r->qC + miss;
//         r->qG = r->qG + miss;
//         r->qT = r->qT + miss;
	break;
      case 'G':
      case 'g':
	   r->coverage ++;
	r->C++;
        
//         r->qA = r->qA + miss;
        r->qC = r->qC + match;
//         r->qG = r->qG + miss;
//         r->qT = r->qT + miss;
        r->qN = r->qN + miss;
        break;
      case 'C':
      case 'c':
	   r->coverage ++;
	r->G++;
        
//         r->qA = r->qA + miss;
//         r->qC = r->qC + miss;
        r->qG = r->qG + match;
//         r->qT = r->qT + miss;
        r->qN = r->qN + miss;
	break;
      case 'A':
      case 'a':
	   r->coverage ++;
	r->T++;
        
//         r->qA = r->qA + miss;
//         r->qC = r->qC + miss;
//         r->qG = r->qG + miss;
        r->qT = r->qT + match;
        r->qN = r->qN + miss;        
	break;
      case 'N':
      case 'X':
      case '-':
      case 'W':
      case 'w':
      case 'S':
      case 's':
      case 'M':
      case 'm':
      case 'K':
      case 'k':
      case 'R':
      case 'r':
      case 'Y':
      case 'y':
      case 'B':
      case 'b':
      case 'd':
      case 'D':
      case 'H':
      case 'h':
      case 'V':
      case 'v':
        r->coverage ++;
	r->N++;
        r->qN = r->qN + 1.0;
	break;
      default:
	fprintf(stderr, "seg ='%s'",seq);
	assert(0);//A sequence should not contain anything but DNA characters... 
		  // TODO if this is a limiting factor, simply skip the entry but keep i small
		  // Also, make sure the same is true for qualities...
        break;
      }
  }
  
  
  
}

void qualityToResult(resultsVector *rv, char * seq, unsigned int pos, unsigned int length)
{
  
  unsigned int i,j;
  
  unsigned int value;
  
  result *r;
  for(i = 0; i < length; i++)
  {
    if((int)seq[i] < 33)
    {
      fprintf(stderr, "%s\n",seq);
      fprintf(stderr, "seq[%u]='%c'\n",i, seq[i]);
    assert((int)seq[i] >= 33);
    }
  value = (int)seq[i] - 33;
  
  r=&(rv->results[pos + i]);
  

  if(value + 1 > r->indivErrorLength)
  {
    r->indivError = realloc(r->indivError, (value + 1) * sizeof(unsigned int));
    for(j = r->indivErrorLength; j < (value + 1); j++)
    {
      r->indivError[j] = 0;
    }
    r->indivErrorLength = (value +1);
  }
  r->indivError[value]++;
  
  }
  
  
  
}

void revQualityToResult(resultsVector *rv, char * seq, unsigned int pos, unsigned int length)
{
  
  unsigned int i,j;
  
  unsigned int value;
  
  result *r;
  for(i = 0; i < length; i++)
  {
    if((int)seq[length -1 -i] < 33)
    {
      fprintf(stderr, "%s\n",seq);
      fprintf(stderr, "seq[%u]='%c'\n",i, seq[length -1 -i]);
    assert((int)seq[length -1 -i] >= 33);
    }
  value = (int)seq[length -1 -i] - 33;
  
  r=&(rv->results[pos + i]);
  

  if(value + 1 > r->indivErrorLength)
  {
    r->indivError = realloc(r->indivError, (value + 1) * sizeof(unsigned int));
    for(j = r->indivErrorLength; j < (value + 1); j++)
    {
      r->indivError[j] = 0;
    }
    r->indivErrorLength = (value +1);
  }
  r->indivError[value]++;
  
  }
  
  
  
}

void postProcessResults(setting arg, resultsVector *rv)
{
  unsigned int i, j, k, sum, count;
  int median, uQuartile, lQuartile;
  int lower5, upper95;
  double fmedian, fsum;
  int max = -1;
  int min = -1;
  unsigned int second;
  char majorBase;
  char secondBase;
  result *r;
  unsigned int *e;
  
  double sumPerSiteCoverage       = 0;
  double sumQFloorPerSiteCoverage = 0;

	for(i = 0; i < rv->assignedLength; i++)
	{
	  r = &(rv->results[i]);
	  
	  assert( (r->A + r->C + r->G + r->T + r->N) == r->coverage);
  //Find major base and second basefrequency
	  
	  
	  //A
	  max = r->A;
	  majorBase = 'A';
	  
	//C
	 if(max < r->C) 
	 {
	 second = max;
	 max = r->C;
	 majorBase = 'C';
	 secondBase = 'A';
	 }else{
	  second = r->C; 
	 secondBase = 'C';
	 }
	 
	 
	//G
	 if(max < r->G) 
	 {
	 second = max;
	 max = r->G;
         secondBase = majorBase;
	 majorBase = 'G';
	 }else{
	   if(r->G > second)
           {
	  second = r->G; 
          secondBase = 'G';
           }
	 }
	   
	  
	//T
	 if(max < r->T) 
	 {
	 second = max;
         secondBase = majorBase;
	 max = r->T;
	 majorBase = 'T';
	 }else{
	   if(r->T > second)
           {
	  second = r->T;
         secondBase = 'T';
           }
	 }	  
	  
	//N
// 	 if(max < r->N) 
// 	 {
// 	 second = max;
//          secondBase = majorBase;
// 	 max = r->N;
// 	 majorBase = 'N';
// 	 }else{
// 	   if(r->N > second)
//            {
// 	  second = r->N; 
//          secondBase = 'N';
//            }
// 	 }
	 
	 if(max == 0)
         {
         secondBase = 'N';
         majorBase = 'N';    
         }
         if(second == 0)
         {
         secondBase = 'N';    
         }
	 
	 r->majorBase = majorBase;
         r->secondBase_char = secondBase;
	 r->majorBaseCount = max;
	 r->secondBaseCount = second;
	 
	 r->frequencyMajorBase  = (double)max/r->coverage;
	 r->frequencySecondBase = (double)second/r->coverage;
	 
// Find Quality scores
	 min  = -1;
	 max  = -1;
	 sum  =  0;
	 fsum =  0;
	 count=  0;
	 median  = -1;
	 lQuartile  = -1;
	 uQuartile  = -1;
	 lower5 = -1;
	 upper95= -1;
	 
	 fmedian = -1;
	 
	 r->qFloorCoverage = 0;

//	   for(j=0; j < r->indivErrorLength; j++)
//	   {
//	     e = r->indivError;
//	    printf("%u ", e[j]); 
//	   }
//	   printf("\n");
	   for(j=0; j < r->indivErrorLength; j++)
	   {
	     e = r->indivError;//NOTE this is the number of times a score of j has been observed!
	     sum   += j * e[j];
	     fsum  += Q2P(j)  * e[j];
	     count += e[j];
//coverage above the quality threhold	     
	     if(arg.qFloor > 0 && arg.qFloor == j+1)
	     {
	       r->qFloorCoverage = r->coverage - count;
	     }
//Median	     
	     if(count >= (double)r->coverage/2 && median < 0)
	     {
	       if(count > (double)r->coverage/2)
	       {
		 median = j;
		fmedian = Q2P(j);
	       }else{
		 
		 for(k=j+1; k< r->indivErrorLength ;k++)
		 {
		   if(e[k] > 0)
		     break;
 
		 }
		median = floor( (j+k)/2 + 0.5); 
		fmedian= floor( (Q2P(j) + Q2P(k))/2     + 0.5);
	       }
	       
	     }	     
//25% lower quartile	     
	     	     if(count >= (double)r->coverage/4 && lQuartile < 0)
	     {
	       if(count > (double)r->coverage/4)
	       {
		 lQuartile = j;
		}else{
		 
		 for(k=j+1; k< r->indivErrorLength ;k++)
		 {
		   if(e[k] > 0)
		     break;
		 }
		lQuartile = floor( (j+k)/2 + 0.5); 
		}
	       
	     }	
	     
//75% upper quartile	     
	     	     if(count >= 3*(double)r->coverage/4 && uQuartile < 0)
	     {
	       if(count > 3*(double)r->coverage/4)
	       {
		 uQuartile = j;
		}else{
		 
		 for(k=j+1; k< r->indivErrorLength ;k++)
		 {
		   if(e[k] > 0)
		     break;
		 }
		uQuartile = floor( (j+k)/2 + 0.5); 
		}
	       
	     }	
	     
//95% upper 	     
	     	     if(count >= 0.95*r->coverage && upper95 < 0)
	     {
	       if(count > 0.95*r->coverage)
	       {
		 upper95 = j;
		}else{
		 
		 for(k=j+1; k< r->indivErrorLength ;k++)
		 {
		   if(e[k] > 0)
		     break;
		 }
		upper95 = floor( (j+k)/2 + 0.5); 
		}
	       
	     }	
	     
//5% lower 	     
	     	     if(count >= 0.05*r->coverage && lower5 < 0)
	     {
	       if(count > 0.05*r->coverage)
	       {
		 lower5 = j;
		}else{
		 
		 for(k=j+1; k< r->indivErrorLength ;k++)
		 {
		   if(e[k] > 0)
		     break;
		 }
		lower5 = floor( (j+k)/2 + 0.5); 
		}
	       
	     }		     
	     
	     
	     if(e[j] > 0)
	     {
	       max = j;
	     }
	     if(e[j] > 0 && min < 0)
	     {
	       min = j;
	     }
	     
 
	   }
	 
	// assert(sumPerSiteCoverage < ULLONG_MAX - r->coverage);	// If this assertion fail (which it may),
								// we have to code around it by changing the mean on the fly
	 sumPerSiteCoverage += r->coverage;
	 sumQFloorPerSiteCoverage += r->qFloorCoverage;
	 
	 r->meanPhError     = floor((double)sum/r->coverage + 0.5);
	 r->medianPhError   = median;
	 r->lowerPhQuartile = lQuartile;
         r->upperPhQuartile = uQuartile;
	// r->minPhError      = min;
	// r->maxPhError      = max;
	 r->minPhError      = lower5;
	 r->maxPhError      = upper95;
	 
	 r->meanError = fsum/r->coverage;
	 r->medianError = fmedian;
	 r->minError = Q2P(min);
	 r->maxError = Q2P(max);
	}
  
    rv->averageCoverage       = sumPerSiteCoverage/(double)rv->assignedLength;
    rv->averageQFloorCoverage = sumQFloorPerSiteCoverage/(double)rv->assignedLength;
}

// NOTE Old and unused?
void printCSV(FILE *file, resultsVector rv)
{
    unsigned int i;
    result r;

//      fprintf(file, ", As, Cs, Gs, Ns, Ts, coverage, expected_number_of_errors, majorbase_ratio, majorbases, majorsequence, position, probability_of_seq_error, secondbase, secondbase_ratio, indels\n");
     fprintf(file, "position, As, Cs, Gs, Ns, Ts, coverage, expected_number_of_errors, majorbase_ratio, majorbases, majorsequence, probability_of_seq_error, secondbase, secondbase_ratio, indels\n");
     
	for(i = 0; i < rv.assignedLength; i++)
	{
	  r = rv.results[i];
	  
	  
	fprintf(file, "%u, %u, %u, %u, %u, %u, %u, %f, %f, %u, %c, %f, %u, %f, %c, %u\n", 
	       i, r.A, r.C, r.G, r.N, r.T, r.coverage, r.coverage * r.meanError, r.frequencyMajorBase, r.majorBaseCount, r.majorBase, r.meanError, r.secondBaseCount, r.frequencySecondBase, r.secondBase_char, r.numIndels );
	}
  
  
}

void printCSV_quality_aware_bases(FILE *file, resultsVector rv)
{
    unsigned int i;
    result r;

     fprintf(file, "position, As, Cs, Gs, Ns, Ts, coverage, expected_number_of_errors, majorbase_ratio, majorbases, majorsequence, probability_of_seq_error, secondbase, secondbase_ratio, indels, qAs, qCs, qGs, qNs, qTs\n");
     
	for(i = 0; i < rv.assignedLength; i++)
	{
	  r = rv.results[i];
	  
	  
	fprintf(file, "%u, %u, %u, %u, %u, %u, %u, %f, %f, %u, %c, %f, %u, %f, %c, %u, %f, %f, %f, %f, %f\n", 
	       i, r.A, r.C, r.G, r.N, r.T, r.coverage, r.coverage * r.meanError, r.frequencyMajorBase, r.majorBaseCount, r.majorBase, r.meanError, r.secondBaseCount, r.frequencySecondBase, r.secondBase_char, r.numIndels , r.qA, r.qC, r.qG, r.qN, r.qT);
	}
  
  
}


//print as html table element
void pd_c(FILE *file, char* color, const char *format, ...)
{
    va_list args;
    va_start(args,format);
    fprintf(file,"<td>");
    vfprintf(file, format, args);
    fprintf(file,"</td>");
    va_end(args);   
}


void pd(FILE *file, const char *format, ...)
{
  
    va_list args;
    va_start(args,format);
    fprintf(file,"<td>");
//     pd_c(file, NULL, format);
    vfprintf(file, format, args);
    fprintf(file,"</td>");
    va_end(args);   
//     
}

void ph(FILE *file, const char *name, const char *tooltip)
{
  
//     va_list args;
//     va_start(args,format);
    if(tooltip)
    fprintf(file,"<th title=\"%s\">", tooltip);
    else
    fprintf(file,"<th>");
    
    fprintf(file,"%s", name);
//     vfprintf(file, format, args);
    fprintf(file,"</th>");
//     va_end(args);   
//     
}

void printHtml(FILE *file, setting s, resultsVector rv)
{
    unsigned int i;
    result r;
    
    fprintf(file,"<html>\n");
    
    if(s.outFilePrefix)
    {
        char *name;
        name = remove_path_by_strdup(s.outFilePrefix);
        fprintf(file,"<head>\n");
        fprintf(file,"<title>Summary: %s</title>\n", name);
        fprintf(file,"</head>\n");
    }else{
        fprintf(file,"<head>\n");
        fprintf(file,"<title>Summary: viroMapper</title>\n");
        fprintf(file,"</head>\n");
       
    }
    
    fprintf(file,"<body>\n");
    fprintf(file,"<table border=\"1\">\n");
    
    fprintf(file,"<tr>\n");
    //ph prints <th>input1</th> or <th title="input2">input1</th> if input2!=NULL
    ph(file,"Position", "Nucleotide position on the reference genome");
    ph(file,"As", "Number of Adenine bases observed at this position");
    ph(file,"Cs", "Number of Cytosine bases observed at this position");
    ph(file,"Gs", "Number of Guanine bases observed at this position");
    ph(file,"Ts", "Number of Thymin or Uracil bases observed at this position");
    ph(file,"Ns", "Number of undetermined beses at this position");
    ph(file,"coverage", "Total coverage at this site");
//     ph(file,"majorbases", NULL);
    ph(file,"majorsequence", "Most prevalent base at this position");
    ph(file,"majorbase_ratio", "Ratio of occurences of the most prevalent base to the total coverage");
//     ph(file,"position", NULL);
    ph(file,"secondbase", NULL);
    ph(file,"secondbase_ratio", NULL);
    ph(file,"probability_of_seq_error", "Average probability of a sequencing error");
    ph(file,"expected_number_of_errors", "Expected number of sequencing errors at this site.");
//     ph(file,"indels", NULL);
    
    
//     fprintf(file, "<th></th><th>As</th><th>Cs</th><th>Gs</th><th>Ts</th><th>Ns</th><th>coverage</th><th>expected_number_of_errors</th><th>majorbase_ratio</th><th>majorbases</th><th>majorsequence</th><th>position</th><th>probability_of_seq_error</th><th>secondbase</th><th>secondbase_ratio</th><th>indels</th>\n");
    fprintf(file,"</tr>\n");
    
    
	for(i = 0; i < rv.assignedLength; i++)
	{
            
	 r = rv.results[i];
	      
                fprintf(file,"<tr>\n");
                pd(file, "<a id=\"l%d\" name=\"l%d\">%u</a>", i, i, i);
                pd(file, "%u", r.A);
                pd(file, "%u", r.C);
                pd(file, "%u", r.G);
                pd(file, "%u", r.T);
                pd(file, "%u", r.N);
                pd(file, "%u", r.coverage);
                pd(file, "%c", r.majorBase);
                pd(file, "%f", r.frequencyMajorBase);
//                 pd(file, "%u", r.majorBaseCount);
//              pd(file, "%u", i);
                pd(file, "%u", r.secondBaseCount);
                pd(file, "%f", r.frequencySecondBase);
                pd(file, "%f", r.meanError);
                pd(file, "%f", r.coverage * r.meanError);
//                 pd(file, "%u", r.numIndels);
              
	      fprintf(file,"\n</tr>\n");
              
	}
  
  fprintf(file,"</table></body></html>\n");
}



void print_interactive_html_file_js(FILE *file, setting s, resultsVector rv)
{
     unsigned int i, length;
     
     result r;
     
     length = rv.assignedLength;
     
     if(s.outFilePrefix)
        {
        char* temp = remove_path_by_strdup(s.outFilePrefix);
        fprintf(file, "var analysis_name = \"%s\";\n", temp);
        free(temp);
        }
     
      fprintf(file, "var max_line = %u;\n", length);
//      -------------------------------
     fprintf(file, "var As = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.A);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");  

//      -------------------------------
     fprintf(file, "var Cs = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.C);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");  
     
//      -------------------------------
     fprintf(file, "var Gs = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.G);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");       
     
     fprintf(file, "var Ns = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.N);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
//     --------------------------------------
          fprintf(file, "var Ts = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.T);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");      
     
     
//     --------------------------------------
          fprintf(file, "var coverage = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.coverage);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");      
     
     
     
     //     --------------------------------------
          fprintf(file, "var expected_number_of_errors = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        if( isnan(r.coverage * r.meanError) )
        {
        fprintf(file, "undefined");    
        }else{
        fprintf(file, "%f", r.coverage * r.meanError);
        }
        
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     

          //     --------------------------------------
          fprintf(file, "var majorbase_ratio = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        if( isnan(r.frequencyMajorBase) )
        {
        fprintf(file, "undefined");    
        }else{
        fprintf(file, "%f", r.frequencyMajorBase);
        }
        
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
     
     
               //     --------------------------------------
          fprintf(file, "var majorsequence = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          

        fprintf(file, "\"%c\"", r.majorBase);
        
        
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
     
          
     
               //     --------------------------------------
          fprintf(file, "var secondsequence = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          

        fprintf(file, "\"%c\"", r.secondBase_char);
        
        
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
     
          //     --------------------------------------
          fprintf(file, "var probability_of_seq_error = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        if( isnan(r.frequencyMajorBase) )
        {
        fprintf(file, "undefined");    
        }else{
        fprintf(file, "%f", r.meanError);
        }
        
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
     
     
               //     --------------------------------------
          fprintf(file, "var secondbase_ratio = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        if( isnan(r.frequencyMajorBase) )
        {
        fprintf(file, "undefined");    
        }else{
        fprintf(file, "%f", r.frequencySecondBase);
        }
        
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
     
     
 //     --------------------------------------
          fprintf(file, "var secondbase = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.secondBaseCount);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
              
//     --------------------------------------
          fprintf(file, "var min = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.minPhError);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
          
//     --------------------------------------
          fprintf(file, "var Q1 = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.lowerPhQuartile);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
          
//     --------------------------------------
          fprintf(file, "var Median = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.medianPhError);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
     
     //     --------------------------------------
          fprintf(file, "var Q3 = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.upperPhQuartile);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
     //     --------------------------------------
          fprintf(file, "var max = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.maxPhError);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
     //     --------------------------------------
          fprintf(file, "var Mean = [");  
     for(i = 0; i < length; i++)
	{            
        r = rv.results[i];
          
        fprintf(file, "%u", r.meanPhError);
     
            if(i+1 < length)
            {
            fprintf(file, ", ");   
                if( ((i+1) % 100) == 0)
                {
                fprintf(file, "\n");   
                }
            }
        }
     fprintf(file, "];\n");   
     
//    ??  var majorbases = [

     
//    ??  var position = [
//      i
     
     
     
     
     
    	fflush(file);
}




void printGnuplotDat(FILE* file, resultsVector rv)
{
  unsigned int i, sum;
  result r;
  fprintf(file, "#pos \tmin \tQ1 \tMedian \tQ3 \tmax \tMean \tCoverage \tqCoverage \tA \tC \tG \tT\n");
  	for(i = 0; i < rv.assignedLength; i++)
	{
	  r = rv.results[i];
	  
	  sum = r.A + r.C + r.G + r.T;
	 fprintf(file, "%u \t%u \t%u \t%u \t%u \t%u \t%u \t%u \t%u \t%f \t%f \t%f \t%f\n"
	               , i, r.minPhError, r.lowerPhQuartile, r.medianPhError, r.upperPhQuartile, r.maxPhError, r.meanPhError, r.coverage, r.qFloorCoverage, r.A/(double)sum, r.C/(double)sum, r.G/(double)sum, r.T/(double)sum); 
	  
	}
	fflush(file);
}

void printGnuplotDat_quality_aware_bases(FILE* file, resultsVector rv)
{
  unsigned int i, sum;
  double qSum;
  result r;
  fprintf(file, "#pos \tmin \tQ1 \tMedian \tQ3 \tmax \tMean \tCoverage \tqCoverage \tA \tC \tG \tT \tqA \tqC \tqG \tqT\n");
  	for(i = 0; i < rv.assignedLength; i++)
	{
	  r = rv.results[i];
	  
	  sum = r.A + r.C + r.G + r.T;
//        sum and qSum should actually be the same without floating-point errors
          qSum = r.qA + r.qC + r.qG + r.qT;
	 fprintf(file, "%u \t%u \t%u \t%u \t%u \t%u \t%u \t%u \t%u \t%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f\n"
	               , i, r.minPhError, r.lowerPhQuartile, r.medianPhError, r.upperPhQuartile, r.maxPhError, r.meanPhError, r.coverage, r.qFloorCoverage, r.A/(double)sum, r.C/(double)sum, r.G/(double)sum, r.T/(double)sum, r.qA/qSum, r.qC/qSum, r.qG/qSum, r.qT/qSum); 
	  
	}
	fflush(file);
}

void writeConsensusReferenceFile(FILE* file, resultsVector rv, char* name)
{

unsigned int i;
        
fprintf(file, ">%s\n", name);
result r;
    	for(i = 0; i < rv.assignedLength; i++)
	{
        r = rv.results[i];
            if(r.frequencyMajorBase > 0.25)
            {
            fprintf(file, "%c",r.majorBase);
            }else{
            fprintf(file, "N");    
            }
            
        }
        
fprintf(file, "\n");
}

//--Plotting-Functions------------------------------------------------
void plotQualityProfile(setting arg, int length, char* tempFileName)
{
unsigned int howMany = arg.numPlots;
unsigned int every;


if(howMany < length)
every = floor((double)length / howMany);
else
every = 1;



FILE * gnuplotPipe = popen("gnuplot -persistent", "w");

fprintf(gnuplotPipe,"set bars 2.0\n");
fprintf(gnuplotPipe,"set style fill solid 0.25 border\n");
fprintf(gnuplotPipe,"set yrange[0:]\n");
fprintf(gnuplotPipe,"set boxwidth %.2f\n", (double)every*2.0/3.0);
fprintf(gnuplotPipe,"set grid ytics lc rgb \"#bbbbbb\" lw 1 lt 0\n");
if(!arg.concealResults)
{
//fprintf(gnuplotPipe, "set term wxt title 'Quality Profiles'\n");
fprintf(gnuplotPipe, "set term x11 title 'Quality Profiles'\n");
//fprintf(gnuplotPipe, "plot '%s' every %u using 1:3:2:6:5 with candlesticks title 'Range [0.05:0.95] and Quartiles' whiskerbars, '%s' every %u using 1:4:4:4:4 with candlesticks lt -1 title 'Median', '%s' every %u u 1:7 smooth csplines title 'Mean'\n"
//		   , tempFileName, every, tempFileName, every, tempFileName, every);
fprintf(gnuplotPipe, "plot '%s' every %u using 1:3:2:6:5 with candlesticks title 'Range [0.05:0.95] and Quartiles' whiskerbars, '%s' every %u using 1:4:4:4:4 with candlesticks lt -1 title 'Median', '%s' u 1:7 smooth csplines title 'Mean'\n"
		   , tempFileName, every, tempFileName, every, tempFileName);
fflush(gnuplotPipe);
}

if(arg.outFilePrefix)
{
char *tempPNGName;
tempPNGName = (char*)calloc(strlen(arg.outFilePrefix)+16, sizeof(char));
sprintf(tempPNGName, "%s.qualityProfile", arg.outFilePrefix);

fprintf(gnuplotPipe,"set title 'Quality Profiles'\n");
fprintf(gnuplotPipe,"set terminal png\n");
fprintf(gnuplotPipe,"set output '%s.png'\n", tempPNGName);
//fprintf(gnuplotPipe, "plot '%s' every %u using 1:3:2:6:5 with candlesticks title 'Range [0.05:0.95] and Quartiles' whiskerbars, '%s' every %u using 1:4:4:4:4 with candlesticks lt -1 title 'Median', '%s' every %u u 1:7 smooth csplines title 'Mean'\n"
//		   , tempFileName, every, tempFileName, every, tempFileName, every);
fprintf(gnuplotPipe, "plot '%s' every %u using 1:3:2:6:5 with candlesticks title 'Range [0.05:0.95] and Quartiles' whiskerbars, '%s' every %u using 1:4:4:4:4 with candlesticks lt -1 title 'Median', '%s' u 1:7 smooth csplines title 'Mean'\n"
		   , tempFileName, every, tempFileName, every, tempFileName);

fflush(gnuplotPipe);

print_selective("Quality profiles output as:\t \"%s.png\"\n",tempPNGName);
free(tempPNGName);  
}

pclose(gnuplotPipe);
  
}


void plotCoverage(setting arg, char* tempFileName)
{


FILE * gnuplotPipe = popen("gnuplot -persistent", "w");

fprintf(gnuplotPipe,"set yrange[0:]\n");

if(!arg.concealResults)
{
fprintf(gnuplotPipe, "set term x11 title 'Coverage of Sites'\n");
fprintf(gnuplotPipe, "plot '%s' using 1:8 with lines title 'Coverage'\n", tempFileName);
if(arg.qFloor >= 0){
fprintf(gnuplotPipe, "replot '%s' using 1:9 with lines title 'Quality limited Coverage Q >= %u'\n", tempFileName, arg.qFloor);  
}
fflush(gnuplotPipe);
}

if(arg.outFilePrefix)
{
char *tempPNGName;
tempPNGName = (char*)calloc(strlen(arg.outFilePrefix)+17, sizeof(char));
sprintf(tempPNGName, "%s.coverageProfile", arg.outFilePrefix);

fprintf(gnuplotPipe,"set title 'Coverage of Sites'\n");
fprintf(gnuplotPipe,"set terminal png\n");
fprintf(gnuplotPipe,"set output '%s.png'\n", tempPNGName);
if(arg.qFloor < 0){
fprintf(gnuplotPipe, "plot '%s' using 1:8 with lines title 'Coverage'\n", tempFileName);
}else{
fprintf(gnuplotPipe, "plot '%s' using 1:8 with lines title 'Coverage', '%s' using 1:9 with lines title 'Quality limited Coverage Q >= %u'\n", tempFileName, tempFileName, arg.qFloor);  
}
fflush(gnuplotPipe);

print_selective("Coverage profile output as:\t \"%s.png\"\n",tempPNGName);
free(tempPNGName);  
}


pclose(gnuplotPipe);
  
}




void plotBaseFrequencies(setting arg, char* tempFileName)
{


FILE * gnuplotPipe = popen("gnuplot -persistent", "w");

fprintf(gnuplotPipe,"set yrange[0:]\n");

if(!arg.concealResults)
{
fprintf(gnuplotPipe, "set term x11 title 'Relative Base Frequencies'\n");
fprintf(gnuplotPipe, "plot '%s' using 1:10 with lines title 'A'\n", tempFileName);
fprintf(gnuplotPipe, "replot '%s' using 1:11 with lines title 'C'\n", tempFileName);
fprintf(gnuplotPipe, "replot '%s' using 1:12 with lines title 'G'\n", tempFileName);
fprintf(gnuplotPipe, "replot '%s' using 1:13 with lines title 'T'\n", tempFileName);
fflush(gnuplotPipe);
}

if(arg.outFilePrefix)
{
char *tempPNGName;
tempPNGName = (char*)calloc(strlen(arg.outFilePrefix)+22, sizeof(char));
sprintf(tempPNGName, "%s.baseFrequencyProfile", arg.outFilePrefix);

fprintf(gnuplotPipe,"set title 'Relative Base Frequencies'\n");
fprintf(gnuplotPipe,"set terminal png\n");
//fprintf(gnuplotPipe,"set term pngcairo\n");  //NOTE this produces nicer results but is less portable...
fprintf(gnuplotPipe,"set output '%s.png'\n", tempPNGName);

fprintf(gnuplotPipe, "plot '%s' using 1:10 with lines title 'A', '%s' using 1:11 with lines title 'C',  '%s' using 1:12 with lines title 'G', '%s' using 1:13 with lines title 'T'\n", tempFileName, tempFileName, tempFileName, tempFileName);
fflush(gnuplotPipe);

print_selective("Base frequencies output as:\t \"%s.png\"\n",tempPNGName);
free(tempPNGName);  
}


pclose(gnuplotPipe);
  
}


void plotMajorBase(setting arg, char* tempFileName)
{


FILE * gnuplotPipe = popen("gnuplot -persistent", "w");

fprintf(gnuplotPipe,"set yrange[0.86:]\n");

if(!arg.concealResults)
{
fprintf(gnuplotPipe, "set term x11 title 'Relative Base Frequencies'\n");
fprintf(gnuplotPipe, "plot '%s' using 1:($10) with dots title 'A'\n", tempFileName);
fprintf(gnuplotPipe, "replot '%s' using 1:($11) with dots title 'C'\n", tempFileName);
fprintf(gnuplotPipe, "replot '%s' using 1:($12) with dots title 'G'\n", tempFileName);
fprintf(gnuplotPipe, "replot '%s' using 1:($13) with dots title 'T'\n", tempFileName);
fflush(gnuplotPipe);
}

if(arg.outFilePrefix)
{
char *tempPNGName;
tempPNGName = (char*)calloc(strlen(arg.outFilePrefix)+22, sizeof(char));
sprintf(tempPNGName, "%s.baseFrequencyProfile", arg.outFilePrefix);

fprintf(gnuplotPipe,"set title 'Relative Base Frequencies'\n");
fprintf(gnuplotPipe,"set terminal png\n");
//fprintf(gnuplotPipe,"set term pngcairo\n");  //NOTE this produces nicer results but is less portable...
fprintf(gnuplotPipe,"set output '%s.png'\n", tempPNGName);

fprintf(gnuplotPipe, "plot '%s' using 1:($10) with dots title 'A', '%s' using 1:($11) with dots title 'C',  '%s' using 1:($12) with dots title 'G', '%s' using 1:($13) with dots title 'T'\n", tempFileName, tempFileName, tempFileName, tempFileName);
fflush(gnuplotPipe);

print_selective("Base frequencies output as:\t \"%s.png\"\n",tempPNGName);
free(tempPNGName);  
}


pclose(gnuplotPipe);
  
}


void printAvgCoverage(setting arg, resultsVector rv)
{
int l, diff = 0;
if(arg.qFloor > 0)
{
  if(rv.averageQFloorCoverage == 0)
diff = floor(log(rv.averageCoverage)/log(10) + 1) - ( 1 + floor(log(arg.qFloor)/log(10) + 1) );
else
  diff = floor(log(rv.averageCoverage)/log(10) + 1) - ( floor(log(rv.averageQFloorCoverage)/log(10) + 1) + floor(log(arg.qFloor)/log(10) + 1) );
}
//printf("avgC %f, avgQC %f, qfllo %d diff %d\n", rv.averageCoverage, rv.averageQFloorCoverage, arg.qFloor, diff);
//exit(1);

print_selective("\nAverage coverage per site:    ");

	for(l = 0; l > diff; l--)
	printf(" ");

print_selective("%.1f\n", rv.averageCoverage);

if(arg.qFloor > 0)
{
  print_selective("Average coverage with Q >= %d:  ", arg.qFloor);

	for(l = 0; l < diff; l++)
	printf(" ");

print_selective("%.1f\n", rv.averageQFloorCoverage);
}
print_selective("\n");
}

void printSamLine(int position, int position_offset, char* referenceName, char* name, char* seq, char* qual, int length,  unsigned int isComplemented, unsigned int mapping_quality, Placement* placement)
{
    
    unsigned int flag = 0;
    
    //QNAME
    char* t_name = strdup(name);
    cut_at_first_space(t_name);
    printf("%s\t",t_name);
    free(t_name);
    
    //FLAG
    if(isComplemented)
        flag = (flag | 16);
    
    if(position < 0)
        flag = (flag | 4);
    
    printf("%u\t", flag);
    
    //RNAME
    if(position < 0)
        printf("*\t");
    else
        printf("%s\t", referenceName);
    
    //POS
    if(position < 0)  
        printf("0\t");    
    else
    {
        if(placement)
        assert(position == placement->start_leftmost);
        
//      position_offset is used to translate from the placement on the concatenated sequences to the placement on the actual reference
        printf("%d\t", (position+1-position_offset));
    }
    
    //MAPQ TODO
    if(position < 0)
    {
    printf("255\t");
    }else{
    printf("%u\t", mapping_quality);    
    }
    //CIGAR TODO clipping
    if(position < 0)
    printf("**\t");
    else
    {
    printf("%dM\t", length);
    }
    //RNEXT
    printf("0\t");
    
    //PNEXT
    printf("*\t");
    
    //TLEN
    printf("*\t");
    
    //SEQ
    printf("%.*s\t", length, seq);
    
    //qualities
    printf("%.*s\n", length, qual);

    fflush(stdout);
}

void printSamHeader(char* nameReference, unsigned int referenceLength)
{
    printf("@HD\tVN:1.0\tSO:unsorted\n");
    printf("@SQ\tSN:%s\tLN:%u\n", nameReference, referenceLength);
}

void printSamHeader_additional_reference(char* nameReference, unsigned int referenceLength)
{
    printf("@SQ\tSN:%s\tLN:%u\n", nameReference, referenceLength);
}

void writeToSam(FILE* file, int position, char*name, char* seq)
{
    
    
    
}


// TODO istead of a loop for each position, a look-up table may be a cleaner solution... 
char* get_reference_name_by_position(globalVariables g, int pos)
{
    if(pos<0)
        return NULL;
    
    int i;
    
    for(i = 0; i < g.number_of_references; i++)
        {
            if( pos < (g.references_individual_start[i] + g.references_individual_lengths[i])  )
            {
            return g.reference_names[i];
            }
        }
       
    
    return NULL;
}

// TODO istead of a loop for each position, a look-up table may be a cleaner solution... 
int get_reference_offset_by_position(globalVariables g, int pos)
{
    if(pos<0)
        return 0;
    
    int i;
    
    for(i = 0; i < g.number_of_references; i++)
        {
            if( pos < (g.references_individual_start[i] + g.references_individual_lengths[i])  )
            {
            return g.references_individual_start[i];
            }
        }
       
    
    return 0;
}