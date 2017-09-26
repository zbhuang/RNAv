#ifndef _COMMON_H_
#define _COMMON_H_

#include "constdef.h"
#include "timecalc.h"

#include <string.h>	//gcc3.4

int isRNABase(char nucleotide);
int isRNABase(int nucleotide);
int chartonum(char nucleotide);
char numtochar(int base);

bool isCanonicalPair(char basex, char basey);

int isFileName( char *str);

void printPairEmP( double prob[5][5] );
void printPairEmP( int pairprob[5][5] );
void printPairHz( int pairprob[5][5], double prob[5][5]);
void printBaseEmP( double prob[4] );
void printSingleEmProb( double baseprob[4] );
int obtainSequence(  char * filename, char* &bufseq );

//calculate normalization of two matrix
void normalizeMatrix( double orgm[5][5], double priorm[5][5], double h);

int isEven(int iVal);

#endif
