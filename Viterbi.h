#ifndef _VITERBI_H_
#define _VITERBI_H_

#include "constdef.h"
#include "common.h"
#include "DataLoader.h"
#include "HMMBuilder.h"

class VTBCandidate    //Viterbi candidates
{
private:
     void freeAllocMemory( );
     void initAllMembers( );

public:
    int pos[3];    //position in 3D table:  seqlen, numstates, D/I/M
    int loc[2];    //0 beg ========== 1 end
     double prob;    //log probability
     int cadno;      //-1, only one;
     char *finseq;   //final sequence with gap
     char *finstr;   //final structure
     int finlen;
public:
    VTBCandidate( );
    ~VTBCandidate( );
    int hasData( );
    double getLogOdds(  );
    double getLog( );
     void initVTBCandidates( );
    void printAlgnedSeq( );
};

class VTBSearch
{
private:
    int *searchSeq;
    int seqlen;    // the length of search sequence
    int numstates;
  	double *** Vt;  //dynamic program table:             seqlen * numstates * 3
  	char *** Tb;    //trace back table storing pointers: seqlen * numstates * 3
    VTBCandidate *tops;
    int numtops;
    
    //deal with null model (allowing some insertions)
    int allowNumIns;
    
private:
    int obtainSequence(  char * filename );
    void freeSearchSeq(  );
    int malloc3DTables(  );
  	void free3DTables(  );
    void init3DTables(  );
    void initTopCands(  );
    int update3DandSeq( int newlen,  int newnumstate );
    // Pick the maximum value from a, b and c
    double	max(double a, double b, double c);
    double	max(double a, double b);
    char	maxstate (double a, double b,double c);
    char	maxstate (double a, double b);
    
    int getTopnMaxProbs(  );
    int getGlobalProbs(  );
    double insertCurPos( int idi, int idj, int idk,  double cprob );
    int isInRemSeq( int *seq, int size,  int num);
    int VTBTraceBack(  VTBCandidate &cad );
    void VTBimplement( Loop *curLoop );
    void output3Dtables( );
    
public:
    VTBSearch(  );
    VTBSearch( int nummax );
    ~VTBSearch( );
    
    VTBCandidate * getCandHead( );
    
    int setAllowedInsNum( int num );
    int VTBSearchString( int *instr, int newlen, Loop *curloop );    
    double VTBSearchMax( Loop *curloop, char *instr, int start, int &end, char **retpath=NULL, char **retSeq=NULL);//interface
    int VTBSearchString( char *instr, int newlen, Loop *curloop );  //read from char* array
    int VTBSearchFile( char *filename,  Loop *curloop ); //read from files, once

    //****** print out *******
    void printSearchSeq( );
    void printResultSeq( );
};

#endif
