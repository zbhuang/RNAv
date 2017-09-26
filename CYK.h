#ifndef _CYK_H_
#define _CYK_H_

#include "constdef.h"
#include "common.h"
#include "DataLoader.h"
#include "ScfgBuilder.h"

typedef struct Triple
{
    int prodtype;
    char self;    //For pair:  P: x -,   Q: - x,   W: - -,  S: x x
    char pre;
    char state; //right nonterminal
    int id;    //right nonterminal
} Triple;

class CYKCandidate    //CYK candidates
{
private:
     void freeAllocMemory( );
     void initAllMembers( );

public:
     int pos[3];     //position in 3D tables
     int loc[2][2];  //0 ls..le========== 1 rs..re
     int preg[2][2]; //pair region
     int numins[2][2]; //0 ls..le========== 1 rs..re
     double prob;    //log probability
     double penalty; //penalty for loop
     double indvpen[3]; //
     int distance;
     int cadno;      //-1, only one;
     char *finseq;   //final sequence with gap
     char *finstr;   //final structure
     int finlen;
public:
    CYKCandidate( );
    ~CYKCandidate( );
    CYKCandidate & operator=(  CYKCandidate & org );
    int hasData( );
    void calcRandomLog( double* bgLogProb );
    double getLogOdds(  );
    double getLog(  );
    void initCYKCandidates( );
    void printAlgnedSeq( );
};


class CYKSearch
{
private:
    int *searchSeq;
    int seqlen;    // the length of search sequence
    int numntm;    // the # of nonterminals
    double ***Prob;
    Triple ***Gram;
    double ***P;   //auxiliary Pb 
    Triple ***G;   //auxiliary Ga
    bool **rec;   // record the candidate for merging candidates
    
    bool candidmerge;
    bool candshortlen;
    int  allowshift;
    double scoredroprate;
    
    CYKCandidate *tops;
    int numtops;
    
    //SCFG
    double  offsetE[2];   //distance of stem end to pasta end. 0: left 1: right
    double  offsetSD[2]; 
    double  armE[2];      //0: left 1: right
    double  armSD[2];
    double  midE;         //middle loop region
    double  midSD;
    
    //condition for searching
    double coestd;
    int allowvar;
    int maxWid;
    int midLen[2];
    int d1; //index of sequence postion
    int d2; //index of width
    int exof;
    
    //length penalty paramters
    double plowerbound;  //if > lowerbound, we add penalty, else we don't
    double pcoeff;      //adjust the penalty level
    
private:
  	int obtainSequence(  char *filename );
    int isInRemSeq( int *seq, int size,  int num);
    
    void storePosIndex( );
    void shiftPosIndex( int num=1 );
    
    void resetRec( );
  	void initTopCands( );
  	int malloc3DTables( );
  	void free3DTables( );
    void init3DTables( );
    int update3DandSeq( int newlen,  int newnumntm );

    int getTopnMaxProbs( int base );
    int getMaximumProbs( CYKCandidate &cad );
    int compareSimilarCands( CYKCandidate * vault, CYKCandidate &cad, int currentNum );
    double insertCurPos( int base, int idi,  int idj, int idk, double cprob );
    double obtainRuleProb( int begid, int endid, int size, GramRule *prod);
    char computeMaxforPair(  GramRule *prod,  int ps,  int wd,   int dx,  int dy, double *curprob);
    char computeMaxforTMNPair(  GramRule *prod,  int j, int dx,  int dy,  double *curprob);
    
    int setBoundaryConds( Stem *curstem, int sequlen );
    double calcLengthPenalty( double dis, double sd );
    int CYKTraceBack(  CYKCandidate &cad );
    
    int CYKimplement( GramRule *prod,  int size,  NontermProp *ntms );
    int CYKOnTwoRegion( GramRule *prod,  int size,  NontermProp *ntms, int midReg[] );
    
    int checkShiftedSeq( int *instr, int newlen, int shiftnum=1);
    int CYKShift( GramRule *prod,  int size,  NontermProp *ntms, int shiftnum=1 );
    int obtainCandidate(  );
    int getLocalMaxIndex( CYKCandidate * vault, bool hasget[ ], int num );
    int rankCandidates( CYKCandidate * vault );
    
    
public:
    CYKSearch(  );
    CYKSearch( int topsize, double stdtimes, int extravar=0 );
    ~CYKSearch( );
    void setParameter( int candsize, double stdtimes, int extravar=0 );
    void setCandMergParameter( bool mergevalid, bool shortlen, int shiftnum, double droprate=PERCENTSCOREDROP );
    void setLengthPenaltyParameters( double lbound, double adjustk );
    
    int getTopCandNum( );
    CYKCandidate * getCandHead( );
    
    int CYKSearchString( char *instr, int newlen, Stem *curstem, int haverun=0, int shiftnum=1 );  //read from char* array
    int CYKSearchString( int  *instr, int newlen, Stem *curstem, int haverun=0, int shiftnum=1 );
    int CYKSearchTwoRegion( int *instr, int newlen, int midreg[2], Stem *curstem, int regNo  );
    int CYKSearchFile( char *filename,  Stem *curstem ); //read from files, for test only
    
    void printSearchSeq( );
    void printResultSeq( );
};

#endif
