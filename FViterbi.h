#pragma warning( disable: 4786 )

#ifndef _FVITERBI_H_
#define _FVITERBI_H_

#include "constdef.h"
#include "common.h"
#include "DataLoader.h"
#include "HMMBuilder.h"

class FVTBCandidate    //Viterbi candidates
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
    FVTBCandidate( );
    ~FVTBCandidate( );
    FVTBCandidate & operator=(  FVTBCandidate &org );

    void initFVTBCandidates( );    
    int hasData( );
    double getScore( );
    void getPosition( int pos[2] );
    void setAbsPosition( int baseidx );
    char * getAlgSequence( );
    char * getAlgStructure( );
    int getBothShortSeqs( char ** algseq, char ** algstu );
    void printAlgnedSeq( );
};

class HmmFilterResult
{
public:
    int expNum;
    int actNum;
    FVTBCandidate * hits;
    
public:
    HmmFilterResult( int topN );
    ~HmmFilterResult( );
    
    void sortHits( );
    int getActualNum( );
    FVTBCandidate * getResultArray(  );
    void setActualNum( int realNum );
    void printCands( );
};

class FVTBSearch
{
private:
    int *searchSeq;
    int seqlen;    // the length of search sequence
    int numstates;
  	double *** Vt;  //dynamic program table:             seqlen * numstates * 3
  	char *** Tb;    //trace back table storing pointers: seqlen * numstates * 3
  	double *** V;  //assistant Vt
  	char *** T;    //assistant Tb
    FVTBCandidate *tops;
    int numtops;
    
    //for linear computation
    int haverun;
    double newThreshold;
    FVTBCandidate prehit;  //for linear computation
    int preBegLoc;

    HmmFilterResult * hfr;
    
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

	//detect memory leak problem in rnatops zbhuang 20090815
	void free3DTablesAssist(  );

    void storePosIndex(  );
    void shiftPosIndex( int shiftnum );
    int checkShiftedSeq( int *instr, int newlen, int shiftnum);
    
    // Pick the maximum value from a, b and c
    double	max(double a, double b, double c);
    double	max(double a, double b);
    char	maxstate (double a, double b,double c);
    char	maxstate (double a, double b);
    
    int getTopnMaxProbs(  );
    int getNewProbs(  );
    double insertCurPos( int idi, int idj, int idk,  double cprob );
    int isInRemSeq( int *seq, int size,  int num);
    int VTBTraceBack(  FVTBCandidate &cad );
    void VTBimplement( Loop *curLoop );
    void VTBshift( Loop *curLoop, int shiftnum );
    
    void output3Dtables( char stateTable );
    int obtainCandidate( );
    
public:
    FVTBSearch( );
    FVTBSearch( int nummax );
    ~FVTBSearch( );
    
    FVTBCandidate * getCandHead( );
    
    void setHitThreshold( double probT );
    int setAllowedInsNum( int num );
    int FVTBSearchString( int *instr, int newlen, Loop *curloop );    
    int FVTBSearchString( char *instr, int newlen, Loop *curloop );  //read from char* array
    FVTBCandidate * SearchOneNoMerge( Loop *curloop, char *instr );
    FVTBCandidate * SearchOneMerge( Loop *curloop, char *instr );
    int SearchTopNMerge( Loop *curloop, char *instr, int subwsize, int topN );
    HmmFilterResult * getSearchResult(  );

    //****** print out *******
    void printSearchSeq( );
    void printResultSeq( );
};

#endif
