#ifndef _CANDSCANNER_H_
#define _CANDSCANNER_H_

#include "constdef.h"
#include "common.h"
#include "ScfgBuilder.h"
#include "CYK.h"
#include "HMMBuilder.h"
#include "Viterbi.h"

typedef struct MaxElems
{
    int uid;
    double prob;
} MaxElems;

typedef struct CandContain
{
    char type;
    char idchar;   //ID of char
    string idstr;
    int UID;       //ID of stem or loop
    int num;
    
    int hour;
    int minu;
    double secd;
    
    class Candid *pMember;
} CandContain;

//------------------------------------  Candid  -----------------------------------
class Candid  //Candidates
{
public:
    char type; //Loop:L    Stem:S    Null:N
    int CompID;    //ID of Stem or Loop
    int UID;
    union
    {
        int stempos[2][2];  //0 ls..le========== 1 rs..re
        int looppos[2];
    };
    double prob;    //probability
    double penalty; //length penalty
    double indvpen[3]; //
    char *finseq;   //final sequence with gap
    char *finstr;   //final structure
    int finlen;
    Candid *nextCand;
     
public:
    Candid( );
    Candid( char candtype, int sd, int id );  //S stem;  L loop
    ~Candid( );
    Candid( Candid &org );
    Candid & operator=(  Candid & org );
    
    double getProbScore(  );
    void setNullCandid( );
    bool isNullCandid( );
    void copyValue( CYKCandidate &pCand, int bpos);
    void copyValue( VTBCandidate &pCand, int bpos);
    int extractAlgPart( char* &seqlm, char* &alglm, int &lenlm,
                        char* &seqrm, char* &algrm, int &lenrm );
    
    void printAlgnedSeq( );
    void printAlgnedSeqOnly( );
};

//------------------------------------  CandScanner  -----------------------------------
class CandScanner
{
private:
    int *searchSeq;
    int seqlen;
    int topnum;
    int newnum;
    double threshold;
    int shiftnum;
    
    bool allowOffset;
    bool addNullCand;
    
private:
    //topn works:
    void putTopStemCandInArray( CYKCandidate *pResult, int size, int bpos );
    void putTopLoopCandInArray( VTBCandidate *pResult, int size, int bpos );
    int resizeCandArray(  );
    
    //topn doesn't work:
    int putStemCandInChain( Candid* &pTail, CYKCandidate *pResult, int stemid, int size, int bpos, int uid );
    int putLoopCandInChain( Candid* &pTail, VTBCandidate *pResult, int loopid, int size, int bpos, int uid );
    int getTopnCandids( Candid *phead, int compid );
    int getNumOfCandids( Candid *phead );
    
    int isDifferent( VTBCandidate *pVC, int basepos );
    int isDifferent( CYKCandidate *pCC, int basepos );

    void printCandid( Candid *phead, int numcd );
 
public:
    Candid *pCand;
    CandScanner( int topn=5,  double thres=INVALID );
    ~CandScanner(  );
    int setTopn( int topn );
    int setThreshold( double thres );
    int setShiftNum( int num );
    void setOffset( bool enable);
    void addNullCandidate( bool enNullCandidate );
    int getTopnum( );
    int getSeqlen( );
    
    int resetMembers( );
    int setSeqFromFile( char *filename );
    int setSeqFromArray( char *seqArray );
    void freeCandidChain( Candid* & phead );
    int hasData( );
    void printSearchSeq( );

    Candid * scanTopnInTwoRegion( Stem *pstem,  CYKSearch & curcyk, double coeff, int basepos=0);
    //speedup
    Candid * scanTopnInOneRegionSP( Stem *pstem, CYKSearch &curcyk, double coeff, int runned );

    Candid * scanCandidates( Stem * pstem,  CYKSearch& curcyk,  double coef );
    Candid * quickScanCands( Stem * pstem,  CYKSearch& curcyk, double coef, int runned );
};


//------------------------------------  CandCollect  -----------------------------------
class CandCollect
{
private:
    int stemnum;
    int loopnum;
    
    //the times of SD
    double coefStDev;
    
    //speedup
    int *flagrun; //0: fisrt, 1; second...
    int thresnum;
    CYKSearch *cyks;
    
    //set CYK Searching parameters
    bool benOffset;
    bool candidmerge;
    bool candshortlen;
    int  allowshift;
    double scoredroprate;
    
    //set length penalty
    double  plbd;
    double pcoef;
    
    //set the method of scan: two region or one region
    bool bTwoRegion;
    
    //return additional null candidate
    bool nullcandidate;
    
    //for test only
    int tgtseqlen;
    int stemerr[4];
    int hopepos[4];
    int rdmlen;
    
    int hasDataInIdx( char type,  int idx);
    int addToCollect( char LorS, string idstr, int idx, Candid * &phead, int curnum, TimeCalc * timc=NULL );
    void freeInsideCands(  );
    
    void initCYKInstance(int topsize, double coeff, int varvalue );
    
public:
    CandContain *stemgrp;
    CandContain *loopgrp;
    
    CandCollect( );
    CandCollect( int stems,  int loops );
    ~CandCollect( );
    
    void setTwoRegionScan( bool scanMethod=true );
    void setCoeffStandardDev( double coeffsd );
    void setCYKMergeParamter( bool mergevalid, bool shortlen, int shiftnum, double droprate=PERCENTSCOREDROP );
    void setCYKLengthPenalty( double lbound, double adjustk );
    void setNullCandidate( bool enable );
    void enableOffset( bool offsetvalid );
    
    int getAlignmentStr( int stemNo, int candNo, char* &seqlm, char* &alglm, int &lenlm,
                                                 char* &seqrm, char* &algrm, int &lemrm );
    
    void freeAllCandids( );
    void printCandids( char type='B', int displayfmt=1 ); //B: stem and loop;   S: stem;   L:loop
                                                          //0: align part   1:all
    int searchCandids( Stem *pstem, int numsm, int topn, double thres, char * pStr, int basepos=0 );
    int quickStemCands( Stem *pstem, int numsm, int topn,  double thres,  char * pStr);

    //for testing only
    void setStemPosAccuStd(int* posErr,  int* armlen, int halfradmlen ); //0-based
    void checkStemPosAccuracy( Stem *pstem, int* accuracy );

};

#endif

