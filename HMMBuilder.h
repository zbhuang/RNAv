#ifndef _HMMBUILDER_H_
#define _HMMBUILDER_H_

#include "constdef.h"
#include "common.h"
#include "DataLoader.h"
#include "ScfgBuilder.h"

// -------------------------------- class HMMGram  ------------------------------------
class HMMGram
{
public:
    char	 state;  //I:insertion;   M:Matching;   D:Deletion;   Start:B;   End:E
    int     ID;
    int    SubID;
    int    UID;
    double emission[4];
    double trprob[3];

    HMMGram	*nextD;	//Deletion
    HMMGram	*nextI;	//Insertion
    HMMGram	*nextM;	//Matching

    HMMGram( );
    HMMGram( char settype,  int setid,  int setuid,  int subid=0  );
    ~HMMGram( );
    
    void setHMMGram( char settype,  int setid,  int setuid,  int subid  );
    void logGrammarProb( double logfreq[5] );
    int countNumRules( );
};

// ---------------------------------- class Stem ---------------------------------------
class Loop
{
public:
    int UID;
    int Pos[2];    //Start, End
    int looplen;   //the length of loop
    int numseqs;   //the number of training seqs
    char *state;     //states: m(d) or i
    int numstates;   //B0, M1, M2, ..., En (n+1)
    int  **loopseqs;    //this loop part
    int logProb;
    double  bgBaseLogHz[5];
    HMMGram * pLoopGram;

    double  offsetE[2];  //distance of loop end to pasta end
    double  offsetSD[2];
    double  lenE;           //the length of loop E and SD
    double  lenSD;

public:
    Loop( );
    Loop( int seqnum, int lplen );
    virtual ~Loop( );
    int getLoopLength( ); //length > 0
    int isLogProb( );
    int logProbforLoop( );
    void freeHMMGrams( );
    
    int getLeftOffset( double coeff=COESTDEV, int ismax=0 );
    int getRightOffset( double coeff=COESTDEV, int ismax=0 );
    int getLoopSpan( double coeff=COESTDEV, int ismax=0 );
    //print-out functions
    void printLoopGrams( char direct='L' );
  	virtual void printHMMGrams( char direct='L' );
};


// ---------------------------------- class Bulge ---------------------------------------
class Bulge: public Loop
{
public:
    // for a bulge in stem
    class ScfgGram	*nextI;	//Scfg Grammer Rules:  trprob[0]
    class ScfgGram	*nextS;	//Scfg Grammer Rules:  trprob[1]
    Bulge *nextLp;
    double trprob[2];
    
    //for DP algorithm:  reserved
    NontermProp *nonterm;  
    int numntm;   // the # of nonterminal
    GramRule * prod;
    int numprod;  // the # of productions

private:
    //for them in stems
    void revLoopSeqs( );
    void genSingleHMMProd( HMMGram *curHmm, HMMGram *nextHmm, GramRule *prod, int idx, char diretype );

public:
    Bulge( );
    Bulge( int seqnum, int lplen );
    virtual ~Bulge( );
    int logProbforBulge( double stemBaseHz[5] );
    int countNumAllRules( );
    int hasLinkVState( );
    int countNumAllNtms( );
    
    void arrayNonterminal( NontermProp *pntms );

    char getFirstGramState( );
    int getFirstGramUID( );
    char getLastGramState( );
    int getLastGramUID( );
    int serialProducts( GramRule * prod,  int hidx, char diretype );
  	virtual void printHMMGrams( char direct='L' );
};

// ------------------------------- class HMMBuilder ----------------------------------
class HMMBuilder
{
private:
    Loop	* allloops;
  	int      numloops;
    int    **sequences;
  	char	* pasta;
    int      seqlength;
    int      numseqs;      //the number of training seqs
    double  pseudocount;
    double  pseudocntgap;  //not used

    int   setTrainingData( DataLoader &trdata );
    void  seperateLoop(  );

    char  determColomnState( Loop *curloop, int curpos );
    int   obtainNextM( Loop *curloop, int curidx, int &dist );
    void  markState(Loop *curloop, int idx, int dist );
    void  countTransState( HMMGram *preHM, Loop *curloop, int idx, int dist, int **Tri );
    void  calcComTransProb( HMMGram *preM, HMMGram *preI, HMMGram *preD, int **Trans );
    void  calcEndTransProb( HMMGram *preM, HMMGram *preI, HMMGram *preD,  int **Trans );
    void  calcEmissionProb( HMMGram *curHG, Loop *curloop, int pos );
    void  calcInsertEmProb( HMMGram *curHG, Loop *curloop, int idx, int dist );
    int   generateHMM( Loop *curloop, int uid=0, int subid=0 );
    
public:
    HMMBuilder( );
    ~HMMBuilder( );

    void setPseudocnt( double pdcnt=PSEUDOCNT, double pcgap=PSDCNTGAP );
    int buildHMM( DataLoader &trdata );
    Loop * getAllLoops( int &num );
    //interface to SCFG
    Bulge * buildBulge( int **sequences, int seqnum, int posbeg, int len, char LorR, int &gramuid, int subid );

    void freeAllLoops( );
    //print-out functions
    void printAllGrams( int trda=0 );
};

#endif
