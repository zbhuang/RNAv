#ifndef _SCFGBUILDER_H_
#define _SCFGBUILDER_H_

#include "constdef.h"
#include "common.h"
#include "DataLoader.h"
#include "HMMBuilder.h"

// -------------- struct ---------------
typedef struct BulgeInfo
{
    char LorR;
    int beg;
    int end;
} BulgeInfo;

// -------------------------------- class ScfgGram  ------------------------------------
class ScfgGram
{
public:
    char state; //stem: R:right_insert; L:left_insert; S:Match_Deletion; B:Start; E: End;
                //V: virtual state for left bulge --> right bulge
                //loop: M:match;   D:deletion  I:insertion
    int  ID;    //readable ID
    int  SubID; //sub readable ID
    int  UID;   //universial ID
    double emprobS[5][5];
    double emprobI[2][4];
    double trprob[5];        //0--LH,  1--LI, 2--S&M   3--RI,  4--RH
     
    class Bulge *nextLH;     //0
    ScfgGram *nextLI;        //1
    ScfgGram *nextS;         //2
    ScfgGram *nextRI;        //3
    class Bulge *nextRH;     //4

private:
    void log1DProb( double mat[],  int  );
    void log2DProb( double* mat[], int, int  );
    
public:
    ScfgGram( );
    ScfgGram(char settype,  int setid,  int setuid, int subid=0 );
    ~ScfgGram( );
    
    void logGrammarProb( double logfreq[5], double priorMatrix[5][5], double h );
    int countNumRules(  );
    int hasLinkVState(  );
   
};

// ---------------------------------- class Stem ---------------------------------------
class Stem
{
public:
    int UID;
    char idchar;   //id of stem
    string idstr;
    int Pos[4];    //LStart, LEnd, RStart, REnd;
    int stemlen;   //the length of stem
    int numseqs;   //the number of training seqs
    int distance;  //length between two arms: number of residues
    char *state[2];     //two states: s or i(insertion)
    char *armpasta[2];	//two parts of stem pasta, left part and right part
    int  **armseq[2];    //two parts of stem sequence, left part and right part

    int  **armorg[2];    //***Mar-20-2009 orginal sequences 
    int numorgs;         //***Mar-20-2009
    vector<bool> goodstem;
    
    double  bgBaseLogHz[5];
    ScfgGram * pStemGram;

    //column pairs freq
    int ***colpairNum;
    int sumcolpair[5][5];
    double stemStrength;
    
    //prior parameters
    double priorPairFreq[5][5];
    double priorcoef;
    //for CYK
    NontermProp *nonterm;
    int numntm;   // the # of nonterminal
    GramRule * prod;
    int numprod;  // the # of productions

    StemLocation *regInfo;

private:
    int logProb;
    //good stem detection
    bool isGoodStem( int idx );
    
    //log Prob
    void logProbforStem( );
    void calcColumnPairFreq( int colno, int Lpos, int Rpos );
    void calcStemStrength( );
    
    //Array all of the nonts and rules ( related functions )
    int countNumAllNtms( );
    int countNumAllRules( );
    void initNontermArray( );
    void initGramProd( );
    void arrayNonterminal( );
    
    //Generate all of the rules in an array
    int genHMMtoHMM( int idx,  Bulge *curHG,   Bulge *nexHG,   char dirtype);
    int genHMMtoHMMviaV( int idx, Bulge *curHG, Bulge *nexHG, char dirtype, double vprob);
    
    int genHMMtoSCFG( int idx,  Bulge *curHG,   ScfgGram *curSG,  char dirtype);
    int genHMMtoSCFGviaV( int idx, Bulge *curHG, ScfgGram *curSG, char dirtype, double vprob);
    
    int genSCFGtoHMM( int idx,  ScfgGram *curSG,  Bulge *curHG,  char dirtype);
    int genSCFGtoHMMviaV( int idx,  ScfgGram *curSG,  Bulge *curHG, char dirtype, double vprob);
    
    int genSCFGtoSCFG( int idx,  ScfgGram *curSG,  ScfgGram *nextSG);
    int genSCFGtoSCFGviaV( int idx,  ScfgGram *curSG,  ScfgGram *nextSG, double vprob );

    int genHMMtoVProd( int idx,  Bulge *curBg,  ScfgGram *curVt,  char dirtype);
    int genSCFGtoVProd( int idx,  ScfgGram *curSG,  ScfgGram *curVt );
    
    int serialProducts( );

    //allocate memory
    void freeFreqMem(  );
    void allocateFreqMem(  );
public:
    Stem( );
    ~Stem( );
    void freeScfgGrams( );
    int isProdsArray( );
    
    void setPriorPairFreq( double priorMatrix[5][5], double h );
    void pickConservedOnes( );
    //strength
    double getStemStrength( );
    //# of sequences involved into modeling
    double getPortionGoodStem( );
    
    //region
    int getStemSpan( double coeff=COESTDEV, int ismax=1 );
    int getMidlen( double coeff=COESTDEV, int ismax=0 );
    int getBothArmRegion( int no, int region[4], double coeff  );
    int getNumPairRegion(  );
    
    //print-out functions
   	void printPairFreq( bool isParallel=false);
   	void printTrainInfo( );
    void printOrgSeqs( );
  	void printScfgGram( );
   	void printNonterms( );
    void printAllRules( );

    //count the number of pair freq
    void countPairFreq(   );
    
  	//generate Productions in an array for CYK search
    int generateProdsArray( );
};

// ------------------------------- class BuildScfg ----------------------------------
class ScfgBuilder
{
private:
    Stem	*allstems;
  	 int    numstems;
    char	*stemchid;
    int    **sequences;
   	char	*pasta;
    int    seqlength;
    int    numseqs;      //the number of training seqs
    int    numorgs;
    double  pseudocount;
    double  pseudocntgap;
    
    //dataLoader class
    DataLoader *dl;
    
    int sumPairFreq[5][5];
    double priorPairFreq[5][5];
    double priorcoef;
    char * priorfn;
    //training file processing functions
   	int setTrainingData( DataLoader &trdata );
   	void seperateStem( );
  	
   	//probability calculation 
    char determPairState(Stem *curstem, int lidx, int ridx );
    int obtainNextS( Stem *curstem, int begl, int begr, int &distl, int &distr);
   	void markState( Stem *curstem, int begl, int begr, int distl, int distr);
  	 void countBulge( int *bulge, Stem *pstem, int begl,  int begr,  int distl,  int distr );
    int judgeBulge( int *bulge,    BulgeInfo* &leftB,    BulgeInfo* &rightB,
                             Stem *pstem,   int begl,  int begr,  int distl,  int distr );
    
    void calcIBaseEmProbRight( ScfgGram *pSG, Stem *pstem, int beg, int dist );
    void calcIBaseEmProbLeft( ScfgGram *pSG, Stem *pstem, int beg, int dist );
    void calcITransPRight( Stem *pstem,  int beg,  int dist,  double *trprb);
    void calcITransPLeft( Stem *pstem,  int beg,  int dist,  double *trprb);
    void calcPairEmProb( ScfgGram *pSG, Stem *pstem, int lpos, int rpos);
    void calcSTransP( ScfgGram *pSG,  ScfgGram *pVT, Stem *curstem,
                      int *bulge,   BulgeInfo *leftB,  BulgeInfo *rightB,
                      int begl,  int begr,  int distl,  int distr );
    void postprecssBothBulge( ScfgGram *preG,  ScfgGram *curV );
    void calcLoopTransPLeft(  Stem *pstem, int beg, int dist, double &lpi, double &lplp );
    void calcLoopTransPRight( Stem *pstem, int beg, int dist, double &lpi, double &lplp );
    int modelRightArm(  int numbulge,    BulgeInfo *rightB,  Stem *curstem,
                        ScfgGram *preG,  ScfgGram *curG,     int idr, int disr, int idx, int uid );
    int modelLeftArm( int numbulge,    BulgeInfo *leftB,  Stem *curstem,
                      ScfgGram *preG,  ScfgGram *curG,    int idl, int disl, int idx, int uid );
    //key upper-level function
  	 void generateSCFG( Stem *curstem );
  	 void countPairInformation(  );
    int putRNAPairHz( );
    
public:
	
    ScfgBuilder( );
    ~ScfgBuilder( );
    
    Stem * getAllStems( int &num );
    void setPseudocnt( double pdcnt=PSEUDOCNT, double pcgap=PSDCNTGAP );
    void setPriorCoef( double h );
    void setPriorFilePath( char * fn );
    int buildSCFG( DataLoader &trdata );
    //get the maximum strength
    int getMaxStrengthStemNo( );
    //free all
    void freeAllStems( );
    
    //print-out functions
    void printSeqFile( );
   	void printAllGrams( int trinfo=0 );
   	void printStemGram( char stemch, int trinfo=0 );
   	void printPairFreq( char stemch=' ' ); //' ': print all
    void printSumPairHz( bool isParallel=true );
    void printPriorFreqMatrix( );
};
#endif

