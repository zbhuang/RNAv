#ifndef _DATALOADER_H_
#define _DATALOADER_H_

#include "constdef.h"
#include "common.h"

class StemLocation
{
public:
    char charid;
    string strid;
    int  Pos[4];
    int  stemlen;
    int  distance;     //between two arms in PASTA
    double  * offsetE[2];    //distance of stem end to pasta end 0: left 1: right
    double  * offsetSD[2];   //xxxxxAAAA......AAAAxxxxx
    int numOffset;
    double  woffsetE[2];
    double  woffsetSD[2];
    double  armE[2];   //for the length of two arms 0: left 1: right
    double  armSD[2];
    double  midE;      //for the length of middle region
    double  midSD;

private:
    int getArmRegionAt( int idx, int reg[4], double coeff );
    int getArmRegionWhole( int reg[4], double coeff );
    
public:
    StemLocation( );
    ~StemLocation( );
    StemLocation & operator=(  StemLocation & org );
    
    void releaseEandSDArray( );
    int offsetAllE( double distance );
    int getStemRegion( int idx, int reg[4], double coeff );
    int obtainMostEnds( char LorR, double coeff );
    
    string getStemStrID( );
};

class LoopLocation
{
public:
    int  id;
    int  Pos[2];
    int  looplen;  //in PASTA
	double  offsetE[2];  //distance of loop end to pasta end
    double  offsetSD[2];
    double  lenE;           //the length of loop E and SD
    double  lenSD;
public:
    LoopLocation( );
    ~LoopLocation( );
};

class DataLoader
{
private:
	char *	filename;
	vector <string> orgpasta;
	vector <string> orgcommns;
    vector <string> orgtrains;
    int     hasCommentLines;
	int     seqlength;
	int     numtrainseqs;
	int     numPastaLines;
    double  splitThreshold;
    
    int  basePos;  //base position to calculate the offsets
    int  memoryAid;
    
    //extension variables
    int filterPos[2];
    int extensionLength[2];  //0: left part,  1: right part
   
    int readOriginaldata( vector <string> & tsdata );
    int inputMyData( vector <string> ts );
    void copyOrigSeqs( vector <string> ts );
    int countSeqLength( vector <string> ts );
    int countPastaLines( vector <string> ts );
    int countTrainingSeqNum( vector <string> ts );
    int checkSeqLength( );
    int copyPastaLines( );
    int copyTrainingSeqs( );
    int detectCommentLines( vector <string> ts );
    int seperateCommentAndTrainingSeq( vector <string> ts );
    string removeSpaceEnd( char * str  );
    int isEmptyLine( char * str  );
    int is2ndPastaLineFormat( char curbase );
  
    int getStemNum( );
    int getEndLoopPos( int &nextBeg );
    int locateStemPosition( );
    int locateLoopPosition( );
    void analyzePasta( );
    void calcLoopStemStat( );
    void offsetSumFromLeft( );  //sum of E's and SD's from left side
    //void calcOffsetStatOvr( );  //sum of overall E and SD from both sides
    //void offsetCalFromBoth(  ); //deprecated
    //void offsetCalFromLeft(  );
    
    void calcEandSD( int begpos, int endpos,  double & eval, double & sdval, int *subset=NULL, int thegrp=0 );
    void calcSpecialEandSD( int begpos, int endpos, double & eval, double & sdval, int *subset=NULL, int thegrp=0 );
    void sumEandSD( int beg, int end, double &sumE, double &sumSD );

    void mergePastaLine( );
    //filter
    int calcExtensionLength(  );

    //group training seqs
    int makeAllEpositive( );
    void groupSeqsOnFilter(  );
    void groupTrainSeq(  );
    int SplitTrainingSeq( int endpos, int *mark, int thegrp, int numgrps  );
    int getSubseqLen( int seqno, int begpos, int endpos );
    bool needSplit( double curE, double curSD );
    int calcGrpOffset( int stemno, int * leftMark, int* rightMark  );
    void printESDEveryStem(  );
    void printVectorStrs( vector <string> ts );

public:
	string *  stemchid;
	string *  mgpasta;
	char *  pasta[2];
	int  **	trainingseq;
	double  BaseHz[4];      //A, C, G, U

	StemLocation	*stemary;
	int	numstems;
	LoopLocation	*loopary;
	int	numloops;
	
	DataLoader( );
	~DataLoader( );
	
	int inputData( char * infileName );
	int inputData( vector <string> ts );
	
	void countBaseFreq( );
	void setBaseFreq( double HzA, double HzC, double HzG,  double HzU);
	void setSplitThreshold( double cutoff );
	void copyLogBgBaseHz( double * baseArray, int basenum );
	//filter setting
	void setFilterEndPosition( int filterBeg, int filterEnd, int clusteringBasePos=-1 );

	int getSeqLength( );
	int getNumOfTrainSeq( );
	int getScanWinLength( double coeff = 3, char sumAll='+' );
	int getMiniWinLength( double coeff = 3 );
    int getBothExtLength( int & leftExt, int & rightExt );
    int getSearchRangeBeforeHit( int begHitPos, int endHitPos, double coeff, int & posBeg, int & posEnd );
    
    int getPastaLines( );
	
	int ** getPointerOfSeqs( );
	char * getPointerOfPasta( );
    int getCertainTrainSeq( int idx, char *&seqbuf );
	void printBaseFreq( );
	void printPasta( );
	void printTrainingSeq( int idx=-1 );
	void printStemLoop( char which='B' );  //  S: stem L:loop B: both
	void printAll( );
	void printStringAll( );

    //test
    void getStemPosInSeq( int seqno,  int stemno,  int pos[4] );
    int getSeqBegPos( char *tgtseq, int seqno );

	int **	getPTrainingSeqs();
	char *	getPPasta();
	int		getNumstem();
	int		getNumloop();
	StemLocation *	getStemarray();
	LoopLocation *	getLooparray();

};

#endif

