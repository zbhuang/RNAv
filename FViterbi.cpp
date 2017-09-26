
#include "FViterbi.h"

//********************************************************************************
//*******************    implementation of VTBCandidate class      ***************
//********************************************************************************
FVTBCandidate::FVTBCandidate( )
{
    #ifdef DEBUG_DEV
        printf("FVTBCandidate::FVTBCandidate\n");
    #endif
    initAllMembers( );
}

FVTBCandidate::~FVTBCandidate( )
{
    #ifdef DEBUG_DEV
        printf("FVTBCandidate::~FVTBCandidate\n");
    #endif
    freeAllocMemory( );
}

void FVTBCandidate::initAllMembers( )
{
    #ifdef DEBUG_DEV
        printf("FVTBCandidate::initAllMembers\n");
    #endif
    int i=0;
    
    pos[0]=-1;
    pos[1]=-1;
    pos[2]=-1;
    for(i=0; i<2; i++ )
        loc[i]=INVALID;
    
     prob=INVLDPROB;     //probability
     cadno=-1;      //-1, only one;
     finseq=NULL;   //final sequence with gap
     finstr=NULL;   //final structure
     finlen=0;
}

void FVTBCandidate::initFVTBCandidates( )
{
    #ifdef DEBUG_DEV
        printf("FVTBCandidate::initFVTBCandidates\n");
    #endif
    freeAllocMemory( );
    initAllMembers( );
}

void FVTBCandidate::freeAllocMemory( )
{
    #ifdef DEBUG_DEV
        printf("FVTBCandidate::freeAllocMemory\n");
    #endif
    if( finseq!=NULL )
    {
        delete [ ] finseq;
        finseq=NULL;
    }
    if( finstr!=NULL )
    {
        delete [ ] finstr;
        finstr=NULL;
    }
}

int FVTBCandidate::hasData( )
{
    #ifdef DEBUG_DEV
        printf("FVTBCandidate::hasData\n");
    #endif
    if( loc[0]==INVALID )
        return 0;
    else
        return 1;
}

double FVTBCandidate::getScore(  )
{
    return prob;
}

void FVTBCandidate::getPosition( int pos[2] )
{
    pos[0]=loc[0];
    pos[1]=loc[1];
}

void FVTBCandidate::setAbsPosition( int baseidx )
{
    loc[0]=loc[0] + baseidx;
    loc[1]=loc[1] + baseidx;
}

int FVTBCandidate::getBothShortSeqs( char ** algseq, char ** algstu )
{
    int i=0, j=0, actLen=0, beg=-1;
    char *curalgseq=NULL;
    char *curalgstu=NULL;
    
    if( hasData( ) ) {
        for( i=0; i<finlen; i++ )
        {
            if( finstr[i]!='.' ) {
                actLen++;
                if( beg==-1 )
                    beg=i;
            }
        }
        curalgseq = new char [ actLen+1 ];
        curalgstu = new char [ actLen+1 ];
        for( i=beg, j=0; j<actLen; i++, j++ ) {
            curalgseq[j]=finseq[i];
            curalgstu[j]=finstr[i];
        }
        curalgseq[j]='\0';
        curalgstu[j]='\0';
    }
    *algseq = curalgseq;
    *algstu = curalgstu;
    return 1;
}

char * FVTBCandidate::getAlgSequence(  )
{
    char * algseq = NULL;
    
    if( hasData( ) ) {
        algseq = new char [ finlen+1 ];
        strcpy( algseq, finseq );
    }
    else
        algseq = NULL;
    return algseq;
}

char * FVTBCandidate::getAlgStructure( )
{
    char * algsut = NULL;
    
    if( hasData( ) ) {
        algsut = new char [ finlen+1 ];
        strcpy( algsut, finstr );
    }
    else
        algsut = NULL;
    return algsut;
}

void FVTBCandidate::printAlgnedSeq( )
{
    #ifdef DEBUG_DEV
        printf("FVTBCandidate::printAlgnedSeq\n");
    #endif
    if( cadno != -1 )
    {
        printf("No.%d:  ", cadno);
    }
    printf("Max Log-Prob: %10.6f at [%d, %d]\n", prob, loc[0], loc[1]);
    //printf("Result:\n");
    printf("   ");
    printf("%s\n", finseq );
    printf("   ");
    printf("%s\n\n", finstr );
}

    int pos[3];    //position in 3D table:  seqlen, numstates, D/I/M
    int loc[2];    //0 beg ========== 1 end
    double prob;    //log probability
    int cadno;      //-1, only one;
    char *finseq;   //final sequence with gap
    char *finstr;   //final structure
    int finlen;

FVTBCandidate & FVTBCandidate::operator=(  FVTBCandidate &org )
{
    if (this == &org)
        return *this;

    pos[0]=org.pos[0];
    pos[1]=org.pos[1];
    pos[2]=org.pos[2];

    loc[0]=org.loc[0];
    loc[1]=org.loc[1];

 	prob=org.prob;
 	cadno=org.cadno;

 	finlen=org.finlen;
    if( finseq!=NULL )
    {
        delete [ ] finseq;
    	finseq=NULL;
    }
    if( finstr!=NULL )
    {
    	delete [ ] finstr;
    	finstr=NULL;
    }
 	if( finlen>0)
    {
 	    finseq=new char[finlen+1];
 	    strcpy( finseq, org.finseq);
 	    finstr=new char[finlen+1];
 	    strcpy( finstr, org.finstr);
 	} 	  
 	return *this;
}

//********************************************************************************
//*******************    implementation of HmmFilterResult class         ***************
//********************************************************************************
HmmFilterResult::HmmFilterResult( int topN )
{
    #ifdef DEBUG_DEV
        printf("HmmFilterResult::HmmFilterResult\n");
    #endif
    expNum = topN;
    actNum = 0;
    hits = new FVTBCandidate[ expNum ] ;
}

HmmFilterResult::~HmmFilterResult(  )
{
    #ifdef DEBUG_DEV
        printf("HmmFilterResult::~HmmFilterResult\n");
    #endif
    if( hits!=NULL ) {
        delete [ ] hits;
        hits=NULL;
    }
}

int HmmFilterResult::getActualNum( )
{
    return actNum;
}

void HmmFilterResult::setActualNum( int realNum )
{
    actNum = realNum;
}

FVTBCandidate * HmmFilterResult::getResultArray(  )
{
    return hits;
}

void HmmFilterResult::printCands(  )
{
    #ifdef DEBUG_DEV
        printf("HmmFilterResult::printCands\n");
    #endif
    int i=0;
    int algpos[2];
    algpos[0]=-1;
    algpos[1]=-1;
    double probv=0.0;
    if( actNum==0) {
        printf("No hits found\n");
        return;
    }
    
    for(i=0; i<actNum; i++ ) {
        hits[i].getPosition( algpos );
        probv = hits[i].getScore( );
        printf("prob=%f, position:[%d - %d]\n", probv, algpos[0], algpos[1]);
        printf("  %s\n  %s\n\n", hits[i].finseq, hits[i].finstr);
    }
}

void HmmFilterResult::sortHits( )
{
    #ifdef DEBUG_DEV
        printf("HmmFilterResult::sortHits\n");
    #endif
    //printCands(  );
    if( actNum==0 || hits==NULL)
        return;
    //sort the final hits
    int i=0, j=0, k=0, stiNum=0;
    int addindex=-10;
    int *sti = new int [ actNum ];
    for( i=0; i<actNum; i++ )
        sti[i]=-1;
    
    FVTBCandidate * realCand = new FVTBCandidate [ actNum ];
    for( i=0; i<actNum; i++ ) {
        //printf("xxx stiNum=%d/actNum=%d xxx\n", stiNum, actNum );
        addindex=-10;
        if( i==0 ) {
            stiNum++;
            sti[0]=0;
        } else {
            for( j=0; j<stiNum; j++ ) {
                if( hits[ sti[j] ].getScore( ) < hits[i].getScore( )) {
                    addindex=j;
                    break;
                }
            }
            if( addindex != -10 ) {
                for( k=stiNum-1; k>=addindex; k-- )
                    sti[k+1]=sti[k];
                sti[ addindex ]=i;
            } else {
                sti[ stiNum ]=i;
            }
            stiNum++;
        }
    }
    //copy
    for( i=0; i<actNum; i++ ) {
        realCand[i]=hits[ sti[i] ];
        //printf(" %d(%f)  ", sti[i], hits[ sti[i] ].getScore( ) );
    }
    //printf("\n");
    delete [ ] hits;
    hits = realCand;
    
    if( sti!=NULL ) {
        delete [ ] sti;
        sti=NULL;
    }
    return;
}


//********************************************************************************
//*******************    implementation of FVTBSearch class         ***************
//********************************************************************************
FVTBSearch::FVTBSearch(  )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::FVTBSearch\n");
    #endif
    searchSeq=NULL;
    seqlen=0;
    numstates=0;
  	Vt=NULL;
  	Tb=NULL;
    numtops=1;
    allowNumIns=0;
    tops=new FVTBCandidate[ numtops ];
    
    haverun=0;
    newThreshold=0.0;
    preBegLoc=-1;
    hfr=NULL;
}

FVTBSearch::FVTBSearch( int nummax )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::FVTBSearch( parameter )\n");
    #endif
    searchSeq=NULL;
    seqlen=0;
    numstates=0;
  	Vt=NULL;
  	Tb=NULL;
  	V=NULL;
  	T=NULL;
  	haverun=0;
    numtops=nummax;
    allowNumIns=0;
    tops=new FVTBCandidate[ numtops ];
    
    haverun=0;
    newThreshold=0.0;
    preBegLoc=-1;
    hfr=NULL;
}

FVTBSearch::~FVTBSearch(  )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::~FVTBSearch\n");
    #endif
    freeSearchSeq( );
    free3DTables( );

	//detect memory leak problem in rnatops zbhuang 20090815
	this->free3DTablesAssist();

    if(tops!=NULL)
    {
        delete [ ] tops;
        tops=NULL;
    }
    if( hfr != NULL )
    {
        delete hfr;
        hfr=NULL;
    }
}

void FVTBSearch::initTopCands( )
{
    int i=0;
    for(i=0; i<numtops; i++ )
    {
        tops[i].initFVTBCandidates( );
    }
}

int FVTBSearch::malloc3DTables(  )
{
    #ifdef DEBUG_DEV
       printf("FVTBSearch::malloc3DTables\n");
    #endif
    
    int i=0, j=0;
    //assistant
    V = new double** [ seqlen+2  ];
    T = new char** [ seqlen+2  ];
    
    //***** 1D  [seq][ ][ ]
    Vt = new double** [ seqlen+2  ];
    if( Vt==NULL )
        return 0;
    
    Tb = new char** [ seqlen+2  ];
    if( Tb==NULL )
        return 0;
    
    //***** 2D  [ ][ state ][ ]
    for( i=0; i<seqlen+2; i++ )
    {
        Vt[i]=new double* [ numstates ];
        if( Vt[i]==NULL )
            return 0;
        Tb[i]=new char* [ numstates ];
        if( Tb[i]==NULL )
            return 0;
    }
    
    //***** 3D  [ ][ ][ M/I/D ]
    for( i=0; i<seqlen+2; i++ )
    {
        for( j=0; j<numstates; j++ )
        {
            Vt[i][j]=new double [ 3 ];
            if( Vt[i][j]==NULL )
                return 0;
            Tb[i][j]=new char [ 3 ];
            if( Tb[i][j]==NULL )
                return 0;
        }
    }
    
    storePosIndex( );
    init3DTables ( );
    
    return 1;
}

void FVTBSearch::storePosIndex( )
{
    int i=0;
    for( i=0; i<seqlen+2; i++ )
    {
        V[i]=Vt[i];
        T[i]=Tb[i];
    }
}

void FVTBSearch::freeSearchSeq(  )
{
    #ifdef DEBUG_DEV
       printf("FVTBSearch::freeSearchSeq\n");
    #endif
    if( searchSeq!=NULL )
    {
        delete [ ] searchSeq;
        searchSeq=NULL;
    }
}

void FVTBSearch::free3DTables(  )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::free3DTables\n");
    #endif
    int i=0, j=0;
    //for Vt
    if( Vt!=NULL )
    {
        for( i=0; i<seqlen+2; i++ )
        {
            if( Vt[i]!=NULL )
            {
                for( j=0; j<numstates; j++ )
                {
                    if( Vt[i][j]!=NULL )
                    {
                        delete [ ] Vt[i][j];
                        Vt[i][j]=NULL;
                    }
                }
                delete [ ] Vt[i];
                Vt[i]=NULL;
            }
        }
        delete [ ] Vt;
        Vt=NULL;
    }

    //for Tb
    if( Tb!=NULL )
    {
        for( i=0; i<seqlen+2; i++ )
        {
            if( Tb[i]!=NULL )
            {
                for( j=0; j<numstates; j++ )
                {
                    if( Tb[i][j]!=NULL )
                    {
                        delete [ ] Tb[i][j];
                        Tb[i][j]=NULL;
                    }
                }
                delete [ ] Tb[i];
                Tb[i]=NULL;
            }
        }
        delete [ ] Tb;
        Tb=NULL;
    }
}

//detect memory leak problem in rnatops zbhuang 20090815
void FVTBSearch::free3DTablesAssist()
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::free3DTables\n");
    #endif
	//int i=0, j=0;
	//assistant
    //for V
    if( V!=NULL )
    {
        delete [ ] V;
        V=NULL;
    }

    //for T
    if( T!=NULL )
    {
        delete [ ] T;
        T=NULL;
    }
}
void FVTBSearch::init3DTables(  )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::init3DTables\n");
    #endif
    
    int i=0, j=0, k=0;
    for(i=0; i<seqlen+2; i++)
    {
        for(j=0; j<numstates; j++)
        {
            for(k=0; k<3; k++)
            {
                Vt[i][j][k]=INVLDPROB;
                Tb[i][j][k]=' ';
            }
        }
    }
}

int FVTBSearch::obtainSequence(  char * filename )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::obtainSequence\n");
    #endif
    int i=0, templen=0;
    char ch;
    FILE * fp=NULL;

    //----Open the search file
    fp = fopen(filename, "r");
    if( fp == NULL )
    {
        printf(" Fail to open %s\n", filename);
        return 0;
    }

    //get the size of str: len
    fseek( fp, 0, SEEK_SET);
    while ( !feof(fp) )
    {
        ch=fgetc (fp);
        if( isRNABase(ch) )
            templen++;
    }
    
    //allocate the memory
    searchSeq= new int[templen];
    //get the str
    fseek( fp, 0, SEEK_SET);
    while ( !feof(fp) )
    {
        ch=fgetc (fp);
        if( isRNABase(ch) )
            searchSeq[i++]=chartonum(ch);
    }
    if( i!=templen )
    {
        printf("Wrong with length in target file.\n");
        return 0;
    }
    fclose(fp);
    return templen;
}

int	FVTBSearch::update3DandSeq( int newlen,  int newnumstate )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::update3DandSeq\n");
    #endif
    int change=0;
    
    if( newlen!= seqlen || newnumstate!= numstates )
    {
        free3DTables( );
        change = 1;
    }

    if( newlen!= seqlen )
    {
        seqlen = newlen;
        if( searchSeq!=NULL )
        {
            delete [ ] searchSeq;
            searchSeq=NULL;
        }
        searchSeq=new int[seqlen];
    }
    
    if( newnumstate!= numstates )
    {
        numstates = newnumstate;
    }
    
    if( change )
    {
        if( malloc3DTables( )==0 )
        {
            printf("Fail to allocate memory in malloc3DTables\n");
            return 0;
        }
    }
    return 1;
}

int FVTBSearch::FVTBSearchString( char *instr, int newlen, Loop *curloop )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::FVTBSearchString\n");
    #endif
    
    int i=0, retval=0;
    int *tempstr = NULL;
    
    if( newlen < 0 )
        return 0;
    
    //copy searchSeq
    tempstr = new int[newlen];
    for( i=0; i<newlen; i++ )
         tempstr[i]=chartonum( instr[i] );
  
    retval = FVTBSearchString( tempstr, newlen, curloop );
    
    delete [ ] tempstr;
    return retval;
}

int FVTBSearch::FVTBSearchString( int *instr, int newlen, Loop *curloop )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::FVTBSearchString\n");
    #endif
    
    if( !curloop->isLogProb( ) )
        curloop->logProbforLoop( );
        
    if( curloop->getLoopLength( )<=0 )
        return 0;    
    
    int i=0, status=0;
    int shiftnum=1;
    int newnumstate=0;
    
    newnumstate = curloop->numstates;
    
    //check shifted nts
    status=checkShiftedSeq( instr, newlen, shiftnum);
    if( status==1)
        haverun=1;
    else
        haverun=0;
    
    //1.check the dimension of 3D
    update3DandSeq( newlen, newnumstate );
    
    //2.copy searchSeq
    for( i=0; i<seqlen; i++ )
    {
        searchSeq[i]= instr[i];
    }
    
    //3.do VTB profile
    if( haverun==0 ) {
        VTBimplement( curloop );
    }
    else {
        shiftPosIndex( shiftnum );
        VTBshift( curloop, shiftnum );
    }

    //printResultSeq( );
    return 1;
}

int FVTBSearch::checkShiftedSeq( int *instr, int newlen, int shiftnum)
{
    int i=0;
    int retval=1;
    if( seqlen!=newlen ) {
        retval=0;
    }
    else {
        for( i=shiftnum; i<seqlen; i++ ) {
             if( instr[i-shiftnum] !=searchSeq[i]) {
                 retval=0;
                 break;
             }
        }
    }
    return retval;
}

void FVTBSearch::shiftPosIndex( int shiftnum )
{
    int i=0, j=0, k=0;
    
    //reset   [seq][state][ D=0/I=1/M=2 ]
    for(i=0; i<numstates; i++)
    {
        for( j=0; j<3; j++ )
        {
            V[seqlen+1][i][j]=INVLDPROB;
            T[seqlen+1][i][j]=' ';
            V[0][i][j]=INVLDPROB;
            T[0][i][j]=' ';
        }
    }
    for( i=0; i<seqlen+2; i++ )
    {
        k=(i+shiftnum)%(seqlen+2);
        V[i]=V[k];
        T[i]=T[k];
    }
}

FVTBCandidate * FVTBSearch::getCandHead( )
{
    return  tops;
}

int FVTBSearch::setAllowedInsNum( int num )
{
    allowNumIns = num;
    return num;
}

void FVTBSearch::setHitThreshold( double probT )
{
    newThreshold = probT;
}

void FVTBSearch::printSearchSeq( )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::printSearchSeq\n");
    #endif
    int i=0;
    printf("\n======================= Search Sequence =========================\n");
    printf("Target Sequences of length %d:\n", seqlen );
    printf("   ");
    for(i=0; i<seqlen; i++)
    {
        printf("%c", numtochar(searchSeq[i]));
    }
    printf("\n\n");
}

void FVTBSearch::printResultSeq( )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::printResultSeq\n");
    #endif
    int i=0;
    for( i=0; i<numtops; i++ )
    {
        if( tops[i].hasData( ) )
        {
            tops[i].printAlgnedSeq( );
        }
    }
}

void FVTBSearch::output3Dtables( char stateTable )  //table= M, I, D, A
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::output3Dtables\n");
    #endif
    int i=0, j=0;
    if( stateTable=='M' || stateTable=='A') {
        printf("Match table( x:seq=%d+1; y:states=%d )\n", seqlen, numstates );
        for( i=0; i<numstates; i++ )
        {
            if( i==0 )
                printf("B0:  ");
            else if(i==numstates-1 )
                printf("E%d:  ", i);
            else
                printf("M%d:  ", i);
            for( j=0; j<seqlen+2; j++ )
            {
                if( Vt[j][i][ HMM_M ] > INVLDPROB )
                    printf(" %10.6f(%c) ", Vt[j][i][ HMM_M ], Tb[j][i][ HMM_M ]);
                else
                    printf(" %10s(%c) ", "-INF", Tb[j][i][ HMM_M ]);
            }
            printf("\n");
        }
        printf("\n");
    }
    if( stateTable=='I' || stateTable=='A') {
        printf("Insertion table( x:seq=%d+1; y:states=%d )\n", seqlen, numstates );
        for( i=0; i<numstates; i++ )
        {
            if( i==0 )
                printf("B0:  ");
            else if(i==numstates-1 )
                printf("E%d:  ", i);
            else
                printf("M%d:  ", i);
            for( j=0; j<seqlen+2; j++ )
            {
                if( Vt[j][i][ HMM_I ] > INVLDPROB )
                    printf(" %10.6f(%c) ", Vt[j][i][ HMM_I ], Tb[j][i][ HMM_I ]);
                else
                    printf(" %10s(%c) ", "-INF", Tb[j][i][ HMM_I ]);
            }
            printf("\n");
        }
        printf("\n");
    }
    if( stateTable=='D' || stateTable=='A') {
        printf("Deletion table( x:seq=%d+1; y:states=%d )\n", seqlen, numstates );
        for( i=0; i<numstates; i++ )
        {
            if( i==0 )
                printf("B0:  ");
            else if(i==numstates-1 )
                printf("E%d:  ", i);
            else
                printf("M%d:  ", i);
                
            for( j=0; j<seqlen+2; j++ )
            {
                if( Vt[j][i][ HMM_D ] > INVLDPROB )
                    printf(" %10.6f(%c) ", Vt[j][i][ HMM_D ], Tb[j][i][ HMM_D ]);
                else
                    printf(" %10s(%c) ", "-INF", Tb[j][i][ HMM_D ]);
            }
            printf("\n");
        }
        printf("\n");
    }
}


// Pick the maximum value from a and b
double FVTBSearch::max(double fromI, double fromM )
{
    double temp;
	
    if (fromI > fromM )
        temp = fromI;
    else
        temp = fromM;
    return temp;
}

// Pick the maximum value from a, b and c
double FVTBSearch::max(double fromD, double fromI, double fromM )
{
    double temp;
	
    if ( fromD > fromI )
        temp=fromD;
    else
        temp=fromI;
    
    if ( temp > fromM )
        return temp;
    else
        return fromM;
}

// Give the type of the max one
char FVTBSearch::maxstate (double fromI, double fromM )
{
    char type;
    if (fromI > fromM )
        type = 'I';
    else
        type = 'M';
  	return type;
}

// Give the type of the max one
char FVTBSearch::maxstate (double fromD, double fromI, double fromM )
{
    char type;
    double  temp=0.0;
    if ( fromD > fromI )
    {
        temp = fromD;
        type = 'D';
    }
    else
    {
        temp = fromI;
        type = 'I';
    }
    if ( temp > fromM )
        return type;
    else
        return 'M';
}

//insert Top N number of positions
double FVTBSearch::insertCurPos( int idi, int idj, int idk,  double cprob )
{
    int i=0, j=0;
    for( i=0; i<numtops; i++ )
    {
        if( cprob > tops[i].prob )
        {
            for( j=numtops-1; j>i; j--)
            {
                tops[j].pos[0]=tops[j-1].pos[0];
                tops[j].pos[1]=tops[j-1].pos[1];
                tops[j].pos[2]=tops[j-1].pos[2];
                tops[j].prob=tops[j-1].prob;
            }
            tops[i].pos[0]=idi;
            tops[i].pos[1]=idj;
            tops[i].pos[2]=idk;
            tops[i].prob=cprob;
            break;
        }
    }
    return tops[ numtops-1 ].prob;
}

int FVTBSearch::getNewProbs(  )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::getNewProbs\n");
    #endif
    
    int i=0, k=0, num=0;
    int floor=numstates-1;
    double thres=newThreshold, curprob=INVLDPROB;
    
    //[seq][state][ D=0/I=1/M=2 ]
    i=seqlen+1;  
    for(k=0; k<3; k++ )
    {
        curprob = V[i][floor][k];
        if( curprob >= thres && curprob!=INVLDPROB)
        {
            thres=insertCurPos( i, floor,  k,  curprob );
      	}
    }
    //count the final number of positions
    for(i=0; i<numtops; i++ )
    {
        if( tops[i].pos[0]!=-1)
            num++;
    }
    return num;
}

// get the n candidates with max probability
int FVTBSearch::getTopnMaxProbs(  )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::getTopnMaxProbs\n");
    #endif
    
    int i=0, k=0, num=0;
    int floor=numstates-1;
    double thres=newThreshold, curprob=INVLDPROB;
    
    for(i=0; i<seqlen+2; i++)    //[seq][state][ D=0/I=1/M=2 ]
    {
        for(k=0; k<3; k++ )
        {
            curprob = V[i][floor][k];
            if( curprob >= thres && curprob!=INVLDPROB)
            {
                thres=insertCurPos( i, floor,  k,  curprob );
            }
      	}
    }
    //count the final number of positions
    for(i=0; i<numtops; i++ )
    {
        if( tops[i].pos[0]!=-1)
            num++;
    }
    return num;
}

int FVTBSearch::isInRemSeq( int *seq, int size,  int num)
{
    int i=0, found=0;
    for( i=0; i<size; i++ )
    {
        if( num == seq[i] )
        {
            seq[i]=INVALID;
            found++;
        }
    }
    return found;
}

int FVTBSearch::VTBTraceBack( FVTBCandidate &cad )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::VTBTraceBack\n");
    #endif

    int i=0, j=0, k=0;
    int alld=0, numdel=0;
    int newlen=0; 
    int templen=5*seqlen+numstates;
    int *rem=new int[ templen ];
    int size=0;
    int tbok=1;
    char *align=new char[ templen+1 ];
    
    char curstate=' ';
    
    for( i=0; i<templen; i++)
    {
        align[i]='.';
        rem[i]=INVALID;
    }
    align[i]='\0';
    
    //***** [seq][state][ D=0/I=1/M=2 ]
    i=cad.pos[0];
    j=cad.pos[1];
    k=cad.pos[2];
    
    curstate = T[i][j][k];
    i--;
    j--;
    cad.loc[1]=i-1;
    cad.loc[0]=i-1;
    while( j>0 )
    {
        if( i<1 )
        {
            tbok=0;
            break;
        }
        cad.loc[0]=i-1;
        //printf(" j=%d, i=%d, state=%c\n", j, i, curstate);
        switch( curstate )
        {
            case 'D':  //0
                      curstate = T[i][j][HMM_D];
                      rem[ size++ ]=i-1;
                      j--;
                      break;
            case 'I':  //1
                      curstate = T[i][j][HMM_I];
                      align[i-1] = AlgHmmIns;
                      i--;
                      break;
            case 'M':  //2
                      curstate = T[i][j][HMM_M];
                      align[i-1] = AlgHmmMatch;
                      i--;
                      j--;
                      break;
            case ' ': 
                      break;
            default:
                      printf(" Unknown nonterminal:\"%c\" id: %d!\n", curstate, k);
                      return 0;
        }
    }
    if( tbok )
    {
        //***********process the output string
        newlen=seqlen + size;
        cad.finlen=newlen;
        cad.finstr= new char[newlen+1];
        cad.finseq= new char[newlen+1];
        //for( i=0; i<size; i++ )
        //    printf("%d  ", rem[i] );
        //printf("\n" );
        for( i=0, k=0; i<newlen; i++)
        {
            numdel=isInRemSeq( rem, size,  i-1-alld);
            if( numdel )
            {
                //printf("k=%d, numdel=%d\n", k, numdel);
                for(j=0; j<numdel; j++ )
                {
                    cad.finseq[i+j]='-';
                    cad.finstr[i+j]=AlgHmmDel;
                }
                i+=numdel-1;
                alld+=numdel;
            }
            else
            {
                cad.finseq[i]= numtochar(searchSeq[k]);
                cad.finstr[i]= align[k];
                k++;
            }
        }
        cad.finseq[i]='\0';
        cad.finstr[i]='\0';
    }
    if( rem!=NULL )
        delete [ ] rem;
    if( align!=NULL )
        delete [ ] align;
    return tbok;
}

void FVTBSearch::VTBimplement( Loop *curLoop )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::VTBimplement\n");
    #endif
    int i=0;
    int nt=0, st=0;
    int residue=-1;
    double	fromM=0, fromI=0, fromD=0;
    
    HMMGram *preM = NULL,	*preI = NULL,	*preD = NULL;
    HMMGram *curM = NULL,	*curI = NULL,	*curD = NULL;    
    HMMGram *beg = curLoop->pLoopGram;
    
    //****** initialization stage ****** [ seq ][ state ][ D=0/I=1/M=2 ]
    initTopCands( );
    init3DTables( );
    Vt[0][0][ HMM_M ]=0;
    
//    // initialize the first column in D
//    Vt[0][1][ HMM_D ] = Vt[0][0][ HMM_M ] + beg->trprob[ HMM_D ];
//    Tb[0][1][ HMM_D ] = 'M';
//    curD=beg->nextD;
//    for (i=2; i<numstates; i++)//  2-----numstates-1
//    {
//       	Vt[0][i][ HMM_D ] = Vt[0][i-1][ HMM_D ] + curD->trprob[ HMM_D ];
//       	Tb[0][i][ HMM_D ] = 'D';
//        preD = curD;
//       	curD=curD->nextD;
//    }
    
    if( seqlen == 0)
    {
        //M  seqlen+1,  numstates
        if( preD!=NULL )
        {
            fromD = Vt[seqlen][numstates-2][ HMM_D ] + preD->trprob[ HMM_M ];
            Vt[ seqlen+1 ][numstates-1][ HMM_M ] = fromD;
            Tb[ seqlen+1 ][numstates-1][ HMM_M ] = 'D';
        }
        else
        {
            fromM = Vt[seqlen][numstates-2][ HMM_M ] + beg->trprob[ HMM_M ];
            Vt[ seqlen+1 ][numstates-1][ HMM_M ] = fromM;
            Tb[ seqlen+1 ][numstates-1][ HMM_M ] = 'M';
        }
    }
    else   //seqlen>0
    {
        // initialize the first row in I    ****** [ seq ][ state ][ D=0/I=1/M=2 ]
        curI=beg->nextI;

        Vt[0][0][ HMM_M ] =0;
        Tb[0][0][ HMM_M ] = ' ';
        for (i=1; i<seqlen+1; i++)
        {
           	Vt[i][0][ HMM_M ] = 0;
           	Tb[i][0][ HMM_M ] = ' ';
        }

        //****** main loop stage ****** [seq][state][ D=0/I=1/M=2 ]
        if( beg->nextM->state != 'E' )
        {
            for ( nt=1; nt<seqlen+1; nt++ )
            {
                preM = beg;
                preI = beg->nextI;
                curM = beg->nextM;
                curD = beg->nextD;
                curI = curM->nextI;
            
                residue = searchSeq[nt-1]-1;
                st=1;
                // M--match
                fromM = Vt[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
    	        fromI = Vt[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];             
    	        Vt[nt][st][ HMM_M ] = max( fromI, fromM ) + curM->emission[ residue ];
    	        Tb[nt][st][ HMM_M ] = maxstate( fromI, fromM );

                // I--insertion
                fromM = Vt[nt-1][st][ HMM_M ] + curM->trprob[ HMM_I ];
                fromI = Vt[nt-1][st][ HMM_I ] + curI->trprob[ HMM_I ];
                fromD = Vt[nt-1][st][ HMM_D ] + curD->trprob[ HMM_I ];
                Vt[nt][st][ HMM_I ] = max( fromD, fromI, fromM ) + curI->emission[ residue ];
                Tb[nt][st][ HMM_I ] = maxstate( fromD, fromI, fromM );
        
                // D--deletion
                fromM = Vt[nt][st-1][ HMM_M ] + preM->trprob[ HMM_D ];
                fromI = Vt[nt][st-1][ HMM_I ] + preI->trprob[ HMM_D ];
                Vt[nt][st][ HMM_D ] = max( fromD, fromI, fromM );
                Tb[nt][st][ HMM_D ] = maxstate( fromI, fromM );

    	        preM = curM;
    	        preD = curD;
    	        preI = curI;
    	        curM = preM->nextM;
    	        curD = preM->nextD;
                curI = curM->nextI;
      	
                for ( st=2; st<numstates-1;  st++ )
                {
                    // M--match
                    fromM = Vt[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
                    fromI = Vt[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];
                    fromD = Vt[nt-1][st-1][ HMM_D ] + preD->trprob[ HMM_M ];
                    Vt[nt][st][ HMM_M ] = max(fromD, fromI, fromM ) + curM->emission[ residue ];
                    Tb[nt][st][ HMM_M ] = maxstate(fromD, fromI, fromM);

                    // I--insertion
                    fromM = Vt[nt-1][st][ HMM_M ] + curM->trprob[ HMM_I ];
                    fromI = Vt[nt-1][st][ HMM_I ] + curI->trprob[ HMM_I ];
                    fromD = Vt[nt-1][st][ HMM_D ] + curD->trprob[ HMM_I ];
                    Vt[nt][st][ HMM_I ] = max(fromD, fromI, fromM ) + curI->emission[ residue ];
                    Tb[nt][st][ HMM_I ] = maxstate(fromD, fromI, fromM);
        
                    // D--deletion
                    fromM = Vt[nt][st-1][ HMM_M ] + preM->trprob[ HMM_D ];
                    fromI = Vt[nt][st-1][ HMM_I ] + preI->trprob[ HMM_D ];
                    fromD = Vt[nt][st-1][ HMM_D ] + preD->trprob[ HMM_D ];
                    Vt[nt][st][ HMM_D ] = max(fromD, fromI, fromM );
                    Tb[nt][st][ HMM_D ] = maxstate(fromD, fromI, fromM);

                    preM = curM;
                    preD = curD;
                    preI = curI;
                   
                    curM = curM->nextM;
                    curD = preM->nextD;
                    curI = curM->nextI;
                }
        
                st=numstates-1;
                // E--end state
                fromM = Vt[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
    	        fromI = Vt[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];
    	        fromD = Vt[nt-1][st-1][ HMM_D ] + preD->trprob[ HMM_M ];
    	        Vt[nt][st][ HMM_M ] = max(fromD, fromI, fromM );
                Tb[nt][st][ HMM_M ] = maxstate(fromD, fromI, fromM);
            }
            //last column & last row
            //M  seqlen+1,  numstates,
            fromM = Vt[seqlen][numstates-2][ HMM_M ] + preM->trprob[ HMM_M ];
            fromI = Vt[seqlen][numstates-2][ HMM_I ] + preI->trprob[ HMM_M ];
            fromD = Vt[seqlen][numstates-2][ HMM_D ] + preD->trprob[ HMM_M ];
            Vt[ seqlen+1 ][numstates-1][ HMM_M ] = max(fromD, fromI, fromM );
            Tb[ seqlen+1 ][numstates-1][ HMM_M ] = maxstate(fromD, fromI, fromM);
        }
        else
        {
            preM = beg;
            preI = beg->nextI;
            // curM = beg->nextM;
            // E--end state
            fromM = Vt[seqlen-1][0][ HMM_M ] + preM->trprob[ HMM_M ];
            fromI = Vt[seqlen-1][0][ HMM_I ] + preI->trprob[ HMM_M ];
            Vt[seqlen][1][ HMM_M ] = max( fromI, fromM );
            Tb[seqlen][1][ HMM_M ] = maxstate( fromI, fromM);
        
            //last column
            //M  seqlen+1,  numstates,
            fromM = Vt[seqlen][numstates-2][ HMM_M ] + preM->trprob[ HMM_M ];
            fromI = Vt[seqlen][numstates-2][ HMM_I ] + preI->trprob[ HMM_M ];
            Vt[ seqlen+1 ][numstates-1][ HMM_M ] = max(fromI, fromM );
            Tb[ seqlen+1 ][numstates-1][ HMM_M ] = maxstate(fromI, fromM);
        }
    }
    obtainCandidate(  );
}

void FVTBSearch::VTBshift( Loop *curLoop, int shiftnum )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::VTBshift\n");
    #endif
    int nt=0, st=0;
    int residue=-1;
    double	fromM=0, fromI=0, fromD=0;
    
    HMMGram *preM = NULL,	*preI = NULL,	*preD = NULL;
    HMMGram *curM = NULL,	*curI = NULL,	*curD = NULL;    
    HMMGram *beg = curLoop->pLoopGram;
    
    //****** initialization stage ****** [ seq ][ state ][ D=0/I=1/M=2 ]
    initTopCands( );
    V[0][0][ HMM_M ]=0;
    
    if( seqlen == 0)
    {
        printf("The length of sequence is 0, please check it\n" );
        exit(0);
    }
    else   //seqlen>0
    {
        // initialize the first row in I    ****** [ seq ][ state ][ D=0/I=1/M=2 ]
        curI=beg->nextI;

        V[seqlen][0][ HMM_M ] = 0;
        T[seqlen][0][ HMM_M ] = ' ';

        //****** main loop stage ****** [seq][state][ D=0/I=1/M=2 ]
        if( beg->nextM->state != 'E' )
        {
            for ( nt=seqlen+1-shiftnum; nt<seqlen+1; nt++ )
            {
                preM = beg;
                preI = beg->nextI;
                curM = beg->nextM;
                curD = beg->nextD;
                curI = curM->nextI;

                residue = searchSeq[nt-1]-1;
                st=1;
                // M--match
                fromM = V[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
    	        fromI = V[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];             
    	        V[nt][st][ HMM_M ] = max( fromI, fromM ) + curM->emission[ residue ];
    	        T[nt][st][ HMM_M ] = maxstate( fromI, fromM );

                // I--insertion
                fromM = V[nt-1][st][ HMM_M ] + curM->trprob[ HMM_I ];
                fromI = V[nt-1][st][ HMM_I ] + curI->trprob[ HMM_I ];
                fromD = V[nt-1][st][ HMM_D ] + curD->trprob[ HMM_I ];
                V[nt][st][ HMM_I ] = max( fromD, fromI, fromM ) + curI->emission[ residue ];
                T[nt][st][ HMM_I ] = maxstate( fromD, fromI, fromM );
        
                // D--deletion
                fromM = V[nt][st-1][ HMM_M ] + preM->trprob[ HMM_D ];
                fromI = V[nt][st-1][ HMM_I ] + preI->trprob[ HMM_D ];
                V[nt][st][ HMM_D ] = max( fromD, fromI, fromM );
                T[nt][st][ HMM_D ] = maxstate( fromI, fromM );

    	        preM = curM;
    	        preD = curD;
    	        preI = curI;
    	        curM = preM->nextM;
    	        curD = preM->nextD;
                curI = curM->nextI;
      	
                for ( st=2; st<numstates-1;  st++ )
                {
                    // M--match
                    fromM = V[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
                    fromI = V[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];
                    fromD = V[nt-1][st-1][ HMM_D ] + preD->trprob[ HMM_M ];
                    V[nt][st][ HMM_M ] = max(fromD, fromI, fromM ) + curM->emission[ residue ];
                    T[nt][st][ HMM_M ] = maxstate(fromD, fromI, fromM);

                    // I--insertion
                    fromM = V[nt-1][st][ HMM_M ] + curM->trprob[ HMM_I ];
                    fromI = V[nt-1][st][ HMM_I ] + curI->trprob[ HMM_I ];
                    fromD = V[nt-1][st][ HMM_D ] + curD->trprob[ HMM_I ];
                    V[nt][st][ HMM_I ] = max(fromD, fromI, fromM ) + curI->emission[ residue ];
                    T[nt][st][ HMM_I ] = maxstate(fromD, fromI, fromM);
        
                    // D--deletion
                    fromM = V[nt][st-1][ HMM_M ] + preM->trprob[ HMM_D ];
                    fromI = V[nt][st-1][ HMM_I ] + preI->trprob[ HMM_D ];
                    fromD = V[nt][st-1][ HMM_D ] + preD->trprob[ HMM_D ];
                    V[nt][st][ HMM_D ] = max(fromD, fromI, fromM );
                    T[nt][st][ HMM_D ] = maxstate(fromD, fromI, fromM);

                    preM = curM;
                    preD = curD;
                    preI = curI;
                   
                    curM = curM->nextM;
                    curD = preM->nextD;
                    curI = curM->nextI;
                }
        
                st=numstates-1;
                // E--end state
                fromM = V[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
    	        fromI = V[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];
    	        fromD = V[nt-1][st-1][ HMM_D ] + preD->trprob[ HMM_M ];
    	        V[nt][st][ HMM_M ] = max(fromD, fromI, fromM );
                T[nt][st][ HMM_M ] = maxstate(fromD, fromI, fromM);
            }
            //last column & last row
            //M  seqlen+1,  numstates,
            fromM = V[seqlen][numstates-2][ HMM_M ] + preM->trprob[ HMM_M ];
            fromI = V[seqlen][numstates-2][ HMM_I ] + preI->trprob[ HMM_M ];
            fromD = V[seqlen][numstates-2][ HMM_D ] + preD->trprob[ HMM_M ];
            V[ seqlen+1 ][numstates-1][ HMM_M ] = max(fromD, fromI, fromM );
            T[ seqlen+1 ][numstates-1][ HMM_M ] = maxstate(fromD, fromI, fromM);
        }
        else
        {
            preM = beg;
            preI = beg->nextI;
            // curM = beg->nextM;
            // E--end state
            fromM = V[seqlen-1][0][ HMM_M ] + preM->trprob[ HMM_M ];
            fromI = V[seqlen-1][0][ HMM_I ] + preI->trprob[ HMM_M ];
            V[seqlen][1][ HMM_M ] = max( fromI, fromM );
            T[seqlen][1][ HMM_M ] = maxstate( fromI, fromM);
        
            //last column
            //M  seqlen+1,  numstates,
            fromM = V[seqlen][numstates-2][ HMM_M ] + preM->trprob[ HMM_M ];
            fromI = V[seqlen][numstates-2][ HMM_I ] + preI->trprob[ HMM_M ];
            V[ seqlen+1 ][numstates-1][ HMM_M ] = max(fromI, fromM );
            T[ seqlen+1 ][numstates-1][ HMM_M ] = maxstate(fromI, fromM);
        }
    }
    obtainCandidate(  );
}

int FVTBSearch::obtainCandidate(  )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::obtainCandidate\n");
    #endif
    
    int i=0, nowtops=0, status=0;
    
    //output3Dtables( 'M' );
    //****** find the maximal prob ****** [seq][state][ D=0/I=1/M=2 ]
    if( haverun==0 )
        nowtops=getTopnMaxProbs(  );  //get maximum prob
    else
        nowtops=getNewProbs(  );     //global alignment
    for(i=0; i<nowtops; i++)
    {
        tops[i].cadno=i+1;
        status = VTBTraceBack( tops[ i ] );
        if( status==0 )
            tops[i].initFVTBCandidates( );
    }
    return 1;
}

HmmFilterResult * FVTBSearch::getSearchResult(  )
{
    return hfr;
}

int FVTBSearch::SearchTopNMerge( Loop *curloop, char *instr, int subwsize, int topN )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::SearchTopNMerge\n");
    #endif
    if( hfr!=NULL )
        delete hfr;
    hfr = new HmmFilterResult( topN );
    
    int i=0, j=0, wholelen=0;
    int curHitNum=0, minHitIdx=0;
    double probv=0.0, minscore=10000;
    FVTBCandidate *hits = hfr->hits;
    wholelen=strlen( instr );
    if( wholelen <= 0 )
        return 0;
    char *segmStr = new char[ subwsize+1 ];
    FVTBCandidate * fcd=NULL;
    //printf("wholelen=%d,  subwsize=%d,   topN:%d\n", wholelen, subwsize,  topN );

    for( i=0; i<=wholelen-subwsize; i++ )
    {
        //printf("   ... subidx: %d ...\n", i);
        //1. search sequence
        for( j=0; j<subwsize; j++ )
            segmStr[j]=instr[i+j];
        segmStr[j]='\0';
        
        //2. search sequence
        fcd = SearchOneMerge( curloop, segmStr );
        if( fcd == NULL ) {
            continue;
        }
        else {
            fcd->setAbsPosition( i );
            probv = fcd->getScore( );
            //printf("score: %f\n", probv);
            if( curHitNum<topN ){
                hits[curHitNum] = *fcd;
                if( probv < minscore ) {
                    minscore = probv;
                    minHitIdx = curHitNum;
                }
                curHitNum++;
            }
            else {
                if( probv > minscore ) {
                    hits[minHitIdx] = *fcd;
                    minscore = 10000;
                    for( j=0; j<topN; j++ ) {
                        probv = hits[j].getScore( );
                        if( probv < minscore ) {
                            minscore = probv;
                            minHitIdx = j;
                        }
                    }
                }
            }
        }
        if( fcd!=NULL )
            delete [ ] fcd;
    }
    delete [ ] segmStr;
    hfr->setActualNum( curHitNum );
    //hfr->printCands(  );
    hfr->sortHits( );
    return curHitNum;
}

FVTBCandidate * FVTBSearch::SearchOneMerge( Loop *curloop, char *instr )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::SearchOneMerge\n");
    #endif
    int i=0, curlen=0;
    int *subseq=NULL;
    FVTBCandidate *curhit=NULL;
    int reportHit=0;
    
    //process the passing parameters
    curlen=strlen( instr );
    if( curlen <= 0)
        curlen=0;
    subseq= new int[ curlen ];
    for( i=0; i<curlen; i++ )
        subseq[i]= chartonum( instr[i] );
    
    //Loop model is empty
    if( curloop->getLoopLength( )<=0 )
    {
        printf("Null loop model can not be applied as filter\n");
        exit(0);
    }
    else  //Loop model is not empty
    {
       //printf("  ~~~SearchOneMerge\n");
       FVTBSearchString( subseq, curlen, curloop );
        
        if( haverun==0 ) //the first frame
        {
            reportHit=0;
            prehit = tops[0];
            preBegLoc = tops[0].loc[0];
        }
        else
        {
            if( tops->hasData( ) )
            {
                if(tops[0].loc[0]==preBegLoc-1)
                {
                    if( tops[0].prob >prehit.prob )
                        prehit = tops[0];
                    preBegLoc = tops[0].loc[0];
                }
                else
                {
                    reportHit=1;
                }
            }
        }               
        if( reportHit )
        {   
            if( prehit.hasData( ) ) {
                curhit = new FVTBCandidate [1];
                curhit[0] = prehit;
            }
            else {
                curhit = NULL;
            }
            //update prehit
            prehit = tops[0];
            preBegLoc = tops[0].loc[0];
        }
        else
        {
            curhit = NULL;
        }
    }
    delete [ ] subseq;
    return curhit;
}

FVTBCandidate * FVTBSearch::SearchOneNoMerge( Loop *curloop, char *instr )
{
    #ifdef DEBUG_DEV
        printf("FVTBSearch::SearchOneNoMerge\n");
    #endif
    int i=0, curlen=0;
    int *subseq=NULL;
    FVTBCandidate *curhit=NULL;
    
    //process the passing parameters
    curlen=strlen( instr );
    if( curlen <= 0)
        curlen=0;

    subseq= new int[ curlen ];
    for( i=0; i<curlen; i++ )
        subseq[i]= chartonum( instr[i] );
    
    //Loop model is empty
    if( curloop->getLoopLength( )<=0 )
    {
        printf("Null loop model can not be applied as filter\n");
        exit(0);
    }
    else  //Loop model is not empty
    {
        FVTBSearchString( subseq, curlen, curloop );
        if( tops[0].hasData( ) )
        {
            curhit = new FVTBCandidate [1];
            curhit[0] = tops[0];
        }
        else
        {
            curhit = NULL;
        }
    }
    delete [ ] subseq;
    return curhit;
}

