#pragma warning(disable: 4786)

#include "HMMBuilder.h"

//********************************************************************************
//*******************    implementation of HMMGram class         ****************
//********************************************************************************
HMMGram::HMMGram( )
{
    #ifdef DEBUG_DEV
        printf("HMMGram::HMMGram\n");
    #endif
    int i=0;
    
    state=' ';  //I:insertion; M:Matching; D:Deletion; Start:B; End:E
    ID=-1;
    SubID=-1;
    UID=-1;
    for(i=0; i<4; i++ )
        emission[i]=ZEROPROB;
    for(i=0; i<3; i++ )
        trprob[i]=ZEROPROB;
    nextD=NULL;
    nextI=NULL;
    nextM=NULL;
}

HMMGram::HMMGram( char settype,  int setid,  int setuid,  int subid  )
{
    #ifdef DEBUG_DEV
        printf("HMMGram::HMMGram( parameter )\n");
    #endif
    
    int i=0;
    
    state=settype;  //I:insertion; M:Matching; D:Deletion; Start:B; End:E
    ID=setid;
    SubID=subid;
    UID=setuid;
    for(i=0; i<4; i++ )
        emission[i]=ZEROPROB;
    for(i=0; i<3; i++ )
        trprob[i]=ZEROPROB;
    nextD=NULL;
    nextI=NULL;
    nextM=NULL;
}

HMMGram::~HMMGram( )
{
    nextD=NULL;
    nextI=NULL;
    nextM=NULL;
}

void HMMGram::logGrammarProb( double logfreq[5] )
{
    int i=0;
    for(i=0; i<4; i++ )
    {
        if( emission[i] != ZEROPROB )
        {
            if( LOGODDS )
                emission[i]=log( emission[i] )-logfreq[i+1];
            else
                emission[i]=log( emission[i] );
        }
        else
        {
            emission[i]=INVLDPROB;
        }
    }
    for(i=0; i<3; i++ )
    {
        if( trprob[i] != ZEROPROB )
            trprob[i]=log( trprob[i] );
        else
            trprob[i]=INVLDPROB;
    }
}

int HMMGram::countNumRules(  )
{
    #ifdef DEBUG_DEV
        printf("HMMGram::countNumRules\n");
    #endif
    int num=0;
    if( nextD != NULL )
    	 num++;
    if( nextI != NULL )
    	 num++;
    if( nextM != NULL )
    	 num++;
    return num;
}

//********************************************************************************
//*******************    implementation of Loop class         ****************
//********************************************************************************
Loop::Loop( )
{
    #ifdef DEBUG_DEV
       printf("Loop::Loop\n");
    #endif
    int i=0;
    
    UID=-1;
    Pos[0]=-1;
    Pos[1]=-1;
    looplen=0;
    numseqs=0;
    state=NULL;
    numstates=0;
    loopseqs=NULL;
    pLoopGram=NULL;
    logProb=0;
    for( i=0; i<2; i++ )
    {
        offsetE[i]=0.0;
        offsetSD[i]=0.0;
    }
    lenE=0.0;
    lenSD=0.0;
    for( i=0; i<5; i++ )
    {
        bgBaseLogHz[i]=0.0;
    }
}

Loop::Loop( int seqnum, int lplen )
{
    #ifdef DEBUG_DEV
       printf("Loop::Loop( parameter )\n");
    #endif
    int i=0;
    
    UID=-1;
    Pos[0]=-1;
    Pos[1]=-1;
    looplen=lplen;
    numseqs=seqnum;
    state=NULL;
    numstates=0;
    loopseqs=NULL;
    pLoopGram=NULL;
    for( i=0; i<2; i++ )
    {
        offsetE[i]=0.0;
        offsetSD[i]=0.0;
    }
    lenE=0.0;
    lenSD=0.0;
    logProb=0;
    for( i=0; i<5; i++ )
    {
        bgBaseLogHz[i]=0.0;
    }
}

Loop::~Loop( )
{
    #ifdef DEBUG_DEV
       printf("Loop::~Loop\n");
    #endif
    int i=0;
    
    freeHMMGrams( );
    
    if( state!=NULL )
    {
        delete [ ] state;
        state=NULL;
    }
    if( loopseqs != NULL )
    {
        for(i=0; i<numseqs; i++ )
        {
            if( loopseqs[i]!=NULL )
            {
                delete [ ] loopseqs[i];
                loopseqs[i]=NULL;
            }
        }
        delete [ ] loopseqs;
        loopseqs=NULL;
    }
}

int Loop::getLoopLength( )
{
    #ifdef DEBUG_DEV
       printf("Loop::getLoopLength\n");
    #endif
     return looplen;
}

int Loop::getLoopSpan( double coeff, int ismax )
{
    int retval=0;
    double sum=0.0;
    if( ismax )
        sum=lenE + coeff*lenSD;
    else
        sum=lenE - coeff*lenSD;
    sum = sum>0 ? sum : 0;
    retval=(int)sum;
    return retval;
}

int Loop::getLeftOffset( double coeff, int ismax )
{
    int retval=0;
    if( ismax )
        retval = (int)( offsetE[0] + coeff*offsetSD[0]);
    else
        retval = (int)( offsetE[0] - coeff*offsetSD[0]);
    retval = retval>0 ? retval : 0;
    return retval;    
}

int Loop::getRightOffset( double coeff, int ismax )
{
    int retval=0;
    if( ismax )
        retval = (int)( offsetE[1] + coeff*offsetSD[1]);
    else
        retval = (int)( offsetE[1] - coeff*offsetSD[1]);
    retval = retval>0 ? retval : 0;
    return retval;    
}

void Loop::freeHMMGrams( )
{
    #ifdef DEBUG_DEV
       printf("Loop::freeHMMGrams\n");
    #endif
    
    HMMGram *curM=NULL, *curI=NULL;
    HMMGram *nexD=NULL, *nexM=NULL;
    
    curM=pLoopGram;
    while( curM!=NULL)
    {
        nexM=curM->nextM;
        nexD=curM->nextD;
        curI=curM->nextI;
        
        if( nexD!=NULL )
        {
            delete nexD;
            nexD=NULL;
        }
        if( curI!=NULL )
        {
            delete curI;
            curI=NULL;
        }
        delete curM;
        curM=nexM;
    }
}

int Loop::isLogProb( )
{
    return logProb;
}

int Loop::logProbforLoop( )
{
    #ifdef DEBUG_DEV
       printf("Loop::logProbforLoop\n");
    #endif
    
    HMMGram *curM=NULL, *curI=NULL;
    HMMGram *nexD=NULL;
    
    if( isLogProb( ) )
    {
        return 1;
    }
    
    curM=pLoopGram;
    while( curM!=NULL)
    {
        curM->logGrammarProb( bgBaseLogHz );
        nexD=curM->nextD;
        if( nexD!=NULL )
        {
            nexD->logGrammarProb( bgBaseLogHz );
        }
        curI=curM->nextI;
        if( curI!=NULL )
        {
            curI->logGrammarProb( bgBaseLogHz );
        }
        curM=curM->nextM;
    }
    
    logProb=1;
    return 1;
}


void Loop::printLoopGrams( char direct )
{
    #ifdef DEBUG_DEV
       printf("Loop::printLoopGrams\n");
    #endif
    char sym[2];
    HMMGram *curM=NULL, *nexM=NULL;
    HMMGram *curD=NULL, *nexD=NULL;
    HMMGram *curI=NULL;
    
    curM = pLoopGram;
    nexM = pLoopGram->nextM;
    while( nexM!=NULL )
    {
        curI=curM->nextI;
        nexD=curM->nextD;
        
        if( direct=='L'||direct==' ')
        {
            if( curM->state=='B')
                sym[0]=' ';
            else
                sym[0]='x';
            sym[1]=' ';
        }
        else
        {
            sym[0]=' ';
            if( curM->state=='B')
                sym[1]=' ';
            else
                sym[1]='x';
        }
        
        //M-->M
        printf("--------------------- %d -> %d ---------------------\n", curM->ID, nexM->ID );
         printf( "[%c%d]->%c[%c%d]%c    P=%f\n",  
                           curM->state,   curM->ID,
                 sym[0],   nexM->state,   nexM->ID,   sym[1],    curM->trprob[ HMM_M ] );
        if( curM->state!='B' )
            printSingleEmProb( curM->emission );
        
        //M-->I
        if(  curI!=NULL )
        {
             printf( "[%c%d]->%c[%c%d]%c    P=%f\n", 
                           curM->state,    curM->ID,
                 sym[0],   curI->state,    curI->ID,   sym[1],    curM->trprob[ HMM_I ] );
            if( curM->state!='B' )
                printSingleEmProb( curM->emission );
        }
        //M-->I
        if(  nexD!=NULL )
        {
             printf( "[%c%d]->%c[%c%d]%c    P=%f\n",
                          curM->state,    curM->ID,
                 sym[0],  nexD->state,    nexD->ID,   sym[1],   curM->trprob[ HMM_D ] );
            if( curM->state!='B' )
                printSingleEmProb( curM->emission );
        }
        
        if( curI!=NULL )
        {
            //  I-->M
            printf( "[%c%d]->%c[%c%d]%c     P=%f\n",
                              curI->state,     curI->ID, 
                    sym[0],   nexM->state,     nexM->ID,   sym[1],  curI->trprob[ HMM_M ] );
            printSingleEmProb( curI->emission );
            
            //  I-->I
            if(  curI!=NULL )
            {
                 printf( "[%c%d]->%c[%c%d]%c     P=%f\n",
                              curI->state,     curI->ID, 
                    sym[0],   curI->state,     curI->ID,   sym[1],  curI->trprob[ HMM_I ] );
                printSingleEmProb( curI->emission );
            }
	
            //  I-->D
            if( nexD!=NULL)
            {
                printf( "[%c%d]->%c[%c%d]%c     P=%f\n",
                              curI->state,     curI->ID, 
                    sym[0],   nexD->state,     nexD->ID,   sym[1],  curI->trprob[ HMM_D ] );
                printSingleEmProb( curI->emission );
            }
        }

        if( curD!=NULL )
        {
            //  D-->M
            printf( "[%c%d]-> [%c%d]%c     P=%f\n",
                              curD->state,   curD->ID, 
                    nexM->state,   nexM->ID,  sym[1],   curD->trprob[ HMM_M ] );
            
            //  D-->I
            if( curI!=NULL )
            {
                 printf( "[%c%d]-> [%c%d]%c     P=%f\n",
                              curD->state,     curD->ID, 
                    curI->state,     curI->ID,  sym[1],   curD->trprob[ HMM_I ] );
            }
	
            //  D-->D
            if( nexD!=NULL)
            {
                printf( "[%c%d]-> [%c%d]%c     P=%f\n",
                              curD->state,     curD->ID, 
                    nexD->state,     nexD->ID,  sym[1],   curD->trprob[ HMM_D ] );
            }
        }
        curD=curM->nextD;
        curM=nexM;
        nexM=curM->nextM;
    }
}

void Loop::printHMMGrams( char direct )
{
    if( isLogProb( ) )
        printf("=============== Base-e Log Probability ===============\n");
    else
        printf("===================   Probability   =================\n");
    printLoopGrams( direct );
}

//********************************************************************************
//********************    implementation of Bulge class        *******************
//********************************************************************************
Bulge::Bulge( )
{
    #ifdef DEBUG_DEV
       printf("Bulge::Bulge\n");
    #endif
    nextLp=NULL;
    nextI=NULL; 
    nextS=NULL;
    trprob[0]=0.0;
    trprob[1]=0.0;
    
    nonterm=NULL;
    numntm=0;
    prod=NULL;
    numprod=0;
}

Bulge::Bulge( int seqnum, int lplen ):Loop( seqnum, lplen )
{
    #ifdef DEBUG_DEV
       printf("Bulge::Bulge( parameter )\n");
    #endif

    nextLp=NULL;
    nextI=NULL; 
    nextS=NULL;
    trprob[0]=0.0;
    trprob[1]=0.0;
    
    nonterm=NULL;
    numntm=0;
    prod=NULL;
    numprod=0;
}

Bulge::~Bulge( )
{
    #ifdef DEBUG_DEV
       printf("Bulge::~Bulge\n");
    #endif

    if( nonterm!=NULL )
    {
        delete [ ] nonterm;
        nonterm=NULL;
    }
    
    if( prod!=NULL )
    {
        delete [ ] prod;
        prod=NULL;
    }
}

int Bulge::logProbforBulge( double stemBaseHz[5] )
{
    #ifdef DEBUG_DEV
       printf("Bulge::logProbforBulge\n");
    #endif
    int i=0;
    
    if( LOGODDS )
    {
        for( i=0; i<5; i++ )
            bgBaseLogHz[i]=stemBaseHz[i];
    }
    logProbforLoop( );
    
    if( trprob[0]!= 0 )
        trprob[0]=log( trprob[0]);
    else
        trprob[0]=INVLDPROB;
        
    if( trprob[1]!= 0 )
        trprob[1]=log( trprob[1]);
    else
        trprob[1]=INVLDPROB;
    
    logProb=1;
    return 1;
}

int Bulge::hasLinkVState( )
{
    #ifdef DEBUG_DEV
       printf("Bulge::hasLinkVState\n");
    #endif
    if( nextS!=NULL )
        return 1;
    else
        return 0;    
}

int Bulge::countNumAllRules( )
{
    #ifdef DEBUG_DEV
       printf("Bulge::countNumAllRules\n");
    #endif
    int num=0;
    HMMGram *curM=NULL, *nxD=NULL, *curI=NULL;
    
    curM=pLoopGram;
    while( curM!=NULL)
    {
        num += curM->countNumRules(  );
        nxD=curM->nextD;
        if( nxD!=NULL )
        {
            num += nxD->countNumRules(  );
        }
        curI=curM->nextI;
        if( curI!=NULL )
        {
            num += curI->countNumRules(  );
        }
        curM=curM->nextM;
    }
    if( nextI != NULL)
        num++;
    if( nextS != NULL)
        num++;
    if( nextLp != NULL)
        num++;
    return num;
}

int Bulge::countNumAllNtms( )
{
    #ifdef DEBUG_DEV
       printf("Bulge::countNumAllNtms\n");
    #endif
    int num=0;
    HMMGram *curM=NULL, *nexD=NULL, *curI=NULL;
    
    curM=pLoopGram;
    while( curM!=NULL)
    {
        num++;
        
        nexD=curM->nextD;
        if( nexD!=NULL )
            num++;

        curI=curM->nextI;
        if( curI!=NULL )
            num++;

        curM=curM->nextM;
    }
    return num;
}

char Bulge::getFirstGramState(   )
{
    return pLoopGram->state;
}

int Bulge::getFirstGramUID(   )
{
    return pLoopGram->UID;
}

char Bulge::getLastGramState(  )
{
    HMMGram *curM=NULL, *preM=NULL;
    
    curM=pLoopGram;
    preM=curM;
    while( curM!=NULL)
    {
        preM=curM;
        curM=curM->nextM;
    }
    return preM->state;
}

int Bulge::getLastGramUID(  )
{
    HMMGram *curM=NULL, *preM=NULL;
    
    curM=pLoopGram;
    preM=curM;
    while( curM!=NULL)
    {
        preM=curM;
        curM=curM->nextM;
    }
    return preM->UID;
}

void Bulge::revLoopSeqs(  )
{
    #ifdef DEBUG_DEV
        printf("Bulge::revLoopSeqs\n");
    #endif

    int i=0, j=0;
    int *temprev=NULL;
    
    temprev=new int[looplen];
    //revert the training seqs
    for( i=0; i<numseqs; i++ )
    {
        for( j=0; j<looplen; j++ )
        {
            temprev[j]=loopseqs[i][looplen-j-1];
        }
        for( j=0; j<looplen; j++ )
        {
            loopseqs[i][j]=temprev[j];
        } 
    }
    //revert the state seq
    for( j=0; j<looplen; j++ )
    {
         temprev[j]=state[looplen-j-1];
    }
    for( j=0; j<looplen; j++ )
    {
        state[j]=temprev[j];
    }
    delete [ ] temprev;
}

void Bulge::arrayNonterminal( NontermProp *pntms )
{
    #ifdef DEBUG_DEV
        printf("Bulge::arrayNonterminal\n");
    #endif
    
    int idx=0;
    HMMGram *curM=NULL, *nexD=NULL, *curI=NULL;
    
    curM=pLoopGram;
    while( curM!=NULL)
    {
        idx=curM->UID;
        pntms[idx].uid=idx;
        pntms[idx].state=curM->state;
        pntms[idx].id=curM->ID;
        
        nexD=curM->nextD;
        if( nexD!=NULL )
        {
            idx=nexD->UID;
            pntms[idx].uid=idx;
            pntms[idx].state=nexD->state;
            pntms[idx].id=nexD->ID;
        }
        
        curI=curM->nextI;
        if( curI!=NULL )
        {
            idx=curI->UID;
            pntms[idx].uid=idx;
            pntms[idx].state=curI->state;
            pntms[idx].id=curI->ID;
        }
        curM=curM->nextM;
    }
}

void Bulge::genSingleHMMProd( HMMGram *curHmm,  HMMGram *nextHmm, 
                            GramRule * prod,  int idx,   char diretype )
{
    int i=0;
    int transit=0;
    int category=0;
    if( diretype=='L' )
    {
        if( curHmm->state=='B' )
    	       category = KEEP_POS;
	       else if( curHmm->state=='M' )
    	       category = LEFT_INS;
        else if( curHmm->state=='I' )
            category = LEFT_INS;
        else if( curHmm->state=='D' )
    	       category = KEEP_POS;
    }
    else
    {
        if( curHmm->state=='B' )
    	       category = KEEP_POS;
	       else if( curHmm->state=='M' )
    	       category = RIGHT_INS;
        else if( curHmm->state=='I' )
            category = RIGHT_INS;
        else if( curHmm->state=='D' )
    	       category = KEEP_POS;
    }

    if( nextHmm->state=='M' || nextHmm->state=='E' )
        transit=HMM_M;
    else if( nextHmm->state=='I')
        transit=HMM_I;
    else if( nextHmm->state=='D')
        transit=HMM_D;

    prod[idx].category=category;
    prod[idx].direct=diretype;
    prod[idx].ltype=curHmm->state;
    prod[idx].lid=curHmm->UID;
    prod[idx].rtype=nextHmm->state;
    prod[idx].rid=nextHmm->UID;
    
    //calculate probability
    if( curHmm->state=='B' )
        prod[idx].singlep=curHmm->trprob[ transit ];
    else if( curHmm->state=='D')
        prod[idx].singlep=curHmm->trprob[ transit ];
    else
    {
        for( i=0; i<4; i++ )
            prod[idx].basep[i] = curHmm->trprob[ transit ] + curHmm->emission[i];
    }
}

int Bulge::serialProducts( GramRule * prod,  int hidx, char diretype )
{
    #ifdef DEBUG_DEV
        printf("Bulge::serialProducts\n");
    #endif
    
    HMMGram *curM=NULL,  *curI=NULL,  *curD=NULL;
    HMMGram *nexM=NULL,  *nexD=NULL;
    
    curM=pLoopGram;
    nexM=curM->nextM;

    while( nexM!=NULL )
    {
        nexD=curM->nextD;
        curI=curM->nextI;
                               //***** M
        genSingleHMMProd( curM, nexM, prod, hidx++, diretype );	 //M -> xM    
        if( curM->nextD!=NULL )      //M -> xD
            genSingleHMMProd( curM, nexD, prod, hidx++, diretype );
        if( curM->nextI!=NULL )    	//M -> xI
             genSingleHMMProd( curM, curI, prod, hidx++, diretype );

        if( curI!=NULL )    	 //***** I
        {
            genSingleHMMProd( curI, curI, prod, hidx++, diretype );  //I -> xI
            if( curI->nextD!=NULL )          //I -> xD
                genSingleHMMProd( curI, nexD, prod, hidx++, diretype );
            if( curI->nextM!=NULL )          //I -> xM
                genSingleHMMProd( curI, nexM, prod, hidx++, diretype );
        }

        if( curD!=NULL )    	 //***** D
        {
            if( curD->nextD!=NULL )          //D -> xD
                genSingleHMMProd( curD, nexD, prod, hidx++, diretype );
            if( curD->nextI!=NULL )          //D -> xI
                genSingleHMMProd( curD, curI, prod, hidx++, diretype );
            if( curD->nextM!=NULL )          //D -> xM
                genSingleHMMProd( curD, nexM, prod, hidx++, diretype );
        }
        curD=nexD;
        curM=nexM;
        nexM=curM->nextM;
    }
    return hidx;
}


void Bulge::printHMMGrams( char direct )
{
    #ifdef DEBUG_DEV
       printf("Bulge::printHMMGrams\n");
    #endif
    printf("********* Bulge HMM Gram: **********\n");
    printLoopGrams( direct );
    printf("********* End of Bulge HMM Gram **********\n\n");
}


//********************************************************************************
//*******************    implementation of HMMBuilder class         ****************
//********************************************************************************
HMMBuilder::HMMBuilder( )
{
    #ifdef DEBUG_DEV
       printf("HMMBuilder::HMMBuilder\n");
    #endif	
    allloops=NULL;
  	numloops=0;
    sequences=NULL;
  	pasta=NULL;
    seqlength=0;
    numseqs=0;
    pseudocount = PSEUDOCNT;
    pseudocntgap = PSDCNTGAP;
}

HMMBuilder::~HMMBuilder( )
{
    #ifdef DEBUG_DEV
       printf("HMMBuilder::~HMMBuilder\n");
    #endif
    int i=0;
    if( sequences!=NULL )
    {
        for(i=0; i<numseqs; i++ )
        {
            if( sequences[i]!=NULL)
            {
                delete [ ] sequences[i];
                sequences[i]=NULL;
            }
        }
        delete [ ] sequences;
        sequences=NULL;
    }
    if( pasta!=NULL)
    {
        delete [ ] pasta;
        pasta=NULL;
  	}
}

void HMMBuilder::freeAllLoops( )
{
    #ifdef DEBUG_DEV
       printf("HMMBuilder::freeAllLoops\n");
    #endif
    if( allloops!=NULL )
    {
        delete [ ] allloops;
        allloops=NULL;
    }
}

Loop * HMMBuilder::getAllLoops( int &num  )
{
    #ifdef DEBUG_DEV
       printf("HMMBuilder::getAllLoops\n");
    #endif
    num=numloops;
    return allloops;
}

void HMMBuilder::setPseudocnt( double pdcnt, double pcgap )
{
    pseudocount = pdcnt;
    pseudocntgap = pcgap;
}

int HMMBuilder::setTrainingData( DataLoader &trdata )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::setTrainingData\n");
    #endif
    int i=0, j=0;
    if( trdata.numloops<=0 )
    {
        numloops=0;
        return 0;
    }
    numseqs = trdata.getNumOfTrainSeq( );
    seqlength = trdata.getSeqLength( );
    
    //copy pasta
    pasta=new char[ seqlength+1 ];
    for(i=0; i<seqlength; i++ )
    {
        pasta[i] = trdata.pasta[0][i];
    }
    pasta[i]='\0';
    
    //copy training sequences
    sequences = new int* [numseqs];
    for(i=0; i<numseqs; i++ )
    {
        sequences[i] = new int [ seqlength ];
        for( j=0; j<seqlength; j++ )
        {
            sequences[i][j]=trdata.trainingseq[i][j];
        }
    }
    
    //import loop
    numloops=trdata.numloops;
    allloops=new Loop[numloops];
    for(i=0; i<numloops; i++)
    {
        allloops[i].UID = trdata.loopary[i].id;
        allloops[i].numseqs = numseqs;
        allloops[i].Pos[0] = trdata.loopary[i].Pos[0];
        allloops[i].Pos[1] = trdata.loopary[i].Pos[1];
        allloops[i].looplen = trdata.loopary[i].looplen;
        
        allloops[i].lenE = trdata.loopary[i].lenE;
        allloops[i].lenSD = trdata.loopary[i].lenSD;
        for( j=0; j<2; j++ )
        {
            allloops[i].offsetE[j] = trdata.loopary[i].offsetE[j];
            allloops[i].offsetSD[j] = trdata.loopary[i].offsetSD[j];
        }
        trdata.copyLogBgBaseHz( allloops[i].bgBaseLogHz, 5 );
    }
    return 1;
}

void HMMBuilder::seperateLoop(  )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::seperateLoop\n");
    #endif
    int i=0, j=0, k=0;
    int beg=0, end=0, len=0;
    
    for( i=0; i<numloops; i++ )
    {
        len=allloops[i].looplen;
        beg=allloops[i].Pos[0];
        end=allloops[i].Pos[1];
        allloops[i].loopseqs = new int*[ numseqs ];
        for( j=0; j<numseqs; j++ )
        {
            allloops[i].loopseqs[j]= new int[len];
            for( k=beg; k<=end; k++ )
            {
                allloops[i].loopseqs[j][k-beg] = sequences[j][k];
            }
        }
    }
}

char HMMBuilder::determColomnState( Loop *curloop, int curpos )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::determColomnState\n");
    #endif
    
    int gaps=0, bases=0;
    int i=0;
    int seqnum=curloop->numseqs;
    char type=' ';
    
    for( i=0; i<seqnum; i++ )
    {
        //printf("%d", curloop->loopseqs[i][curpos] );
        if( curloop->loopseqs[i][curpos]==0 )
            gaps++;
    }
    bases = seqnum - gaps;
    if( bases>=gaps)
        type='M';
    else
        type='I';
    return type;
}
    
int HMMBuilder::obtainNextM( Loop *curloop, int curidx, int &dist )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::obtainNextM\n");
    #endif
    
    int status=0;
    int beg=curidx;
    char type=' ';
    
    while( 1 )
    {
        if( curidx < curloop->looplen )
        {
            type=determColomnState( curloop,  curidx );
            if( type=='M')
            {
                status=1;
                break;
            }
            else
                curidx++;
        }
        else
        {
            status=0;
            break;
        }
    }
    dist=curidx-beg;
    return status;
}

void HMMBuilder::markState(Loop *curloop, int idx, int dist )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::markState\n");
    #endif
    int i=0;
    int len=idx+dist;
    if( len < curloop->looplen )
        curloop->state[ len ]='m';
    for( i=0; i<dist; i++ )
    {
        curloop->state[idx+i]='i';
    }
}

void HMMBuilder::countTransState( HMMGram *preHM, Loop *curloop, int idx, int dist, int **Tri )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::countTransState\n");
    #endif
    int i=0, j=0;
    int rowins=0;
    int seqnum=curloop->numseqs;
    int *precol=NULL;
    int *curcol=NULL;

    precol=new int[seqnum];
    curcol=new int[seqnum];

    for( i=0; i<3; i++ ) //0:D   1:I   2:M
    {
        Tri[0][i]=0;  //0: D
        Tri[1][i]=0;  //1: I
        Tri[2][i]=0;  //2: M
    }
    
    for( i=0; i<seqnum; i++ )
    {
        if( preHM->state=='B')
            precol[i]=5;
        else
            precol[i]=curloop->loopseqs[i][idx-1];
            
        if( idx+dist==curloop->looplen ) //'E'
            curcol[i]=5;
        else
            curcol[i]=curloop->loopseqs[i][idx+dist];
    }

    for( i=0; i<seqnum; i++ )
    {
        rowins=0;
        for( j=0; j<dist; j++ )
        {
            if( rowins && curloop->loopseqs[i][idx+j]!=0 )
            {
                 Tri[1][1]++;
            }
            if( curloop->loopseqs[i][idx+j] != 0)
            {
                 rowins=1;
            }
        }
        if( rowins )
        {
            if( precol[i]!=0 )
               Tri[2][1]++;  //M-I
            else
               Tri[0][1]++;  //D-I
            
            if( curcol[i]!=0 )
               Tri[1][2]++;  //I-M
            else
               Tri[1][0]++;  //I-D
        }
        else
        {
            if( precol[i]!=0)  //M
            	 if( curcol[i]!=0 )
            	      Tri[2][2]++;  //M-M
            	 else
            	      Tri[2][0]++;  //M-D
            else   //D
            	 if( curcol[i]!=0 )
            	      Tri[0][2]++;  //D-M
            	 else
            	      Tri[0][0]++;  //D-D
        }
    }
    delete [ ] precol;
    precol=NULL;
    delete [ ] curcol;
    curcol=NULL;
} 

void HMMBuilder::calcComTransProb( HMMGram *preM, HMMGram *preI, HMMGram *preD,  int **Trans )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::calcComTransProb\n");
    #endif
    double pdocnt=pseudocount;
    double sum=0.0;
    
    if( preM!=NULL)
    {
        sum = Trans[2][0]+Trans[2][1]+Trans[2][2]+3*pdocnt;  //0:D   1:I   2:M
        preM->trprob[HMM_M]= (Trans[2][2] + pdocnt)/sum;
        preM->trprob[HMM_I]= (Trans[2][1] + pdocnt)/sum;
        preM->trprob[HMM_D]= (Trans[2][0] + pdocnt)/sum;
        #ifdef DEBUG_DEV
            printf("preM: M: %6.4f  I: %6.4f  D: %6.4f\n",
               preM->trprob[HMM_M], preM->trprob[HMM_I], preM->trprob[HMM_D] );
        #endif
    }
    if( preI!=NULL)
    {
        sum = Trans[1][0]+Trans[1][1]+Trans[1][2]+3*pdocnt;  //0:D   1:I   2:M
        preI->trprob[HMM_M]= (Trans[1][2] + pdocnt)/sum;
        preI->trprob[HMM_I]= (Trans[1][1] + pdocnt)/sum;
        preI->trprob[HMM_D]= (Trans[1][0] + pdocnt)/sum;
        #ifdef DEBUG_DEV
            printf("preI: M: %6.4f  I: %6.4f  D: %6.4f\n",
               preI->trprob[HMM_M], preI->trprob[HMM_I], preI->trprob[HMM_D] );
        #endif
    }
    if( preD!=NULL)
    {
        sum = Trans[0][0]+Trans[0][1]+Trans[0][2]+3*pdocnt;  //0:D   1:I   2:M
        preD->trprob[HMM_M]= (Trans[0][2] + pdocnt)/sum;
        preD->trprob[HMM_I]= (Trans[0][1] + pdocnt)/sum;
        preD->trprob[HMM_D]= (Trans[0][0] + pdocnt)/sum;
        #ifdef DEBUG_DEV
            printf("preD: M: %6.4f  I: %6.4f  D: %6.4f\n",
               preD->trprob[HMM_M], preD->trprob[HMM_I], preD->trprob[HMM_D] );
        #endif
    }
}

void HMMBuilder::calcEndTransProb( HMMGram *preM, HMMGram *preI, HMMGram *preD,  int **Trans )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::calcEndTransProb\n");
    #endif
    double pdocnt=pseudocount;
    double sum=0;
    
    if( preM!=NULL)
    {
        sum = Trans[2][1]+Trans[2][2]+2*pdocnt;  //0:D   1:I   2:M
        preM->trprob[HMM_M]= (Trans[2][2] + pdocnt)/sum;
        preM->trprob[HMM_I]= (Trans[2][1] + pdocnt)/sum;
    }
    if( preI!=NULL)
    {
        sum = Trans[1][1]+Trans[1][2]+2*pdocnt;  //0:D   1:I   2:M
        preI->trprob[HMM_M]= (Trans[1][2] + pdocnt)/sum;
        preI->trprob[HMM_I]= (Trans[1][1] + pdocnt)/sum;
    }
    if( preD!=NULL)
    {
        sum = Trans[0][1]+Trans[0][2]+2*pdocnt;  //0:D   1:I   2:M
        preD->trprob[HMM_M]= (Trans[0][2] + pdocnt)/sum;
        preD->trprob[HMM_I]= (Trans[0][1] + pdocnt)/sum;
    }
}

void HMMBuilder::calcEmissionProb( HMMGram *curHG, Loop *curloop, int pos )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::calcEmissionProb\n");
    #endif
    int i=0;
    int seqnum=curloop->numseqs;
    double pdocnt=pseudocount;
    double base[5]={0.0, 0.0, 0.0, 0.0, 0.0};
    double sum=0.0;
    
    for( i=0; i<seqnum; i++ )
    {
        base[ curloop->loopseqs[i][pos] ]++;
    }
    sum = base[1]+base[2]+base[3]+base[4];
    for( i=1; i<5; i++ )
    {
        curHG->emission[i-1]= (base[i] + pdocnt)/(sum + 4*pdocnt);
    }    
}

void HMMBuilder::calcInsertEmProb( HMMGram *curHG, Loop *curloop, int idx, int dist )
{
    int i=0, j=0;
    int seqnum=curloop->numseqs;
    double pdocnt=pseudocount;
    double base[5]={0.0, 0.0, 0.0, 0.0, 0.0};
    double sum=0.0;
    
    for(i=0; i<seqnum; i++ )
    {
        for(j=0; j<dist; j++ )
        {
            base[ curloop->loopseqs[i][idx+j] ]++;
        }
    }
    sum=base[1]+base[2]+base[3]+base[4];
    for( i=1; i<5; i++ )
    {
        curHG->emission[i-1]= (base[i] + pdocnt)/(sum + 4*pdocnt);
    }
}

int HMMBuilder::generateHMM( Loop *curloop, int uid, int subid )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::generateHMM\n");
    #endif
    
    int i=0, pi=0;
    int status=0;
    int **Trans=NULL;
    int len=curloop->looplen;
    int idx=0, dist=0;

    Trans=new int*[3];
    for( i=0; i<3; i++ )
    {
        Trans[i]=new int [3];
        Trans[i][0]=0;
        Trans[i][1]=0;
        Trans[i][2]=0;
    }
    HMMGram *preM=NULL, *preD=NULL, *preI=NULL;
    HMMGram *curM=NULL, *curD=NULL;

    curloop->state = new char[len+1];
    memset( curloop->state, '\0', len+1 );

    preM = new HMMGram('B', 0, uid++, subid );
    preI = new HMMGram('I', 0, uid++, subid );
    curloop->pLoopGram=preM;
    while( 1 )
    {
        #ifdef DEBUG_DEV
            printf("\n");
        #endif
        status=obtainNextM( curloop, idx, dist );
        markState(curloop, idx, dist );
        #ifdef DEBUG_DEV
            printf("idx=%d, dist=%d, i=%d, status=%d\n", idx, dist, pi+1, status);
        #endif
        if( status!=0)
            curM = new HMMGram( 'M',  pi+1, uid++, subid);
        else
        {
            curM = new HMMGram( 'E',  pi+1, uid++, subid);
            curloop->numstates=pi+2;
        }
        
        if( status!=0 )
            curD = new HMMGram( 'D',  pi+1, uid++, subid);
        else
            curD = NULL;

        countTransState( preM, curloop, idx, dist, Trans );
        #ifdef DEBUG_DEV
        int j=0;
        printf("----------->pos=%d,  dis=%d, curpos=%d\n", idx, dist, idx+dist);
        for(i=0; i<3; i++ )
        {
            for(j=0; j<3; j++ )
                 printf(" %2d ", Trans[i][j]);
            printf("\n");
        }
        #endif
        if( curM->state != 'E' )
            calcEmissionProb( curM, curloop, idx+dist );
        calcInsertEmProb( preI, curloop, idx, dist  );
        if( preM!=NULL)
        {
            preM->nextM=curM;
            preM->nextI=preI;
            preM->nextD=curD;
        }
        if( preI!=NULL )
        {
            preI->nextM=curM;
            preI->nextI=preI;
            preI->nextD=curD;
        } 
        if( preD!=NULL )
        {
            preD->nextM=curM;
            preD->nextI=preI;
            preD->nextD=curD;
        }
        if( status!=0)
            calcComTransProb( preM, preI, preD, Trans );
        else
            calcEndTransProb( preM, preI, preD, Trans );
        
        if( status!=0)
            preI = new HMMGram( 'I',  pi+1, uid++, subid);
        else
            break;
        preD = curD;
        preM = curM;
        curM = NULL;
        curD = NULL;
        idx=idx+dist+1;
        pi++;
    }

    for( i=0; i<3; i++ )
    {
        delete [ ] Trans[i];
        Trans[i]=NULL;
    }
    delete [ ] Trans;
    Trans=NULL;
     
    return uid;
}

Bulge * HMMBuilder::buildBulge( int **sequences, int seqnum,   int posbeg,
                               int len,   char LorR,  int &gramuid,  int subid )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::buildBulge\n");
    #endif
    int i=0, j=0;

    Bulge *ploop=new Bulge( seqnum, len );
    //set the training seqs
    ploop->loopseqs=new int* [ seqnum ];
    for( i=0; i<seqnum; i++ )
    {
        ploop->loopseqs[i]=new int [ len ];
        for( j=0; j<len; j++ )
        {
            if( LorR=='L')
        	    ploop->loopseqs[i][j]=sequences[i][posbeg+j];
            else
                ploop->loopseqs[i][j]=sequences[i][posbeg-j];
        }
    }
    
    //************ print *************
    gramuid = generateHMM( ploop, gramuid, subid );

    return ploop;
}

int HMMBuilder::buildHMM( DataLoader &trdata )
{
    #ifdef DEBUG_DEV
        printf("HMMBuilder::buildHMM\n");
    #endif
    int status=0;
    int i=0;
    
    status=setTrainingData( trdata ); //---- Load training pasta
    if( !status )
        return 0;
    //---- Seperate the loop from sequences
    seperateLoop(  );
    
    //---- Generate HMM for every stem
    for(i=0; i<numloops; i++ )
    {
        generateHMM( &allloops[i] );
    }
    //---- Print out the grammar
    return 1;
}

void HMMBuilder::printAllGrams( int trda )
{
    #ifdef DEBUG_DEV
    printf("HMMBuilder::printAllGrams\n");
    #endif
    
    int i=0, j=0;
    int no=0;
    int len=0;
    int numofseqs=0;
    Loop *curloop=NULL;

    for( no=0; no<numloops; no++)
    {
        if( trda != 0 )
        {
    	    curloop=&allloops[no];
    	    len=curloop->looplen;

            numofseqs=curloop->numseqs;
        
            //print sequences
            printf("\n========================================================\n" );
            printf("====================== LOOP No. %-2d======================\n", no+1);
            printf("---> Training Loops(Length %d):\n", len );
    	    printf(" State:  %s\n", curloop->state );
        	for(i=0; i<numofseqs; i++)
        	{
        	   printf("    %3d  ", i+1);
               for( j=0; j<len; j++)
               {
                   printf("%c", numtochar(curloop->loopseqs[i][j]));
               }
               printf("\n");
        	}
    	}
    	else
    	{
    	    printf("\n====================== LOOP No. %-2d======================", no);
    	}        
    	//print grammar
    	printf("\n---> HMM Grammar:\n");
        allloops[no].printHMMGrams( );
    }
}
