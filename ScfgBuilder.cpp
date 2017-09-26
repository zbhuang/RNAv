#pragma warning(disable: 4786)

#include "ScfgBuilder.h"

//********************************************************************************
//*******************    implementation of ScfgGram class         ****************
//********************************************************************************
ScfgGram::ScfgGram( )
{
    #ifdef DEBUG_DEV
       printf("ScfgGram::ScfgGram\n");
    #endif
    int i=0, j=0;
    state=' ';
    ID=-1;
    SubID=-1;
    UID=-1;

    nextLH=NULL;
    nextLI=NULL;
    nextS=NULL;
    nextRI=NULL;
    nextRH=NULL;
    
    for(i=0; i<5; i++ )
        for(j=0; j<5; j++ )
        {
            emprobS[i][j]=ZEROPROB;
        }
    for(i=0; i<2; i++ )
        for(j=0; j<4; j++ )
        {
            emprobI[i][j]=ZEROPROB;
        }
    for(i=0; i<5; i++ )
        trprob[i]=ZEROPROB;
}

ScfgGram::ScfgGram( char settype,  int setid,  int setuid, int subid )
{
    #ifdef DEBUG_DEV
        printf("ScfgGram::ScfgGram( parameter )\n");
    #endif
    int i=0, j=0;
    state=settype;
    ID=setid;
    SubID=subid;
    UID=setuid;
    nextLH=NULL;
    nextLI=NULL;
    nextS=NULL;
    nextRI=NULL;
    nextRH=NULL;
    for(i=0; i<5; i++ )
        for(j=0; j<5; j++ )
            emprobS[i][j]=ZEROPROB;
    for(i=0; i<2; i++ )
        for(j=0; j<4; j++ )
            emprobI[i][j]=ZEROPROB;
    for(i=0; i<5; i++ )
        trprob[i]=ZEROPROB;
}

ScfgGram::~ScfgGram( )
{
    nextLH=NULL;
    nextLI=NULL;
    nextS=NULL;
    nextRI=NULL;
    nextRH=NULL;
}

void ScfgGram::log1DProb( double prob[],  int size )
{
    #ifdef DEBUG_DEV
    printf("ScfgGram::log1DProb\n");
    #endif
    int i=0;
    for(i=0; i<size; i++ )
    {
        if( prob[i] != INVLDPROB )
        {
             prob[i]=log( prob[i] );
        }
    }
}

void ScfgGram::log2DProb( double *prob[], int row, int col )
{
    #ifdef DEBUG_DEV
    printf("ScfgGram::log2DProb\n");
    #endif
    int i=0, j=0;
    for(i=0; i<row; i++ )
    {
        for(j=0; j<col; j++)
        {
            if( prob[i][j]!= INVLDPROB )
            {
                prob[i][j]=log( prob[i][j] );
            }
        }
    }
}

void ScfgGram::logGrammarProb( double logfreq[5], double priorMatrix[5][5], double h )
{
    int i=0, j=0;
    //pair emission    
    if( state == 'S' )
    {
        //(1-h)*Mt + h*Mp  considering the prior matrix
        normalizeMatrix( emprobS, priorMatrix, h );
        
        //calculate log/logodds
        for(i=0; i<5; i++ )
        {
            for(j=0; j<5; j++)
            {
                if( emprobS[i][j]!= ZEROPROB )
                {
                    if( LOGODDS )
                        emprobS[i][j]=log( emprobS[i][j] )-logfreq[i]-logfreq[j];//0,1,2,3,4
                    else
                        emprobS[i][j]=log( emprobS[i][j] );
                }
                else
                {
                    emprobS[i][j]=INVLDPROB;
                }
            }
        }
    }
    else  if( state == 'R' || state == 'L')  //single emission 
    {
        for(i=0; i<2; i++ )
        {
             for(j=0; j<4; j++)
            {
                if( emprobI[i][j]!= ZEROPROB )
                {
                    if( LOGODDS )
                        emprobI[i][j]=log( emprobI[i][j] )-logfreq[j+1]; //0,1,2,3
                    else
                        emprobI[i][j]=log( emprobI[i][j] );
                }
                else
                {
                    emprobS[i][j]=INVLDPROB;
                }
            }
        }
    }
    
    for(i=0; i<SCFG_NUM_NEXTP; i++ )
    {
        if( trprob[i] != ZEROPROB )
        {
             trprob[i]=log( trprob[i] );
        }
        else
        {
            trprob[i]=INVLDPROB;
        }
    }
}

int ScfgGram::hasLinkVState( )
{
    if(  nextS!=NULL ) 
        return 1;
    else
        return 0;
}

int ScfgGram::countNumRules( )
{
    int num=0;
    if( nextLH!=NULL )
        num++;
    if( nextLI!=NULL )
        num++;
    if( nextRI!=NULL )
        num++;
    if( nextRH!=NULL )
        num++;
    if( nextS!=NULL )
        num++;
    return num;
}

//********************************************************************************
//********************        implementation of Stem class       *****************
//********************************************************************************
Stem::Stem( )
{
    #ifdef DEBUG_DEV
        printf("Stem::Stem\n");
    #endif
    
    int i=0, j=0;
    idchar=' ';
    for(i=0; i<4; i++ )
        Pos[i]=-1;
    stemlen=-1;
    numseqs=-1;
    distance=-1;
    logProb=0;
    numntm=0;
    numprod=0;
    
    for(i=0; i<2; i++ )
    {
        state[i]=NULL;
        armpasta[i]=NULL;
        armseq[i]=NULL;
    }
    pStemGram=NULL;
    nonterm=NULL;
    prod=NULL;
    for( i=0; i<5; i++ )
    {
        bgBaseLogHz[i]=0.0;
    }
    for( i=0; i<5; i++ )
    {
        for( j=0; j<5; j++ )
        {
            priorPairFreq[i][j]=0.0;
            sumcolpair[i][j]=0;
        }
    }
    colpairNum=NULL;
    stemStrength=0.0;
    regInfo=NULL;
}

Stem::~Stem( )
{
    #ifdef DEBUG_DEV
       printf("Stem::~Stem\n");
    #endif
    
    int j=0, k=0;
    
    freeScfgGrams( );
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
    for(j=0; j<2; j++ )
    {
        if( armpasta[j]!=NULL)  {
            delete [ ] armpasta[j];
            armpasta[j]=NULL;
        }
        if( state[j]!=NULL) {
            delete [ ] state[j];
            state[j]=NULL;
        }
        //the seqs involved in modeling
        for( k=0; k<numseqs; k++ ) {
            if( armseq[j][k]!=NULL) {
                delete [ ] armseq[j][k];
                armseq[j][k]=NULL;
            }
        }
        if( armseq[j]!=NULL )  {
            delete [ ] armseq[j];
            armseq[j]=NULL;
        }
        // the original sequences
        for( k=0; k<numorgs; k++ ) {
            if( armorg[j][k]!=NULL) {
                delete [ ] armorg[j][k];
                armorg[j][k]=NULL;
            }
        }
        if( armorg[j]!=NULL )  {
            delete [ ] armorg[j];
            armorg[j]=NULL;
        }
    }
    freeFreqMem( );
}

void Stem::freeScfgGrams( )
{
    #ifdef DEBUG_DEV
       printf("Stem::freeScfgGrams\n");
    #endif

    ScfgGram *curSG=pStemGram, *nexSG=NULL;
    ScfgGram *lIG=NULL, *rIG=NULL;
    Bulge *lBg=NULL, *rBg=NULL, *nxBg=NULL;

    while( curSG!=NULL )
    {
        lBg=curSG->nextLH;
        lIG=curSG->nextLI;
        nexSG=curSG->nextS;
        rIG=curSG->nextRI;
        rBg=curSG->nextRH;
        
        delete curSG;
        curSG=NULL;
        if( lIG!=NULL )
        {
            delete lIG;
            lIG=NULL;
            if( lBg!=NULL )
            {   
                while( 1 )
                {
                     lIG=lBg->nextI;
                     nxBg=lBg->nextLp;
                     delete lBg;
                     lBg=NULL;
                     if( lIG!=NULL )
                     {
                        delete lIG;
                        lIG=NULL;
                    }
                    if( nxBg!=NULL )
                        lBg=nxBg;
                    else
                        break;
                }
          	}
        }
        //*********
        if( rIG!=NULL )
        {
            delete rIG;
            rIG=NULL;
            if( rBg!=NULL )
            {   
                while( 1 )
                {
                     rIG=rBg->nextI;
                     nxBg=rBg->nextLp;
                  	delete rBg;
                     rBg=NULL;
                     
                     if( rIG!=NULL )
                     {
                        delete rIG;
                        rIG=NULL;
                    }
                    if( nxBg!=NULL )
                        rBg=nxBg;
                    else
                        break;
                }
          	}
        }
        curSG=nexSG;
    }
}

//--------------- log prob ----------------
void Stem::logProbforStem( )
{
    #ifdef DEBUG_DEV
        printf("Stem::logProbforStem\n");
    #endif

    ScfgGram *curSG=pStemGram;
    ScfgGram *curLIG=NULL, *curRIG=NULL;
    Bulge *curBg=NULL;

    while( curSG!=NULL )
    {
        curSG->logGrammarProb( bgBaseLogHz, priorPairFreq, priorcoef );        //S
        //*********
        curLIG=curSG->nextLI;        //Left
        if( curLIG!=NULL )
        {
          	curLIG->logGrammarProb( bgBaseLogHz, priorPairFreq, priorcoef );

            curBg=curSG->nextLH;
            if( curBg!=NULL )
            {   
                while( 1 )
                {
                    if( !curBg->isLogProb( ) )
                        curBg->logProbforBulge( bgBaseLogHz );//log of bulge
                    
                    if( curBg->nextI!=NULL )
                        curBg->nextI->logGrammarProb( bgBaseLogHz, priorPairFreq, priorcoef );
                 
                    if( curBg->nextLp!=NULL )
                        curBg=curBg->nextLp;
                    else
                        break;
                }
          	}
        }
        //*********
        curRIG=curSG->nextRI;
        if( curRIG!=NULL )
        {
          	curRIG->logGrammarProb( bgBaseLogHz, priorPairFreq, priorcoef );

            curBg=curSG->nextRH;
            if( curBg!=NULL )
            {   
                while( 1 )
                {
                    if( !curBg->isLogProb( ))
                        curBg->logProbforBulge( bgBaseLogHz );//log of bulge
                    
                    if( curBg->nextI!=NULL )
                        curBg->nextI->logGrammarProb( bgBaseLogHz, priorPairFreq, priorcoef );
                 
                    if( curBg->nextLp!=NULL )
                        curBg=curBg->nextLp;
                    else
                        break;
                }
          	}
        }
        curSG=curSG->nextS;
    }
}

int Stem::countNumAllRules(  )
{
    #ifdef DEBUG_DEV
        printf("Stem::countNumAllRules\n");
    #endif
    int allrules=0;
    int VLinks[2]={0, 0}; //0: in,  1: out
    Bulge *curBg=NULL;
    ScfgGram *curIns=NULL,  *curVT=NULL, *curSG=pStemGram;
    
    if( curSG==NULL )
    {
        numprod = 0;
        return 1;
    }
    curVT=curSG->nextS;  //S->V
    while( curVT!=NULL )
    {
        allrules += curSG->countNumRules( );
        VLinks[0] += curSG->hasLinkVState( );
        VLinks[1] = curVT->countNumRules( );
        
        curIns=curSG->nextLI;        //S->L
        if( curIns!=NULL )
        {
          	allrules += curIns->countNumRules( );
          	VLinks[0] += curIns->hasLinkVState( );
          	
            curBg=curSG->nextLH;
            if( curBg!=NULL )
            {
                while( 1 )
                {
                    allrules += curBg->countNumAllRules( );
                    VLinks[0] += curBg->hasLinkVState( );
                    
                    if( curBg->nextI!=NULL )
                    {
                        allrules += curBg->nextI->countNumRules( );
                        VLinks[0] += curBg->nextI->hasLinkVState( );
                    }
                    if( curBg->nextLp!=NULL )
                        curBg=curBg->nextLp;
                    else
                        break;
                }
          	}
        }
        
        curIns=curVT->nextRI;    //V->R
        if( curIns!=NULL )
        {
          	allrules += curIns->countNumRules( );
            curBg=curVT->nextRH;
            if( curBg!=NULL )
            {   
                while( 1 )
                {
                    allrules += curBg->countNumAllRules( );
                    if( curBg->nextI!=NULL )
                        allrules += curBg->nextI->countNumRules( );
                    if( curBg->nextLp!=NULL )
                        curBg=curBg->nextLp;
                    else
                        break;
                }
          	}
        }
        //process the number of rules
        allrules += VLinks[1]*VLinks[0] - VLinks[0];
        curSG = curVT->nextS; //V->S
        curVT = curSG->nextS; //S->V
        VLinks[0]=0;
        VLinks[1]=0;
    }

    numprod = allrules;
    return 1;
}    
    
int Stem::countNumAllNtms(  )
{
    #ifdef DEBUG_DEV
        printf("Stem::countNumAllNtms\n");
    #endif
    int allntms=0;
    Bulge *curBg=NULL;
    ScfgGram *curIns=NULL,  *curSG=pStemGram;

    while( curSG!=NULL )
    {
        allntms++;
        curIns=curSG->nextLI;        //Left
        if( curIns!=NULL )
        {
          	allntms++;
            curBg=curSG->nextLH;
            if( curBg!=NULL )
            {   
                while( 1 )
                {
                    allntms += curBg->countNumAllNtms( );
                    if( curBg->nextI!=NULL )
                        allntms++;
                    if( curBg->nextLp!=NULL )
                        curBg=curBg->nextLp;
                    else
                        break;
                }
          	}
        }
        //*********
        curIns=curSG->nextRI;
        if( curIns!=NULL )
        {
          	allntms++;
            curBg=curSG->nextRH;
            if( curBg!=NULL )
            {   
                while( 1 )
                {
                    allntms += curBg->countNumAllNtms( );
                    if( curBg->nextI!=NULL )
                        allntms++;
                    if( curBg->nextLp!=NULL )
                        curBg=curBg->nextLp;
                    else
                        break;
                }
          	}
        }
        curSG=curSG->nextS;
    }
    numntm = allntms;
    return 1;
}

void Stem::arrayNonterminal(  )
{
    #ifdef DEBUG_DEV
        printf("Stem::arrayNonterminal\n");
    #endif

    int idx=0;
    ScfgGram *curSG=NULL, *curIns=NULL;
    Bulge *curBg=NULL;
    
    curSG=pStemGram;
    while( curSG!=NULL )
    {
        idx=curSG->UID;         // S-->S
        nonterm[idx].uid=idx;
        nonterm[idx].state=curSG->state;
        nonterm[idx].id=curSG->ID;
        
        curIns=curSG->nextLI;	 //Left
        if( curIns!=NULL )
        {
          	idx=curIns->UID;
            nonterm[idx].uid=idx;
            nonterm[idx].state=curIns->state;
            nonterm[idx].id=curIns->ID;
            
            curBg=curSG->nextLH;
            if( curBg!=NULL )
            {   
                while( 1 )
                {
                    curBg->arrayNonterminal( nonterm );
                    if( curBg->nextI!=NULL )
                    {
                        idx=curBg->nextI->UID;
                        nonterm[idx].uid=idx;
                        nonterm[idx].state=curBg->nextI->state;
                        nonterm[idx].id=curBg->nextI->ID;
                    }
                    if( curBg->nextLp!=NULL )
                        curBg=curBg->nextLp;
                    else
                        break;
                }
          	}
        }
        curIns=curSG->nextRI;  //Right
        if( curIns!=NULL )
        {
          	idx=curIns->UID;
            nonterm[idx].uid=idx;
            nonterm[idx].state=curIns->state;
            nonterm[idx].id=curIns->ID;
            
            curBg=curSG->nextRH;
            if( curBg!=NULL )
            {   
                while( 1 )
                {
                    curBg->arrayNonterminal( nonterm );
                    if( curBg->nextI!=NULL )
                    {
                        idx=curBg->nextI->UID;
                        nonterm[idx].uid=idx;
                        nonterm[idx].state=curBg->nextI->state;
                        nonterm[idx].id=curBg->nextI->ID;
                    }
                    if( curBg->nextLp!=NULL )
                        curBg=curBg->nextLp;
                    else
                        break;
                }
          	}
        }
        curSG=curSG->nextS;
    }
}

void Stem::initGramProd( )
{
    int i=0, j=0;
    int no=0;
    for( no=0; no<numprod; no++)
    {
        prod[no].category=-1;
        prod[no].ltype=' ';
        prod[no].lid=-1;
        prod[no].rtype=' ';
        prod[no].rid=-1;
        for(i=0; i<5; i++)
        {
            for(j=0; j<5; j++)
            {
                prod[no].pairp[i][j]=INVLDPROB;
            }
        }
    }
}

void Stem::initNontermArray( )
{
    int i=0;
    for( i=0; i<numntm; i++)
    {
        nonterm[i].uid=0;
        nonterm[i].state=' ';
        nonterm[i].id=0;
    }
}

// functions to generate productions for CYK
int Stem::genHMMtoHMM( int idx, Bulge *curHG, Bulge *nexHG, char dirtype )
{
    prod[idx].category=KEEP_POS;
    prod[idx].direct=dirtype;
    prod[idx].ltype=curHG->getLastGramState( );
    prod[idx].lid=curHG->getLastGramUID( );
    prod[idx].rtype=nexHG->getFirstGramState( );
    prod[idx].rid=nexHG->getFirstGramUID( );
    
     prod[idx].singlep = curHG->trprob[ 1 ];
    
    return 1;
}

int Stem::genHMMtoSCFG( int idx, Bulge *curHG, ScfgGram *curSG, char dirtype)
{
    prod[idx].category=KEEP_POS;
    prod[idx].direct=dirtype;
    prod[idx].ltype=curHG->getLastGramState( );
    prod[idx].lid=curHG->getLastGramUID( );
    prod[idx].rtype=curSG->state;
    prod[idx].rid=curSG->UID;
    
    if( curSG->state=='S'|| curSG->state=='V' )
        prod[idx].singlep = curHG->trprob[ 1 ];
    else
        prod[idx].singlep = curHG->trprob[ 0 ];
    
    return 1;
}


int Stem::genSCFGtoHMM( int idx,  ScfgGram *curSG,  Bulge *curHG, char dirtype )
{
    int i=0, j=0;

    prod[idx].direct=dirtype;
    prod[idx].ltype=curSG->state;
    prod[idx].lid=curSG->UID;

    prod[idx].rtype=curHG->getFirstGramState( );
    prod[idx].rid=curHG->getFirstGramUID( );
    
    if( curSG->state=='B' )
    {
        prod[idx].category=KEEP_POS;
        if( dirtype == 'L')
            prod[idx].singlep = curSG->trprob[ SCFG_LFH ];
        else
            prod[idx].singlep = curSG->trprob[ SCFG_RGH ];
    }
    else if( curSG->state=='S' )
    {
        prod[idx].category=PAIR_ADD;
        for(i=0; i<5; i++ )
        {
            for(j=0; j<5; j++ )
            {
                if( dirtype == 'L')
                    prod[idx].pairp[i][j] = curSG->trprob[ SCFG_LFH ]+curSG->emprobS[i][j];
                else
                    prod[idx].pairp[i][j] = curSG->trprob[ SCFG_RGH ]+curSG->emprobS[i][j];
            }
        }
    } 
    else if( curSG->state=='L')
    {
        prod[idx].category=LEFT_INS;
        for(i=0; i<4; i++ )
            prod[idx].basep[i] = curSG->trprob[ SCFG_LFH ] + curSG->emprobI[SCFG_LEFT ][i];
    }
    else if( curSG->state=='R' )
    {
        prod[idx].category=RIGHT_INS;
        for(i=0; i<4; i++ )
            prod[idx].basep[i] = curSG->trprob[ SCFG_RGH ] + curSG->emprobI[SCFG_RIGHT ][i];
    }
    else if( curSG->state=='V' )
    {
        printf("Not allowed V in productions!\n");
        exit( 0 );
    }
    else
    {
        printf("Unknown grammar state in genSCFGtoHMM: %c\n", curSG->state);
        exit( 0 );
    }
    return 1;
}

int Stem::genSCFGtoSCFG( int idx,  ScfgGram *curSG,  ScfgGram *nextSG )
{
    int i=0, j=0, which=-1;
    
    //determine which transition prob should be used
    if( nextSG->state=='L' )
        which=SCFG_LFI;
    else if( nextSG->state=='R' )
        which=SCFG_RGI;
    else
        which=SCFG_MS;

    prod[idx].ltype=curSG->state;
    prod[idx].lid=curSG->UID;
    prod[idx].rtype=nextSG->state;
    prod[idx].rid=nextSG->UID;
    
    if( curSG->state=='B' )
    {
        prod[idx].category = KEEP_POS;
        prod[idx].singlep = curSG->trprob[which];
    }
    else if( curSG->state=='S' )
    {
        if( nextSG->state=='E')
            prod[idx].category = PAIR_TMN;
        else
            prod[idx].category = PAIR_ADD;
        for(i=0; i<5; i++ )
        {
            for(j=0; j<5; j++ )
            {
                 prod[idx].pairp[i][j] = curSG->trprob[which]+curSG->emprobS[i][j];
            }
        }
    }
    else if( curSG->state=='L' )
    {
        if( nextSG->state=='E')
            prod[idx].category = LSINGLE_TMN;
        else
            prod[idx].category = LEFT_INS;
        for(i=0; i<4; i++ )
        {
            prod[idx].basep[i] = curSG->trprob[which] + curSG->emprobI[ SCFG_LEFT ][i];
        }
    }
    else if( curSG->state=='R' )
    {
        if( nextSG->state=='E')
            prod[idx].category = RSINGLE_TMN;
        else
            prod[idx].category = RIGHT_INS;
        for(i=0; i<4; i++ )
        {
            prod[idx].basep[i] = curSG->trprob[which]+curSG->emprobI[ SCFG_RIGHT ][i];
        }
    }
    else if( curSG->state=='V' )
    {
        printf("Not allowed V in productions!\n");
        exit( 0 );
    }
    else
    {
        prod[idx].category = UNKNOWN_NT;
        printf("Unknown nonterminal in genSCFGtoSCFG: %c\n", curSG->state );
        exit( 0 );
    }
    return 1;
}

//functions to convert X-->V-->Y ==> X-->Y
int Stem::genSCFGtoHMMviaV( int idx,  ScfgGram *curSG,  Bulge *curHG, char dirtype, double vprob)
{
    int i=0, j=0;

    prod[idx].direct=dirtype;
    prod[idx].ltype=curSG->state;
    prod[idx].lid=curSG->UID;

    prod[idx].rtype=curHG->getFirstGramState( );
    prod[idx].rid=curHG->getFirstGramUID( );
    
    if( curSG->state=='B' || curSG->state=='S' )
    {
        prod[idx].category=PAIR_ADD;
        for(i=0; i<5; i++ )
            for(j=0; j<5; j++ )
                 prod[idx].pairp[i][j] = curSG->trprob[ SCFG_MS ]+curSG->emprobS[i][j]+ vprob;
    }
    else if( curSG->state=='L')
    {
        prod[idx].category=LEFT_INS;
        for(i=0; i<4; i++ )
            prod[idx].basep[i] = curSG->trprob[ SCFG_MS ] + curSG->emprobI[SCFG_LEFT ][i]+ vprob;
    }
    else if( curSG->state=='R' )
    {            
        printf("Error: Right insertion couldn't be directed to V state at genSCFGtoHMMviaV!\n");
        exit(0);
    }
    else if( curSG->state=='V' )
    {
        printf("Not allowed V in productions!\n");
        exit( 0 );
    }
    else
    {
        printf("Unknown grammar state for genSCFGtoHMM\n");
        exit( 0 );
    }
    return 1;
}

int Stem::genSCFGtoSCFGviaV( int idx,  ScfgGram *curSG,  ScfgGram *nextSG, double vprob )
{
    int i=0, j=0;
    prod[idx].ltype=curSG->state;
    prod[idx].lid=curSG->UID;
    prod[idx].rtype=nextSG->state;
    prod[idx].rid=nextSG->UID;
    
    if( curSG->state=='B' )
    {
        prod[idx].category = KEEP_POS;
        prod[idx].singlep = curSG->trprob[ SCFG_MS ] + vprob;
    }
    else if( curSG->state=='S' )
    {
        if( nextSG->state=='E')
            prod[idx].category = PAIR_TMN;
        else
            prod[idx].category = PAIR_ADD;
        for(i=0; i<5; i++ )
        {
            for(j=0; j<5; j++ )
             {
                 prod[idx].pairp[i][j] = curSG->trprob[ SCFG_MS ]+curSG->emprobS[i][j]+ vprob;
            }
         }
    }
    else if( curSG->state=='L' )
    {
        if( nextSG->state=='E')
            prod[idx].category = LSINGLE_TMN;
        else
            prod[idx].category = LEFT_INS;
        for(i=0; i<4; i++ )
        {
            prod[idx].basep[i] = curSG->trprob[ SCFG_MS ] + curSG->emprobI[ SCFG_LEFT ][i]+ vprob;
         }
    }
    else if( curSG->state=='R' )
    {
        printf("Error: Right insertion couldn't be directed to V state at genSCFGtoSCFGviaV!\n");
        exit(0);
    }
    else if( curSG->state=='V' )
    {
        printf("Not allowed V in productions!\n");
        exit( 0 );
    }
    else
    {
        prod[idx].category = UNKNOWN_NT;
        printf("Unknown nonterminal\n");
        exit( 0 );
    }
    return 1;
}

int Stem::genSCFGtoVProd( int idx,  ScfgGram *curSG,  ScfgGram *curVt )
{
    if( curSG->nextS != NULL && curSG->nextS->state=='V')
    {
        if( curVt->nextLH != NULL )
        {            
            printf("Error: nextLH of V state shouldn't exist!\n");
            exit(0);
        }
        if( curVt->nextLI != NULL )
        {            
            printf("Error: nextLI of V state shouldn't exist!\n");
            exit(0);
        }
        if( curVt->nextS != NULL )
        {
            genSCFGtoSCFGviaV( idx++,  curSG,  curVt->nextS, curVt->trprob[SCFG_MS] );
        }
        if( curVt->nextRI != NULL )
        {
             genSCFGtoSCFGviaV( idx++,  curSG,  curVt->nextRI, curVt->trprob[SCFG_RGI] );
        }
        if( curVt->nextRH != NULL )
        {
            genSCFGtoHMMviaV( idx++, curSG, curVt->nextRH, 'R', curVt->trprob[SCFG_RGH] );
        }
    }
    return idx;
}

int Stem::genHMMtoHMMviaV( int idx, Bulge *curHG, Bulge *nexHG, char dirtype, double vprob)
{
    prod[idx].category=KEEP_POS;
    prod[idx].direct=dirtype;
    prod[idx].ltype=curHG->getLastGramState( );
    prod[idx].lid=curHG->getLastGramUID( );
    prod[idx].rtype=nexHG->getFirstGramState( );
    prod[idx].rid=nexHG->getFirstGramUID( );
     prod[idx].singlep = curHG->trprob[ 1 ] + vprob;
    return 1;
}

int Stem::genHMMtoSCFGviaV( int idx, Bulge *curHG, ScfgGram *curSG, char dirtype, double vprob)
{
    prod[idx].category=KEEP_POS;
    prod[idx].direct=dirtype;
    prod[idx].ltype=curHG->getLastGramState( );
    prod[idx].lid=curHG->getLastGramUID( );
    prod[idx].rtype=curSG->state;
    prod[idx].rid=curSG->UID;
    prod[idx].singlep = curHG->trprob[ 1 ]+ vprob;
    return 1;
}

int Stem::genHMMtoVProd( int idx,  Bulge *curBg,  ScfgGram *curVt,  char dirtype)
{
    if( curBg->nextS != NULL && curBg->nextS->state=='V')
    {
        if( curVt->nextLH != NULL )
        {            
            printf("Error: nextLH of V state shouldn't exist!\n");
            exit(0);
        }
        if( curVt->nextLI != NULL )
        {            
            printf("Error: nextLI of V state shouldn't exist!\n");
            exit(0);
        }
        if( curVt->nextS != NULL )
        {
            genHMMtoSCFGviaV( idx++, curBg,  curVt->nextS,  dirtype, curVt->trprob[ SCFG_MS ] );
        }
        if( curVt->nextRI != NULL )
        {
            genHMMtoSCFGviaV( idx++, curBg,  curVt->nextRI,  dirtype, curVt->trprob[ SCFG_RGI ] );
        }
        if( curVt->nextRH != NULL )
        {
            genHMMtoHMMviaV( idx++, curBg, curVt->nextRH, dirtype, curVt->trprob[SCFG_RGH] );
        }
    }
    return idx;
}

int Stem::serialProducts( )
{
    #ifdef DEBUG_DEV
        printf("Stem::serialProducts\n");
    #endif
    
    int idx=0;
    ScfgGram *curSG=pStemGram;
    ScfgGram * curVt=NULL, *nexSG=NULL, *curLIG=NULL, *curRIG=NULL;
    Bulge *curBg=NULL;
    
    curVt=curSG->nextS;
    while( curVt!=NULL )
    {
        nexSG=curVt->nextS;
        //S-->V
        idx=genSCFGtoVProd( idx, curSG, curVt );
        
        curLIG=curSG->nextLI;
        if( curLIG!=NULL )
        {
            //S-->L (left insertion)
            genSCFGtoSCFG( idx++, curSG, curLIG);
          	
            curBg=curSG->nextLH;
            if( curBg!=NULL )
            {   
                //S-->Bulge
                genSCFGtoHMM( idx++, curSG, curBg, 'L' );
                while( 1 )
                {
                    genSCFGtoSCFG( idx++, curLIG, curLIG); //L-->xL
                    genSCFGtoHMM( idx++, curLIG, curBg, 'L');  //L-->Bulge
                    idx = curBg->serialProducts( prod, idx, 'L' );  //Bulge
                    curLIG=curBg->nextI;
                    genHMMtoSCFG( idx++, curBg, curLIG, 'L');   //Bulge-->xL
                    if( curBg->nextLp!=NULL )
                        genHMMtoHMM( idx++, curBg, curBg->nextLp, 'L'); //Bulge-->Bulge
                    else
                    {                      
                        idx=genHMMtoVProd( idx, curBg, curVt, 'L');  //Bulge-->V
                        break;
                    }
                    curBg=curBg->nextLp;
                }
          	}
          	//  L-->xL
            genSCFGtoSCFG( idx++, curLIG, curLIG);
            //  L-->xV
            idx=genSCFGtoVProd( idx, curLIG, curVt );
        }
        
        //***SCFG ******* S->Right
        curRIG=curVt->nextRI;
        if( curRIG!=NULL )
        {
            //S-->R
            curBg=curVt->nextRH;
            if( curBg!=NULL )
            {
                //S-->Bulge
                while( 1 )
                {
                    //R-->xR
                    genSCFGtoSCFG( idx++, curRIG, curRIG);
                    //R-->Bulge
                    genSCFGtoHMM( idx++, curRIG, curBg, 'R');
                    idx = curBg->serialProducts( prod, idx, 'R');	 //Bulge
                    curRIG=curBg->nextI;
                    genHMMtoSCFG( idx++, curBg, curRIG, 'R');     //Bulge-->xL
                    
                    if( curBg->nextLp!=NULL )
                        genHMMtoHMM( idx++, curBg, curBg->nextLp, 'R'); //Bulge-->Bulge
                    else
                    {
                        genHMMtoSCFG( idx++, curBg, nexSG, 'R');  //Bulge-->S
                        break;
                    }
                    curBg=curBg->nextLp;
                }
          	}
          	//  R-->xR
            genSCFGtoSCFG( idx++, curRIG, curRIG);
            //  R-->xS
            genSCFGtoSCFG( idx++, curRIG, nexSG);
        }
        curSG=nexSG;
        curVt=curSG->nextS;
    }

    if( idx!=numprod)
    {
        printf(" Wrong to calculate the number of productions: %d!=%d.\n", idx, numprod);
        exit( 0 );
    }
    else
        return 1;
}

int Stem::isProdsArray(  )
{
    if( logProb==1 && prod!=NULL )
        return 1;
    else
        return 0;
}

/////////
int Stem::getNumPairRegion(  )
{
    return regInfo->numOffset;
}

int Stem::getBothArmRegion( int no, int region[4], double coeff  )
{
    int i=0, status=0;
    status = regInfo->getStemRegion( no, region, coeff );
    //for( i=0; i<4; i++ )
    //    region[i]=(region[i]<0) ? 0 : region[i];
    return status;
}

int Stem::getMidlen( double coeff, int ismax )
{
    int retval=0;
    if( ismax )
        retval=(int)( regInfo->midE + coeff*regInfo->midSD );
    else
        retval=(int)( regInfo->midE - coeff*regInfo->midSD );
    retval = retval>0 ? retval : 0;
    return  retval;
}

int Stem::getStemSpan( double coeff, int ismax )
{
    int retval=0;
    double sum=0.0;
    if( ismax )
        sum = regInfo->armE[0] + regInfo->armE[1] + regInfo->midE + 
              coeff*(regInfo->armSD[0] + regInfo->armSD[1] + regInfo->midSD);
    else
        sum = regInfo->armE[0] + regInfo->armE[1] + regInfo->midE - 
              coeff*(regInfo->armSD[0] + regInfo->armSD[1] + regInfo->midSD);
    sum = sum>0 ? sum : 0;
    retval=(int)ceil(sum);
    return retval;
}

void Stem::setPriorPairFreq( double priorMatrix[5][5], double h )
{
    int i=0, j=0;
    double allpair=0;
    
    //check the value to be valid
    if( h>1 )
        priorcoef=1;
    else if( h<0)
        priorcoef=0;
    else
        priorcoef=h;
    
    if( ONLYCANONICAL )  //just calculate the canonical base pairs
    {
        allpair = priorMatrix[4][1] + priorMatrix[1][4] +
                  priorMatrix[2][3] + priorMatrix[3][2] + 
                  priorMatrix[4][3] + priorMatrix[3][4];
        priorPairFreq[4][1] = priorMatrix[4][1]/allpair;
        priorPairFreq[1][4] = priorMatrix[1][4]/allpair;
        priorPairFreq[4][3] = priorMatrix[4][3]/allpair;
        priorPairFreq[3][4] = priorMatrix[3][4]/allpair;
        priorPairFreq[2][3] = priorMatrix[2][3]/allpair;
        priorPairFreq[3][2] = priorMatrix[3][2]/allpair;
    }
    else    //calculate all of the canonical base pairs
    {
        for(i=0; i<5; i++ )
            for(j=0; j<5; j++ )
            {
                allpair += priorMatrix[i][j];
            }
        for(i=0; i<5; i++ )
            for(j=0; j<5; j++ )
            {
                priorPairFreq[i][j]=priorMatrix[i][j]/allpair;                    
            } 
    }
}

bool Stem::isGoodStem( int idx )
{
    #ifdef DEBUG_DEV
       printf("Stem::isGoodStem\n");
    #endif
    int i=0;
    int bpnum=0;
    int lenl= Pos[1] - Pos[0]+1;
    int lenr= Pos[3] - Pos[2]+1;
    int curl=0, curr=lenr-1;
    
    while( curl<lenl && curr>=0 )  {
        for( i=curl; i<lenl; i++ ) {
            if( armpasta[0][i] != '.' ) {
                curl=i;
                break;
            }
        }        
        for( i=curr; i>=0; i-- ) {
            if( armpasta[1][i] != '.' ) {
                curr=i;
                break;
            }
        }
        char lres = numtochar( armorg[0][idx][curl] );
        char rres = numtochar( armorg[1][idx][curr] );

        if( lres != GAP && rres != GAP ) // not x- or -x or --
        //if( isCanonicalPair( lres,  rres ) ) // canonical pair
            bpnum ++;
        curl++;
        curr--;
    }
    double pertagepair = (double)bpnum / stemlen;
    if( pertagepair >= PERTAGEPAIRVALUE && bpnum >= LEASTNUMPAIR )
        return true;
    else
        return false;
}

void Stem::pickConservedOnes(  )
{
    #ifdef DEBUG_DEV
       printf("Stem::pickConservedOnes\n");
    #endif
    int i=0, j=0, k=0, numGoodStem=0, setNumGoodstem=0;
    int lenl= Pos[1] - Pos[0]+1;
    int lenr= Pos[3] - Pos[2]+1;

    for( i=0; i<numorgs; i++ ) {
        bool curStatus = true;
        curStatus = isGoodStem( i );
        goodstem.push_back( curStatus );
        if( curStatus )
            numGoodStem++;
    }
    double gstemportion = (double)numGoodStem/numorgs;
    if( gstemportion > 0.5 )
    {
        setNumGoodstem = numorgs;
        for( i=0; i<numorgs; i++ )
            goodstem[i]=true;        
    }
    else
    {
        setNumGoodstem = numGoodStem;
    }
    
    //allocate the memory for seqs
    numseqs = setNumGoodstem;
    armseq[0]= new int*[ numseqs ];
    armseq[1]= new int*[ numseqs ];

    for( i=0, j=0; i<numorgs; i++ )
    {
        if( goodstem[i] ) {
            armseq[0][j]= new int[ lenl ];
            for(k=0; k<lenl; k++ )
                armseq[0][j][k]=armorg[0][i][k];
            
            armseq[1][j]= new int[ lenr ];
            for(k=0; k<lenr; k++ )
                armseq[1][j][k]=armorg[1][i][k];
            j++;
        }
    }
}

int Stem::generateProdsArray(  )
{
    #ifdef DEBUG_DEV
       printf("Stem::generateProdsArray\n");
    #endif
    int status=0;
    
    if( pStemGram==NULL )
        return 0;
    
    //1. logarithms of Probabilities
    if( logProb!=1 )
    {
        logProbforStem(  );
        logProb=1;  //set the class flag
    }

    //2. get the number of Nonterminals and Rules
    countNumAllNtms( );
    countNumAllRules( );
    nonterm = new NontermProp[ numntm ];
    initNontermArray( );
    prod = new GramRule[ numprod ];
    initGramProd( );

    //3. Array all nonterminals
    arrayNonterminal(  );

    //4. Array all productions
    status=serialProducts(  );
    if( status ==0 )
    {
        printf("Error: serialize products!\n");
        exit( 0 );
    }
    return 1;
}

void Stem::printNonterms( )
{
    int i=0;
    printf("\n==================== List of Nonterminals: %d =====================\n", numntm);
    for( i=0; i<numntm; i++ )
    {
        printf("%-3d: %c%d\n", nonterm[i].uid, nonterm[i].state, nonterm[i].id);
    }
}

void Stem::printAllRules(  )
{
    int no=0;
    printf("\n");
    printf("================== Production rules in array ===================\n");
    for( no=0; no<numprod; no++ )
    {
        printf("%d: <%d>  ", no+1, prod[no].category );
        if( prod[no].ltype=='B' )
        {
            printf("%c%d->[%c%d] (%lf)\n", 
                    prod[no].ltype, nonterm[ prod[no].lid ].id, 
                    prod[no].rtype, nonterm[ prod[no].rid ].id, prod[no].singlep );
        }
        else if( prod[no].ltype=='S' )
        {
            printf("%c%d->x[%c%d]y\n",
                    prod[no].ltype, nonterm[ prod[no].lid ].id,
                    prod[no].rtype, nonterm[ prod[no].rid ].id );
            printPairEmP( prod[no].pairp );
        }
        else if( prod[no].ltype=='V' )
        {
            printf("%c%d->[%c%d] (%lf)\n",
                    prod[no].ltype, nonterm[ prod[no].lid ].id,
                    prod[no].rtype, nonterm[ prod[no].rid ].id, prod[no].singlep );
        }
        else if( prod[no].ltype=='L' ) 
        {
            printf("%c%d->x[%c%d]\n",
                   prod[no].ltype, nonterm[ prod[no].lid ].id,
                   prod[no].rtype, nonterm[ prod[no].rid ].id );
            printBaseEmP( prod[no].basep );
        }
        else if( prod[no].ltype=='R' ) 
        {
            printf("%c%d->[%c%d]y\n",
                    prod[no].ltype, nonterm[ prod[no].lid ].id,
                    prod[no].rtype, nonterm[ prod[no].rid ].id );
            printBaseEmP( prod[no].basep );
        }
        else if( prod[no].ltype=='E' )
        {
            printf("%c%d->[%c%d] (%lf)\n",
                    prod[no].ltype, nonterm[ prod[no].lid ].id,
                    prod[no].rtype, nonterm[ prod[no].rid ].id, prod[no].singlep );
        }
        else if( prod[no].ltype=='M' )
        {
            printf("%c%d->x[%c%d]\n",
                    prod[no].ltype, nonterm[ prod[no].lid ].id,
                    prod[no].rtype, nonterm[ prod[no].rid ].id );
            printBaseEmP( prod[no].basep );
        }
        else if( prod[no].ltype=='I' ) 
        {
            printf("%c%d->x[%c%d]\n",
                    prod[no].ltype, nonterm[ prod[no].lid ].id,
                    prod[no].rtype, nonterm[ prod[no].rid ].id );
            printBaseEmP( prod[no].basep );
        }
        else if( prod[no].ltype=='D' ) 
        {
            printf("%c%d->[%c%d] (%lf)\n",
                    prod[no].ltype, nonterm[ prod[no].lid ].id,
                    prod[no].rtype, nonterm[ prod[no].rid ].id, prod[no].singlep );
        }
        else
        {
            printf("Unknown nonterminal\n");
        }
    }    
}

void Stem::printScfgGram(  )
{

    char sym[2];
    ScfgGram *curSG=pStemGram;
    ScfgGram *nexSG=NULL, *curLIG=NULL, *curRIG=NULL;
    Bulge *curBg=NULL;
    
    if( logProb )
        printf("=================== Base-e Log Prob =================\n");
    else
        printf("===================   Probability   =================\n");

    nexSG=curSG->nextS;
    while( nexSG!=NULL )
    {
        if( curSG->state=='B'||curSG->state=='V')
        {
            sym[0]=' ';
            sym[1]=' ';
        }
        else
        {
            sym[0]='x';
            sym[1]='y';
        }
        //***SCFG ******* S-->S
        printf("--------------------- %d -> %d ---------------------\n", curSG->ID, nexSG->ID );
         printf( "[%c%d]->%c[%c%d]%c    P=%f\n",  curSG->state, curSG->ID,
                 sym[0], nexSG->state, nexSG->ID, sym[1],     curSG->trprob[ SCFG_MS ] );
        if( curSG->state!='B' && curSG->state!='V')
            printPairEmP( curSG->emprobS );
        
        //***SCFG ******* S->Left
        curLIG=curSG->nextLI;
        if( curLIG!=NULL )
        {
            //  S-->L
          	printf( "[%c%d]->%c[%c%d]%c    P=%f\n", curSG->state, curSG->ID, 
                        sym[0], curLIG->state, curLIG->ID, sym[1],  curSG->trprob[ SCFG_LFI ] );
             if( curSG->state != 'B' && curSG->state!='V')
                printPairEmP( curSG->emprobS );
                
            curBg=curSG->nextLH;
            if( curBg!=NULL )
            {   
                //  S-->Bulge
                printf( "[%c%d]->%c[B%d]%c    P=%f\n", curSG->state, curSG->ID, 
                        sym[0], curSG->ID, sym[1],  curSG->trprob[ SCFG_LFH ] );
                 if( curSG->state != 'B' && curSG->state!='V')
                    printPairEmP( curSG->emprobS );
                    
                while( 1 )
                {
                    //L-->Bulge
                    printf( "[%c%d]->x[B%d]    P=%f\n",
                            curLIG->state, curLIG->ID, 
                            curLIG->ID,  curLIG->trprob[ SCFG_LFH ] );
                    printBaseEmP( curLIG->emprobI[ SCFG_LEFT ] );
                    //L-->xL
                    printf( "[%c%d]->x[%c%d]     P=%f\n",
                            curLIG->state,  curLIG->ID,
                            curLIG->nextLI->state,  curLIG->nextLI->ID,  curLIG->trprob[ SCFG_LFI ] );
                    printBaseEmP( curLIG->emprobI[ SCFG_LEFT ] );

                    curBg->printHMMGrams( 'L' );//Bulge HMM
                    
                    //Bulge-->xL
                    printf( "[E%d]->[%c%d]    P=%f\n", 
                            curLIG->ID,  curBg->nextI->state,
                            curBg->nextI->ID, curBg->trprob[ 0 ] );
                    
                    if( curBg->nextLp != NULL )
                    {
                        //Bulge-->Bulge
                        printf( "[E%d]->[B%d]    P=%f\n",
                                 curLIG->ID,  curLIG->ID,  curBg->trprob[ 1 ] );
                    }
                    else
                    {
                        //Bulge-->S
                        printf( "[E%d]->[%c%d]    P=%f\n", curLIG->ID, 
                                curBg->nextS->state,   curBg->nextS->ID,  curBg->trprob[ 1 ] ); 
                        curLIG=curBg->nextI;
                        break;
                    }
                    curLIG=curBg->nextI;
                    curBg=curBg->nextLp;
                }
          	}
          	//  L-->xL
            printf( "[%c%d]->x[%c%d]     P=%f\n",
                        curLIG->state,     curLIG->ID, 
                        curLIG->nextLI->state,     curLIG->nextLI->ID,
                        curLIG->trprob[ SCFG_LFI ] );
            printBaseEmP( curLIG->emprobI[ SCFG_LEFT ] );

            //  L-->xS
            printf( "[%c%d]->x[%c%d]     P=%f\n",
                        curLIG->state,     curLIG->ID, 
                        curLIG->nextS->state,      curLIG->nextS->ID,
                        curLIG->trprob[ SCFG_MS ] );
            printBaseEmP( curLIG->emprobI[ SCFG_LEFT ] );
  
        }
        //***SCFG ******* S->Right
        curRIG=curSG->nextRI;
        if( curRIG!=NULL )
        {
            //  S-->R
          	printf( "[%c%d]->%c[%c%d]%c    P=%f\n", curSG->state, curSG->ID, 
                    sym[0], curRIG->state, curRIG->ID, sym[1],  curSG->trprob[ SCFG_RGI ] );
             if( curSG->state!='B'&& curSG->state!='V' )
                printPairEmP( curSG->emprobS );
                
            curBg=curSG->nextRH;
            if( curBg!=NULL )
            {   
                //  S-->Bulge
                printf( "[%c%d]->%c[B%d]%c    P=%f\n", curSG->state, curSG->ID, 
                        sym[0], curSG->ID, sym[1],  curSG->trprob[ SCFG_RGH ] );
                 if( curSG->state!='B'&& curSG->state!='V' )
                    printPairEmP( curSG->emprobS );
                    
                while( 1 )
                {
                    //R-->Bulge
                    printf( "[%c%d]->[B%d]x    P=%f\n",
                            curRIG->state, curRIG->ID, 
                            curRIG->ID,  curRIG->trprob[ SCFG_RGH ] );
                    printBaseEmP( curRIG->emprobI[ SCFG_RIGHT ] );
                    //R-->Rx
                    printf( "[%c%d]->[%c%d]x     P=%f\n",
                            curRIG->state,  curRIG->ID,
                            curRIG->nextRI->state,  curRIG->nextRI->ID,  curRIG->trprob[ SCFG_RGI ] );
                    printBaseEmP( curRIG->emprobI[ SCFG_RIGHT ] );

                    curBg->printHMMGrams('R' );//Bulge HMM
                    
                    //Bulge-->R
                    printf( "[E%d]->[%c%d]    P=%f\n", 
                            curRIG->ID,
                            curBg->nextI->state, curBg->nextI->ID, curBg->trprob[ 0 ] );
                    
                    if( curBg->nextLp!=NULL )
                    {
                        //Bulge-->Bulge
                        printf( "[E%d]->[B%d]    P=%f\n",
                                 curRIG->ID,  curRIG->ID,  curBg->trprob[ 1 ] );
                    }
                    else
                    {
                        //Bulge-->S
                        printf( "[E%d]->[%c%d]    P=%f\n", curRIG->ID,
                                curBg->nextS->state,   curBg->nextS->ID,  curBg->trprob[ 1 ] );
                        curRIG=curBg->nextI;
                        break;
                    }
                    curRIG=curBg->nextI;
                    curBg=curBg->nextLp;
                }
          	}
          	//  R-->Rx
            printf( "[%c%d]->[%c%d]x     P=%f\n",
                      curRIG->state,     curRIG->ID,
                      curRIG->nextRI->state,     curRIG->nextRI->ID,
                      curRIG->trprob[ SCFG_RGI ] );
            printBaseEmP( curRIG->emprobI[ SCFG_RIGHT ] );
            //  R-->Sx
            printf( "[%c%d]->[%c%d]x     P=%f\n",
                      curRIG->state,     curRIG->ID, 
                      curRIG->nextS->state,      curRIG->nextS->ID,
                      curRIG->trprob[ SCFG_MS ] );
            printBaseEmP( curRIG->emprobI[ SCFG_RIGHT ] );
        }
        curSG=nexSG;
        nexSG=curSG->nextS;
    }
}

void Stem::freeFreqMem(  )
{
    int i=0, j=0;
    for( i=0; i<stemlen; i++ )
    {
        for( j=0; j<5; j++ )
        {
            if( colpairNum[i][j] != NULL)
            {
                delete [ ] colpairNum[i][j];
                colpairNum[i][j]=NULL;
            }
        }
    }
    for( i=0; i<stemlen; i++ )
    {
        if( colpairNum[i] != NULL)
        {
            delete [ ] colpairNum[i];
            colpairNum[i]=NULL;
        }
    }
    if( colpairNum !=NULL )
    {
        delete [ ] colpairNum;
        colpairNum = NULL;
    }
}
    
void Stem::allocateFreqMem(  )
{
    int i=0, j=0, k=0;
    colpairNum = new int** [stemlen];
    for( i=0; i<stemlen; i++ )
    {
        colpairNum[i] = new int* [5];
        for( j=0; j<5; j++ )
        {
            colpairNum[i][j] = new int [5];
        }
    }
    for( i=0; i<stemlen; i++ )
    {
        for( j=0; j<5; j++ )
        {
            for( k=0; k<5; k++ )
                colpairNum[i][j][k]=0;
        }
    }
}

void Stem::calcColumnPairFreq( int colno, int Lpos, int Rpos )
{
    int i=0;
    int x=0, y=0;
    for( i=0; i<numseqs; i++ )
    {
        x=armseq[0][i][Lpos];
        y=armseq[1][i][Rpos];
        colpairNum[colno][x][y]++;
    }
}

void Stem::countPairFreq(  )
{
    int i=0, j=0, k=0;
    int Lpos=0, Rpos=0;
    int Llen=0, Rlen=0;
    
    //allocate memory
    allocateFreqMem( );
    
    Llen=strlen( armpasta[0] );
    Rlen=strlen( armpasta[1] );
    Rpos=Rlen-1;
    
    for( i=0; i<stemlen; i++ )
    {
        for( j=Lpos; j<Llen; j++ )
        {
            if( toupper(armpasta[0][j])!= toupper(idchar) )
                continue;
            else
            {
                Lpos=j;
                break;
            }
        }
        for( j=Rpos; j>=0; j-- )
        {
            if( toupper(armpasta[1][j])!= toupper(idchar) )
                continue;
            else
            {
                Rpos=j;
                break;
            }
        }
        
        calcColumnPairFreq(  i, Lpos, Rpos );

        for( j=0; j<5; j++ )
            for( k=0; k<5; k++ )
                sumcolpair[j][k] += colpairNum[i][j][k];
        Lpos++;
        Rpos--;
    }
    calcStemStrength( );
}

void Stem::calcStemStrength(  )
{
    #ifdef DEBUG_DEV
        printf("Stem::calcStemStrength\n");
    #endif
    int i=0, j=0, k=0;
    double curstrength=0.0, localmax=0.0;
    for( i=0; i<stemlen; i++ )
    {
        localmax=0.0;
        for( j=0; j<5; j++ ) {
            for( k=0; k<5; k++ ) {
                if( colpairNum[i][j][k] > localmax )
                    localmax = colpairNum[i][j][k];
            }
        }
        curstrength += localmax;
    }
    stemStrength = curstrength;
}

double Stem::getStemStrength( )
{
    return stemStrength;
}

double Stem::getPortionGoodStem( )
{
    return (double)numseqs/numorgs;
}

void Stem::printPairFreq( bool isParallel)
{
    int i=0, j=0, k=0;
    int temphz[5][5];
    double Hzmatrix[5][5];
    
    printf("****************************************************\n" );
    printf("***** Frequency of Base Pair in Stem %s (UID=%2d)*****\n", idstr.data( ), UID);
    printf("Pasta:  %s    %s\n", armpasta[0], armpasta[1] );
    for( i=0; i<stemlen; i++ )
    {
        printf("No.%2d:\n", i+1 );
        
        for( j=0; j<5; j++ )
            for( k=0; k<5; k++ )
            {
                Hzmatrix[j][k]=0.0;
                temphz[j][k]=colpairNum[i][j][k];
            }
        
        for( j=0; j<5; j++ )
            for( k=0; k<5; k++ )
                Hzmatrix[j][k]=(double)colpairNum[i][j][k]/numseqs;
        
        //printf
        if( isParallel )
            printPairHz( temphz, Hzmatrix );
        else
        {
            printf("# of Pairs: \n");
            printPairEmP( temphz );
            printf("Percentage: \n");
            printPairEmP( Hzmatrix );
            printf("\n");
        }
    }
}

void Stem::printTrainInfo( )
{
    int i=0, j=0;
    int lenlarm=0, lenrarm=0;
    
    lenlarm=Pos[1]-Pos[0]+1;
    lenrarm=Pos[3]-Pos[2]+1;
        
    //print sequences
    printf("\n========================================================\n" );
    printf("==================== STEM No.%2d ( %s )====================\n", UID, idstr.data( ));
    printf("---> Training Stems(Distance %d):\n", distance );
    printf(" State:  %s   %s\n", state[0], state[1]);
     printf(" Pasta:  %s   %s\n", armpasta[0], armpasta[1]);
     for(i=0; i<numseqs; i++)
     {
         printf("    %3d  ", i+1);
         for( j=0; j<lenlarm; j++)
         {
             printf("%c", numtochar(armseq[0][i][j]));
         }
         printf("   ");
         for( j=0; j<lenrarm; j++)
         {
             printf("%c", numtochar(armseq[1][i][j]));
         }
         printf("\n");
     }
}

void Stem::printOrgSeqs( )
{
    int i=0, j=0;
    int lenlarm=0, lenrarm=0;
    
    lenlarm=Pos[1]-Pos[0]+1;
    lenrarm=Pos[3]-Pos[2]+1;
        
    //print sequences
    printf("\n========================================================\n" );
    printf("==================== STEM No.%2d ( %s )====================\n", UID, idstr.data( ));
    printf("---> Training Stems(Distance %d):\n", distance );
    printf(" State:  %s   %s\n", state[0], state[1]);
     printf(" Pasta:  %s   %s    Model%:%5.2f\n", armpasta[0], armpasta[1], getPortionGoodStem( ) );
     for(i=0; i<numorgs; i++)
     {
         printf("    %3d  ", i+1);
         for( j=0; j<lenlarm; j++)
         {
             printf("%c", numtochar(armorg[0][i][j]));
         }
         printf("   ");
         for( j=0; j<lenrarm; j++)
         {
             printf("%c", numtochar(armorg[1][i][j]));
         }
         printf("   %s", goodstem[i]?"Good":"Bad");
         printf("\n");
     }
}

//********************************************************************************
//****************        implementation of ScfgBuilder class     ******************
//********************************************************************************
ScfgBuilder::ScfgBuilder( )
{
    #ifdef DEBUG_DEV
    printf("ScfgBuilder::ScfgBuilder\n");
    #endif
    
    int i=0, j=0;
    
    allstems=NULL;
  	numstems=-1;
    sequences=NULL;
    stemchid=NULL;
  	pasta=NULL;
    seqlength=0;
    numseqs=-1;
    pseudocount = PSEUDOCNT;
    pseudocntgap = PSDCNTGAP;
    dl=NULL;
    priorcoef = 0.0;
    priorfn = NULL;
    for(i=0; i<5; i++ )
        for(j=0; j<5; j++ )
            sumPairFreq[i][j]=0;
}

ScfgBuilder::~ScfgBuilder(  )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::~ScfgBuilder\n");
    #endif
    
    int i=0;
    for(i=0; i<numseqs; i++)
    {
        if( sequences[i]!=NULL )
        {
            delete [ ] sequences[i];
            sequences[i]=NULL;
        }
    }
    if( sequences!=NULL )
    {
        delete [ ] sequences;
        sequences=NULL;
    }
    if( pasta!=NULL )
    {
        delete [ ] pasta;
        pasta=NULL;
    }
    if( stemchid!=NULL )
    {
        delete [ ] stemchid;
        stemchid=NULL;
    }
    if( priorfn != NULL)
    {
        delete [ ] priorfn;
        priorfn=NULL;
    }
}

void ScfgBuilder::freeAllStems( )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::freeAllStems\n");
    #endif

    if( allstems!=NULL)
    {
        delete [ ] allstems;
        allstems=NULL;
    }
}

//----------------------
void ScfgBuilder::printSeqFile( )
{
    #ifdef DEBUG_DEV
    printf("ScfgBuilder::printSeqFile\n");
    #endif
    
    int i=0,j=0;
    printf("\n");
    printf("Sequences( length: %d ):\n", seqlength );
    printf("      %s\n", pasta);
    for(i=0; i<numseqs; i++)
    {
        printf("%4d  ", i+1);
        for(j=0; j<seqlength; j++)
        {
            printf("%c", numtochar(sequences[i][j]));
        }
        printf("\n");
    }
}

void ScfgBuilder::printAllGrams( int trinfo )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::printAllGrams\n");
    #endif
    int no=0;
    for( no=0; no<numstems; no++)
    {
        if( trinfo!=0 )
             allstems[no].printTrainInfo( );
         else
    	       printf("====================== STEM No. %-2d======================\n", no);
        printf("\n---> SCFG Grammar:\n");
        allstems[no].printScfgGram( );
    }
}

void ScfgBuilder::printStemGram( char stemch, int trinfo )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::printStemGram\n");
    #endif

    int no=0;
    for( no=0; no<numstems; no++)
    {
        if( toupper(stemch) == toupper(allstems[no].idchar) ) 	
            break;
    }
    if( no>=numstems )
    {
        printf("Not find the stem. Please check the stem char.\n");
        exit(0);
    }
     if( trinfo!=0 )
         allstems[no].printTrainInfo( );
     printf("\n---> SCFG Grammar:\n");
    allstems[no].printScfgGram( );
}

void ScfgBuilder::printPairFreq( char stemch )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::printStemPairFreq\n");
    #endif

    int i=0;
    for( i=0; i<numstems; i++)
    {
        if( stemch!=' ' )
        {
            if( toupper(stemch) == toupper(allstems[i].idchar) )
            {
                allstems[i].printPairFreq(  );
                break;
            }
        }
        else
        {
            allstems[i].printPairFreq(  );
        }
    }
    if( i>=numstems && stemch!=' ')
    {
        printf("Not find the stem. Please check the stem char.\n");
        exit(0);
    }
}

void ScfgBuilder::printSumPairHz( bool isParallel )
{
    int j=0, k=0;
    int totalnum=0;
    int numpair=0;
    int numBase[5];
    double HzMatrix[5][5];
    
    for( j=0; j<5; j++ ){
       numBase[j]=0;
    }
    for( j=0; j<5; j++ )
        for( k=0; k<5; k++ )
        {
            totalnum += sumPairFreq[j][k];
            numBase[j] += sumPairFreq[j][k];
            numBase[k] += sumPairFreq[j][k];
            HzMatrix[j][k]=0.0;
        }
        
    for( j=0; j<5; j++ )
       for( k=0; k<5; k++ )
            HzMatrix[j][k]=(double)sumPairFreq[j][k]/totalnum;
    
    //print
    //base count
    printf("Single Base in Pairs:\n   -:%d  A:%d  C:%d  G:%d  U:%d / All:%d\n", 
            numBase[0], numBase[1], numBase[2], numBase[3], numBase[4], numBase[0]+numBase[1]+numBase[2]+numBase[3]+numBase[4] );
    //check the calculation
    for( j=0; j<numstems; j++ )
        numpair += allstems[j].stemlen;
    printf("Pair = %d(%d)\n", numseqs*numpair, totalnum  );
    
    if( isParallel )
        printPairHz( sumPairFreq, HzMatrix );
    else
    {
        printf("# of Pairs: \n");
        printPairEmP( sumPairFreq );
        printf("Percentage: \n");
        printPairEmP( HzMatrix );
        printf("\n");
    }
}

void ScfgBuilder::setPseudocnt( double pdcnt, double gapcnt  )
{
    pseudocount = pdcnt;
    pseudocntgap = gapcnt;
}

void ScfgBuilder::setPriorCoef( double h )
{
    if( h>1 )
        priorcoef=1;
    else if( h<0)
        priorcoef=0;
    else
        priorcoef=h;
}

Stem * ScfgBuilder::getAllStems( int &num  )
{
    num=numstems;
    return allstems;
}

void ScfgBuilder::seperateStem(  )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::seperateStem\n");
    #endif
    
    int i=0, j=0, k=0;
    int lenl=0, lenr=0;
    int ls=0, le=0, rs=0, re=0;
    int seqsnum = numseqs;

    for( i=0; i<numstems; i++ )
    {
        ls=allstems[i].Pos[0];
        le=allstems[i].Pos[1];
        rs=allstems[i].Pos[2];
        re=allstems[i].Pos[3];
        
        lenl=le-ls+1;
        lenr=re-rs+1;
        
        //allocate the memory for pasta
        allstems[i].armpasta[0]= new char[lenl+1];
        for(j=0; j<lenl; j++)
        {
            allstems[i].armpasta[0][j]=toupper( pasta[ls+j] );
        }
        allstems[i].armpasta[0][j]='\0';

        allstems[i].armpasta[1]= new char[lenr+1];
        for(j=0; j<lenr; j++)
        {
            allstems[i].armpasta[1][j]=toupper( pasta[rs+j] );
        }
        allstems[i].armpasta[1][j]='\0';
        #ifdef DEBUG_DEV
            printf("Left:%s,  Right:%s\n", allstems[i].armpasta[0], allstems[i].armpasta[1]);
        #endif
        //allocate the memory for seqs
        allstems[i].armorg[0]= new int*[ seqsnum ];
        allstems[i].armorg[1]= new int*[ seqsnum ];
        for( j=0; j<seqsnum; j++ )
        {
            allstems[i].armorg[0][j]= new int[ lenl ];
            for(k=0; k<lenl; k++ )
            {
                allstems[i].armorg[0][j][k]=sequences[j][ls+k];
            }
            
            allstems[i].armorg[1][j]= new int[ lenr ];
            for(k=0; k<lenr; k++ )
            {
                allstems[i].armorg[1][j][k]=sequences[j][rs+k];
            }
        }
    }
}

//--------------------------------
char ScfgBuilder::determPairState(Stem *curstem, int lidx, int ridx )
{
    #ifdef DEBUG_DEV
    printf("ScfgBuilder::determPairState\n");
    #endif
    
    int i=0;
    int nongappair=0;
    int curseqnum=curstem->numseqs;

    for( i=0; i<curseqnum; i++ )
    {
        if( curstem->armseq[0][i][lidx]!=0 && curstem->armseq[1][i][ridx]!=0 )
        {
            nongappair++;
        }
    }
    #ifdef DEBUG_DEV
        printf("   BasePair: %d, GapPair: %d in %d seqs.\n", nongappair, curseqnum-nongappair, curseqnum);
    #endif
    if( nongappair>(curseqnum-nongappair))
        return 'S';
    else
        return 'I';
}

//---------------------- determine S or I --------------------
int ScfgBuilder::obtainNextS( Stem *curstem, int begl, int begr, int &distl, int &distr)
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::obtainNextS\n");
    #endif
    
    int status=-1;
    char idchar=curstem->idchar;
    char type=' ';
    char ntl=' ', ntr=' ';
    int lenl=curstem->Pos[1]-curstem->Pos[0]+1;
    int lenr=curstem->Pos[3]-curstem->Pos[2]+1;
    int endl=begl, endr=begr;

    if( endl>=lenl || endr>=lenr)
    {
        return 0;
    }
    
    while(1 )
    {
        ntl=toupper(curstem->armpasta[0][endl]);
        ntr=toupper(curstem->armpasta[1][endr]);
        
        if( idchar==ntl && idchar==ntr )
        {
            type=determPairState(curstem, endl, endr);
            #ifdef DEBUG_DEV
                printf("   =>Type: %c\n", type);
            #endif
            if(type=='S')
            {
                status=1;
                break;
            }
            else
            {
                endl+=1;
                endr-=1;
            }
        }
        else if( idchar!=ntl && idchar==ntr )
        {
            endl+=1;
        }
        else if( idchar==ntl && idchar!=ntr )
        {
            endr-=1;
        }
        else  //ntl!=idchar && ntr!=idchar
        {
            endl+=1;
            endr-=1;
        }
        
        if( endl>=lenl || endr>=lenr)
        {
            status=0;
            break;
        }
    }
    distl=endl-begl;
    distr=begr-endr;
    
    return status;    
}

//calculate the probability of S-->V and V-->next S
void ScfgBuilder::calcSTransP( ScfgGram *pSG,   ScfgGram *pVT,   Stem *curstem,
                               int *bulge,      BulgeInfo *leftB,   BulgeInfo *rightB, 
                               int begl,  int begr,   int distl,  int distr )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::calcSTransP\n");
    #endif
    
    int i=0, j=0;
    int insl=0, insr=0;
    int sl=0, sr=0;
    int cntss=0, cntsl=0, cntsr=0, cntslh=0, cntsrh=0;
    int curseqnum=curstem->numseqs;
    double pdocnt=pseudocount;
    double prob[ SCFG_NUM_NEXTP ];

    for( i=0; i<SCFG_NUM_NEXTP; i++ )
        prob[i]=0.0;
    
    //S->S
    for( i=0; i<curseqnum; i++ )
    {
        //printf("**** i=%d ****\n", i);
        insl=0;
        insr=0;
        sl=0;
        sr=0;
        for( j=0; j<distl; j++ )
        {
            //printf("left:  j=%d / %d :::  pos=%d\n", j, distl, begl+j );
            if( curstem->armseq[ SCFG_LEFT ][i][begl+j]!=0 )
                insl=1;

            if( bulge[0]>0 )
            {
                if( insl && j<leftB[0].beg )
                    sl=1;
            }
            
        }
        if( sl )
            cntslh++;
        if( insl )
            cntsl++;
        
        for( j=0; j<distr; j++ )
        {
            //printf("right: j=%d / %d :::  pos=%d\n", j, distr, begr-j );
            if( curstem->armseq[ SCFG_RIGHT ][i][begr-j]!=0 )
                insr=1;
            if( bulge[1]>0 )
            {
                if(insr && j<rightB[0].beg )
                    sr=1;
            }
        }
        if( sr )
            cntsrh++;
        if( insr )
            cntsr++;
        
        if( !insl && !insr)
            cntss++;
  	}
    //S-->S
    /*prob[SCFG_MS]= (double)(cntss+pdocnt)/(curseqnum+2*pdocnt);
    tempprob=(1-prob[SCFG_MS])*(cntsl+pdocnt)/( cntsr + cntsl + 2*pdocnt );
    if( bulge[0]>0 )
    {
        prob[SCFG_LFI]= tempprob*(cntslh+pdocnt)/(curseqnum+2*pdocnt);
        prob[SCFG_LFH]= tempprob*(curseqnum-cntslh+pdocnt)/(curseqnum+2*pdocnt);
    }
    else
    {
        prob[SCFG_LFI]= tempprob;
        prob[SCFG_LFH]= 0.0;
    }
    tempprob=(1-prob[SCFG_MS])*(cntsr+pdocnt)/( cntsr + cntsl + 2*pdocnt );
    if( bulge[1]>0 )
    {
        prob[SCFG_RGI]= tempprob*(cntsrh+pdocnt)/(curseqnum+2*pdocnt);
        prob[SCFG_RGH]= tempprob*(curseqnum-cntsrh+pdocnt)/(curseqnum+2*pdocnt);
    }
    else
    {
        prob[SCFG_RGI]= tempprob;
        prob[SCFG_RGH]= 0.0;
    }
    for( i=0; i<SCFG_NUM_NEXTP; i++ )
    {
        //printf("prob[%d]=%f\n", i, prob[i] );
        pSG->trprob[i]=prob[i];
    }
    */
    
    //S-->V
    if( bulge[0]>0 )
    {
        pSG->trprob[SCFG_MS]= (double)pdocnt/(curseqnum+3*pdocnt);
        pSG->trprob[SCFG_LFI]= (double)(cntslh+pdocnt)/(curseqnum+3*pdocnt);
        pSG->trprob[SCFG_LFH]= (double)(curseqnum-cntslh+pdocnt)/(curseqnum+3*pdocnt);
    }
    else
    {
        pSG->trprob[SCFG_MS]= (double)(curseqnum-cntsl+pdocnt)/(curseqnum+2*pdocnt);
        pSG->trprob[SCFG_LFI]= (double)(cntsl+pdocnt)/(curseqnum+2*pdocnt);
        pSG->trprob[SCFG_LFH]= 0.0;
    }
    //V-->S
    if( bulge[1]>0 )
    {
        pVT->trprob[SCFG_MS]= (double)pdocnt/(curseqnum+3*pdocnt);
        pVT->trprob[SCFG_RGI]= (double)(cntsrh+pdocnt)/(curseqnum+3*pdocnt);
        pVT->trprob[SCFG_RGH]= (double)(curseqnum-cntsrh+pdocnt)/(curseqnum+3*pdocnt);
    }
    else
    {
        pVT->trprob[SCFG_MS]= (double)(curseqnum-cntsr+pdocnt)/(curseqnum+2*pdocnt);
        pVT->trprob[SCFG_RGI]= (double)(cntsr+pdocnt)/(curseqnum+2*pdocnt);
        pVT->trprob[SCFG_RGH]= 0.0;
    }
}

void ScfgBuilder::postprecssBothBulge( ScfgGram *preG, ScfgGram *curV )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::postprecssBothBulge\n");
    #endif
    
    double ssprob=0.0, stoIns=0.0, stoBge=0.0;
    
    ssprob=preG->trprob[SCFG_MS];
    stoIns=preG->trprob[SCFG_LFI];
    stoBge=preG->trprob[SCFG_LFH];
    preG->trprob[SCFG_LFI]=(1-ssprob)*stoIns/(stoIns + stoBge);
    preG->trprob[SCFG_LFH]=(1-ssprob)*stoBge/(stoIns + stoBge);
    
    stoIns=preG->trprob[SCFG_RGI];
    stoBge=preG->trprob[SCFG_RGH];
    curV->trprob[SCFG_RGI]=stoIns/(stoIns + stoBge);
    curV->trprob[SCFG_RGH]=stoBge/(stoIns + stoBge);
    preG->trprob[SCFG_RGI]=0.0;
    preG->trprob[SCFG_RGH]=0.0;
}

void ScfgBuilder::calcITransPLeft( Stem *pstem,  int beg,  int dist,  double *trprb)
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::calcITransPLeft\n");
    #endif
    
    int i=0, j=0, k=0;
    int isbase=0;
    double pdocnt=pseudocount;
    double iibase=0.0;
    int curseqnum = pstem->numseqs;
    
    for( i=0; i<curseqnum; i++ )    //I->S/H
    {
        for( j=0; j<dist; j++ )
	 	     {
             if( pstem->armseq[ SCFG_LEFT ][i][beg+j]!=0 )
             {
                 isbase++;
        	       break;
            }
        }
    }
    for( i=0; i<curseqnum; i++ )  	//I->I
	 	{
        for( j=0; j<dist-1; j++ )
        {
    	      if( pstem->armseq[ SCFG_LEFT ][i][beg+j]!=0 )
    	      {
    	          for( k=j+1; k<dist; k++ )
    	          {
    	              if( pstem->armseq[ SCFG_LEFT ][i][beg+k]!=0 )
    	              {
    	                  iibase++;
    	                  break;
    	              }
    	          }
    	       }
        }
    }
    trprb[ 0 ]=(isbase+pdocnt)/(iibase+isbase+2*pdocnt); //I->S
    trprb[ 1 ]=(iibase+pdocnt)/(iibase+isbase+2*pdocnt); //I->I
}

void ScfgBuilder::calcITransPRight( Stem *pstem,  int beg,  int dist,  double *trprb)
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::calcITransPRight\n");
    #endif
    
    int i=0, j=0, k=0;
    int isbase=0;
    double pdocnt=pseudocount;
    double iibase=0.0;
    int curseqnum = pstem->numseqs;
    
    for( i=0; i<curseqnum; i++ )    //I->S/H
    {
        for( j=0; j<dist; j++ )
	 	     {
             if( pstem->armseq[ SCFG_LEFT ][i][beg-j]!=0 )
             {
                 isbase++;
        	       break;
            }
        }
    }
    for( i=0; i<curseqnum; i++ )  	//I->I
	 	{
        for( j=0; j<dist-1; j++ )
        {
    	      if( pstem->armseq[ SCFG_LEFT ][i][beg-j]!=0 )
    	      {
    	          for( k=j+1; k<dist; k++ )
    	          {
    	              if( pstem->armseq[ SCFG_LEFT ][i][beg-k]!=0 )
    	              {
    	                  iibase++;
    	                  break;
    	              }
    	          }
    	       }
        }
    }
    trprb[ 0 ]=(isbase+pdocnt)/(iibase+isbase+2*pdocnt); //I->S
    trprb[ 1 ]=(iibase+pdocnt)/(iibase+isbase+2*pdocnt); //I->I
}

void ScfgBuilder::calcPairEmProb( ScfgGram *pSG,  Stem *pstem,  int lpos,  int rpos)
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::calcPairEmProb\n");
    #endif
    
    int i=0, j=0;
    int lidx=0, ridx=0;
    int curseqnum=pstem->numseqs;
    double pdocnt=pseudocount;
    double pdocntgap=pseudocntgap;
    double pair[5][5];

    for( i=0; i<5; i++ )
    {
        for( j=0; j<5; j++ )
        {
            pair[i][j]=0.0;
        }
    }
    //count pairs
    for( i=0; i<curseqnum; i++ )
    {
        lidx=pstem->armseq[ SCFG_LEFT ][i][lpos];
        ridx=pstem->armseq[ SCFG_RIGHT ][i][rpos];
        pair[lidx][ridx]++;
    }
    
    /*
    for( i=0; i<5; i++ )
    {
        for( j=0; j<5; j++ )
        {
            pSG->emprobS[i][j]=(pair[i][j]+pdocnt)/(curseqnum+25*pdocnt);
        }
    }
    */
    for( i=0; i<5; i++ )
    {
        for( j=0; j<5; j++ )
        {
            if( i==0 || j==0 )
                pSG->emprobS[i][j]=(pair[i][j]+pdocntgap)/(curseqnum + 16*pdocnt + 9*pdocntgap);
            else
                pSG->emprobS[i][j]=(pair[i][j]+pdocnt)/(curseqnum + 16*pdocnt + 9*pdocntgap);
        }
    }
}

void ScfgBuilder::calcIBaseEmProbLeft( ScfgGram *pSG, Stem *pstem, int beg,  int dist )
{
    #ifdef DEBUG_DEV
    printf("ScfgBuilder::calcIBaseEmProbLeft\n");
    #endif
    
    int i=0, j=0, idx=0;
    int curseqnum=pstem->numseqs;
    double sum=0;
    double pdocnt=pseudocount;
    double base[5]={ 0.0,  0.0,  0.0,  0.0,  0.0 };

    for( i=0; i<curseqnum; i++ )
        for( j=0; j<dist; j++ )
        {
            idx=pstem->armseq[ SCFG_LEFT ][i][beg+j];
            base[ idx ]++;
        }
    for( i=1; i<5; i++ )
        sum+=base[i];
    for( i=1; i<5; i++ )
        pSG->emprobI[ SCFG_LEFT ][i-1]=(double)( base[i]+pdocnt)/(sum+4*pdocnt);
}

void ScfgBuilder::calcIBaseEmProbRight( ScfgGram *pSG, Stem *pstem, int beg,  int dist )
{
    #ifdef DEBUG_DEV
    printf("ScfgBuilder::calcIBaseEmProbRight\n");
    #endif
    
    int i=0, j=0, idx=0;
    int curseqnum=pstem->numseqs;
    double sum=0;
    double pdocnt=pseudocount;
    double base[5]={ 0.0,  0.0,  0.0,  0.0,  0.0 };

    for( i=0; i<curseqnum; i++ )
        for( j=0; j<dist; j++ )
        {
            idx=pstem->armseq[ SCFG_RIGHT ][i][beg-j];
            base[ idx ]++;
        }
    for( i=1; i<5; i++ )
        sum+=base[i];
    for( i=1; i<5; i++ )
        pSG->emprobI[ SCFG_RIGHT ][i-1]=(double)( base[i]+pdocnt)/(sum+4*pdocnt);
}

void ScfgBuilder::markState( Stem *curstem, int begl, int begr, int distl, int distr)
{
    #ifdef DEBUG_DEV
    printf("ScfgBuilder::markState\n");
    #endif
    
    int i=0;
    int lenlarm=curstem->Pos[1]-curstem->Pos[0]+1;

    if( begl+distl < lenlarm)
    {
        curstem->state[ SCFG_LEFT ][begl+distl]='s';
        curstem->state[ SCFG_RIGHT ][begr-distr]='s';
    }
    for( i=0; i<distl; i++)
    {
        if(curstem->armpasta[ SCFG_LEFT ][begl+i] == curstem->idchar)
        {
            curstem->state[ SCFG_LEFT ][begl+i]='i';
        }
        else
        {
            curstem->state[ SCFG_LEFT ][begl+i]=curstem->armpasta[ SCFG_LEFT ][begl+i];
        }
    }
    for( i=0; i<distr; i++)
    {
        if(curstem->armpasta[ SCFG_RIGHT ][begr-i] == curstem->idchar)
        {
            curstem->state[ SCFG_RIGHT ][begr-i]='i';
        }
        else
        {
            curstem->state[ SCFG_RIGHT ][begr-i]=curstem->armpasta[ SCFG_RIGHT ][begr-i];
        }
    }
}

void ScfgBuilder::countBulge( int *bulge, Stem *pstem, int begl,  int begr,  int distl,  int distr )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::countBulge\n");
    #endif
    
    int i=0;
    int bpos=0;
    for(i=0; i<distl; i++)
    {
        if( (pstem->armpasta[ SCFG_LEFT ][begl+i] != pstem->idchar )&& (!bpos) )
        {
            bulge[0]++;
            bpos=1;
        }
        if( bpos )
        {
            if( pstem->armpasta[ SCFG_LEFT ][begl+i] == pstem->idchar )
            {
                bpos=0;
            }
        }
    }
    bpos=0;
    for(i=0; i<distr; i++)
    {
        if( (pstem->armpasta[ SCFG_RIGHT ][begr-i] != pstem->idchar )&& (!bpos) )
        {
            bulge[1]++;
            bpos=1;
        }
        if( bpos )
        {
            if( pstem->armpasta[ SCFG_RIGHT ][begr-i] == pstem->idchar )
            {
                bpos=0;
            }
        }
    }
}

void ScfgBuilder::calcLoopTransPLeft(  Stem *pstem, int beg, int dist, double &lpi, double &lplp )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::calcLoopTransPLeft\n");
    #endif
    int i=0, j=0;
    double ss=0, si=0;
    int curseqnum=pstem->numseqs;
    double pdocnt=pseudocount;
    
    for( i=0; i<curseqnum; i++ )
    {
        for( j=0; j<dist; j++ )
        {
            if( pstem->armseq[ SCFG_LEFT ][i][beg+j] != 0  )
            {
                si++;
                break;
            }
        }
    }
    ss=curseqnum-si;
    lpi=(si+pdocnt)/(curseqnum+2*pdocnt);
    lplp=(ss+pdocnt)/(curseqnum+2*pdocnt);
}
            
void ScfgBuilder::calcLoopTransPRight(  Stem *pstem, int beg, int dist, double &lpi, double &lplp )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::calcLoopTransPRight\n");
    #endif
    int i=0, j=0;
    double ss=0, si=0;
    int curseqnum=pstem->numseqs;
    double pdocnt=pseudocount;
    
    for( i=0; i<curseqnum; i++ )
    {
        for( j=0; j<dist; j++ )
        {
            if( pstem->armseq[ SCFG_RIGHT ][i][beg-j] != 0  )
            {
                si++;
                break;
            }
        }
    }
    ss=curseqnum-si;
    lpi=(si+pdocnt)/(curseqnum+2*pdocnt);
    lplp=(ss+pdocnt)/(curseqnum+2*pdocnt);
}

//ret value: 0: nothing;   1: left;    2: right;  3: both
// int *bulge: the number of bulge[0] in left arm and bulge[1] in right arm
//      leftB: the position info for left bulges
//     rightB: the position info for right bulges, 
int ScfgBuilder::judgeBulge( int *bulge,    BulgeInfo* &leftB,    BulgeInfo* &rightB,
                             Stem *pstem,   int begl,  int begr,  int distl,  int distr )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::judgeBulge\n");
    #endif
    int i=0, j=0;
    int fpos=0;
    
    countBulge( bulge, pstem,  begl,  begr,  distl,  distr );
    //left arm
    if( bulge[0]>0 )
    {
        leftB = new BulgeInfo[ bulge[0] ];
        for(i=0; i<distl; i++)
        {
            if( (pstem->armpasta[ SCFG_LEFT ][begl+i] != pstem->idchar )&& (!fpos) )
            {
                leftB[j].beg=i;
                leftB[j].LorR='L';
                fpos=1;
            }
            if( fpos )
            {
                if( pstem->armpasta[ SCFG_LEFT ][begl+i] == pstem->idchar  )
                {
                    leftB[j].end=i-1;
                    j++;
                    fpos=0;
                }
                else if( i==distl-1 )
                {
                    leftB[j].end=i;
                    j++;
                    fpos=0;
                }
            }
        }
    }
    //right arm
    j=0;
    fpos=0;
    if( bulge[1]>0 )
    {
        rightB = new BulgeInfo[ bulge[1] ];
        for(i=0; i<distr; i++)
        {
            if( (pstem->armpasta[ SCFG_RIGHT ][begr-i] != pstem->idchar )&& (!fpos) )
            {
                rightB[j].beg=i;
                rightB[j].LorR='R';
                fpos=1;
            }
            if( fpos )
            {
                if( pstem->armpasta[ SCFG_RIGHT ][begr-i] == pstem->idchar )
                {
                    rightB[j].end=i-1;
                     j++;
                    fpos=0;
                }
                else if( i==distr-1 )
                {
                    rightB[j].end=i;
                    j++;
                    fpos=0;
                }
            }
        }
    }

    //for( i=0; i<bulge[0]; i++ )
    //    printf("Left: bulge%d  %d--%d\n", i, leftB[i].beg, leftB[i].end );
    //for( i=0; i<bulge[1]; i++ )
    //    printf("Right: bulge%d  %d--%d\n", i, rightB[i].beg, rightB[i].end );

    if( bulge[0]==0 && bulge[1]==0 )
        return 0;
    else if( bulge[0]>0 && bulge[1]==0 )
        return 1;
    else if( bulge[0]==0 && bulge[1]>0 )
        return 2;
    else
        return 3;
}

int ScfgBuilder::modelLeftArm(int numbulge, BulgeInfo *leftB, Stem *curstem,
                      ScfgGram *preG, ScfgGram *curG,   int idl, int disl, int idx, int uid )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::modelLeftArm\n");
    #endif
    
    int i=0;
    int bulgelen=0;
    int curp=0, dist=0;
    double trprob[2];
    double lpi=0.0, lplp=0.0;
    
    ScfgGram *pIns=NULL;
    Bulge *pLoop=NULL, *preLoop=NULL;
    HMMBuilder pHB;
    if( numbulge == 0)
    {
        pIns = new ScfgGram( 'L',  idx, uid++ );
        preG->nextLI=pIns;
        pIns->nextLI=pIns;
        pIns->nextS=curG;
        calcITransPLeft( curstem, idl, disl, trprob );
        pIns->trprob[SCFG_MS]=trprob[0];
        pIns->trprob[SCFG_LFI]=trprob[1];
        calcIBaseEmProbLeft( pIns, curstem, idl, disl);
    }
    else
    {
        curp=idl;
        for( i=0; i<numbulge; i++ )
        {
            pIns = new ScfgGram( 'L',  idx, uid++, i);
            if( i==0 )
                dist=leftB[i].beg;  // beg:0,1,2,...;
            else
                dist=leftB[i].beg - leftB[i-1].end - 1;
            calcITransPLeft( curstem, curp, dist, trprob );
            calcIBaseEmProbLeft( pIns, curstem, curp, dist);
            pIns->trprob[ SCFG_LFH ]=trprob[0];
            pIns->trprob[ SCFG_LFI ]=trprob[1];
            
            bulgelen=leftB[i].end - leftB[i].beg + 1;
            pLoop = pHB.buildBulge( curstem->armseq[SCFG_LEFT], curstem->numseqs,
                                    idl+leftB[i].beg, bulgelen,  'L',  uid,  i );
            curp=idl+leftB[i].end+1;
            if( i<numbulge-1)
                dist=leftB[i+1].beg-leftB[i].end-1;
            else
                dist=disl-leftB[i].end-1;
            calcLoopTransPLeft( curstem, curp, dist, lpi, lplp);
            pLoop->trprob[0]=lpi;
            pLoop->trprob[1]=lplp;
            
            if( i==0 )
            {
                preG->nextLI=pIns;
                preG->nextLH=pLoop;
            }
            else
            {
                preLoop->nextI=pIns;
                preLoop->nextLp=pLoop;
            }
            pIns->nextLH=pLoop;
            pIns->nextLI=pIns;
            preLoop=pLoop;
        }
    
        pIns = new ScfgGram( 'L', idx, uid++, i);
        curp=idl+leftB[numbulge-1].end+1;
        dist=disl-leftB[ numbulge-1 ].end-1;
        calcITransPLeft( curstem, curp, dist, trprob );
        pIns->trprob[ SCFG_MS ]=trprob[0];
        pIns->trprob[ SCFG_LFI ]=trprob[1];
        calcIBaseEmProbLeft( pIns, curstem, curp, dist);
        
        pLoop->nextI=pIns;
        pLoop->nextS=curG;
        pIns->nextS=curG;
        pIns->nextLI=pIns;
    }
    return uid;
}

int ScfgBuilder::modelRightArm( int numbulge,  BulgeInfo *rightB,  Stem *curstem,
                        ScfgGram *preG, ScfgGram *curG,    int idr, int disr, int idx, int uid )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::modelRightArm\n");
    #endif
    
    int i=0;
    int curp=0, dist=0;
    double trprob[2];
    double lpi=0.0, lplp=0.0;
    int bulgelen=0;
    
    ScfgGram *pIns=NULL;
    Bulge *pLoop=NULL, *preLoop=NULL;
    HMMBuilder pHB;
    
    if( numbulge == 0)
    {
        pIns = new ScfgGram( 'R',  idx, uid++);
        preG->nextRI=pIns;
        pIns->nextRI=pIns;
        pIns->nextS=curG;
        calcITransPRight( curstem, idr, disr, trprob );
        pIns->trprob[ SCFG_MS ]=trprob[0];
        pIns->trprob[ SCFG_RGI ]=trprob[1];
        calcIBaseEmProbRight( pIns, curstem, idr, disr);
    }
    else
    {
        curp=idr;
        for( i=0; i<numbulge; i++ )
        {
            pIns = new ScfgGram( 'R',  idx,  uid++, i);
            if( i==0 )
                dist=rightB[i].beg;  // beg:0,1,2,...;
            else
                dist=rightB[i].beg - rightB[i-1].end - 1;
                
             calcITransPRight( curstem, curp, dist, trprob );
            calcIBaseEmProbRight( pIns, curstem, curp, dist);
            pIns->trprob[ SCFG_RGH ]=trprob[0];
            pIns->trprob[ SCFG_RGI ]=trprob[1];
        
            bulgelen = rightB[i].end - rightB[i].beg + 1;
            pLoop = pHB.buildBulge( curstem->armseq[SCFG_RIGHT], curstem->numseqs,
                                    idr-rightB[i].beg,  bulgelen, 'R',  uid, i );
            curp=idr-rightB[i].end-1;
            if( i<numbulge-1 )
                 dist=rightB[i+1].beg - rightB[i].end - 1;
            else
                 dist=disr-rightB[i].end-1;
            calcLoopTransPRight( curstem, curp, dist, lpi, lplp);
            pLoop->trprob[0]=lpi;   //Bulge...Insert
            pLoop->trprob[1]=lplp;  //Bulge...bulge
            
            if( i==0 )
            {
                preG->nextRI=pIns;
                preG->nextRH=pLoop;
            }
            else
            {
                preLoop->nextI=pIns;
                preLoop->nextLp=pLoop;
            }
            pIns->nextRH=pLoop;
            pIns->nextRI=pIns;
            preLoop=pLoop;
        }
        pIns = new ScfgGram( 'R',  idx, uid++, i);
        curp=idr-rightB[numbulge-1].end-1;
        dist=disr-rightB[ numbulge-1 ].end-1;
        calcITransPRight( curstem, curp, dist, trprob );
        pIns->trprob[ SCFG_MS ]=trprob[0];
        pIns->trprob[ SCFG_RGI ]=trprob[1];
        calcIBaseEmProbRight( pIns, curstem, curp, dist);
    
        pLoop->nextI=pIns;
        pLoop->nextS=curG;
        pIns->nextS=curG;
        pIns->nextRI=pIns;
    }
    return uid;
}


void ScfgBuilder::generateSCFG( Stem *curstem )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::generateSCFG\n");
    #endif
    
    int i=0, j=0,  uid=0;
    int status=0, bulgestatus=0;
    int lenlarm=0, lenrarm=0;
    int idl=-1, idr=-1, disl=0, disr=0;
    int bulge[2]={0, 0};

    BulgeInfo *leftB=NULL,  *rightB=NULL;
    ScfgGram *curG=NULL, *preG=NULL, *preV=NULL;
    
    lenlarm=curstem->Pos[1]-curstem->Pos[0]+1;
    lenrarm=curstem->Pos[3]-curstem->Pos[2]+1;
    
    curstem->state[ SCFG_LEFT ]=new char[ lenlarm+1 ];
    memset( curstem->state[ SCFG_LEFT ], '\0', lenlarm+1 );
    curstem->state[ SCFG_RIGHT ]=new char[ lenrarm+1 ];
    memset( curstem->state[ SCFG_RIGHT ], '\0', lenrarm+1 );
    
    idl=0;            //...idl........idr...
    idr=lenrarm-1;
    
    preG = new ScfgGram('B', 0, uid++ );
    curstem->pStemGram=preG;
    while( 1 )
    {
        //printf("\n\n");
        for( j=0; j<2; j++ )
            bulge[j]=0;
        status=obtainNextS( curstem, idl, idr, disl, disr );
        //printf("*****  Left:%d (%d)  Right:%d (%d) *****\n", idl, disl, idr, disr );
        markState(curstem, idl, idr, disl, disr);
        if( status!=0 )
        {
            #ifdef DEBUG_DEV
                printf("S \n");
            #endif
            curG = new ScfgGram( 'S',  i+1, uid++);
        }
        else
        {
            #ifdef DEBUG_DEV
                printf("E \n");
            #endif
            curG = new ScfgGram( 'E',  i+1, uid++);
        }
        preV=new ScfgGram( 'V',  i, uid++);
        bulgestatus=judgeBulge( bulge, leftB, rightB, curstem, idl, idr, disl, disr );
        calcSTransP( preG, preV, curstem, bulge, leftB, rightB, idl, idr, disl, disr );
        if( curG->state != 'E')
            calcPairEmProb( curG, curstem, idl+disl, idr-disr);
        // ******
        preG->nextS=preV;
        preV->nextS=curG;
        uid=modelLeftArm( bulge[0], leftB, curstem,    preG, preV,    idl, disl, i, uid );
        uid=modelRightArm( bulge[1], rightB, curstem,    preV, curG,    idr, disr, i, uid );
        
        delete [ ] leftB;
        leftB=NULL;
        delete [ ] rightB;
        rightB=NULL;
        
        if( status==0 )
            break;
            
        preG=curG;
        preV=NULL;
        idl=idl+disl+1;
        idr=idr-disr-1;
        disl=0;
        disr=0;
        i++;
    }
}

int ScfgBuilder::setTrainingData( DataLoader &trdata )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::setTrainingData\n");
    #endif
    int i=0, j=0;
    
    if( trdata.numstems==0 )
    {
        numstems=0;
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
    
    //---- Set the IDs of all stems
    numstems=trdata.numstems;
    allstems=new Stem[numstems];

    for(i=0; i<numstems; i++)
    {
        allstems[i].idchar=trdata.stemary[i].charid;
        allstems[i].idstr=trdata.stemary[i].strid;
        allstems[i].UID=i;
        //allstems[i].numseqs=numseqs;
        allstems[i].numorgs=numseqs;
        
        
        allstems[i].Pos[0]=trdata.stemary[i].Pos[0];
        allstems[i].Pos[1]=trdata.stemary[i].Pos[1];
        allstems[i].Pos[2]=trdata.stemary[i].Pos[2];
        allstems[i].Pos[3]=trdata.stemary[i].Pos[3];
        allstems[i].distance=trdata.stemary[i].distance;
        allstems[i].stemlen=trdata.stemary[i].stemlen;
        
        allstems[i].regInfo = &(trdata.stemary[i]);
        trdata.copyLogBgBaseHz( allstems[i].bgBaseLogHz, 5 );
    }
    return 1;
}

void ScfgBuilder::printPriorFreqMatrix( )
{
     printf("Prior Freq Matrix:\n");
     printPairEmP( priorPairFreq );
}

void ScfgBuilder::countPairInformation(  )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::countPairInformation\n");
    #endif
    int i=0, j=0, k=0;
    for(i=0; i<numstems; i++ )
    {
        allstems[i].countPairFreq( );
        //summation
        for(k=0; k<5; k++ )
            for(j=0; j<5; j++ )
                sumPairFreq[k][j] += allstems[i].sumcolpair[k][j];
    }
}

void ScfgBuilder::setPriorFilePath( char * fn )
{
    int len = strlen(fn);
    priorfn = new char[len+1];
    strcpy( priorfn, fn );
}


int ScfgBuilder::putRNAPairHz(  )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::putRNAPairHz\n");
    #endif
    int i=0, j=0;
    double tempval=0.0;
    double tempall=0.0;
    int retval=0;
    FILE * fhz=NULL;
    
    if( priorfn!=NULL )
    {
        //check if the file exists?
        fhz = fopen( priorfn,  "r");
        if( fhz == NULL )
        {
            printf("Please make sure the file exists: %s\n", priorfn);
            exit(0);
        }
        else
        {
            for(i=0; i<5; i++)
            {
                for(j=0; j<5; j++)
     	        {
     	            fscanf( fhz, " %lf ", &tempval);
     	            priorPairFreq[i][j]=tempval;
                }
            }
            retval=1;
        }
        fclose(fhz);
    }
    else
    {
        for(i=0; i<5; i++)
        {
            for(j=0; j<5; j++)
     	        tempall += sumPairFreq[i][j];
     	}
        for(i=0; i<5; i++)
        {
            for(j=0; j<5; j++)
     	        priorPairFreq[i][j]=(double)sumPairFreq[i][j]/tempall;
     	}
        retval=1;
    }
    return retval;
}

int ScfgBuilder::getMaxStrengthStemNo( )
{
    int i=0,  maxIdx=0;
    double curs=0.0, maxs=0.0;
    for( i=0; i<numstems; i++ )
    {
        curs=allstems[i].getStemStrength( );
        printf("StemNo: %2d(%c)    StemLen: %3d    Strength: %6.2f ( %5.2f per column)\n", i, allstems[i].idchar, allstems[i].stemlen, curs, curs/allstems[i].stemlen );
        if(  curs > maxs ) {
            maxIdx=i;
            maxs = curs;
        }
    }
    return maxIdx;
}    

int ScfgBuilder::buildSCFG( DataLoader &trdata )
{
    #ifdef DEBUG_DEV
        printf("ScfgBuilder::buildSCFG\n");
    #endif
    int status=0;
    int i=0;
    
    status=setTrainingData( trdata ); //---- Load training pasta
    if( !status )
        return 0;
    //---- Seperate the stem from sequences
    seperateStem(  );
    
    //count all pair frequency
    countPairInformation(  );
    
    //---- Generate SCFG for every stem
    status=putRNAPairHz( );
        
    for(i=0; i<numstems; i++ )
    {
        allstems[i].setPriorPairFreq( priorPairFreq, priorcoef );
        allstems[i].pickConservedOnes( );
        generateSCFG( &allstems[i] );
    }
    return 1;
}
