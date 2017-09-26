#pragma warning( disable: 4786 )

#include "CYK.h"

//********************************************************************************
//*******************    implementation of CYKCandidate class      ***************
//********************************************************************************
CYKCandidate::CYKCandidate( )
{
    #ifdef DEBUG_DEV
        printf("CYKCandidate::CYKCandidate\n");
    #endif
    initAllMembers( );
}

CYKCandidate::~CYKCandidate( )
{
    #ifdef DEBUG_DEV
        printf("CYKCandidate::~CYKCandidate\n");
    #endif
    freeAllocMemory( );
}

CYKCandidate & CYKCandidate::operator=(  CYKCandidate &org )
{
    int i=0, j=0;
    if (this == &org)
        return *this;

    pos[0]=org.pos[0];
    pos[1]=org.pos[1];
    pos[2]=org.pos[2];
    
    for(i=0; i<2; i++ )
    {
        for(j=0; j<2; j++ )
        {
            loc[i][j]=org.loc[i][j];
            preg[i][j]=org.preg[i][j];
            numins[i][j]=org.numins[i][j];
        }
    }
    prob=org.prob;
    penalty=org.penalty;
    for(i=0; i<3; i++ )
    {
        indvpen[i]=org.indvpen[i];
    }

    distance=org.distance;
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

void CYKCandidate::initAllMembers( )
{
    #ifdef DEBUG_DEV
        printf("CYKCandidate::initAllMembers\n");
    #endif
    int i=0, j=0;
    
    pos[0]=INVALID;
    pos[1]=INVALID;
    pos[2]=INVALID;
    for(i=0; i<2; i++ )
    {
        for(j=0; j<2; j++ )
        {
            loc[i][j]=-1;
            preg[i][j]=-1;
            numins[i][j]=0;
        }
    }
    finlen=-1;
    cadno=-1;
    prob=INVLDPROB;
    for(i=0; i<3; i++ )
    {
        indvpen[i]=INVLDPROB;
    }
    penalty=INVLDPROB;
    finseq=NULL;
    finstr=NULL;
}

void CYKCandidate::initCYKCandidates( )
{
    #ifdef DEBUG_DEV
        printf("CYKCandidate::initCYKCandidates\n");
    #endif
    freeAllocMemory( );
    initAllMembers( );
}

void CYKCandidate::freeAllocMemory( )
{
    #ifdef DEBUG_DEV
        printf("CYKCandidate::freeAllocMemory\n");
    #endif
    if( finseq!=NULL )
    {
        delete [] finseq;
        finseq=NULL;
    }
    if( finstr!=NULL )
    {
        delete [] finstr;
        finstr=NULL;
    }
}

int CYKCandidate::hasData( )
{
    #ifdef DEBUG_DEV
        printf("CYKCandidate::hasData\n");
    #endif
    if( pos[0]==INVALID )
        return 0;
    else
        return 1;
}

double CYKCandidate::getLog(  )
{
    return prob;  
}

void CYKCandidate::printAlgnedSeq( )
{
    #ifdef DEBUG_DEV
        printf("CYKCandidate::printAlgnedSeq\n");
    #endif

    if( cadno != -1 )
    {
        printf("No.%d:  ", cadno);
    }
    printf("Max Log-Prob: %10.6f  at [%d, %d]-[%d, %d]\n",
                prob,  loc[0][0], loc[0][1], loc[1][0],loc[1][1]);
    //printf("Result: \n");
    printf("   ");
    printf("%s\n", finseq );
    printf("   ");
    printf("%s\n\n", finstr );
}

//********************************************************************************
//*******************    implementation of CYKSearch class         ***************
//********************************************************************************
CYKSearch::CYKSearch(  )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKSearch\n");
    #endif
    int i=0;
    
    searchSeq=NULL;
    seqlen=0;
    numntm=0;
    numtops=0;
    
    coestd=COESTDEV;
    allowvar=0;
    maxWid=0;
    midLen[0]=0;
    midLen[1]=0;
    
    exof=2;
    
    //for length penalty
    for(i=0; i<2; i++ )
    {
        offsetE[i]=0.0;
        offsetSD[i]=0.0;
        armE[i]=0.0;
        armSD[i]=0.0;
    }
    midE=0.0;
    midSD=0.0;    

    Prob=NULL;
    Gram=NULL;
    P=NULL;
    G=NULL;
    rec=NULL;
    tops=NULL;

    candidmerge=CYKCANDFILTER;
    candshortlen=PREFERSHORTLEN;
    allowshift=STEMFTALLOW;
    scoredroprate=PERCENTSCOREDROP;
    
    //length penalty
    plowerbound=1.0;
    pcoeff=2.0;
}

CYKSearch::CYKSearch( int candsize, double stdtimes, int extravar )  //variable range and the number of candidates
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKSearch\n");
    #endif
    int i=0;
    
    searchSeq=NULL;
    seqlen=0;
    numntm=0;
    numtops=candsize;
    
    coestd=stdtimes;
    allowvar=extravar;
    maxWid=0;
    midLen[0]=0;
    midLen[1]=0;
    exof=2;
    
    for(i=0; i<2; i++ )
    {
        offsetE[i]=0.0;
        offsetSD[i]=0.0;
        armE[i]=0.0;
        armSD[i]=0.0;
    }
    midE=0.0;
    midSD=0.0;  
    
    Prob=NULL;
    Gram=NULL;
    P=NULL;
    G=NULL;
    rec=NULL;
    tops=new CYKCandidate[ numtops ];
    
    candidmerge=CYKCANDFILTER;
    candshortlen=PREFERSHORTLEN;
    allowshift=STEMFTALLOW;
    scoredroprate=PERCENTSCOREDROP;
    
    //length penalty
    plowerbound=1.0;
    pcoeff=2.0;
}

CYKSearch::~CYKSearch( )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::~CYKSearch\n");
    #endif
    if(tops!=NULL)
    {
       delete[] tops;
       tops=NULL;
    }
    free3DTables( );
    if(searchSeq!=NULL)
    {
        delete [] searchSeq;
        searchSeq=NULL;
    }
}

void CYKSearch::setParameter( int candsize, double stdtimes, int extravar )
{
    numtops=candsize;
    coestd=stdtimes;
    allowvar=extravar;
    tops=new CYKCandidate[ numtops ];
}

//mergevalid:  define whether candidate-merging is used.
//shortlen:    define how many shifting positions are allowed.
//shiftnum:    define whether the shortest stem is selected into top n.
void CYKSearch::setCandMergParameter( bool mergevalid, bool shortlen, int shiftnum, double droprate )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::setCandMergParameter\n");
    #endif
    candidmerge = mergevalid;
    candshortlen = shortlen;
    allowshift = shiftnum;
    scoredroprate = droprate;
}

void CYKSearch::setLengthPenaltyParameters( double lbound, double adjustk )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::setLengthPenaltyParameters\n");
    #endif
    plowerbound=lbound;  //if > lowerbound, we add penalty, else we don't
    pcoeff=adjustk;      //adjust the penalty level
}

int CYKSearch::obtainSequence(  char * filename )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::obtainSequence\n");
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
        printf("There's something wrong with target file.\n");
        return 0;
    }
    fclose(fp);
    return templen;
}

int CYKSearch::getTopCandNum( )
{
    return  numtops;
}

CYKCandidate * CYKSearch::getCandHead( )
{
    return  tops;
}

void CYKSearch::printResultSeq( )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::printResultSeq\n");
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

void CYKSearch::printSearchSeq( )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::printSearchSeq\n");
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

int CYKSearch::setBoundaryConds( Stem *curstem, int sequlen )
{
    maxWid=curstem->getStemSpan( coestd, 1 )+ allowvar;  //get maximum width
    if( sequlen<maxWid )
        maxWid=sequlen;
    
    midLen[0]=curstem->getMidlen( coestd, 0 ) - allowvar;   //get minimum mid loop length
    midLen[0]=midLen[0]>=0 ? midLen[0]: 0;
    	
    midLen[1]=curstem->getMidlen( coestd, 1 ) + allowvar;  //get maximum mid loop length
    midLen[1]=midLen[1]<=maxWid ? midLen[1]: maxWid;
    //check computability
    if( maxWid<=midLen[0] )
    {
        //printf("midlen:[%d~%d] and MaxWidth:%d\n", midLen[0], midLen[1], maxWid );
        //printf("Target sequence is too short. Abort execution in CYK.\n");
        //exit(0);
    }
	return 1;
}

int CYKSearch::CYKSearchFile( char *filename, Stem *curstem )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKSearchFile\n");
    #endif
    
    if( curstem==NULL )
        return 0;
          
    if( !curstem->isProdsArray(  ) )
        curstem->generateProdsArray( );    

    int status=0;
    int numprod=curstem->numprod;

    GramRule * prod = curstem->prod;
    NontermProp *nonterm = curstem->nonterm;
    numntm = curstem->numntm;
    
    seqlen = obtainSequence( filename );
    printSearchSeq( );
    if( seqlen == 0)
        return 0;
    //allocate P and G memory
    setBoundaryConds( curstem, seqlen );
    status = malloc3DTables( );
    if( status==0)
    {
        printf(" error in malloc3D\n");
        return 0;
    }
    CYKimplement( prod, numprod, nonterm );
    //printResultSeq( );
    return 1;
}

int	CYKSearch::update3DandSeq( int newlen,  int newnumntm )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::update3DandSeq\n");
    #endif
    bool change=false;

    if( newlen!= seqlen || newnumntm!= numntm )
    {
        free3DTables( );
        change = true;
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
    
    if( newnumntm!= numntm )
        numntm = newnumntm;
    
    if( change )
    {
        if( malloc3DTables( )==0 )
        {
            //printf("Fail to allocate memory in malloc3DTables\n");
            return 0;
        }
    }
    return 1;
}

int CYKSearch::checkShiftedSeq( int *instr, int newlen, int shiftnum)
{
    int i=0;
    int retval=1;
    if( seqlen!=newlen )
    {
        retval=0;
    }
    else
    {
        for( i=shiftnum; i<d1; i++ )
        {
             if( instr[i-shiftnum] !=searchSeq[i])
             {
                 retval=0;
                 break;
             }
        }
    }
    return retval;    
}

int CYKSearch::CYKSearchString( char *instr, int newlen, Stem *curstem, int firstrun, int shiftnum )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKSearchString(char *)\n");
    #endif
    
    if( newlen<=0 )
        return 0;   

    int i=0, j=0, curlen=0;
    int * bufstr=NULL;

    for( i=0; i<newlen; i++ )
    {
         if( isRNABase( instr[i] ))
             curlen++;
    }
    bufstr = new int[ curlen ];

    for( i=0, j=0; i<newlen; i++ )
    {
        if( isRNABase( instr[i] ) )
             bufstr[j++]=chartonum( instr[i] );
    }
    CYKSearchString( bufstr,  curlen,  curstem,  firstrun, shiftnum );
    return 1;
}

int CYKSearch::CYKSearchString( int *instr, int newlen, Stem *curstem, int haverun, int shiftnum )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKSearchString(int *)\n");
    #endif
    if( newlen<=0 )
        return 0;
    if( curstem==NULL )
        return 0;
    if( !curstem->isProdsArray(  ) )
        curstem->generateProdsArray( );
        
    int i=0, status = 0;
    int numprod=0, newnumntm=0;

    GramRule * prod = curstem->prod;
    NontermProp *nonterm = curstem->nonterm;
    numprod=curstem->numprod;
    newnumntm=curstem->numntm;
    
    //copy avg and sd
    for(i=0; i<2; i++ )
    {
        offsetE[i]=curstem->regInfo->woffsetE[i];
        offsetSD[i]=curstem->regInfo->woffsetSD[i];
        armE[i]=curstem->regInfo->armE[i];
        armSD[i]=curstem->regInfo->armSD[i];
    }
    midE = curstem->regInfo->midE;
    midSD = curstem->regInfo->midSD;
    
    //check speedup sequence
    if( haverun==1 )
    {
        status=checkShiftedSeq( instr, newlen, shiftnum);
        if( status==0 )
        {
            printf("Error: sequence shifting %d residues\n", shiftnum);
            haverun=0;
        }
    }
    //1.check the dimension of 3D
    setBoundaryConds( curstem, newlen );
    status = update3DandSeq( newlen, newnumntm );
    if( status == 0)
        return 0;
    //2.copy searchSeq
    for( i=0; i<seqlen; i++ )
    {
         searchSeq[i] = instr[i];
    }
    //3.do CYK profile
    if( haverun==0)
    {
        CYKimplement( prod, numprod, nonterm );
    }
    else
    {
        shiftPosIndex( shiftnum );
        CYKShift( prod, numprod, nonterm, shiftnum );
    }
    return 1;
}

int CYKSearch::CYKSearchTwoRegion( int *instr, int newlen, int midreg[2], Stem *curstem, int regNo )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKSearchTwoRegion(int *)\n");
    #endif
    if( newlen<=0 )   return 0;
    if( curstem==NULL )    return 0;
    if( !curstem->isProdsArray(  ) )
        curstem->generateProdsArray( );
        
    int i=0;
    int status=0;
    int numprod=0, newnumntm=0;
    int midreglen=0;

    GramRule * prod = curstem->prod;
    NontermProp *nonterm = curstem->nonterm;
    numprod=curstem->numprod;
    newnumntm=curstem->numntm;
    
    //copy avg and sd
    for(i=0; i<2; i++ )
    {
        offsetE[i]=curstem->regInfo->offsetE[i][ regNo ];
        offsetSD[i]=curstem->regInfo->offsetSD[i][ regNo ];
        armE[i]=curstem->regInfo->armE[i];
        armSD[i]=curstem->regInfo->armSD[i];
    }
    midE = curstem->regInfo->midE;
    midSD = curstem->regInfo->midSD;

    //1.check boundary
    setBoundaryConds( curstem, newlen );
    midreglen = midreg[1]-midreg[0]+1;
    //printf("m0: %d;      m1: %d ---- midRegLen:%d\n",midreg[0], midreg[1], midreglen );
    if( midreglen > midLen[0] )
        midLen[0]= midreglen;
    if( midreglen > midLen[1] )
    {
        //printf("Skip: given mid-loop length is beyond the maximum allowed length\n");
        initTopCands(  );
        return 0;
    }

    //2.check the dimension of 3D
    status  = update3DandSeq( newlen, newnumntm );
    if( status == 0)
        return 0;
        
    //3.copy searchSeq
    for( i=0; i<seqlen; i++ ) {
        searchSeq[i] = instr[i];
    }

    //4.do CYK profile
    /*printf("mid len: %d-%d\n", midLen[0], midLen[1]);
    printf("mid reg: %d-%d   ", midreg[0], midreg[1]);
    if( midreg[0]< midreg[1] ) {
        for( i=midreg[0]; i<=midreg[1]; i++ )
            printf( "%c", numtochar( instr[i]  ));
        printf( "+\n" );
    }
    else  {
        for( i=midreg[1]; i<=midreg[0]; i++ )
            printf( "%c", numtochar( instr[i]  ));
        printf( "-\n" );
    }
    */
    CYKOnTwoRegion( prod, numprod, nonterm, midreg );
    return 1;
}

int CYKSearch::malloc3DTables( )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::malloc3DTables\n");
    #endif
    int i=0, j=0;
    //set the 1D and 2D's dimension
    d1=seqlen;
    d2=maxWid-midLen[0]-1+exof;
    if(d1<=0 || d2<=0)
    {
        return 0;
    }
    //assistant registers
    P=new double** [ d1 ];
    G=new Triple** [ d1 ];
    
    //***** merging candidate
    rec=new bool* [ d1 ];
    for(i=0; i<d1; i++)
    {
        rec[i] = new bool[ d2 ];
    }
    
    //***** 1D  [x][ d1 ][ ]
    Prob=new double** [ d1 ];
    if( Prob==NULL)
    {
        printf("Fail to allocate memory for double Prob[][][]\n");
        return 0;
    }
    Gram=new Triple** [ d1 ];
    if( Gram==NULL)
    {
        printf("Fail to allocate memory for Triple* Gram[][][]\n");
        return 0;
    }
    //***** 2D  [ ][x][ ]
    for( i=0; i<d1; i++ )
    {
        Prob[i]=new double*[ d2 ];
        if( Prob[i]==NULL)
        {
            printf("Fail to allocate memory for double Prob[][][]\n");
            return 0;
        }
        Gram[i]=new Triple* [ d2 ];
        if( Gram[i]==NULL)
        {
            printf("Fail to allocate memory for Triple* Gram[][][]\n");
            return 0;
        }
    }
    //***** 3D  [ ][ ][x]
    for( i=0; i<d1; i++ )
    {
        for( j=0; j<d2; j++ )
        {
            Prob[i][j]= new double [ numntm ];
            if( Prob[i][j]==NULL)
            {
                printf("Fail to allocate memory for double Prob[][][]\n");
                 return 0;
            }
            Gram[i][j]=new Triple [ numntm ];
            if( Gram[i][j]==NULL)
            {
                printf("Fail to allocate memory for Triple* Gram[][][]\n");
                return 0;
            }
        }
    }
    storePosIndex( );
    init3DTables ( );
    resetRec( );
    return 1;
}

void CYKSearch::free3DTables( )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::free3DTables\n");
    #endif
    
    int i=0, j=0;
    if( rec != NULL )
    {
        for(i=0; i<d1; i++)
        {
            if( rec[i] != NULL )
            {
                delete [ ] rec[i];
                rec[i]=NULL;
            }
        }
        delete [ ] rec;
        rec = NULL;
    }

    //free p, G
    if( P!=NULL )
    {
        delete [ ] P;
        P=NULL;
    }
    if( G!=NULL )
    {
        delete [ ] G;
        G=NULL;
    }
    //free Prob
    if( Prob!=NULL)
    { 
        for(i=0; i<d1; i++ )
        {
            if( Prob[i]!=NULL)
            {
                for(j=0; j<d2; j++ )
                {
                    if( Prob[i][j]!=NULL)
                    {
                        delete [ ] Prob[i][j];
                        Prob[i][j]=NULL;
                    }
                }
                delete [ ] Prob[i];
                Prob[i]=NULL;
            }
        }
        delete [ ] Prob;
        Prob=NULL;
    }
    //free Gram
    if( Gram!=NULL)
    { 
        for(i=0; i<d1; i++ )
        {
            if( Gram[i]!=NULL)
            {
                for(j=0; j<d2; j++ )
                {
                    if( Gram[i][j]!=NULL)
                    {
                        delete [ ] Gram[i][j];
                        Gram[i][j]=NULL;
                    }
                }
                delete [ ] Gram[i];
                Gram[i]=NULL;
            }
        }
        delete [ ] Gram;
        Gram=NULL;
    }
}

void CYKSearch::initTopCands( )
{
    int i=0;
    for(i=0; i<numtops; i++ )
    {
        tops[i].initCYKCandidates( );
    }
}

void CYKSearch::storePosIndex( )
{
    int i=0;
    for( i=0; i<d1; i++ )
    {
        P[i]=Prob[i];
        G[i]=Gram[i];
    }
}

void CYKSearch::shiftPosIndex( int num )
{
    int i=0, j=0, k=0;
    
    //reset    
    for(i=0; i<num; i++)
    {
        for(j=0; j<d2; j++)
        {
            for(k=0; k<numntm; k++)
            {
                P[i][j][k]=INVLDPROB;
            }
        }
    }
    for( i=0; i<d1; i++ )
    {
        k=(i+num)%d1;
        P[i]=P[k];
        G[i]=G[k];
    }
}

//initialization
void CYKSearch::init3DTables(  )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::init3DTables\n");
    #endif
    
    int i=0, j=0, k=0;
    
    for(i=0; i<d1; i++)
    {
        for(j=0; j<d2; j++)
        {
            for(k=0; k<numntm; k++)
            {
                P[i][j][k]=INVLDPROB;
            }
        }
    }
}

void CYKSearch::resetRec(  )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::init3DTables\n");
    #endif
    
    int i=0, j=0;
    
    for(i=0; i<d1; i++)
    {
        for(j=0; j<d2; j++)
        {
            rec[i][j] = false;
        }
    }
}

//--------------------- CYK Algorithm ---------------------------
int CYKSearch::isInRemSeq( int *seq, int size,  int num)
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

int CYKSearch::CYKTraceBack( CYKCandidate &cad  )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKTraceBack\n");
    #endif

    int i=0, w=0, k=0, j=0, alld=0, found=0;
    int rm=0, lm=0, newlen=0;
    int side=-1;
    bool needpostp=false;
    int cp=midLen[0]+1-exof;
    int *rrem=new int[seqlen*3+50];    // ]]]  make it large enough to contain many GAP
    int *lrem=new int[seqlen*3+50];    // [[[
    char *align=new char[seqlen*3+300];
    char curtype=' ';
    Triple *curTriple=NULL;
    
    for( i=0; i<seqlen; i++)
    {
        align[i]='.';
        rrem[i]=INVALID;
        lrem[i]=INVALID;
    }
    align[i]='\0';
    
    i=cad.pos[0];
    w=cad.pos[1];
    k=cad.pos[2];
    cad.loc[0][0]=i;
    cad.loc[1][1]=i+w+cp;
    curTriple = &G[i][w][k];
    curtype = G[i][w][k].self;
    while( curtype!=' ' )
    { 
        curtype=G[i][w][k].self;
        //printf("[%d][%d][%d]: %c   score:%f\n", i, w, k, curtype, P[i][w][k] );
        side=G[i][w][k].prodtype;

        cad.loc[0][1]=i;
        cad.loc[1][0]=i+w+cp;
        switch( curtype )
        {
            case 'B': break;
            case 'E': break;
            case 'V': break;
            case 'L':
                      align[i]= AlgScfgLeftIns;
                      i=i+1;
                      w=w-1;
                      break;
            case 'R':
                      align[i+w+cp]=AlgScfgRightIns;
                      w=w-1;
                      break;
            //==== Bulge HMM ====
            case 'M':
                      if( side == LEFT_INS )
                      {
                      	 align[i]=AlgHmmMatch;
                         i=i+1;
                      }
                      else
                      {
                      	 align[i+w+cp]=AlgHmmMatch;
                      }
                      w=w-1;
                      break;
            case 'I':
                      align[i]=AlgHmmIns;
                      if( side == LEFT_INS )
                      {
                      	 align[i]=AlgHmmIns;
                         i=i+1;
                      }
                      else
                      {
                      	 align[i+w+cp]=AlgHmmIns;
                      }
                      w=w-1;
                      break;
            case 'D':
                      if( side == LEFT_INS )
                      {
                          lrem[ lm++ ]=-i;   //"negative" means HMM
                      }
                      else
                      {
                          rrem[ rm++ ]=-i-w-cp;
                      }
                      break;
            //==== pair S ====
            case 'S':
                      align[i]=AlgScfgLeftArm;
                      align[i+w+cp]=AlgScfgRightArm;
                      i=i+1;
                      w=w-2;
                      break;
            case 'P':
                      align[i]=AlgScfgLeftArm;     //x....-
                      rrem[ rm++ ]=i+w+cp;
                      i=i+1;
                      w=w-1;
                      break;
            case 'Q':
                      align[i+w+cp]=AlgScfgRightArm;   //-....x
                      lrem[ lm++ ]=i;
                      w=w-1;
                      break;
            case 'W':
                      lrem[ lm++ ]=i;
                      rrem[ rm++ ]=i+w+cp;
                      break;
            //==== end of Pair S ====
            //==== Terminal pair S ====
            case 's':
                      align[i]=AlgScfgLeftArm;
                      align[i+w+cp]=AlgScfgRightArm;
                      i=i+1;
                      w=w-2;
                      break;
            case 'p':
                      cad.loc[1][0]=i+w+cp+1;
                      align[i]=AlgScfgLeftArm;     //x....-
                      rrem[ rm++ ]=i+w+cp;
                      i=i+1;
                      w=w-1;
                      break;
            case 'q':
                      cad.loc[0][1]=i-1;
                      align[i+w+cp]=AlgScfgRightArm;   //-....x
                      lrem[ lm++ ]=i;
                      w=w-1;
                      break;
            case 'w':
                      if( i>i+w+cp )
                          needpostp=true;
                      cad.loc[0][1]=i-1;
                      cad.loc[1][0]=i+w+cp+1;
                      lrem[ lm++ ]=i;
                      rrem[ rm++ ]=i+w+cp;
                      break;
            //==== end ====
            default:
                      printf(" Unknown nonterminal: %c id: %d!\n", curtype, k);
                      return 0;
        }
        if( (curTriple->prodtype == PAIR_TMN )||
            (curTriple->prodtype == SINGLE_TMN )||
            (curTriple->prodtype == LSINGLE_TMN )||
            (curTriple->prodtype == RSINGLE_TMN ) )
        {
            if ( curtype == 'L' )  //adjust the index for special cases: xxxL...xxx, xxx...Rxxx
				cad.loc[1][0] ++;
			if ( curtype == 'R' )
				cad.loc[0][1] --;
            break;
        }
        k = curTriple->id;
        curtype=curTriple->state;
        curTriple = &G[i][w][k];
    }
    //printf("xxx\n");
    //***********calculate the pair region
    for(i=0; i<seqlen; i++ )
    {
        //count the number of insertion at most left side
        if( align[i]==AlgScfgLeftIns && cad.preg[0][0]==-1 )
            cad.numins[0][0]++;
        if( align[i]==AlgScfgLeftArm && cad.preg[0][0]==-1 )
            cad.preg[0][0]=i;
        if( align[i]==AlgScfgRightIns && cad.preg[1][0]==-1 )
            cad.numins[1][0]++;
        if( align[i]==AlgScfgRightArm && cad.preg[1][0]==-1 )
            cad.preg[1][0]=i;

        //count the number of insertion at most right side
        if( align[seqlen-i-1]==AlgScfgLeftIns && cad.preg[0][1]==-1 )
            cad.numins[0][1]++;
        if( align[seqlen-i-1]==AlgScfgLeftArm && cad.preg[0][1]==-1 )
            cad.preg[0][1]=seqlen-i-1;
        if( align[seqlen-i-1]==AlgScfgRightIns && cad.preg[1][1]==-1 )
            cad.numins[1][1]++;
        if( align[seqlen-i-1]==AlgScfgRightArm && cad.preg[1][1]==-1 )
            cad.preg[1][1]=seqlen-i-1;
    }
    //calculate distance
    cad.distance = cad.loc[1][0]-cad.loc[0][1]-1;
    
    //calculate the penalty of variant length
    if( pcoeff!=0 )
    {
        double sd=0.0, dis=0.0;
        //left part
        dis = abs(cad.loc[0][0] - coestd*offsetSD[0]);
        sd = offsetSD[0];
        cad.indvpen[0]=calcLengthPenalty( dis, sd );
        
        //middle part
        dis = abs(cad.distance-midE);
        sd = midSD;
        cad.indvpen[1]=calcLengthPenalty( dis, sd );
        
        cad.penalty = cad.indvpen[0] + cad.indvpen[1];
    }
    else
    {
        cad.penalty = 0;
    }
    
    //printf("Traceback:\n%s\n", align);
    //printf("%d - %d ... %d - %d\n", cad.preg[0][0], cad.preg[0][1], cad.preg[1][0], cad.preg[1][1] );
    //***********process the output string
    newlen=seqlen + lm + rm;
    //cout<<"seqlen:"<<seqlen<<"newlen:"<<newlen<<" lm:"<<lm<<" rm:"<<rm<<endl;
    cad.finlen=newlen;
    cad.finstr= new char[newlen+1];
    cad.finseq= new char[newlen+1];
    for( i=0, k=0; i<newlen; i++)
    {
        found=isInRemSeq( rrem, rm,  i-1-alld);
        if( found )  //right arm
        {
            for( j=0; j<found; j++ )
            {
                cad.finseq[i+j] = '-';
                cad.finstr[i+j] = AlgScfgRightArm;
            }
            i+=found-1;
            alld+=found;
            continue;
        }
        found=isInRemSeq( rrem, rm,  -(i-1-alld));
        if( found )
        {
            for( j=0; j<found; j++ )
            {
                cad.finseq[i+j] = '-';
                cad.finstr[i+j] = AlgHmmDel;
            }
            i+=found-1;
            alld+=found;
            continue;
        }
        found=isInRemSeq( lrem, lm,  i-alld);
        if( found )  //left arm
        {
            for( j=0; j<found; j++ )
            {
                cad.finseq[i+j] = '-';
                cad.finstr[i+j] = AlgScfgLeftArm;
            }
            i+=found-1;
            alld+=found;
            continue;
        }
        found=isInRemSeq( lrem, lm,  -(i-alld));
        if( found )
        {
            for( j=0; j<found; j++ )
            {
                cad.finseq[i+j] = '-';
                cad.finstr[i+j] = AlgHmmDel;
            }
            i+=found-1;
            alld+=found;
            continue;
        }
        cad.finseq[i]= numtochar(searchSeq[k]);
        cad.finstr[i]= align[k];
        k++;
    }

    delete [ ] rrem;
    delete [ ] lrem;
    delete [ ] align;
    
    cad.finseq[i]='\0';
    cad.finstr[i]='\0';
    //postprocessing
    if( needpostp )
    {
        for(i=0; i<cad.finlen-1; i++ )
        {
            if( cad.finstr[i]==AlgScfgRightArm &&
                cad.finstr[i+1]==AlgScfgLeftArm && 
                cad.finseq[i]=='-' && 
                cad.finseq[i+1]=='-')
            {
                cad.finstr[i]=AlgScfgLeftArm;
                cad.finstr[i+1]=AlgScfgRightArm;
                break;
            }
        }                
    }
    return 1;
}

double CYKSearch::calcLengthPenalty( double dis, double sd )
{
    double k=0.0, penalty=0.0;
    
    if( sd ==0 && SDZEROINFINITY )
    {
        if( dis > SDZEROALLOWBOUND )
            penalty = INVLDPROB;
        else
            penalty = 0;
    }
    else
    {
        if( sd==0 )
            sd = 1;
        k = dis/sd;
        penalty = k > plowerbound ? log( 1.0/(k*k*pcoeff) ):0;
    }
    return penalty;
}

// get Probability in GramRule
double CYKSearch::obtainRuleProb( int begid, int endid, int size, GramRule *prod)
{
    #ifdef _DEBUG_DEV
    printf("CYKSearch::obtainRuleProb\n");
    #endif
    
    int i=0;
    double retprob=INVLDPROB;

    for(i=0; i<size; i++ )
    {
        if( prod[i].lid==begid && prod[i].rid==endid )
        {
            retprob=prod[i].singlep;
            break;
        }
    }
    return retprob;
}

//insert Top N number of positions
double CYKSearch::insertCurPos( int base, int idi,  int idj, int idk, double cprob )
{
    int i=0, j=0;
    int lowbd=base, upperbd=base+numtops-1;
    for( i=lowbd; i<=upperbd; i++ )
    {
        if( cprob > tops[i].prob )
        {
            for( j=upperbd; j>i; j--)
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
    return tops[ upperbd ].prob;
}
    
// get another max Probability
// base: the number of candidates which have been already get.
int CYKSearch::getTopnMaxProbs( int base )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::getTopnMaxProbs\n");
    #endif
    
    int i=0, j=0, d=0, ps=0;
    double thres=INVLDPROB, curprob=INVLDPROB;
    int bfloor=BEGUID;
    //printf("d2=%d, d1=%d\n", d2, d1);
    for( d=exof; d<d2; d++ )
    {
        for( ps=0; ps<d1-d+exof; ps++ )
        {
            curprob=P[ps][d][bfloor];
            if( curprob >= thres && curprob!=INVLDPROB)
            {
                thres=insertCurPos( base, ps, d, bfloor, curprob );
            }
        }
    }
    //count the final number of positions
    for(i=base; i<base+numtops; i++ )
    {
        if( tops[i].pos[0]!=-1)
            j++;
        //printf("max prob: %f\n", tops[i].prob);
    }
    return j;
}

int CYKSearch::getMaximumProbs( CYKCandidate &cad )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::getMaximumProbs\n");
    #endif
    
    int found=0, d=0, ps=0;
    double maxprob=INVLDPROB, curprob=INVLDPROB;
    int bfloor=BEGUID;
    
    //find the maximum prob
    for( d=exof; d<d2; d++ )
    {
        for( ps=0; ps<d1-d+exof; ps++ )
        {
            if( !rec[ ps ][ d ] )
            {
                curprob=P[ps][d][bfloor];
                if( curprob > maxprob && curprob!=INVLDPROB)
                {
                    found = 1;
                    cad.pos[0]=ps;
                    cad.pos[1]=d;
                    cad.pos[2]=bfloor;
                    cad.prob=curprob;
                    maxprob=curprob;
                }
            }
        }
    }
    if( found==1 )
    {
        rec[ cad.pos[0] ][ cad.pos[1] ]=true;
        return 1;
    }
    else
        return 0;    
}

// return: >=0: replace  -1: discard  -100: not similar, add this
int CYKSearch::compareSimilarCands( CYKCandidate * vault, CYKCandidate &cad, int currentNum )
{
    #ifdef _DEBUG_DEV
        printf("CYKSearch::compareSimilarCands\n");
    #endif
    int i=0, found=-100;
    int similarpos=0;
    bool shorter=false;
    int curallins=0, topallins=0;

    for( i=0; i<currentNum; i++ )
    {
        similarpos = 0;
        shorter = false;
        found = -100;
        if( vault[i].hasData( ) )
        {
            //For condition 2: compair the position of pair region
            if( abs( vault[i].preg[0][0]-cad.preg[0][0]) <= allowshift )
                similarpos++;
            if( abs( vault[i].preg[0][1]-cad.preg[0][1]) <= allowshift )
                similarpos++;
            if( abs( vault[i].preg[1][0]-cad.preg[1][0]) <= allowshift )
                similarpos++;
            if( abs( vault[i].preg[1][1]-cad.preg[1][1]) <= allowshift )
                similarpos++;
            //For condition 3: compair the number of side insertion
            curallins = cad.numins[0][0] + cad.numins[0][1] + cad.numins[1][0] + cad.numins[1][1];
            topallins = vault[i].numins[0][0] + vault[i].numins[0][1] + vault[i].numins[1][0] + vault[i].numins[1][1];
            //printf("No.%d: %d-%d:%d-%d,    temp: %d-%d:%d-%d\n", i, 
            //    vault[i].numins[0][0], vault[i].numins[0][1], vault[i].numins[1][0], vault[i].numins[1][1],
            //    cad.numins[0][0], cad.numins[0][1], cad.numins[1][0], cad.numins[1][1]);
            if( curallins < topallins )
                shorter=true;
        }

        if( similarpos >=3 )  //if similar
        {
            if( shorter && candshortlen)  //replace ith one with current one
                found=i;
            else
                found=-1;  //discard the current one
        }
        if( found>=-1)
            break;
    }
    return found;  // >=0: replace  -1: discard  -100: not similar, add this
}

char CYKSearch::computeMaxforTMNPair(  GramRule *prod,  int j, int dx,  int dy,  double *curprob)
{
    #ifdef _DEBUG_DEV
        printf("CYKSearch::computeMaxforTMNPair\n");
    #endif
    char type=' ';
    if( j==0 )
    {
        *curprob = prod->pairp[0][0];
        type = 'w';
    }
    else if(j==1)
    {
        if( prod->pairp[dx][0]>=prod->pairp[0][dy])
        {
            *curprob = prod->pairp[dx][0];
            type = 'p';
        }
        else
        {
            *curprob = prod->pairp[0][dy];
            type = 'q';
        }
    }
    else if( j<=midLen[1]-midLen[0]+exof )  // X.....(dist+1).....X
    {
        *curprob = prod->pairp[dx][dy];
        type = 's';
    }
    return type;
}

char CYKSearch::computeMaxforPair(  GramRule *prod,  int ps,  int d,   int dx,  int dy,  double *curprob)
{
    #ifdef _DEBUG_DEV
        printf("CYKSearch::computeMaxforPair\n");
    #endif
    
    int i=0, j=0;
    char type=' ';
    double maxprob=INVLDPROB;
    double tmpprob[4];
    
    tmpprob[0] = prod->pairp[dx][dy] + P[ps+1][d-2][ prod->rid ]; //S
    tmpprob[1] = prod->pairp[dx][0] + P[ps+1][d-1][ prod->rid ];  //P
    tmpprob[2] = prod->pairp[0][dy] + P[ps][d-1][ prod->rid ];    //Q
    tmpprob[3] = prod->pairp[0][0] + P[ps][d][ prod->rid ];       //W
    for( i=0; i<4; i++ )
    {
        if( tmpprob[i] > maxprob )
        {
            maxprob=tmpprob[i];
            j=i;
        }
    }
    if( j==0 )
        type='S';
    else if( j==1 )
        type='P';  //Prob: X...-
    else if( j==2 )
        type='Q';  //Prob: -...X
    else // j==3
        type='W';  //Prob: -...-
    *curprob=maxprob;
    return type;
}

int CYKSearch::CYKimplement( GramRule *prod,  int size,  NontermProp *ntms )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKimplement\n");
    #endif
    int i=0, j=0;
    int ps=0, wd=0;
    int dx=0, dy=0;
    char ltype=' ';
    double curprob=INVLDPROB;
    
    //assistant debug variables
    //double dbgprob=0.0;
    //int preps=0, prewd=0;
    
    //**** initialization ****
    initTopCands(  );
    init3DTables(  );
    resetRec(  );


    //**** main loop **** 0(-1), 1(0), 2(1)...
    for( j=0; j<d2; j++ )  //width: the length of subsequence
    {
        if( j<2)
            wd = midLen[0]+1;
        else
            wd = midLen[0]+1+j-exof;
        for( ps=0; ps<d1-wd; ps++ ) //position: beginning position of subsequence
        {
            for(i=size-1; i>=0; i-- )
            {
                curprob = INVLDPROB;
                ltype = prod[i].ltype;
                if( prod[i].category == PAIR_TMN )
                {
                    dx = searchSeq[ps];
                    dy = searchSeq[ps+wd];
                    ltype = computeMaxforTMNPair( &prod[i],  j, dx, dy,  &curprob);
                }
                else if( prod[i].category == PAIR_ADD )
                {
                    if( j>=exof )  // X.....(dist+1).....X
                    {
                        dx = searchSeq[ps];
                        dy = searchSeq[ps+wd];
                        ltype = computeMaxforPair( &prod[i], ps, j, dx, dy, &curprob);
                    }
                }
                else if( prod[i].category == LSINGLE_TMN )
                {
                    dx = searchSeq[ps]-1;
                    if( j>=exof && j<=midLen[1]-midLen[0]+exof )
                        curprob = prod[i].basep[dx];
                }
                else if( prod[i].category == RSINGLE_TMN )
                {
                    dy = searchSeq[ps+wd]-1;
                    if( j>=exof && j<=midLen[1]-midLen[0]+exof )
                        curprob = prod[i].basep[dy];
                }
                else if( prod[i].category == LEFT_INS )
                {
                    dx = searchSeq[ps]-1;
                    if( j>=exof )
                        curprob = prod[i].basep[dx] + P[ps+1][j-1][ prod[i].rid ];
                }
                else if( prod[i].category == RIGHT_INS )
                {
                    dy = searchSeq[ps+wd]-1;
                    if( j>=exof )
                        curprob = prod[i].basep[dy] + P[ps][j-1][ prod[i].rid ];
                }
                else if( prod[i].category == KEEP_POS )
                {
                    curprob = prod[i].singlep + P[ps][j][ prod[i].rid ];
                }
                else
                {
                    continue;
                }
                if( curprob > P[ps][j][prod[i].lid] )
                {
                    //printf("[%d][%d][%d]: ltype: %c, rid=%d,  %f-->%f \n", ps, j, prod[i].lid,
                    //        ltype, prod[i].rid, P[ps][j][prod[i].lid], curprob);
                    P[ps][j][prod[i].lid]=curprob;
                    G[ps][j][prod[i].lid].prodtype=prod[i].category;
                    G[ps][j][prod[i].lid].self=ltype;
                    G[ps][j][prod[i].lid].state=prod[i].rtype;
                    G[ps][j][prod[i].lid].id=prod[i].rid;
            	}
            }
        }
    }
    
    obtainCandidate( );
    return 1;
}

int CYKSearch::CYKShift( GramRule *prod,  int size,  NontermProp *ntms, int shiftednum )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKShift\n");
    #endif
    int i=0, j=0, k=0;
    int ps=0, wd=0;
    int dx=0, dy=0;
    char ltype=' ';
    double curprob=INVLDPROB;
    
    //**** initialization ****
    initTopCands(  );
    resetRec( );
    
    for( k=0; k<shiftednum; k++ )
    {
        ps=seqlen-midLen[0]-2-k;
        for( j=0; j<d2; j++ )
        {
            //actual width of a pair 3  4 5 6  7  ==> 3+3+1=7
            if( j<exof)
                wd = midLen[0]+1;
            else
                wd = midLen[0]+1+j-exof;
            
            for(i=size-1; i>=0; i-- )
            {
                curprob = INVLDPROB;
                ltype = prod[i].ltype;
                if( prod[i].category == PAIR_TMN )
                {
                    dx = searchSeq[ps];
                    dy = searchSeq[ps+wd];
                    ltype = computeMaxforTMNPair( &prod[i],  j, dx, dy,  &curprob);
                }
                else if( prod[i].category == PAIR_ADD )
                {
                    if( j>=exof )  // X.....(dist+1).....X
                    {
                        dx = searchSeq[ps];
                        dy = searchSeq[ps+wd];
                        ltype = computeMaxforPair( &prod[i], ps, j, dx, dy, &curprob);
                        //printf("j=%d,(wd=%d, pos=%d, id=%d)%c-%c, type=%c, prob=%f-->org:%f\n",j, wd, ps, prod[i].lid, numtochar(dx), numtochar(dy), ltype, curprob, P[ps][j][prod[i].lid]);
                    }
                }
                else if( prod[i].category == LSINGLE_TMN )
                {
                    dx = searchSeq[ps]-1;
                    if( j>=exof && j<=midLen[1]-midLen[0]+exof )
                        curprob = prod[i].basep[dx];
                }
                else if( prod[i].category == RSINGLE_TMN )
                {
                    dy = searchSeq[ps+wd]-1;
                    if( j>=exof && j<=midLen[1]-midLen[0]+exof )
                        curprob = prod[i].basep[dy];
                }
                else if( prod[i].category == LEFT_INS )
                {
                    dx = searchSeq[ps]-1;
                    if( j>=exof )
                        curprob = prod[i].basep[dx] + P[ps+1][j-1][ prod[i].rid ];
                }
                else if( prod[i].category == RIGHT_INS )
                {
                    dy = searchSeq[ps+wd]-1;
                    if( j>=exof )
                        curprob = prod[i].basep[dy] + P[ps][j-1][ prod[i].rid ];
                }
                else if( prod[i].category == KEEP_POS )
                {
                    curprob = prod[i].singlep + P[ps][j][ prod[i].rid ];
                }
                else
                {
                    continue;
                }
                if( curprob > P[ps][j][prod[i].lid] )
                {
                    P[ps][j][prod[i].lid]=curprob;
                    G[ps][j][prod[i].lid].prodtype=prod[i].category;
                    G[ps][j][prod[i].lid].self=ltype;
                    G[ps][j][prod[i].lid].state=prod[i].rtype;
                    G[ps][j][prod[i].lid].id=prod[i].rid;
            	}
            }
            if( j>=exof )
                ps--;
        }
    }
    
    obtainCandidate( );
    return 1;
}

//two region based CYK algorithm
int CYKSearch::CYKOnTwoRegion( GramRule *prod,  int size,  NontermProp *ntms, int midReg[] )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::CYKOnTwoRegion\n");
    #endif
    int i=0, j=0;
    int ps=0, wd=0;
    int begps=0, endps=0, tempval=0;
    int dx=0, dy=0;
    char ltype=' ';
    double curprob=INVLDPROB;
    
    //assistant debug variables
    //double dbgprob=0.0;
    //int preps=0, prewd=0;
    
    //**** initialization ****
    initTopCands(  );
    init3DTables(  );
    resetRec(  );

    //**** main loop **** 0(-1), 1(0), 2(1)...
    for( j=0; j<d2; j++ )  //width: the length of subsequence
    {
        if( j<2)
            wd = midLen[0]+1;
        else
            wd = midLen[0]+1+j-exof;
            
        //calculate the left starting position
        tempval=midReg[1]+1-wd;
        begps = tempval>0?tempval:0;

        tempval= d1-wd;  //actually seqlen-1-(wd-1)
        endps = tempval<midReg[0]?tempval:midReg[0];
                    
        //printf("j:%d, wd: %d -- ps: (%d ~ %d)\n",j, wd, begps, endps-1);
        for( ps=begps; ps<endps; ps++ ) //position: beginning position of subsequence
        {
            for(i=size-1; i>=0; i-- )
            {
                curprob = INVLDPROB;
                ltype = prod[i].ltype;
                if( prod[i].category == PAIR_TMN )
                {
                    dx = searchSeq[ps];
                    dy = searchSeq[ps+wd];
                    ltype = computeMaxforTMNPair( &prod[i],  j, dx, dy,  &curprob);
                }
                else if( prod[i].category == PAIR_ADD )
                {
                    if( j>=exof )  // X.....(dist+1).....X
                    {
                        dx = searchSeq[ps];
                        dy = searchSeq[ps+wd];
                        ltype = computeMaxforPair( &prod[i], ps, j, dx, dy, &curprob);
                    }
                }
                else if( prod[i].category == LSINGLE_TMN )
                {
                    dx = searchSeq[ps]-1;
                    if( j>=exof && j<=midLen[1]-midLen[0]+exof )
                        curprob = prod[i].basep[dx];
                }
                else if( prod[i].category == RSINGLE_TMN )
                {
                    dy = searchSeq[ps+wd]-1;
                    if( j>=exof && j<=midLen[1]-midLen[0]+exof )
                        curprob = prod[i].basep[dy];
                }
                else if( prod[i].category == LEFT_INS )
                {
                    dx = searchSeq[ps]-1;
                    if( j>=exof )
                        curprob = prod[i].basep[dx] + P[ps+1][j-1][ prod[i].rid ];
                }
                else if( prod[i].category == RIGHT_INS )
                {
                    dy = searchSeq[ps+wd]-1;
                    if( j>=exof )
                        curprob = prod[i].basep[dy] + P[ps][j-1][ prod[i].rid ];
                }
                else if( prod[i].category == KEEP_POS )
                {
                    curprob = prod[i].singlep + P[ps][j][ prod[i].rid ];
                }
                else
                {
                    continue;
                }
                if( curprob > P[ps][j][prod[i].lid] )
                {
                    //printf("[%d][%d][%d]: ltype: %c, rid=%d,  %f-->%f \n", ps, j, prod[i].lid,
                    //        ltype, prod[i].rid, P[ps][j][prod[i].lid], curprob);
                    P[ps][j][prod[i].lid]=curprob;
                    G[ps][j][prod[i].lid].prodtype=prod[i].category;
                    G[ps][j][prod[i].lid].self=ltype;
                    G[ps][j][prod[i].lid].state=prod[i].rtype;
                    G[ps][j][prod[i].lid].id=prod[i].rid;
            	}
            }
        }
    }
    
    obtainCandidate( );
    return 1;
}


int CYKSearch::obtainCandidate(  )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::obtainCandidate\n");
    #endif
    int i=0, j=0, status=0;
    int cadnum=0, nowtops=0;
    bool needmore=true;

    if ( candidmerge )
    {
        //********* traceback eliminating similar candidates ********
        int times=0, currentnum=0;
        CYKCandidate temp;
        CYKCandidate *vault=new CYKCandidate[ NUMCHECKCANDID ];
        CYKCandidate *temptop=new CYKCandidate[ numtops ];
        bool hasget[ NUMCHECKCANDID ];
        while( needmore )
        {
            currentnum=0;
            times=0;
            // initialization
            for( i=0; i<NUMCHECKCANDID; i++ )
            {
                hasget[ i ]=false;
                vault[i].initCYKCandidates( );
            }
            // 1: get these many candidates & put them in array
            for( i=0; i<NUMCHECKCANDID; i++)
            {
                hasget[ i ]=false;
                temp.initCYKCandidates( );
                status = getMaximumProbs( temp );
                if( status==0 )
                {
                    needmore=false;
                    break;
                }
                CYKTraceBack( temp  );
                vault[i] = temp;
                currentnum++;
            }
            // 2: choose the maximum candidate with length penalty
            while( nowtops < numtops && times < currentnum )
            {
                times++;
                j = getLocalMaxIndex( vault, hasget, currentnum );
                if( j==-1 )
                    break;
                hasget[j]=true;
                //**** similar check ****
                status = compareSimilarCands( temptop, vault[j], nowtops );
                if ( status > -1)    //replace the one
                    temptop[ status ] = vault[j];
                else if( status < -1 )  //not similar, add it to topn
                {
                    temptop[ nowtops ] = vault[j];
                    nowtops++;
                }
                //status == -1:  similar -> do nothing
            }
            if( nowtops==numtops )
                needmore=false;
        }
        //3: rank the candidates
        rankCandidates( temptop );
        delete [ ] vault;
        delete [ ] temptop;
    }
    else
    {
        //******************* traceback orginal version*********************
        nowtops=getTopnMaxProbs( cadnum );
        for(j=0; j<nowtops; j++)
        {
            tops[cadnum+j].cadno=j+1;
            CYKTraceBack( tops[cadnum+j] );
            tops[cadnum+j].distance=tops[cadnum+j].loc[1][0]-tops[cadnum+j].loc[0][1]-1;
        }
        cadnum+=numtops;
    }
    return 1; 
}

int CYKSearch::getLocalMaxIndex( CYKCandidate * vault, bool hasget[ ], int num )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::getLocalMaxIndex\n");
    #endif
    int i=0, j=-1;
    double maxscore=INVLDPROB;
    for( i=0; i<num; i++ )
    {
        if( !hasget[ i ] )
        {
            if( maxscore < vault[i].prob + vault[i].penalty )
            {
                j=i;
                maxscore = vault[i].prob + vault[i].penalty;
            }
        }
    }
    return j;
}

int CYKSearch::rankCandidates( CYKCandidate * unsortedtops )
{
    #ifdef DEBUG_DEV
        printf("CYKSearch::rankCandidates\n");
    #endif
    int i=0, j=0, k=0;
    bool *flag = new bool[ numtops ];
    
    for( i=0; i<numtops; i++ )
        flag[i]=false;

    for( i=0; i<numtops; i++ )
    {
        j = getLocalMaxIndex( unsortedtops, flag, numtops );
        if( j==-1)
            break;
        flag[ j ]= true;
        tops[ k++ ]= unsortedtops[j];
        tops[ k-1 ].prob = tops[ k-1 ].prob + tops[ k-1 ].penalty;
    }
    delete [] flag;
    return 1;
}


