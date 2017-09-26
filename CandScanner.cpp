#pragma warning( disable: 4786 )

#include "CandScanner.h"

//********************************************************************************
//*******************    implementation of Candidate class      ******************
//********************************************************************************
Candid::Candid( )
{
    type=' ';
    UID=INVALID;
    CompID=INVALID;
    stempos[0][0]=INVALID;
    stempos[0][1]=INVALID;
    stempos[1][0]=INVALID;
    stempos[1][1]=INVALID;
    prob=INVLDPROB;    //probability
    penalty=INVLDPROB;
    indvpen[0]=INVLDPROB;
    indvpen[1]=INVLDPROB;
    indvpen[2]=INVLDPROB;
    finseq=NULL;
    finstr=NULL;
    finlen=INVALID;
     
    nextCand=NULL;
}

Candid::Candid( char candtype, int sd, int id )
{
    type=candtype;
    UID=id;
    CompID=sd;
    if( candtype=='S' )
    {
        stempos[0][0]=INVALID;
        stempos[0][1]=INVALID;
        stempos[1][0]=INVALID;
        stempos[1][1]=INVALID;
    }
    else if( candtype=='L' )
    {
        looppos[0]=INVALID;
        looppos[1]=INVALID;
    }
    prob=INVLDPROB;
    penalty=INVLDPROB;
    indvpen[0]=INVLDPROB;
    indvpen[1]=INVLDPROB;
    indvpen[2]=INVLDPROB;
    finseq=NULL;
    finstr=NULL;
    finlen=INVALID;
    nextCand=NULL;
}

Candid::~Candid( )
{
    if(  finseq!=NULL )
    {
        delete [ ] finseq;
        finseq=NULL;
    }
    if(  finstr!=NULL )
    {
        delete [ ] finstr;
        finstr=NULL;
    }
}

Candid::Candid( Candid &org )
{
    type=org.type; //Loop:L    Stem:S
    CompID=org.CompID;    //ID of Stem or Loop
    if( type=='L' )
    {
        looppos[0]=org.looppos[0];
        looppos[1]=org.looppos[1];
    }
    else
    {
        stempos[0][0]=org.stempos[0][0];
        stempos[0][1]=org.stempos[0][1];
        stempos[1][0]=org.stempos[1][0];
        stempos[1][1]=org.stempos[1][1];
    }
    prob=org.prob;
    penalty=org.penalty;
    indvpen[0]=org.indvpen[0];
    indvpen[1]=org.indvpen[1];
    indvpen[2]=org.indvpen[2];
    nextCand=org.nextCand;
     
    finlen=org.finlen;
     
    finseq=new char[finlen+1];
    strcpy( finseq, org.finseq);
    
    finstr=new char[finlen+1];
    strcpy( finstr, org.finstr);
}

Candid & Candid::operator=(  Candid &org )
{
    if (this == &org)
        return *this;

    type=org.type;
    CompID=org.CompID;
    if( type=='L' )
    {
        looppos[0]=org.looppos[0];
        looppos[1]=org.looppos[1];
    }
    else
    {
        stempos[0][0]=org.stempos[0][0];
        stempos[0][1]=org.stempos[0][1];
        stempos[1][0]=org.stempos[1][0];
        stempos[1][1]=org.stempos[1][1];
    }
    prob=org.prob;
    penalty=org.penalty;
    indvpen[0]=org.indvpen[0];
    indvpen[1]=org.indvpen[1];
    indvpen[2]=org.indvpen[2];
    nextCand=org.nextCand;
     
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

bool Candid::isNullCandid( )
{
    if( type == 'N' )
        return true;
    else
        return false;
}

void Candid::setNullCandid( )
{
    type='N';
    UID=INVALID;
    CompID=INVALID;
    stempos[0][0]=INVALID;
    stempos[0][1]=INVALID;
    stempos[1][0]=INVALID;
    stempos[1][1]=INVALID;
    prob=0;    //probability
    penalty=0;
    indvpen[0]=0;
    indvpen[1]=0;
    indvpen[2]=0;
    finseq=NULL;
    finstr=NULL;
    finlen=INVALID;
     
    nextCand=NULL;
}

void Candid::copyValue( CYKCandidate &pCand, int bpos)
{
    int length=0;

    stempos[0][0]=pCand.loc[0][0] + bpos;
    stempos[0][1]=pCand.loc[0][1] + bpos;
    stempos[1][0]=pCand.loc[1][0] + bpos;
    stempos[1][1]=pCand.loc[1][1] + bpos;
    prob=pCand.prob;
    penalty=pCand.penalty;
    indvpen[0]=pCand.indvpen[0];
    indvpen[1]=pCand.indvpen[1];
    indvpen[2]=pCand.indvpen[2];
    
    length=pCand.finlen;
    finlen=length;
    if( finseq!=NULL )
    {
        delete [ ] finseq;
        finseq=NULL;
    }
    finseq=new char[ length+1 ];
    strcpy( finseq,  pCand.finseq );
    
    if( finstr!=NULL )
    {
        delete [ ] finstr;
        finstr=NULL;
    }
    finstr=new char[ length+1 ];
    strcpy( finstr,  pCand.finstr );
}

void Candid::copyValue( VTBCandidate &pCand, int bpos)
{
    int length=0;
    looppos[0]=pCand.loc[0] + bpos;
    looppos[1]=pCand.loc[1] + bpos;
    prob=pCand.prob;
    
    length=pCand.finlen;
    finlen=length;
    if( finseq!=NULL )
    {
        delete [ ] finseq;
        finseq=NULL;
    }
    finseq=new char[ length+1 ];
    strcpy(finseq, pCand.finseq );
    if( finstr!=NULL )
    {
        delete [ ] finstr;
        finstr=NULL;
    }
    finstr=new char[ length+1 ];
    strcpy( finstr, pCand.finstr );
}

double Candid::getProbScore(  )
{
    double orgscore=0.0;
    orgscore = prob - penalty;  //current score = prob+penalty ==> prob = current score - penalty
    return orgscore;
}

int Candid::extractAlgPart( char* &seqlm, char* &alglm, int &lenlm,
                            char* &seqrm, char* &algrm, int &lenrm )
{
    int i=0, j=0;
    int alglen[2]={0,0};
    int localpos[4]={-100, -100, -100, -100};
    //position
    for( i=0; i<finlen; i++)
    {
        if( (finstr[i]==AlgScfgLeftArm ||finstr[i]== AlgScfgLeftIns ) && localpos[0]==-100)
             localpos[0]=i;
        if( (finstr[i]==AlgScfgRightArm ||finstr[i]==AlgScfgRightIns) && localpos[2]==-100 )
             localpos[2]=i;
    }
    for( i=finlen-1; i>=0; i--)
	{
        if( (finstr[i]==AlgScfgRightArm||finstr[i]==AlgScfgRightIns ) && localpos[3]==-100 )
             localpos[3]=i;
        if( (finstr[i]==AlgScfgLeftArm||finstr[i]==AlgScfgLeftIns ) && localpos[1]==-100 )
             localpos[1]=i;
    }
    //length
    if( localpos[0]==-100 && localpos[1]==-100 )
        alglen[0]=0;
    else
        alglen[0]=localpos[1]-localpos[0]+1;
    
    if( localpos[2]==-100 && localpos[3]==-100 )
        alglen[1]=0;
    else
        alglen[1]=localpos[3]-localpos[2]+1;

    //cout<< "localpos0:"<<localpos[0]<< "    localpos1:"<<localpos[1]<<"    alglen0:"<<alglen[0] << endl;
    //cout<< "localpos2:"<<localpos[2]<< "    localpos3:"<<localpos[3]<<"    alglen1:"<<alglen[1] << endl;
    //copy string
    seqlm=new char[alglen[0]+1];
    alglm=new char[alglen[0]+1];
    for( i=0, j=0; i<alglen[0]; i++, j++ )
    {
         seqlm[j]=finseq[localpos[0]+i ];
         alglm[j]=finstr[localpos[0]+i];
    }
    seqlm[j]='\0';
    alglm[j]='\0';

    seqrm=new char[alglen[1]+1];
    algrm=new char[alglen[1]+1];
    for( i=0, j=0; i<alglen[1]; i++, j++ )
    {
         seqrm[j]=finseq[localpos[2]+i];
         algrm[j]=finstr[localpos[2]+i];
    }
    seqrm[j]='\0';
    algrm[j]='\0';
    
    lenlm=alglen[0];
    lenrm=alglen[1];
    return 1;
}

void Candid::printAlgnedSeq(  )
{
    #ifdef DEBUG_DEV
        printf("Candid::printAlgnedSeq\n");
    #endif
    if( type=='S' )
    {
        printf("Score:%10.6f  at [%d, %d]-[%d, %d]  penalty %6.4f(%6.4f, %6.4f)\n",
                prob, stempos[0][0], stempos[0][1], stempos[1][0],stempos[1][1], penalty,
                indvpen[0], indvpen[1]);
    }
    else if( type=='L' )
    {
        printf("Score:%10.6f  at [%d, %d]\n",
                prob, looppos[0], looppos[1]);
    }
    //printf("Result: \n");
    printf("   ");
    printf("%s\n", finseq );
    printf("   ");
    printf("%s\n\n", finstr );
}

void Candid::printAlgnedSeqOnly(  )
{
    #ifdef DEBUG_DEV
        printf("Candid::printAlgnedSeqOnly\n");
    #endif

    int lenlm=0, lenrm=0;
    char *seqlm=NULL, *seqrm=NULL, *alglm=NULL, *algrm=NULL;
    
    if( type!='S' )
        return;
    printf("Score:%10.6f   at [%d, %d]-[%d, %d]\n",
                prob,  stempos[0][0], stempos[0][1], stempos[1][0],stempos[1][1]);
    //printf("Result: \n");
    //determine the local positions
    extractAlgPart( seqlm, alglm, lenlm,  seqrm, algrm, lenrm );
    //output
    printf("   ");
    printf("%s ..... %s\n", seqlm, seqrm );
    printf("   ");
    printf("%s ..... %s\n", alglm, algrm );
    printf("\n");
    
    if( seqlm!=NULL )
    {
        delete [ ] seqlm;
        seqlm=NULL;
    }
    if( seqrm!=NULL )
    {
        delete [ ] seqrm;
        seqrm=NULL;
    }
    if( alglm!=NULL )
    {
        delete [ ] alglm;
        alglm=NULL;
    }
    if( algrm!=NULL )
    {
        delete [ ] algrm;
        algrm=NULL;
    }
}

//********************************************************************************
//*******************    implementation of CandScanner class      ****************
//********************************************************************************
CandScanner::CandScanner( int topn,  double thres )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::CandScanner( parameter )\n");
    #endif
    searchSeq=NULL;
    topnum=topn;
    newnum=INVALID;
    shiftnum=1;
    threshold=thres;
    pCand=NULL;
    
    allowOffset = ENOFFSET;
    addNullCand = false;
}

CandScanner::~CandScanner(  )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::~CandScanner\n");
    #endif
    if( searchSeq!=NULL )
    {
		//detect memory leak problem in rnatops zbhuang 20090815
        //delete searchSeq;
		delete [] searchSeq;
        searchSeq=NULL;
    }
}

int CandScanner::setTopn( int topn )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::setTopNumber\n");
    #endif
    topnum=topn;
    return 1;
}

int CandScanner::setThreshold( double thres )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::setThreshold\n");
    #endif
    threshold=thres;
    return 1;
}

int CandScanner::setShiftNum( int num )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::setThreshold\n");
    #endif
    shiftnum=num;
    return 1;
}

void CandScanner::setOffset( bool enable)
{
    #ifdef DEBUG_DEV
        printf("CandScanner::setOffset\n");
    #endif
    allowOffset=enable;
}

void CandScanner::addNullCandidate( bool enNullCandidate )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::addNullCandidate\n");
    #endif
    addNullCand=enNullCandidate;
}

int CandScanner::setSeqFromFile( char *filename )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::setSeqFromFile\n");
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
    seqlen=templen;
    return templen;
}

void CandScanner::printSearchSeq( )
{
    int i=0;
    printf("The sequence to search:\n");
    printf("   ");
    for( i=0; i<seqlen; i++ )
    {
        printf("%c", numtochar( searchSeq[i] ) );
    }
    printf("\n");
}

int CandScanner::setSeqFromArray( char *seqArray )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::setSeqFromArray\n");
    #endif
    int i=0, templen=0, len=0;
    char ch;

    //----calculate the length of seq
    templen = strlen(seqArray);
    for(i=0; i<templen; i++ )
    {
        ch=seqArray[i];
        if( isRNABase( ch ) )
            len++;
    }
    
    searchSeq= new int[len];
    for(i=0, len=0; i<templen; i++ )
    {
        ch=seqArray[i];
        if( isRNABase( ch ) )
            searchSeq[len++]=chartonum(ch);;
    }
    seqlen=len;
    return len;
}

int CandScanner::isDifferent( CYKCandidate *pCC, int basepos )
{

    Candid *curCand=pCand;
    while( curCand!=NULL )
    {
        if(  (curCand->prob == pCC->prob)
        	 &&(curCand->stempos[0][0] == pCC->loc[0][0] + basepos)
           &&(curCand->stempos[0][1] == pCC->loc[0][1] + basepos)
           &&(curCand->stempos[1][0] == pCC->loc[1][0] + basepos)
           &&(curCand->stempos[1][1] == pCC->loc[1][1] + basepos) )
        {
            return 0;
        }
         curCand=curCand->nextCand;
     }
    return 1;    
}

int CandScanner::isDifferent( VTBCandidate *pVC, int basepos )
{

    Candid *curCand=pCand;
    while( curCand!=NULL )
    {
        if(  (curCand->prob == pVC->prob)
        	 &&(curCand->looppos[0] == pVC->loc[0] + basepos)
           &&(curCand->looppos[1] == pVC->loc[1] + basepos) )
        {
            return 0;
        }
         curCand=curCand->nextCand;
     }
    return 1;    
}

int CandScanner::putStemCandInChain( Candid* &pTail,  CYKCandidate *pResult, 
                                     int stemid,  int size,   int bpos,  int uid )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::putStemCandInChain\n");
    #endif

    int i=0, j=0;
    Candid *curEC=NULL, *preEC=NULL;
    
    for( i=0; i<size; i++ )
    {
        if( pResult[i].prob <= threshold )
            continue;
        if( pResult[i].hasData( ) && isDifferent(&pResult[i], bpos) )
        {
            curEC=new Candid( 'S', stemid, uid++ );
 	           curEC->copyValue( pResult[i], bpos );
            
            preEC=curEC;
            if( pCand==NULL )
                 pCand=curEC;
            else
                pTail->nextCand=curEC;
            break;
        }
    }
    for( j=i+1; j<size; j++ )
    {
        if( pResult[j].prob <= threshold )
            continue;
        if( pResult[j].hasData( ) && isDifferent(&pResult[j], bpos) )
        {
            curEC=new Candid( 'S', stemid, uid++ );
 	           curEC->copyValue( pResult[j], bpos );
            
    	       preEC->nextCand=curEC;
    	       preEC=curEC;
    	 	}
    }
    if( curEC!=NULL )
        pTail=curEC;
    return uid;
}

int CandScanner::putLoopCandInChain( Candid* &pTail,  VTBCandidate *pResult, 
                                     int loopid,  int size,   int bpos,  int uid )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::putLoopCandInChain\n");
    #endif

    int i=0, j=0;
    Candid *curEC=NULL, *preEC=NULL;
    
    for( i=0; i<size; i++ )
    {
        if( pResult[i].prob <= threshold )
            continue;
        if( pResult[i].hasData( )&& isDifferent(&pResult[i], bpos) )
        {
            curEC=new Candid( 'L', loopid, uid++ );
 	           curEC->copyValue( pResult[i], bpos );
            
            preEC=curEC;
            if( pCand==NULL )
                 pCand=curEC;
            else
                pTail->nextCand=curEC;
            break;
        }
    }
    for( j=i+1; j<size; j++ )
    {
        if( pResult[j].prob <= threshold )
            continue;
        if( pResult[j].hasData( ) && isDifferent(&pResult[j], bpos) )
        {
            curEC=new Candid( 'L', loopid, uid++ );
 	           curEC->copyValue( pResult[j], bpos );
            
    	       preEC->nextCand=curEC;
    	       preEC=curEC;
    	 	}
    }
    if( curEC!=NULL )
        pTail=curEC;
    return uid;
}

int CandScanner::getNumOfCandids( Candid * phead )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::getNumOfCandids\n");
    #endif
    int num=0;
    Candid *curCand=phead;
    while( curCand!=NULL )
    {
        num++;
         curCand=curCand->nextCand;
     }
    return num;
}

int CandScanner::getTopnCandids( Candid * phead, int compid  )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::getTopnCandids\n");
    #endif
    int i=0, j=0;

    double cutoff=INVALID;
    MaxElems *topelem=new MaxElems[ newnum ];
    
    for( i=0; i<newnum; i++ )
    {
        topelem[i].prob=INVALID;
        topelem[i].uid=INVALID;
    }
    
    Candid *curCand=phead;
    while( curCand!=NULL )
    {
        if( curCand->prob > cutoff )
        {
            for( i=0; i<newnum; i++ )
            {
                if( curCand->prob > topelem[i].prob )
                {
                    //insert this one here
                    for( j=newnum-1; j>i; j-- )
                    {
                        topelem[j].uid=topelem[j-1].uid;
                        topelem[j].prob=topelem[j-1].prob;
                    }
                    topelem[i].uid=curCand->UID;
                    topelem[i].prob=curCand->prob;
                    break;
                }
            }
        }
        cutoff=topelem[newnum-1].prob;
        curCand=curCand->nextCand;
    }
    
    //get the real ones of class Candid 
    curCand=phead;
    while( curCand!=NULL )
    {
        if( curCand->prob >= cutoff )
        {
            for( i=0; i<newnum; i++ )
            {
                if( curCand->UID==topelem[i].uid )
                {
                     pCand[i]=*curCand;
                     pCand[i].CompID=compid;
                     pCand[i].UID=i;
                     pCand[i].nextCand=NULL;
                    break;
                }
          	}
        }
        curCand=curCand->nextCand;
    }
    //printCandid( pCand, topnum );
    delete [ ] topelem;
    return 1;
}

void CandScanner::freeCandidChain( Candid* & phead )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::freeCandidChain\n");
    #endif
    Candid *curCand=phead, *nxCand=NULL;
    while( curCand!=NULL )
    {
        nxCand=curCand->nextCand;
        delete curCand;
         curCand=nxCand;
     }
     phead=NULL;
}

int CandScanner::getTopnum( )
{
    return newnum;    
}

int CandScanner::getSeqlen( )
{
    return seqlen;    
}

int CandScanner::hasData( )
{
    if( pCand!=NULL && pCand->prob!=INVLDPROB )
        return 1;
    else
        return 0;
}

int CandScanner::resetMembers( )
{
    pCand=NULL;
    newnum=0;
    return 1;
}

void CandScanner::printCandid( Candid *phead, int numcd )
{
    int i=0;
    Candid *curCand=phead;
    
    printf("**************Debug output*****************\n" );
    if( phead->nextCand==NULL )
    {
        for( i=0; i<numcd; i++ )
        {	
            printf("=====>No.%d/%d\n",i+1, numcd );
            phead[i].printAlgnedSeq( );
         }
     }
     else
     {
        while( curCand!=NULL )
        {	
            printf("==>No.%d/%d\n", ++i, numcd );
            curCand->printAlgnedSeq( );
             curCand=curCand->nextCand;
         }
     }
}

int CandScanner::resizeCandArray(  )
{
    int i=0, oneNull=0;
    Candid * phead=NULL;
    newnum=0;
    
    for( i=0; i<topnum; i++ )
    {
        if( pCand[i].prob > threshold )
            newnum++;
    }
    if ( addNullCand ) {
        oneNull=1;
    }
    phead=pCand;
    pCand=NULL;
    newnum = newnum>topnum?topnum:newnum;
    
    pCand=new Candid[ newnum + oneNull];
    if ( addNullCand ) {
        pCand[ newnum ].setNullCandid( );
    }

    for( i=0; i<newnum; i++ )
    {
        pCand[i]=phead[i];
        pCand[i].UID=i;
    }
    newnum+=oneNull; //add one more for null candidate

    delete [ ] phead;
    return newnum;
}

void CandScanner::putTopStemCandInArray( CYKCandidate *pResult, int size, int bpos )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::putTopStemCandInArray\n");
    #endif
    int i=0, j=0, k=0;

    
    for( i=0; i<size; i++ )
    {
        if( pResult[i].hasData( ) )
        {
            if( pResult[i].prob <= threshold )
                continue;
            if( pResult[i].prob <= pCand[topnum-1].prob )
                continue;
            for( j=0; j<topnum; j++ )
            {
                  //(pCand[j].prob ==pResult[i].prob)&&
                 if(  (pCand[j].stempos[0][0] == pResult[i].loc[0][0] + bpos)
                   &&(pCand[j].stempos[0][1] == pResult[i].loc[0][1] + bpos)
                   &&(pCand[j].stempos[1][0] == pResult[i].loc[1][0] + bpos)
                   &&(pCand[j].stempos[1][1] == pResult[i].loc[1][1] + bpos) )
                    break;
                
                if( pResult[i].prob > pCand[j].prob )
                {
                    for( k=topnum-1; k>j; k-- )
                        pCand[k]=pCand[k-1];
                    pCand[j].copyValue( pResult[i], bpos );
                    break;
                }
            }
        }
    }
}

void CandScanner::putTopLoopCandInArray( VTBCandidate *pResult, int size, int bpos )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::putTopLoopCandInArray\n");
    #endif
    int i=0, j=0, k=0;
    
    for( i=0; i<size; i++ )
    {
        if( pResult[i].hasData( ) && isDifferent(&pResult[i], bpos) )
        {
            if( pResult[i].prob <= threshold )
                continue;
            if( pResult[i].prob <= pCand[topnum-1].prob )
                continue;
            for( j=0; j<topnum; j++ )
            {
                if(  (pCand[j].prob == pResult[i].prob)
                   &&(pCand[j].looppos[0] == pResult[i].loc[0] + bpos)
                   &&(pCand[j].looppos[1] == pResult[i].loc[1] + bpos) )
                    break;
                if( pResult[i].prob > pCand[j].prob )
                {
                    for( k=topnum-1; k>j; k-- )
                        pCand[k]=pCand[k-1];
                    pCand[j].copyValue( pResult[i], bpos );
                    break;
                }
            }
        }
    }
}

Candid * CandScanner::scanTopnInOneRegionSP( Stem *pstem, CYKSearch &curcyk, double coeff, int runned )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::scanTopnInOneRegionSP( Stem )\n");
    #endif
    int i=0;
    int offset[4]={0, 0, 0, 0};
    int scanlen=0;
    int *seqseg=NULL;
    
    if( pstem==NULL )
        return NULL;
    if( topnum <= 0 )
    {
        printf("Please set the top\n");
        exit(0);
    }
    CYKCandidate *pResult=NULL;

    pCand=new Candid[ topnum ];
    for( i=0; i<topnum; i++ )
    {
        pCand[i].type='S';
        pCand[i].CompID=pstem->UID;
        pCand[i].UID=i;
    }
    
    if( allowOffset )  //disable or enable Offset
    {
        pstem->getBothArmRegion( -1,  offset,  coeff  );
        scanlen = offset[3]-offset[0]+1;
    }
    else
    {
        offset[0]=0;
        offset[3]=0;
        scanlen = seqlen;
    }

    seqseg=new int[ scanlen ]; 
    for( i=0; i<scanlen; i++ )  //the number of windows
    {
        seqseg[i]=searchSeq[ offset[0] + i ];
    }
    curcyk.CYKSearchString( seqseg, scanlen, pstem, runned, shiftnum );

    pResult = curcyk.getCandHead( );
    putTopStemCandInArray( pResult, topnum, offset[0] );
    resizeCandArray(  );
    
    if( seqseg!=NULL )
    {
        delete [ ] seqseg;
        seqseg=NULL;
    }
    return pCand;
}

//topinWin: the number of top n candidates in a scanning window
Candid * CandScanner::scanTopnInTwoRegion( Stem *pstem,  CYKSearch & curcyk, double coeff, int basepos )
{
    #ifdef DEBUG_DEV
        printf("CandScanner::scanTopnInTwoRegion( Stem )\n");
    #endif
    int i=0, j=0, numreg=0, scanlen=0, status=0;
    int *seqseg=new int[ seqlen ];
    
    if( pstem==NULL )
        return NULL;
    if( topnum <= 0 )
    {
        printf("Please set the topn\n");
        exit(0);
    }

    CYKCandidate *pResult=NULL;

    pCand=new Candid[ topnum ];
    for( i=0; i<topnum; i++ )
    {
        pCand[i].type='S';
        pCand[i].CompID=pstem->UID;
        pCand[i].UID=i;
    }
    
    //**** beg ****
    int region[4]={0, 0, 0, 0};
    int midreg[2]={0, 0};
    int offset[2]={0, 0};
    numreg = pstem->getNumPairRegion(  );
    for( j=0; j<numreg; j++ )
    {
        pstem->getBothArmRegion( j, region, coeff );
        for( i=0; i<4; i++ )
            region[i] += basepos;

        offset[0]= region[0]<0 ? 0 : region[0];
        if( offset[0]>seqlen-2 )
            continue;
        
        midreg[0]= region[1]<=seqlen-1 ? region[1] : seqlen-1;
        if( midreg[0] < 0 )
            continue;

        midreg[1]= region[2]<0 ? 0 : region[2];
        if( midreg[1]>seqlen-1 )
            continue;
        
        offset[1]=region[3]<seqlen-1 ? region[3] : seqlen-1;
        if( offset[1] < 0 )
            continue;
//		printf("Stem %c, region %d, [%d-%d][%d-%d] (Region:[%d-%d][%d-%d])\n", pstem->idchar, j, offset[0], midreg[0], midreg[1], offset[1], region[0], region[1], region[2], region[3] );
        /*for(i=offset[0]; i<=midreg[0]; i++ )   printf("%c", numtochar(searchSeq[i]) );
        printf("   ");
        for(i=midreg[1]; i<=offset[1]; i++ )   printf("%c", numtochar(searchSeq[i]) );
        */
        if( allowOffset )  //disable or enable Offset
        {
            midreg[0] = midreg[0] - offset[0]+1;
            midreg[1] = midreg[1] - offset[0]-1;
            scanlen = offset[1]-offset[0]+1;
        }
        else
        {
            offset[0]=0;
            offset[1]=0;
            scanlen = seqlen;
        }
        
        for( i=0; i<scanlen; i++ )  {
            seqseg[i]=searchSeq[ offset[0] + i ];
            //printf("%c", numtochar(seqseg[i]) );
        }
        //printf("\n");

        status = curcyk.CYKSearchTwoRegion( seqseg, scanlen, midreg, pstem, j );
        if( status == 0 )
            continue;
        //**** end ****
        //curcyk.printResultSeq( );
        pResult = curcyk.getCandHead( );
        putTopStemCandInArray( pResult, topnum, offset[0] );
    }
    resizeCandArray(  );

    if( seqseg!=NULL )
    {
        delete [ ] seqseg;
        seqseg=NULL;
    }
    return pCand;
}

//deprecated
Candid * CandScanner::quickScanCands( Stem * pstem,  CYKSearch& curcyk, double coef, int runned )
{
    Candid * pcd=NULL;
    if( pstem==NULL )
        return NULL;
    if( topnum <= 0 )
    {
        printf("Please set the top\n");
        exit(0);
    }
    pcd = scanTopnInOneRegionSP( pstem, curcyk, coef, runned  );
    return pcd;
}

//deprecated
Candid * CandScanner::scanCandidates( Stem * pstem,  CYKSearch& curcyk,  double coef )
{
    Candid * pcd=NULL;
    if( pstem==NULL )
        return NULL;
    if( topnum <= 0 )
    {
        printf("Please set the top\n");
        exit(0);
    }
    pcd = scanTopnInTwoRegion( pstem,  curcyk, coef  );
    return pcd;
}

//********************************************************************************
//*******************    implementation of CandCollect class      ****************
//********************************************************************************
CandCollect::CandCollect( )
{
    int i=0;
    stemnum=INVALID;
    loopnum=INVALID;
    coefStDev=COESTDEV;
    stemgrp=NULL;
    loopgrp=NULL;
    flagrun=NULL;
    cyks=NULL;
    for( i=0; i<4; i++ )
    {
        stemerr[i]=0;
        hopepos[i]=0;
    }
    rdmlen=-1;
    
    benOffset = ENOFFSET;
    candidmerge=CYKCANDFILTER;
    candshortlen=PREFERSHORTLEN;
    allowshift=STEMFTALLOW;
    scoredroprate=PERCENTSCOREDROP;
    
    //length penalty
    plbd=1.0;
    pcoef=2.0;
    //set method for scanning
    bTwoRegion = true;
    nullcandidate = false;
}

CandCollect::CandCollect( int stems,  int loops )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::CandCollect\n");
    #endif
    int i=0;
    stemnum=stems;
    loopnum=loops;
    coefStDev=COESTDEV;
    flagrun=NULL;
    cyks=NULL;
    if( stems>0 )
        stemgrp=new CandContain[stemnum];
    else
        stemgrp=NULL;
        
    if( loops>0 )
        loopgrp=new CandContain[loopnum];
    else
        loopgrp=NULL;
        
    for(i=0; i<stemnum; i++ )
    {
        stemgrp[i].type=' ';
        stemgrp[i].idchar=' ';
        stemgrp[i].UID=INVALID;
        stemgrp[i].num=INVALID;
        stemgrp[i].pMember=NULL;
        stemgrp[i].hour=0;
        stemgrp[i].minu=0;
        stemgrp[i].secd=0.0;
    }
    for(i=0; i<loopnum; i++ )
    {
        loopgrp[i].type=' ';
        loopgrp[i].idchar=' ';
        loopgrp[i].UID=INVALID;
        loopgrp[i].num=INVALID;
        loopgrp[i].pMember=NULL;
        loopgrp[i].hour=0;
        loopgrp[i].minu=0;
        loopgrp[i].secd=0.0;
    }
    for( i=0; i<4; i++ )
    {
        stemerr[i]=0;
        hopepos[i]=0;
    }
    rdmlen=-1;
    
    benOffset = ENOFFSET;
    candidmerge=CYKCANDFILTER;
    candshortlen=PREFERSHORTLEN;
    allowshift=STEMFTALLOW;
    scoredroprate=PERCENTSCOREDROP;
    
    //length penalty
    plbd=1.0;
    pcoef=2.0;
    
    //set method for scanning
    bTwoRegion = true;
    nullcandidate = false;
}

CandCollect::~CandCollect( )
{
    if( stemgrp!=NULL )
    {
        delete [ ] stemgrp;
        stemgrp=NULL;
    }
    if( loopgrp!=NULL )
    {
        delete [ ] loopgrp;
        loopgrp=NULL;
    }
    if( flagrun!=NULL )
    {
        delete [ ] flagrun;
        flagrun=NULL;
    }
    if( cyks!=NULL)
    {
        delete [ ] cyks;
        cyks=NULL;
    }
}

void CandCollect::setTwoRegionScan( bool scanMethod )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::setCoeffStandardDev\n");
    #endif
    bTwoRegion = scanMethod;
}

void CandCollect::setCoeffStandardDev( double coeffsd )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::setCoeffStandardDev\n");
    #endif
    coefStDev = coeffsd;
}

int CandCollect::hasDataInIdx( char type,  int idx)
{
    #ifdef DEBUG_DEV
        printf("CandCollect::hasDataInIdx\n");
    #endif
    int retval=-1;
    if( type=='S')
    {
        if( stemgrp[idx].UID==INVALID )
            retval= 0;
        else
            retval= 1;
    }
    else if( type=='L')
    {
        if( loopgrp[idx].UID==INVALID )
            retval= 0;
        else
            retval= 1;
    }
    
    return retval;
}

//B: stem and loop;   S: stem;   L:loop
void CandCollect::printCandids( char type, int displayfmt )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::printCandids\n");
    #endif
    int i=0, j=0;
    if( type=='S' || type=='B' )
    {
        for(i=0; i<stemnum; i++ )
        {
            //printf("================== Candidates of Stem %d ( %d:%d:%4.2f ) ===============\n", 
            //        i, stemgrp[i].hour, stemgrp[i].minu, stemgrp[i].secd);
            printf("================== Candidates of Stem %s (%d)  ==================\n", stemgrp[i].idstr.data( ), i );
            if( stemgrp[i].pMember==NULL )
                printf("  Skip without scanning.\n");
            else
            {
                if( stemgrp[i].num == 0)
                    printf("  No candidate\n");
                else
                {
                    for(j=0; j<stemgrp[i].num; j++ )
                    {
                        printf("No.%d   ", j+1 );
                        if( displayfmt==1 )
                            stemgrp[i].pMember[j].printAlgnedSeq( );
                        else
                            stemgrp[i].pMember[j].printAlgnedSeqOnly( );
                    }
                }
            }
        }
    }
    
    if( type=='L' || type=='B' )
    {
        for(i=0; i<loopnum; i++ )
        {
            //printf("================== Candidates of Loop %d ( %d:%d:%4.2f ) ===============\n", 
            //        i, loopgrp[i].hour, loopgrp[i].minu, loopgrp[i].secd);
            printf("================== Candidates of Loop %d ==================\n",  i);
            if( loopgrp[i].pMember==NULL )
                printf("  Skip without scanning.\n");
            else
            {
                if( loopgrp[i].num == 0)
                    printf("  No candidate\n");
                else
                {
                    for(j=0; j<loopgrp[i].num; j++ )
                    {
                        printf("No.%d\n", j+1 );
                        loopgrp[i].pMember[j].printAlgnedSeq( );
                    }
                }
            }
        }
    }
}

void CandCollect::freeAllCandids( )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::freeAllCandids\n");
    #endif

    freeInsideCands(  );
    if( stemgrp!=NULL )
    {
        delete [ ] stemgrp;
        stemgrp=NULL;
    }
    
    if( loopgrp!=NULL )
    {
        delete [ ] loopgrp;
        loopgrp=NULL;
    }
}

void CandCollect::freeInsideCands(  )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::freeInsideCands\n");
    #endif
    int i=0;
    if( stemgrp!=NULL )
    {
        for(i=0; i<stemnum; i++ )
        {
            if( stemgrp[i].pMember!=NULL )
            {
                delete [ ] stemgrp[i].pMember;
                stemgrp[i].pMember=NULL;
            }
            stemgrp[i].type=' ';
            stemgrp[i].idchar=' ';
            stemgrp[i].UID=INVALID;
            stemgrp[i].num=INVALID;
            stemgrp[i].hour=0;
            stemgrp[i].minu=0;
            stemgrp[i].secd=0.0;
        }
    }
    if( loopgrp!=NULL )
    {
        for(i=0; i<loopnum; i++ )
        {
            if( loopgrp[i].pMember!=NULL )
            {
                delete [ ] loopgrp[i].pMember;
                loopgrp[i].pMember=NULL;
            }
            loopgrp[i].type=' ';
            loopgrp[i].idchar=' ';
            loopgrp[i].UID=INVALID;
            loopgrp[i].num=INVALID;
            loopgrp[i].hour=0;
            loopgrp[i].minu=0;
            loopgrp[i].secd=0.0;
        }
    }
}

int CandCollect::addToCollect( char LorS, string idstr, int idx, Candid * &phead, int curnum, TimeCalc* timc )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::addToCollect\n");
    #endif
    
    if( stemnum==INVALID && loopnum==INVALID )
    {
        printf("Please initialze the instance of CandCollect\n");
        exit( 0 );
    }

    if( hasDataInIdx( LorS, idx ) )
    {
        printf("Error: The element at %d[%c] already has real data in addToCollect.\n", idx, phead->type);
        return 0;
    }
    else
    {
        if( LorS=='S' )
        {
            stemgrp[idx].pMember=phead;
            stemgrp[idx].type='S';
            stemgrp[idx].idchar=idstr.at(0);
            stemgrp[idx].idstr=idstr;
            stemgrp[idx].UID=idx;
            stemgrp[idx].num=curnum;
            if( timc != NULL )
                timc->getTime( stemgrp[idx].hour, stemgrp[idx].minu, stemgrp[idx].secd);
        }
        else if( LorS=='L' )
        {
            loopgrp[idx].pMember=phead;
            loopgrp[idx].type='L';
            //loopgrp[idx].idchar=idchar;
            loopgrp[idx].UID=idx;
            loopgrp[idx].num=curnum;
            if( timc!=NULL )
                timc->getTime( loopgrp[idx].hour, loopgrp[idx].minu, loopgrp[idx].secd);
        }
        else
        {
            printf("Unknown type for Stem or Loop\n");
        }
    }
    return 1;
}

void CandCollect::setStemPosAccuStd(int* posErr,  int* armlen, int halfradmlen ) //0-based
{
    #ifdef DEBUG_DEV
        printf("CandCollect::setStemPosAccuStd\n");
    #endif
    int i=0;
    for( i=0; i<4; i++)
        stemerr[i]=posErr[i];

    hopepos[0]=halfradmlen;
    hopepos[1]=halfradmlen+armlen[0]-1;
    hopepos[2]=tgtseqlen-halfradmlen-armlen[1];
    hopepos[3]=tgtseqlen-halfradmlen-1;
    
    rdmlen=halfradmlen;
}

void CandCollect::checkStemPosAccuracy( Stem *pstem, int* accuracy )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::checkStemPosAccuracy\n");
    #endif
    
    int i=0, j=0;
    
    if( stemnum!=1 )
        return;
    for( i=0; i<stemnum; i++ )
    {
        if( stemgrp[i].pMember==NULL )
            continue;
        for( j=0; j<stemgrp[i].num; j++ )
        {
            if(    (stemgrp[i].pMember[j].stempos[0][0]<= hopepos[0] + stemerr[0])
                && (stemgrp[i].pMember[j].stempos[0][0]>= hopepos[0] - stemerr[0])
                
                && (stemgrp[i].pMember[j].stempos[0][1]<= hopepos[1] + stemerr[1])
                && (stemgrp[i].pMember[j].stempos[0][1]>= hopepos[1] - stemerr[1])
                
                && (stemgrp[i].pMember[j].stempos[1][0]<= hopepos[2] + stemerr[2])
                && (stemgrp[i].pMember[j].stempos[1][0]>= hopepos[2] - stemerr[2])
                
                && (stemgrp[i].pMember[j].stempos[1][1]<= hopepos[3] + stemerr[3])
                && (stemgrp[i].pMember[j].stempos[1][1]>= hopepos[3] - stemerr[3]) )
             {
                   accuracy[j]++;
                   break;
             }
        }
    }
}

//pos: int pos[4]
int CandCollect::getAlignmentStr( int stemNo, int candNo,
                                  char* &seqlm,   char* &alglm,   int &lenlm,
                                  char* &seqrm,   char* &algrm,   int &lemrm )
{
    Candid *pcand=NULL;
    
    if( stemNo >= stemnum || stemNo < 0 )
    {
        printf("The index of stem is out of range: 0-%d.\n", stemnum-1);
        return -1;
    }
    if( candNo >= stemgrp[ stemNo ].num  || candNo < 0 )
    {
        printf("The index of candid for stem %c is out of range: 0-%d.\n", stemgrp[stemNo].idchar, stemgrp[stemNo].num-1 );
        return -1;
    }
    pcand = &(stemgrp[stemNo].pMember[candNo]);
    pcand->extractAlgPart( seqlm, alglm, lenlm,   seqrm, algrm, lemrm );
    return 1;
}

//parameter-setting function
void CandCollect::initCYKInstance(int topsize, double coeff, int varvalue )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::initCYKInstance\n");
    #endif
    int i=0;
    int deftop=0;

    if( flagrun==NULL )
    {
        flagrun=new int[ stemnum ];
        for( i=0; i<stemnum; i++ )
            flagrun[i]=0;
    }
    
    if( cyks==NULL )
    {
        deftop=topsize>0?topsize:CANDIDNUM;
        cyks = new CYKSearch[ stemnum ];
        for( i=0; i<stemnum; i++ )
        {
            cyks[i].setParameter( deftop, coeff, varvalue );
            cyks[i].setCandMergParameter( candidmerge, candshortlen, allowshift, scoredroprate);
            cyks[i].setLengthPenaltyParameters( plbd, pcoef );
        }
    }
}

void CandCollect::enableOffset( bool offsetvalid )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::enableOffset\n");
    #endif
    benOffset = offsetvalid;
}

//set the candidate-merging parameters
//mergevalid: whether candidate-merging is used.
//shortlen:   whether the shortest stem is selected into top n. 
//shiftnum:   how many shifting positions are allowed
//droprate:   reserved( just set it to be -1 )
void CandCollect::setCYKMergeParamter(bool mergevalid, bool shortlen, int shiftnum, double droprate )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::setCYKMergeParamter\n");
    #endif
    
    candidmerge = mergevalid;
    candshortlen = shortlen;
    allowshift = shiftnum;
    scoredroprate = droprate;
}

void CandCollect::setCYKLengthPenalty( double lbound, double adjustk )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::setCYKLengthPenalty\n");
    #endif
    //length penalty
    plbd = lbound;
    pcoef = adjustk;
}

void CandCollect::setNullCandidate( bool enableNull )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::setNullCandidate\n");
    #endif
    nullcandidate = enableNull;
}

//deprecated
int CandCollect::quickStemCands( Stem *pstem, int numsm, int topn,  double thres,  char * pStr)
{
    #ifdef DEBUG_DEV
        printf("CandCollect::quickStemCands\n");
    #endif
    
    int i=0;
    TimeCalc durtime;
    Candid *pCD=NULL;
    CandScanner scan( topn, thres );
    scan.setOffset( benOffset );
    if( isFileName( pStr )  )
        scan.setSeqFromFile( pStr );
    else
        scan.setSeqFromArray( pStr );
    freeInsideCands(  );
    
    //initialize the speedup variables
    if( flagrun==NULL )
        initCYKInstance( topn,  coefStDev, 0 );

    for( i=0; i<numsm; i++ )
    {
        durtime.resetTime( );
        durtime.countBegin( );
        pCD=scan.scanTopnInOneRegionSP( &pstem[i], cyks[i], coefStDev, flagrun[i] );
        flagrun[i]=1;
        durtime.countEnd( );
        addToCollect( 'S', pstem[i].idstr, i, pCD, scan.getTopnum( ), &durtime );
    }
    return 1;
}

//topn: select out the best n number of candidates
int CandCollect::searchCandids( Stem *pstem, int numsm,  int topn,  double thres,  char * pStr, int basepos )
{
    #ifdef DEBUG_DEV
        printf("CandCollect::searchCandids\n");
    #endif
    
    int i=0;
    TimeCalc durtime;
    Candid *pCD=NULL;
    CandScanner scan( topn, thres );
    scan.setOffset( benOffset );
    if( isFileName( pStr )  )
        scan.setSeqFromFile( pStr );
    else
        scan.setSeqFromArray( pStr );
    freeInsideCands(  );
    
    //initialize the speedup variables
    if( flagrun==NULL )
        initCYKInstance( topn,  coefStDev, 0 );

    for( i=0; i<numsm; i++ )
    {
        durtime.resetTime( );
        durtime.countBegin( );
        //cout<<endl <<"STEM: "<<pstem[i].idstr<<" "<<i<< endl;
		//if( i==65 )
		//	cout<<"Now here"<<endl;
        if( bTwoRegion ) {  //two-region method, but couldn't speedup
            //scan.addNullCandidate( nullcandidate );
            pCD=scan.scanTopnInTwoRegion( &pstem[i], cyks[i], coefStDev, basepos  );
        }
        else  { //one-region method, but can speedup
            pCD=scan.scanTopnInOneRegionSP( &pstem[i], cyks[i], coefStDev, flagrun[i] );
        }
        flagrun[i]=1;
        durtime.countEnd( );
        addToCollect( 'S', pstem[i].idstr, i, pCD, scan.getTopnum( ), &durtime );
        //cout<<"Finish "<<endl;
    }
    return 1;
}

