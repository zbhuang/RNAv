#pragma warning(disable: 4786)

#include "DataLoader.h"

//********************************************************************************
//*******************    implementation of StemLocation class       ****************
//********************************************************************************
StemLocation::StemLocation( )
{
    int i=0;
    
    charid=' ';
    for(i=0; i<4; i++)
        Pos[i]=-1;
    stemlen = 0;
    distance = 0;     //between two arms in PASTA
    numOffset = 0;
    for(i=0; i<2; i++) 
    {
        offsetE[i] = NULL;    //distance of stem end to pasta end 0: left 1: right
        offsetSD[i] = NULL;   //xxxxxAAAA......AAAAxxxxx
        armE[i]=0.0;   //for the length of two arms 0: left 1: right
        armSD[i]=0.0;
    }
    midE=0.0;      //for the length of middle region
    midSD=0.0;
}

StemLocation::~StemLocation( )
{
    releaseEandSDArray( );
}

void StemLocation::releaseEandSDArray( )
{
    int i=0;
    numOffset = 0;
    for (i=0; i<2; i++)
    {
        if( offsetE[i] != NULL)
        {
            delete [ ] offsetE[i];
            offsetE[i] = NULL;
        }
        if( offsetSD[i] != NULL)
        {
            delete [ ] offsetSD[i];
            offsetSD[i] = NULL;
        }
    }
}

StemLocation & StemLocation::operator=(  StemLocation &org )
{
    int i=0, j=0;
    
    charid=org.charid;
    strid=org.strid;
    for(i=0; i<4; i++)
        Pos[i]=org.Pos[i];
    stemlen = org.stemlen;
    distance = org.distance;     //between two arms in PASTA
    numOffset = org.numOffset;
    for(i=0; i<2; i++) 
    {
        armE[i]=org.armE[i];   //for the length of two arms 0: left 1: right
        armSD[i]=org.armSD[i];
        woffsetE[i]=org.woffsetE[i];
        woffsetSD[i]=org.woffsetSD[i];
    }
    midE=org.midE;      //for the length of middle region
    midSD=org.midSD;
    //release the old memory and allocate the new
    for(i=0; i<2; i++) 
    {
        if( offsetE[i]==NULL )
        {
            delete [ ] offsetE[i];
            offsetE[i]= NULL;
        }
        if( offsetSD[i]==NULL )
        {
            delete [ ] offsetSD[i];
            offsetSD[i]= NULL;
        }
        offsetE[i] = new double[ numOffset ];    //distance of stem end to pasta end 0: left 1: right
        offsetSD[i] = new double[ numOffset ];   //xxxxxAAAA......AAAAxxxxx
    }
    for(i=0; i<2; i++) 
    {
        for(j=0; j<numOffset; j++)
        {
            offsetE[i][j] = org.offsetE[i][j];    //distance of stem end to pasta end 0: left 1: right
            offsetSD[i][j] = org.offsetSD[i][j];   //xxxxxAAAA......AAAAxxxxx
        }
    }
    return *this;
}

int StemLocation::getArmRegionAt( int idx, int reg[4], double coeff )
{
    int i=0, retval=0;
    double dbeg=0, dend=0;
    if( idx>=numOffset || idx<0) {
        printf("Make sure the index within the range: [0 - %d]\n", numOffset-1);
        for( i=0; i<4; i++ )
            reg[i]=-1;
        retval=-1;
    }
    else {
        //left arm
        dbeg = offsetE[0][idx] - coeff*offsetSD[0][idx];
        dend = offsetE[0][idx] + armE[0] + coeff*(offsetSD[0][idx]+armSD[0]);
        reg[0]=(int)floor( dbeg );
        reg[1]=(int)ceil( dend );
        
        //right arm
        dbeg = offsetE[1][idx] - coeff*offsetSD[1][idx];
        dend = offsetE[1][idx] + armE[1] + coeff*(offsetSD[1][idx]+armSD[1]);
        reg[2]=(int)floor( dbeg );
        reg[3]=(int)ceil( dend );
        
        retval=1;
    }
    return retval;
}

int StemLocation::getArmRegionWhole( int reg[4], double coeff )
{
    double dbeg=0, dend=0;
    //left arm
    dbeg = woffsetE[0] - coeff*woffsetSD[0];
    dend = woffsetE[0] + armE[0] + coeff*(woffsetSD[0] + armSD[0]);
    reg[0]=(int)floor( dbeg );
    reg[1]=(int)ceil( dend );
       
    //right arm
    dbeg = woffsetE[1] - coeff*woffsetSD[1];
    dend = woffsetE[1] + armE[1] + coeff*(woffsetSD[1] + armSD[1]);
    reg[2]=(int)floor( dbeg );
    reg[3]=(int)ceil( dend );

    return 1;
}

int StemLocation::getStemRegion( int idx, int reg[4], double coeff )
{
    int status=0;
    if( idx<0 )
        status=getArmRegionWhole(  reg, coeff );
    else
        status=getArmRegionAt( idx, reg, coeff );
    return status;
}

int StemLocation::offsetAllE( double distance )
{
    int i=0;
    for(i=0; i<numOffset; i++ ) {
        if( offsetE[0][i] < 0 )
            offsetE[0][i] += distance;
        else
            offsetE[0][i] += distance-1;
            
        if( offsetE[1][i] < 0 )
            offsetE[1][i] += distance;
        else
            offsetE[1][i] += distance-1;
    }
    woffsetE[0] += distance;
    woffsetE[1] += distance-1;
    return 1;
}

int StemLocation::obtainMostEnds( char LorR, double coeff )
{
    int i=0, retDist=0;
    double endDist=0.0, curDist=0.0;
    if( LorR=='L' ) {
        for(i=0; i<numOffset; i++ ) {
            curDist = offsetE[0][i] - coeff*offsetSD[0][i];
            if( curDist < endDist )
                endDist = curDist;
        }
    }
    else {
        for(i=0; i<numOffset; i++ ) {
            curDist = offsetE[1][i] + armE[1] + coeff*(offsetSD[1][i]+armSD[1]);
            if( curDist > endDist )
                endDist = curDist;
        }
    }
    retDist = (int)ceil( abs( endDist ) );
    return retDist;
}

string StemLocation::getStemStrID( )
{
    return strid;
}

//********************************************************************************
//*******************    implementation of DataLoader class       ****************
//********************************************************************************

LoopLocation::LoopLocation( )
{
    int i=0;
    
    id=-1;
    for(i=0; i<2; i++)
        Pos[i]=-1;
    looplen = 0;
    for(i=0; i<2; i++) 
    {
        offsetE[i] = 0.0;    //distance of stem end to pasta end 0: left 1: right
        offsetSD[i] = 0.0;   
    }
    lenE=0.0;      //for the length of middle region
    lenSD=0.0;
}

LoopLocation::~LoopLocation( )
{
    //int i=0;
}
//********************************************************************************
//*******************    implementation of DataLoader class       ****************
//********************************************************************************
DataLoader::DataLoader(  )
{
    int i=0;
    stemchid=NULL;

    for( i=0; i<2; i++ )
        pasta[i]=NULL;
    mgpasta=NULL;
    trainingseq=NULL;
    seqlength=0;
    numtrainseqs=0;

    splitThreshold = 13.0;
    basePos=0;
    filename=NULL;
    hasCommentLines=0;
    stemary=NULL;
    numstems=0;
    loopary=NULL;
    numloops=0;
    for( i=0; i<4; i++ )
        BaseHz[i]=-1.0;
        
    filterPos[0]=0;
    filterPos[1]=0;
    extensionLength[0]=0;
    extensionLength[1]=0;

    //detect memory leak problem in rnatops (linux-valgrind) zbhuang 20090815
    this->numPastaLines = 0;
}

DataLoader::~DataLoader( )
{
    int i=0;
    
    if( filename!=NULL )
    {
        delete [ ] filename;
        filename=NULL;
    }
    
    for( i=0; i<numPastaLines; i++ )
    {
        if( pasta[i]!=NULL )
        {        
            delete [ ] pasta[i];
            pasta[i]=NULL;
        }
    }
        
    if( stemchid!=NULL )
    {
        delete [ ] stemchid;
        stemchid=NULL;
    }

    if( mgpasta!=NULL )
    {
        delete [ ] mgpasta;
        mgpasta=NULL;
    }
  
    if( trainingseq!=NULL )
    {
        for( i=0; i<numtrainseqs;  i++ )
        {
            delete[ ] trainingseq[i];
            trainingseq[i]=NULL;
        }
        delete [ ] trainingseq;
        trainingseq=NULL;
    }
    
    if( stemary!=NULL )
    {
        delete [ ] stemary;
        stemary=NULL;
    }
    if( loopary!=NULL )
    {
        delete [ ] loopary;
        loopary=NULL;
    }
}

int DataLoader::getSeqLength( )
{
    return seqlength;  
}

int DataLoader::getNumOfTrainSeq( )
{
    return numtrainseqs;
}

int ** DataLoader::getPointerOfSeqs( )
{
    return trainingseq;
}

char * DataLoader::getPointerOfPasta( )
{
    return pasta[0];
}

//get the size of str: len
int DataLoader::countSeqLength( vector <string> ts )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::countSeqLength\n");
    #endif
    int length = (int)ts[0].length();
    return length;
}

int DataLoader::is2ndPastaLineFormat( char curbase )
{
    int result=0;
    switch( curbase )
    {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case '.':
                result=1;
                break;
        default:
            	result=0;
            	break;
    }
    return result;
}

int DataLoader::countPastaLines( vector <string> ts )
{
    //check pastaline
    int numLines=2, j=0;
    char base=' ';
    
    if( (int)ts[1].length( ) != seqlength )
        numLines = 1;
    else {
    	for (j=0; j<seqlength; j++) {
    	    base=ts[1][j];
            if( is2ndPastaLineFormat( base ) == 0 ) {
                 numLines=1;
                 break;
            }
        }
    }
    return numLines;
}

int DataLoader::checkSeqLength( )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::checkSeqLength\n");
    #endif
    int i=0, templen=0;
    int size = orgtrains.size( );

	for (i=0; i<size; i++) {
        templen = (int)orgtrains[i].length( );
//        for(int j=0; j<templen; j++ )
//            printf("%c", numtochar( orgtrains[i][j] ) );
//        printf("\n");
	    if( templen != seqlength ) {
            printf("The length of Seq No.%d is %d, which doesn't equal to pasta line:%d\n", i, templen, seqlength);
            return 0;
        }
	}

	return 1;
}

void DataLoader::mergePastaLine(  )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::mergePastaLine\n");
    #endif
    
    int i=0, j=0;
    char ch[10];
    
    mgpasta = new string [ seqlength ];
    for( i=0; i<seqlength; i++ )
    {
        for( j=0; j<numPastaLines; j++ )
            ch[j]=toupper( pasta[j][i] );
        ch[j]='\0';
        mgpasta[i].assign( ch );
    }
}

int DataLoader::readOriginaldata( vector <string> & tsdata )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::readOriginaldata(  )\n");
    #endif
    int i=0, len=0;
    char * str=new char[ 50000 ];
    string rmstr;
    memset (str,'\0',50000 );

    //Open the training data file
    FILE * pfdata = fopen(filename, "r");
    if( pfdata == NULL ) {
        printf(" Fail to open %s\n", filename);
        exit( 0 );
    }
    
    fseek( pfdata, 0, SEEK_SET);
    fgets(str, 50000, pfdata);
    while( !feof( pfdata ) )
    {
        len=strlen(str);
        for(i=0; i<len; i++ ) {
            if( str[i]=='\r' || str[i]=='\n' )
                str[i]='\0';
        }
        if( isEmptyLine( str )==1 )
            continue;
        else {
            rmstr = removeSpaceEnd( str );
            tsdata.push_back( rmstr );
        }
        memset (str,'\0',50000 );
        fgets(str, 50000, pfdata);
    }
    if( isEmptyLine( str )!=1 ) {
        rmstr = removeSpaceEnd( str );
        tsdata.push_back( rmstr );
    }
/*    
     cout<<"xxxxxxxxxxxxxxxxx"<<endl;
     for( i=0; i<(int)tsdata.size( ); i++ )
         cout<<tsdata[i]<<endl;
     cout<<"xxxxxxxxxxxxxxxxx"<<endl;
*/     
    delete [ ] str;
    fclose( pfdata );
    return 1;
}

string DataLoader::removeSpaceEnd( char * str  )
{
    int i=0, j=0;
    int thelen = strlen( str );
    //remove the tail space
    for(i=thelen-1; i>=0; i--) {
        if( str[i]==' ' )
            str[i]='\0';
        else
            break;
    }
    //remove the heading space
    char * retbegpos = NULL;
    for(i=0, j=0; i<thelen; i++ ) {
        if( str[i]==' ' )
            continue;
        else {
            retbegpos = &( str[i] );
            break;
        }
    }
    string retval( retbegpos );
    // cout<<"remspace: "<< retval<<endl;
    return retval;
}

int DataLoader::isEmptyLine( char * str  )
{
    int i=0;
    int actlen=0;
    int orglen=strlen( str );
    for(i=0; i<orglen; i++ ) {
        if( str[i]!=' ' && str[i]!='\r' )
            actlen++;
    }
    if( actlen<=0 )
        return 1;
    else
        return 0;
}

int DataLoader::detectCommentLines( vector <string> ts )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::detectCommentLines\n");
    #endif
    int i=0, hasComment=0;
    int size = (int)ts.size( );
    for( i=numPastaLines; i<size; i++ ) {
        if( ts[i][0]=='>' ){
            hasComment=1;
            break; 
        }
        if( hasComment )
            break;
    }
    return hasComment;
}

int DataLoader::countTrainingSeqNum( vector <string> ts )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::countTrainingSeqNum\n");
    #endif
    int numts=0;
    int size = (int)ts.size( );

    if( hasCommentLines ) {
        numts = (int)(size-numPastaLines)/2;
    }
    else
        numts = size-numPastaLines;
    return numts;
}

int DataLoader::seperateCommentAndTrainingSeq( vector <string> ts )
{
    int i=0;
    // 1 ****** pasta line
    for( i=0; i<numPastaLines; i++ ) {
        orgpasta.push_back( ts[i] );
    }
    // 2 ****** comment lines
    for( i=numPastaLines; i<(int)ts.size( ); i++ )
    {
        if( ts[i][0]=='>')
            orgcommns.push_back( ts[ i ] );
        else
            orgtrains.push_back( ts[ i ] );
    }
    return 1;
}

int DataLoader::copyPastaLines( )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::copyPastaLines( )\n");
    #endif
    int i=0, j=0;
    //read pastaline
    for( i=0; i<numPastaLines; i++ ) {
        pasta[i] =new char[seqlength+1];
        if( pasta[i]==NULL ) {
            printf("Fail to allocate memory in DataLoader::copyPastaLines\n");
            return 0;
        }
        else {
        	for (j=0; j<seqlength; j++) {
                pasta[i][j]=orgpasta[i][j];
            }
            pasta[i][j]='\0';
        }
    }
    return 1;
}

int DataLoader::copyTrainingSeqs( )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::copyTrainingSeqs( )\n");
    #endif
    int j=0, k=0;
    trainingseq = new int* [ numtrainseqs ];
    if( trainingseq==0) {
        printf("Fail to allocate memory!\n");
        return 0;
    }
    
    //read training seqs
    for( k=0; k<numtrainseqs; k++ ) {
        trainingseq[k]= new int [seqlength];
        for (j=0; j<seqlength; j++) {
            trainingseq[k][j] = chartonum( orgtrains[k][j] );
        }
    }
    return 1;
}

void DataLoader::copyOrigSeqs( vector <string> ts )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::copyOrigSeqs( vector <string> )\n");
    #endif
    orgtrains=ts;
}

int DataLoader::inputData( char * infileName )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::inputData( char * )\n");
    #endif
    int namelen=0, status=0;
    vector <string> tsf;
    
    //Set filename
    namelen=strlen( infileName );
    filename = new char[namelen+1];
    strcpy( filename, infileName);

    //read training seqs
    readOriginaldata( tsf );
    //call vector InputData function
    status = inputData( tsf );

    return status;
}

int DataLoader::inputData( vector <string> ts )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::inputData( vector <string> )\n");
    #endif
    int status=0;
    
    //check length for hfstr
    if ( (int)ts.size()<= 1 ) {
        printf("There is NO training seqs\n");
        return 0;
    }
    
    seqlength = countSeqLength( ts );
    numPastaLines = countPastaLines( ts );
    //cout<< "numPastaLines: "<<numPastaLines <<endl;
    hasCommentLines = detectCommentLines( ts );
    //cout<< "has Commentline? "<<hasCommentLines <<endl;
    numtrainseqs = countTrainingSeqNum( ts );
    //cout<< "numtrainseqs: "<<numtrainseqs <<endl;
    if( seqlength==0 || numtrainseqs==0 ) {
        printf("Error: the length of some training seq=0 or the number of training seqs=0!\n");
        exit( 0 );
    }
    seperateCommentAndTrainingSeq( ts );
    //cout<<"numPastaLines:"<<numPastaLines << "  --- size: "<<orgpasta.size( )<<endl;
    //cout<<"numtrainseqs:"<<numtrainseqs << "  --- size: "<<orgtrains.size( )<<endl;
    //printStringAll( );
    status=checkSeqLength( );
    if( status==0 )
        exit( 0 );

    copyPastaLines( );
    copyTrainingSeqs( );
	
    mergePastaLine( );
    analyzePasta( );
    
    calcLoopStemStat( );
    groupSeqsOnFilter( );
    
    calcExtensionLength(  );
    
    return 1;
}

void DataLoader::printVectorStrs( vector <string> ts )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::printVectorStrs( )\n");
    #endif
    int size = (int) ts.size( );
    int i=0;
    cout<<"# of seqs: "<<size<<endl;
    for( i=0; i<size; i++ )
        cout<<ts[i]<<endl;
    cout<<endl;
}

void DataLoader::countBaseFreq( )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::countBaseFreq\n");
    #endif
    int i=0, j=0, allbasenum=0;
    int numbase[5]={0,0,0,0,0};

    for(i=0; i<numtrainseqs; i++ )
    {
        for(j=0; j<seqlength; j++)
        {
            numbase[ trainingseq[i][j] ]++;
        }
    }
    allbasenum = numbase[1]+numbase[2]+numbase[3]+numbase[4];
    for(i=1; i<5; i++ )
    {
        BaseHz[i-1]=( double )numbase[i]/allbasenum;
    }  
}

void DataLoader::setBaseFreq( double HzA, double HzC, double HzG,  double HzU)
{
    #ifdef DEBUG_DEV
        printf("DataLoader::setBaseFreq\n");
    #endif
    BaseHz[0]=HzA;
    BaseHz[1]=HzC;
    BaseHz[2]=HzG;
    BaseHz[3]=HzU;
}

void DataLoader::setSplitThreshold( double cutoff )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::setSplitThreshold\n");
    #endif
    splitThreshold = cutoff;
}
//if BaseHz is valid, copy log into the target array
void DataLoader::copyLogBgBaseHz( double * baseArray, int basenum )
{
    int i=0;
    //calculate the freq of ACUG basd on training seqs
    if( BaseHz[0]==-1.0)
         countBaseFreq( );
    
    //according to different basenum: 4 or 5
    if( basenum==4 )
    {
        for( i=0; i<4; i++ )
            baseArray[i]=log( BaseHz[i] );
    }
    else if( basenum==5 )
    {
        baseArray[0]=0.0;
        for( i=0; i<4; i++ ) 
            baseArray[i+1]=log( BaseHz[i] );
    }
    else
        printf("Error: copyLogBgBaseHz, basenum should be 4/5, but yours is %d\n", basenum);
}

int DataLoader::getStemNum(  )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::getStemNum\n");
    #endif
    
    int i=0, j=0, repeat=0;
    string *mark=NULL;
    int size=0;
    char pastalp[10];
    memset( pastalp, '\0', 10);
    for(i=0; i<numPastaLines; i++ )
        pastalp[i]=PASTALOOP;

    mark = new string[ seqlength ];
    
    for( i=0; i<seqlength; i++ )
    {
        repeat=0;
        if( mgpasta[i].compare( pastalp )==0 )
            continue;
        for(j=0; j<size; j++ )
        {
            if( mark[j].compare( mgpasta[i])==0 )
            {
                repeat=1;
                break;
            }
        }
        if( !repeat )
        {
            mark[size].assign( mgpasta[i] );
            size++;
        }
    }
    stemchid= new string[size];
    for( i=0; i<size; i++ )
    {
        stemchid[i]=mark[i];
    }
    //printf("%s   %s\n", mark, stemchid);
    delete [ ] mark;
    mark=NULL;

    return size;
}

int DataLoader::getEndLoopPos( int &nextBeg )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::getEndLoopPos\n");
    #endif
    int i=0;
    int curPos=nextBeg;
    int temp=seqlength-1;
    
    if( nextBeg >= seqlength )
        return seqlength-1;
    
    for( i=0; i<numstems; i++ )
    {
        //printf("%d: nextBeg=%d\n",i, nextBeg);
        
        if( ( curPos <= stemary[i].Pos[0] )&&( temp > stemary[i].Pos[0] ) )
        {
            temp=stemary[i].Pos[0];
            nextBeg=stemary[i].Pos[1]+1;
        }

        if( ( curPos <= stemary[i].Pos[2] )&&( temp > stemary[i].Pos[2] ) )
        {
            temp=stemary[i].Pos[2];
            nextBeg=stemary[i].Pos[3]+1;
        }
    }
    if( nextBeg == curPos )
        return seqlength-1;
    else
        return temp-1;
}

int DataLoader::locateLoopPosition(  )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::locateLoopPosition\n");
    #endif
    int i=0, lowbd=0;
    for( i=0; i<numloops; i++ )
    {
        //printf("loop: %d lowbd=%d\n", i, lowbd);
        loopary[i].Pos[0]=lowbd;
        loopary[i].Pos[1]=getEndLoopPos( lowbd );
        loopary[i].looplen=loopary[i].Pos[1] - loopary[i].Pos[0] + 1;
    }
    return 1;
}

//beginning position, ending position, distance
int DataLoader::locateStemPosition( )
{
    #ifdef DEBUG_DEV
       printf("DataLoader::locateStemPosition\n");
    #endif
    int i=0, j=0;
    int leftone=0, rightone=0;
    int ls=0, le=0, rs=0, re=0;
    int stemno=0, rlen=0, llen=0;
    string strid;

    for( stemno=0; stemno<numstems; stemno++ )
    {
        ls=-100;
        le=-100;
        rs=-100;
        re=-100;
        llen=0;
        rlen=0;
        strid=stemary[stemno].strid;
        i=0;
        j=seqlength-1;
        while( i<=j)
        {
            if( strid.compare( mgpasta[i] )==0 )
            {
                if( ls==-100 )
                    ls=i;
                leftone =1;
                llen++;
                le=i;
            }
            if( leftone )
            {
                while( !rightone && (j>i) )
                {
                    if( strid.compare( mgpasta[j] )==0 )
                    {
                        if( re==-100 )
                            re=j;
                        rightone =1;
                        rlen++;
                        rs=j;
                    }
                     j--;
                }
            }
            leftone=0;
            rightone=0;
            i++;
        }
        stemary[stemno].Pos[0]=ls;
        stemary[stemno].Pos[1]=le;
        stemary[stemno].Pos[2]=rs;
        stemary[stemno].Pos[3]=re;
        stemary[stemno].distance=rs-le-1;
        if( llen!=rlen)
        {
            printf("The length of both arms in stem \"%s\" doesn't equal.\n", strid.data( ) );
            exit(1);
        }
        else
        {
            stemary[stemno].stemlen=llen;
        }
    }
    #ifdef DEBUG_DEV
        for(i=0; i<numstems; i++)
        {
            printf("stem %d: %3d--%3d, %3d--%3d (Len: %d)\n",  i,
                       stemary[i].Pos[0],   stemary[i].Pos[1],
                       stemary[i].Pos[2],   stemary[i].Pos[3], stemary[i].stemlen );
        }
    #endif
    return 1;
}


int DataLoader::getSeqBegPos( char *tgtseq, int seqno )
{
    int i=0, j=0;
    int numbase=0;
    int lengtgt=0, idxpos=-1;
    char thebase=' ';
    char *trainstr = NULL;
    bool found=false;
    
    for(i=0; i<seqlength; i++ )
    {
       if( trainingseq[ seqno ][ i ] != 0 )
            numbase++;
    }
    trainstr = new char[numbase+1];
    for(i=0; i<seqlength; i++ )
    {
        thebase = trainingseq[ seqno ][ i ];
        if(  thebase != 0 )
            trainstr[j++]=numtochar( thebase );
    }
    trainstr[numbase]='\0';
    
    lengtgt=strlen( tgtseq );
    for( i=0; i<lengtgt; i++ )
    {
        if( tgtseq[i]==trainstr[0] )
        {
            for(j=0; j<numbase; j++ )
            {
                if( tgtseq[i+j]!=trainstr[j])
                    break;
                else
                {
                    idxpos=i;
                    found=true;
                }
            }
        }
        if( found )
            break;
    }
    delete [ ] trainstr;
    return idxpos;
}

void DataLoader::getStemPosInSeq( int seqno,  int stemno,  int pos[4])
{
    int i=0;
    int numbase=0;
    int sp[4];
    for(i=0; i<4; i++ )
    {
        sp[i]=stemary[stemno].Pos[i];
        pos[i]=-100;
    }
    for(i=0; i<seqlength; i++ )
    {
        if( pos[0]==-100 && sp[0]==i )
        {
            pos[0]=numbase;
        }
        else if( pos[1]==-100 && sp[1]==i )
        {
            if( trainingseq[ seqno ][ i ] != 0 )
                pos[1]=numbase;
            else
                pos[1]=numbase-1;
        }
        else if( pos[2]==-100 && sp[2]==i )
        {
            pos[2]=numbase;
        }
        else if( pos[3]==-100 && sp[3]==i )
        {
            if( trainingseq[ seqno ][ i ] != 0 )
                pos[3]=numbase;
            else
                pos[3]=numbase-1;
        }
        
        if( trainingseq[ seqno ][ i ] != 0 )
            numbase++;
    }
}

//no nested stem in another stem's arm
void DataLoader::analyzePasta(  )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::analyzePasta\n");
    #endif
    int i=0;
    
    //for stems
    numstems = getStemNum(  );
    //printf("number: %d\n", numstems);
    stemary = new StemLocation[ numstems ];
    //printf("number: %d\n", numstems);
    for(i=0; i<numstems; i++)
    {
        stemary[i].strid=stemchid[i];
        stemary[i].charid=stemchid[i].at(0);
        //printf("%d--%c\n", i, stemary[i].charid);
    }
    locateStemPosition( );
    
    //for loops
    numloops = 2*numstems+1;
    loopary = new LoopLocation[ numloops ];
    for(i=0; i<numloops; i++)
    {
        loopary[i].id = i;
    }
    locateLoopPosition( );
}

//subset==0 is the involving seq
void DataLoader::calcEandSD( int begpos, int endpos, double & eval, double & sdval, int *subset, int thegrp )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::calcEandSD\n");
    #endif
    int j=0, k=0;
    int numofseq=0;
    double variance=0.0, sumlen=0.0;
    int *len = new int[ numtrainseqs ];
    bool *setind = new bool[ numtrainseqs ];
    
    //check the begpos and endpos
    if( begpos<0  ) {
        printf("Error in DataLoader::calcEandSD:  begpos<0 \n");
        exit(0);
    }
    
    if( subset==NULL )
        thegrp = 0;
    //set the indicator of subset for calculating the E and SD
    for( j=0; j<numtrainseqs; j++ )
    {
        len[j]=0;
        if( subset==NULL )
            setind[j] = true;
        else
            setind[j] = (subset[j]==thegrp) ? true : false ;
    }
    //initialization
    eval=0.0;
    sdval=0.0;
    for( j=0; j<numtrainseqs; j++ )
    {
        if( setind[j] ) //is in the set
        {
            numofseq++;
            for( k=begpos; k<=endpos; k++ ) {
                if( trainingseq[j][k] != 0 )
                    len[j]++;
            }
            sumlen += len[j];
        }
    }
    //calculate E
    if( numofseq<=0 ) {
        eval=0.0;
        sdval=0.0;
    }
	else if( numofseq==1 ) {
        eval=sumlen;
        sdval=0.0;
	}
    else {
        eval=sumlen/numofseq;
        //calculate SD
        for( j=0; j<numtrainseqs; j++ ) {
            if( setind[j] ) //is in the set
                variance+=( len[j]-eval )*( len[j]-eval );
        }
        sdval=sqrt( variance/(numofseq-1) );
    }
    delete [ ] len;
    delete [ ] setind;
}

//get the statistical length of scanning window
//sumAll: default '+': sum of all avgs and lens;  '-': overall avg and len in training seqs
int DataLoader::getScanWinLength( double coeff, char sumAll )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::getScanWinLength\n");
    #endif
    int i=0;
    double sumE=0.0, sumSD=0.0;
    int overlen=0;
    
    if( sumAll=='-' )
    {
        calcEandSD( 0, seqlength-1, sumE, sumSD );
        //printf("Single calc: E=%f,  SD=%f\n", sumE, sumSD );
    }
    else
    {    
        for(i=0; i<numstems; i++ )
        {
            sumE += stemary[i].armE[0] + stemary[i].armE[1];
            sumSD += stemary[i].armSD[0] + stemary[i].armSD[1];
        }
        for(i=0; i<numloops; i++ )
        {
            sumE += loopary[i].lenE;
            sumSD += loopary[i].lenSD;
        }
        //printf("SumAll calc: E=%f,  SD=%f\n", sumE, sumSD );
    }
    overlen = (int) ceil( sumE + coeff*sumSD );
    //printf("Window size = %d\n", overlen );
    return overlen;
}

//calculate minimal length of scanning window
int DataLoader::getMiniWinLength( double coeff )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::getMiniWinLength\n");
    #endif
    int i=0, j=0, curlen=0;
    int minlen=seqlength;
    //get the minimum length: minimum length of sequence
    for(i=0; i<numtrainseqs; i++ ) {
        curlen=0;
        for(j=0; j<seqlength; j++ ) {
            if( trainingseq[i][j]!=0 )
                curlen++;            
        }
        if( curlen < minlen )
            minlen = curlen;
    }
    
    //get the minimum length: AVG - 3*SD
    //double sumE=0.0, sumSD=0.0;
    //calcEandSD( 0, seqlength-1, sumE, sumSD );
    //minlen = ceil( sumE - coeff*sumSD );
    return minlen;
}

//begpos: the beginning position of conserved region, 0-based index
//endpos: the end position of conserved region, 0-based index
//hitLen: the length of hit
//leftExt: left extension length
//leftExt: right extension length
//return: overall length of extended sequence
int DataLoader::getBothExtLength( int & leftExt, int & rightExt )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::getBothExtLength\n");
    #endif
    leftExt = extensionLength[0];
    rightExt = extensionLength[1];
    return leftExt+rightExt;
}

//begpos: the beginning position of conserved region, 0-based index
//endpos: the end position of conserved region, 0-based index
int DataLoader::getSearchRangeBeforeHit( int begHitPos, int endHitPos, double coeff, int & posBeg, int & posEnd )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::getSearchRangeBeforeHit\n");
    #endif
    double leftAvg=0.0, leftSD=0.0;
    
    //check the position is valid
    if( endHitPos < 0  || begHitPos >= seqlength || begHitPos > endHitPos ) {
        printf("Please check the conserved region positions for extension length calculation\n");
        exit(0);
    }

    begHitPos = ( begHitPos < 0 ) ? 0 : begHitPos;
    //left part before hit
    calcEandSD( 0, begHitPos-1, leftAvg, leftSD );
    posBeg = ( int ) ceil( leftAvg + coeff*leftSD );
    posEnd = ( int ) floor( leftAvg - coeff*leftSD );
    
    return 1;    
}

int DataLoader::getPastaLines(  )
{
    return numPastaLines;
}

void DataLoader::calcLoopStemStat(  )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::calcLoopStemStat\n");
    #endif
    int i=0;
    //************  stem  ***************
    for( i=0; i<numstems; i++ )
    {
        //left arm
        calcEandSD( stemary[i].Pos[0],  stemary[i].Pos[1],  stemary[i].armE[0],  stemary[i].armSD[0] );
        //right arm
        calcEandSD( stemary[i].Pos[2],  stemary[i].Pos[3],  stemary[i].armE[1],  stemary[i].armSD[1] );
        //mid region
        calcEandSD( stemary[i].Pos[1]+1, stemary[i].Pos[2]-1, stemary[i].midE, stemary[i].midSD );
    }
    //************  loop  ***************
    for( i=0; i<numloops; i++ )
    {
        //length for current loop
        if( loopary[i].looplen<=0 )
        {
            loopary[i].lenE=0.0;
            loopary[i].lenSD=0.0;
        }
        else
            calcEandSD( loopary[i].Pos[0], loopary[i].Pos[1], loopary[i].lenE, loopary[i].lenSD );
    }
} 

//**************** functions to group training seqs *********************
void DataLoader::groupTrainSeq(  )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::groupTrainSeq\n");
    #endif
    int i=0, j=0;
    int numleftgrps=0, numrightgrps=0, numgrps=0;
    int * leftMark = new int [numtrainseqs];
    int * rightMark = new int [numtrainseqs];
    //************  stem  ***************
    for( i=0; i<numstems; i++ )
    {
        numleftgrps = 1;
        numrightgrps = 1;
        for(j=0; j<numtrainseqs; j++ )
        {
            leftMark[j] = 1;
            rightMark[j] = 1;
        }
        //printf("Stem %d( %c ):\n", i, stemary[i].charid );
        //for( j=0; j<numtrainseqs; j++ )  printf(" %2d", j+1 );
        //printf("\n  Left Offset\n");
        numleftgrps = SplitTrainingSeq( stemary[i].Pos[0], leftMark, 1, numleftgrps );  //left offset
        //printf("  Right Offset\n");
        numrightgrps = SplitTrainingSeq( stemary[i].Pos[2], rightMark, 1, numrightgrps );  //right offset
        
        /*printf("  FINAL: leftGroup:%d  rightGroup:%d\n", numleftgrps, numrightgrps );
        for( j=0; j<numtrainseqs; j++ )  printf(" %2d", leftMark[j] );
        printf("\n");
        for( j=0; j<numtrainseqs; j++ )  printf(" %2d", rightMark[j] );
        printf("\n\n");
        */
        numgrps = calcGrpOffset( i, leftMark, rightMark );
    }
    //printESDEveryStem( );
	delete [] leftMark;
	delete [] rightMark;
}

void DataLoader::groupSeqsOnFilter( )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::groupSeqsOnFilter( )\n");
    #endif
    for( int i=0; i<numstems; i++ ) {
        calcSpecialEandSD( basePos, stemary[i].Pos[0],  stemary[i].woffsetE[0],  stemary[i].woffsetSD[0] );
        calcSpecialEandSD( basePos, stemary[i].Pos[2],  stemary[i].woffsetE[1],  stemary[i].woffsetSD[1] );
    }
    groupTrainSeq(  );
    //makeAllEpositive(  );
}

void DataLoader::calcSpecialEandSD( int begpos, int endpos,  double & eval, double & sdval, int *subset, int thegrp )
{
    if( begpos==endpos )  {
        eval=0.0;
        sdval=0.0;
    }
    else  {
        if( begpos > endpos ) {
            begpos--;
            if( begpos < 0 ) begpos=0;
            calcEandSD( endpos, begpos,  eval,  sdval, subset, thegrp );
            eval = -eval;
        }
        else {
            begpos++;
            calcEandSD( begpos, endpos,  eval,  sdval, subset, thegrp );
        }
    }
}

int DataLoader::makeAllEpositive( )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::makeAllEpositive\n");
    #endif
    //calculate the offset between filterBeg and basePos
    int FbegBase = basePos - filterPos[0];  //filter Beg to Base point
    cout<<endl<<endl<<endl;
    printStemLoop('S');
    //the amount of nts to adjust
    int baseOffset = extensionLength[0]+ FbegBase+1;
    for( int i=0; i<numstems; i++ )
        stemary[i].offsetAllE( (double)baseOffset );

    cout<<"  basePos: "<<basePos<<"    Offset: "<<baseOffset<<endl;
    cout<<" "<<extensionLength[0]<<"  [ "<<filterPos[0]<<" - "<<filterPos[1]<<" ]  "<<extensionLength[1]<<endl;
    
    return 1;
}

int DataLoader::calcExtensionLength(  )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::calcExtensionLength\n");
    #endif
    int i=0, maxLeft=0, maxRight=0, curLeft=0, curRight=0;
    double leftE=0.0, leftSD=0.0, rightE=0.0, rightSD=0.0;
    double coeff=3.0;
    
    calcEandSD( 0,      filterPos[0]-1,   leftE,    leftSD );
    calcEandSD( filterPos[1]+1, seqlength-1,   rightE,   rightSD );
    extensionLength[0] = (int)ceil(leftE + coeff*leftSD);
    extensionLength[1] = (int)ceil(rightE + coeff*rightSD);
    
    for( i=0; i<numstems; i++ )  {
        curLeft = stemary[i].obtainMostEnds( 'L', coeff );
        curRight = stemary[i].obtainMostEnds( 'R', coeff );
        if( maxLeft < curLeft )
            maxLeft = curLeft;
        if( maxRight < curRight )
            maxRight = curRight;
    }
    maxRight -= filterPos[1]-filterPos[0];
    
    if( extensionLength[0]<maxLeft )
        extensionLength[0] = maxLeft;
    if( extensionLength[1]<maxRight )
        extensionLength[1] = maxRight;
    //cout<<"maxLeft: "<<maxLeft<<"   maxRight: "<<maxRight<<" || wleft: "<<extensionLength[0]<<"   wright: "<<extensionLength[1]<<endl;
    return 1;
}

void DataLoader::setFilterEndPosition( int filterBeg, int filterEnd, int clusteringBasePos )
{
    #ifdef DEBUG_DEV
        printf("DataLoader::setFilterEndPosition\n");
    #endif
    //check the position is valid
    if( filterBeg > filterEnd || filterEnd < 0 ) { //no specified HMM filter
        filterPos[0] = 0;
        filterPos[1] = 0;
        basePos = 0;
    }
    else
    {
        filterPos[0] = ( filterBeg < 0 ) ? 0 : filterBeg;
        filterPos[1] = filterEnd;
        if( clusteringBasePos == -1 )
             basePos=filterPos[0];
        else {
            if( clusteringBasePos<=filterBeg || clusteringBasePos>=filterEnd )
            {
                printf("Make sure the BasePos is within fitlerBeg and filterEnd\n");
                exit(0);
            }
            basePos=clusteringBasePos;
        }
    }
    //groupSeqsOnFilter( clusteringBasePos );
}

bool DataLoader::needSplit( double curE, double curSD )
{
    bool retval = false;
    if( curSD > splitThreshold )
        retval = true;
    return retval;
}

int DataLoader::getSubseqLen( int seqno, int begpos, int endpos )
{
    int i=0, subsize=0;
	int leftEnd=0, rightEnd=0;
	if( begpos < endpos ) {
		leftEnd = begpos;
	    rightEnd = endpos;
	} else {
		leftEnd = endpos;
	    rightEnd = begpos;
	}
    for(i=leftEnd; i<=rightEnd; i++ )
    {
        if( trainingseq[seqno][i] ==1 || trainingseq[seqno][i] ==2 ||
            trainingseq[seqno][i] ==3 || trainingseq[seqno][i] ==4  )
            subsize++;
    }
    return subsize;
}

//thegrp: groupID,   numgrps:# of groups
int DataLoader::SplitTrainingSeq( int endpos, int *mark, int thegrp, int numgrps  )
{
    int i=0;
    double tempE=0.0, tempSD=0.0;
    int curid1=0, curid2=0;
	int leftEnd=0, rightEnd=0;
    if( basePos<endpos ) {
		leftEnd=basePos;
		rightEnd=endpos;
	} else {
		leftEnd=endpos;
		rightEnd=basePos;
	}
    calcEandSD( leftEnd,  rightEnd,  tempE,  tempSD,  mark, thegrp );
    
    /*printf("      Group: %d, E: %.3f, SD: %.3f \n", thegrp, tempE, tempSD );
    for( i=0; i<numtrainseqs; i++ ) {
        if( mark[i] == thegrp )
            printf(" %2d", mark[i] );
        else
            printf("   " );
    }
    */
    if( needSplit( tempE, tempSD ) ) {
        //printf("    need split\n");
        numgrps++;
        curid1 = thegrp;
        curid2 = numgrps;
        for( i=0; i<numtrainseqs; i++ ) {
            if( mark[i]==thegrp ) {
                if( getSubseqLen(i, basePos, endpos) > tempE )
                    mark[i] = curid2;
            }
        }
        numgrps = SplitTrainingSeq( endpos, mark, curid1, numgrps );
        numgrps = SplitTrainingSeq( endpos, mark, curid2, numgrps );
    }
    //printf("\n");
    return numgrps;
}

int DataLoader::calcGrpOffset( int stemno, int * leftMark, int* rightMark  )
{
    int i=0, j=0;
    bool hassame=false;
    int cursum=0;
    int *markgp = new int [ numtrainseqs ];
    int *leftRec = new int [ numtrainseqs ];
    int *rightRec = new int [ numtrainseqs ];
    for(i=0; i<numtrainseqs; i++ ) {
        markgp[i]=-1;
        leftRec[i]=-1;
        rightRec[i]=-1;
    }
    //get different groups
    for(i=0; i<numtrainseqs; i++ ) 
    {
        hassame=false;
        for( j=0; j<cursum; j++ ) {
            if( (leftMark[i] == leftRec[j])&& (rightMark[i] == rightRec[j]))
            {
                hassame=true;
                break;
            }
        }
        if( hassame ) {
            markgp[i]=j;
        }
        else {
            leftRec[ cursum ] = leftMark[i];
            rightRec[ cursum ] = rightMark[i];
            markgp[i]=cursum;
            cursum++;
        }
    }
    /*for( j=0; j<numtrainseqs; j++ )  printf(" %2d", markgp[j] );
    printf("\n\n");*/
    stemary[stemno].numOffset = cursum;
    for( i=0; i<2; i++) {
        stemary[stemno].offsetE[i] = new double[cursum];
        stemary[stemno].offsetSD[i] = new double[cursum];
    }
    
    for( i=0; i<cursum; i++) {
        calcSpecialEandSD( basePos, stemary[stemno].Pos[0],  stemary[stemno].offsetE[0][i],  stemary[stemno].offsetSD[0][i],  markgp, i );
        calcSpecialEandSD( basePos, stemary[stemno].Pos[2],  stemary[stemno].offsetE[1][i],  stemary[stemno].offsetSD[1][i],  markgp, i );
    }

	//detect memory leak problem in rnatops zbhuang 20090815
	if(markgp != NULL) {
		delete [] markgp;
		markgp = NULL;
	}
	if(leftRec != NULL) {
		delete [] leftRec;
		leftRec = NULL;
	}
	if(rightRec != NULL) {
		delete [] rightRec;
		rightRec = NULL;
	}
    return cursum;
}

void DataLoader::printESDEveryStem(  )
{
    int i=0, j=0;
    for( i=0; i<numstems; i++ )
    {
        printf("\nStem%2d(%c), # of groups:%2d\n", i, stemary[i].charid, stemary[i].numOffset );
        for( j=0; j<stemary[i].numOffset; j++ )
        {
            printf("Group %2d:  Left offset: E=%6.2f  SD=%6.2f   Right offset: E=%6.2f  SD=%6.2f \n",
                    j+1 , stemary[i].offsetE[0][j], stemary[i].offsetSD[0][j], 
                   stemary[i].offsetE[1][j], stemary[i].offsetSD[1][j] );
        }
    }
}

void DataLoader::sumEandSD( int beg, int end, double &sumE, double &sumSD )
{
    int i=0;
    
    sumE=0.0;
    sumSD=0.0;
    
    for( i=0; i<numloops; i++ )
    {
        if( loopary[i].Pos[0]>beg && loopary[i].Pos[1]<end )
        {
            sumE += loopary[i].lenE;
            sumSD += loopary[i].lenSD;
        }
    }
    for( i=0; i<numstems; i++ )
    {
        if( stemary[i].Pos[0]>beg && stemary[i].Pos[1]<end )  //check left arm
        {
            sumE += stemary[i].armE[0];
            sumSD += stemary[i].armSD[0];
        }
        if( stemary[i].Pos[2]>beg && stemary[i].Pos[3]<end ) //check right arm
        {
            sumE += stemary[i].armE[1];
            sumSD += stemary[i].armSD[1];
        }
    }
}


void DataLoader::printBaseFreq(  )
{
    int i=0;
    printf("The frequency of residues:\n");
    for(i=1; i<5; i++ )
    {
        printf(" %c: %f  log:%f\n", numtochar(i), BaseHz[i-1], log(BaseHz[i-1]));
    } 
}

void DataLoader::printPasta( )
{
    int i=0; 
    for( i=0; i<numPastaLines; i++ )
        printf("%s\n", pasta[i] );
}

int DataLoader::getCertainTrainSeq( int idx, char *&seqbuf )
{
    int i=0, j=0, len=0;
    if( idx >= numtrainseqs )
    {
        printf("TrainingSequence index is beyond the max. number: %d-1\n", numtrainseqs );
        exit(0);
    }
    else
    {
        for( i=0; i<seqlength; i++ )
        {
            if( trainingseq[ idx ][i]!= 0 )
                len++;
        }
        seqbuf = new char[ len+1 ];
        seqbuf[len] = '\0';
        for( i=0, j=0; i<seqlength; i++ )
        {
            if( trainingseq[ idx ][i]!=0 )
                seqbuf[j++] = numtochar( trainingseq[idx][i] );
        }
        return len;
    }
}

void DataLoader::printTrainingSeq( int idx )
{
    int i=0, j=0;
    if( idx < 0 )
    {
        for( i=0; i<numtrainseqs; i++ )
        {
            for( j=0; j<seqlength; j++ )
            {
                printf("%c", numtochar( trainingseq[i][j] ) );
            }
            printf("\n" );
        }
    }
    else if( idx >= numtrainseqs )
    {
        printf("TrainingSequence index is beyond the max. number: %d-1\n", numtrainseqs );
        return;

    }
    else
    {
        for( j=0; j<seqlength; j++ )
        {
            printf("%c", numtochar( trainingseq[ idx ][j] ) );
        }
        printf("\n" );
    }
    printf("\n" );
}

//S: stem L:loop B: both
void DataLoader::printStemLoop( char which )
{
    int i=0, j=0;
    //print Stem
    printf("Base Position is: %d\n", basePos);
    if( which=='S' || which=='B' )
    {
        printf("Stem position:\n");
        printf("       Left       Right   Stem   LOffset      LArm       MidLoop      RArm        ROffset\n");
        printf("  No Beg   End  Beg   End  Len   E     SD    E    SD     E    SD     E     SD     E     SD\n");
        for( i=0; i<numstems; i++ )
        {
            for( j=0; j<stemary[i].numOffset; j++ )
            {
                printf("\n");
                if( j == 0 )
                {
                    printf("%s %2d %3d---%3d  %3d---%3d %3d  %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    stemary[i].strid.data( ), i+1,
                    stemary[i].Pos[0],  stemary[i].Pos[1],       stemary[i].Pos[2], stemary[i].Pos[3], 
                    stemary[i].stemlen, stemary[i].offsetE[0][j],   stemary[i].offsetSD[0][j],
                    stemary[i].armE[0], stemary[i].armSD[0],     stemary[i].midE, stemary[i].midSD,
                    stemary[i].armE[1], stemary[i].armSD[1],     stemary[i].offsetE[1][j], stemary[i].offsetSD[1][j]);
                }
                else
                {
                    printf("%c %2c %3c   %3c  %3c   %3c %3c  %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                         ' ', ' ', ' ',  ' ', ' ', ' ', ' ', 
                    stemary[i].offsetE[0][j], stemary[i].offsetSD[0][j],
                    stemary[i].armE[0], stemary[i].armSD[0],   stemary[i].midE, stemary[i].midSD,
                    stemary[i].armE[1], stemary[i].armSD[1],   stemary[i].offsetE[1][j], stemary[i].offsetSD[1][j]);
                }
            }
        }
        printf("\n");
    }
    //print Loop
    if( which=='L' || which=='B' )
    {
        printf("Loop position:\n");
        printf("           Loop          LOffset           Length          ROffset\n");
        printf("#No    Beg     End       E     SD         E     SD         E     SD\n");
        for( i=0; i<numloops; i++ )
        {
            printf("%3d    %3d----%3d    %6.2f  %6.2f   %6.2f  %6.2f   %6.2f  %6.2f\n",
                    i+1,   loopary[i].Pos[0],  loopary[i].Pos[1],
                    loopary[i].offsetE[0], loopary[i].offsetSD[0],  loopary[i].lenE,  loopary[i].lenSD,
                    loopary[i].offsetE[1], loopary[i].offsetSD[1] );
        }
        printf("\n");
    }
}

void DataLoader::printStringAll( )
{
    int i=0;
    printf("PastaLine ( %d lines):\n", numPastaLines );
    for( i=0; i<numPastaLines; i++ )
        cout<<orgpasta[i]<<endl;
    printf("Comment ( %d lines ):\n", numtrainseqs);
    for( i=0; i<(int)orgcommns.size( ); i++ )
        cout<<orgcommns[i]<<endl;
    printf("Sequence Lines:\n");
    for( i=0; i<(int)orgtrains.size( ); i++ )
        cout<<orgtrains[i]<<endl;
}
  
void DataLoader::printAll( )
{
    printf("Training sequences:\n" );
    printPasta( );
    printTrainingSeq( );
    printStemLoop( );
}


//added by zbhuang 20070701
int ** DataLoader::getPTrainingSeqs()
{
  return this->trainingseq;
}

char * DataLoader::getPPasta()
{
  return this->pasta[0];
}

int  DataLoader::getNumstem()
{
  return this->numstems;
}

int  DataLoader::getNumloop()
{
  return this->numloops;
}

StemLocation * DataLoader::getStemarray()
{
  return this->stemary;
}

LoopLocation * DataLoader::getLooparray()
{
  return this->loopary;
}

