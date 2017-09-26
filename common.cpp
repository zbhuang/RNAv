#include "common.h"

//-------------------------------
int isRNABase(char nucleotide)
{
    int result;
    nucleotide = toupper( nucleotide);
    switch( nucleotide )
    {
        //case GAP:
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'U':
               result=1;
               break;
        default:
            	 //printf("There is invalid residue: %c\n", nucleotide);
            	 result=0;
    }
    return result;
}

int isRNABase(int nucleotide)
{
    int result;
    nucleotide = toupper( nucleotide);
    switch( nucleotide )
    {
        //case GAP:
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
               result=1;
               break;
        default:
            	 //printf("There is invalid residue: %c\n", nucleotide);
            	 result=0;
    }
    return result;
}

int chartonum(char nucleotide)
{
    int result;
    nucleotide = toupper( nucleotide);
    switch( nucleotide )
    {
        case GAP:  result=0;
                   break;
        case 'A':  result=1;
                   break;
        case 'C':  result=2;
                   break;
        case 'G':  result=3;
                   break;
        case 'T':
        case 'U':  result=4;
                   break;
        default:   
            	 printf("There is invalid residue: %c, which is treated as gap.\n", nucleotide);
            	 result=0;
    }
    return result;
}

char numtochar(int base)
{
    char result;
    switch( base )
    {
        case 0:  result=GAP;
                   break;
        case 1:  result='A';
                   break;
        case 2:  result='C';
                   break;
        case 3:  result='G';
                   break;
        case 4:  result='U';
                   break;
        default: printf("There is invalid residue: %d, which is treated as gap\n", base);
                 result=' ';
    }
    return result;
}

//-------------------------------
bool isCanonicalPair(char basex, char basey)
{
    char bex=' ';
    char bey=' ';
    bool result=false;
    
    bex = toupper( basex);
    bey = toupper( basey); 
    
    if( bex=='U' && bey=='A' )
        result=true;
    else if( bex=='A' && bey=='U' )
        result=true;
    else if( bex=='G' && bey=='C' )
        result=true;
    else if( bex=='C' && bey=='G' )
        result=true;
    else if( bex=='U' && bey=='G' )
        result=true;
    else if( bex=='G' && bey=='U' )
        result=true;
    else
        result=false;
    return result;
}

int isFileName( char *str)
{
    int i=0, len=0;
    int flag=0;
    char base=' ';
    
    len=strlen( str );
    for( i=0; i<len; i++ )
    {
        base=toupper( str[i] );
        if( base=='\r' || base=='\n')
            break;
        if( base!='A' && base!='C' && base!='G' && base!='U' && base!='T' && base!='-' )
        {
            flag=1;
            break;
        }
    }
    return flag;
}

void printPairEmP( double pairprob[5][5] )
{
    int i=0, j=0;
    char base[5]={'-', 'A','C', 'G', 'U'};
    printf("  ");
    for(i=0; i<5; i++ )
    {
        printf("%9c", base[i] );
    }
    printf("\n");
    for(i=0; i<5; i++ )
    {
        printf("%5c", base[i]);
         for(j=0; j<5; j++)
        {
            printf("%9.5f", pairprob[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printPairEmP( int pairprob[5][5] )
{
    int i=0, j=0;
    char base[5]={'-', 'A','C', 'G', 'U'};
    printf("     ");
    for(i=0; i<5; i++ )
    {
        printf("%9c", base[i] );
    }
    printf("\n");
    for(i=0; i<5; i++ )
    {
        printf("%5c", base[i]);
         for(j=0; j<5; j++)
        {
            printf("%9d", pairprob[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printPairHz( int pairprob[5][5], double prob[5][5])
{
    int i=0, j=0;
    char base[5]={'-', 'A','C', 'G', 'U'};
    printf("# of Pairs:                   ");
    printf("        Percentage: \n");
    printf("     ");
    for(i=0; i<5; i++ )
        printf("%5c", base[i] );
    printf("   ");
    for(i=0; i<5; i++ )
        printf("%8c", base[i] );
    printf("\n");
    
    for(i=0; i<5; i++ )
    {
        printf("%5c", base[i]);
        for(j=0; j<5; j++)
            printf("%5d", pairprob[i][j]);
        printf("      ");
        for(j=0; j<5; j++)
            printf("%8.4f", prob[i][j]);
        printf("\n");
    }
    printf("\n");    
}

void printBaseEmP( double baseprob[4] )
{
    int i=0;
    char base[4]={'A','C', 'G', 'U'};
    printf("  ");
    for(i=0; i<4; i++ )
    {
        printf("%9c", base[i] );
    }
    printf("\n");
    printf("%5c", ' ');
    for(i=0; i<4; i++ )
    {
        printf("%9.5f", baseprob[i]);
    }
    printf("\n");
    printf("\n");
}

void printSingleEmProb( double baseprob[4] )
{
    int i=0;
    printf("                    <");
    for(i=0; i<4; i++ )
        printf("%9.5f", baseprob[i]);
    printf(" >\n");
}

//--------------------------------------------------

int obtainSequence(  char * filename, char* &bufseq )
{
    #ifdef DEBUG_DEV
        printf("obtainSequence\n");
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
    bufseq= new char[templen+1];
    bufseq[templen] = '\0';
    //get the str
    fseek( fp, 0, SEEK_SET);
    while ( !feof(fp) )
    {
        ch=fgetc (fp);
        if( isRNABase(ch) )
            bufseq[i++]=ch;
    }
    if( i!=templen )
    {
        printf("There's something wrong with target file.\n");
        return 0;
    }
    fclose(fp);
    return templen;
}

//(1-h)*orgm + h*priorm
void normalizeMatrix( double orgm[5][5], double priorm[5][5], double h)
{
    int i=0, j=0;
    double sumvalue=0.0;
    double ht=0.0, hp=0.0;
    if( h>1 )
        hp=1;
    else if( h<0)
        hp=0;
    else
        hp=h;
    
    ht=1-hp;    
    //print the original matrix
    //printf("Orignal matrix:\n");
    //printPairEmP( orgm );

    //printf("Prior matrix:\n");
    //printPairEmP( priorm );
    
    for(i=0; i<5; i++ )
    {
        for(j=0; j<5; j++ )
        {
            orgm[i][j]=ht*orgm[i][j] + hp*priorm[i][j];
            sumvalue += orgm[i][j];
        }
    }
    
    for(i=0; i<5; i++ )
    {
        for(j=0; j<5; j++ )
        {
            orgm[i][j]=orgm[i][j]/sumvalue;
        }
    }
    //print the new matrix
    //printf("New matrix:\n");
    //printPairEmP( orgm );
    //printf("\n\n");
    
}

int isEven(int iVal)
{
  if (iVal & 1)
    return 0;	//odd
  else
    return 1;	//even
}
