#include "timecalc.h"

//********************************************************************************
//*******************    implementation of Candidate class      ******************
//********************************************************************************
TimeCalc::TimeCalc(  )
{
    beg=0;
    end=0;
    hasbeg=0;
    hasend=0;
    totaltime=0.0;
    
    hour=0;
    minu=0;
    secd=0.0;
}

TimeCalc::~TimeCalc(  )
{
    beg=0;
    end=0;
    hasbeg=0;
    hasend=0;
    totaltime=0.0;
    
    hour=0;
    minu=0;
    secd=0.0;
}

void TimeCalc::resetTime(  )
{
    beg=0;
    end=0;
    hasbeg=0;
    hasend=0;
    totaltime=0.0;
    hour=0;
    minu=0;
    secd=0.0;
}

void TimeCalc::recordTime(  )
{
    if( hasbeg == 0 )
    {
        beg = clock();
        hasbeg = 1;
    }
    else
    {
        end = clock();
        hasend = 1;
        calcDuration(  );
        calcHMSFormat(  );
    }
}

void TimeCalc::countBegin(  )
{
    beg=clock();
    hasbeg = 1;
}

void TimeCalc::countEnd(  )
{
    end=clock();
    hasend = 1;
    calcDuration(  );
    calcHMSFormat(  );
}

int TimeCalc::calcDuration(  )
{
    if( hasbeg==0 )
    {
        printf("Error: don't begin the clock\n");
        return 0;
    }
    if( hasend==0 )
    {
        printf("Error: don't end the clock\n");
        return 0;
    }
    totaltime=(double)(end - beg )/CLOCKS_PER_SEC;
    return 1;
}

int TimeCalc::calcHMSFormat(  )
{
    int status=0;
    int intpart=0, tempm=0;
    
    if( totaltime==0.0 )
    {
        status=calcDuration(  );
        if( !status )
            return 0;
    }
    
    intpart=(int)totaltime;
    tempm = (int)floor( intpart/60.0 );
    
    secd = totaltime - tempm*60;    
    hour = (int)floor( tempm/60.0 );
    minu = tempm%60;

    return 1; 
}

void TimeCalc::getTime( int & hr, int & mn, double & se )
{
    hr = hour;
    mn = minu;
    se = secd;  
}

void TimeCalc::printTime( char * remark )
{
    printf("%s\n", remark );
    if( minu>0 )
         printf("  Total time:  %f s ( %d:%d:%4.2f ) \n", totaltime, hour, minu, secd );
     else
         printf("  Total time:  %f s  \n", totaltime );
}
