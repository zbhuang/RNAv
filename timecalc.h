#ifndef _TIMECALC_H_
#define _TIMECALC_H_

#include "constdef.h"

//------------------------------------  Candid  -----------------------------------
class TimeCalc  //Candidates
{
private:
    clock_t  beg, end;
    int hasbeg, hasend;
    double   totaltime;
    
    int  hour;
    int  minu;
    double  secd;
    
    int calcDuration(  );
    int calcHMSFormat(  );
    
public:
    TimeCalc(  );
    ~TimeCalc(  );
    
    void resetTime(  );
    void countBegin(  );
    void countEnd(  );
    void recordTime(  );
    void getTime( int & hr, int & mn, double & se );
    void printTime( char * remark=NULL );
};

#endif
