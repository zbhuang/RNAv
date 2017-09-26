#ifndef _HEADDEF_H_
#define _HEADDEF_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

//#define DEBUG_DEV
#define BEGUID  0
#define INVLDPROB -10000.0
#define INVALID -10000
#define ZEROPROB 0.0
#define GAP '-'
#define PASTALOOP '.'
#define PSEUDOCNT 0.001
#define PSDCNTGAP 0.001
#define LOGODDS true
#define ENOFFSET true


#define COESTDEV 3.0  //E + cSD
//#define COESTDEV 12.0  //tried by zbhuang 20090615

#define CANDIDNUM 15  //just above threshold

//offset calculation
#define CAL_OFFSET 'L';

//candidates' merging parameter
#define CYKCANDFILTER true
#define STEMFTALLOW 0
#define PREFERSHORTLEN true
#define PERCENTSCOREDROP -1   // <0: don't consider the score

#define NUMCHECKCANDID 30    //for loop region shift penalty

//prior matrix
#define ONLYCANONICAL false

//length penalty
#define SDZEROALLOWBOUND 0
#define SDZEROINFINITY false

//type of production
#define UNKNOWN_NT 0
#define LEFT_INS   1
#define PAIR_ADD   2
#define RIGHT_INS  3
#define BIFUR_NT   4
#define KEEP_POS    5
#define PAIR_TMN   6
#define LSINGLE_TMN 7
#define RSINGLE_TMN 8
#define SINGLE_TMN 9
#define BEGIN_NT   10


//transition prob for model
#define SCFG_NUM_NEXTP  5
#define SCFG_LFH  0
#define SCFG_LFI  1
#define SCFG_MS   2
#define SCFG_RGI  3
#define SCFG_RGH  4

#define SCFG_LEFT  0
#define SCFG_RIGHT 1

//transition prob for model
#define HMM_NUM_NEXTP  3
#define HMM_D    0
#define HMM_I    1
#define HMM_M    2

//symbol for SCFG and HMM in alignment structure
#define AlgScfgLeftIns  'l'
#define AlgScfgRightIns  'r'
#define AlgScfgLeftArm  '['
#define AlgScfgRightArm  ']'

#define AlgHmmIns  'i'
#define AlgHmmDel  'd'
#define AlgHmmMatch  'm'

//for good stem detection
#define PERTAGEPAIRVALUE 0.1
#define LEASTNUMPAIR 1

typedef struct NontermProp
{
     int uid;
     char state;
     int id;
} NontermProp;

typedef struct GramRule     //GramProd
{
    int category;
    char direct;
    
    char ltype;
    int lid;  // UID
     
    char rtype;
    int rid;  // UID
    union
    {
        double pairp[5][5];
        double basep[4];
        double singlep;
    };
} GramRule;

#endif
