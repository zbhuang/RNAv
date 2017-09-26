#ifndef _CONSTANTSINSEARCH_H_
#define _CONSTANTSINSEARCH_H_

//------ rnatops.cpp ------

//searchtype
#define	SEARCH_HMMFILTER			0
#define	SEARCH_SUBSTRUCTUREFILTER		1
#define	SEARCH_WHOLESTRUCTURE			2
#define	SEARCH_STREAMLINE			3

#define SEARCH_FILTERHIT_BASED_WHOLESTRUCTURE	4

//do local structure alignment zbhuang 20090723
//do full combination zbhuang 20090914
#define	SEARCH_LOCALSTRUCTUREALIGNMENT		5	//4

#define K_DEFAULT_VALUE 10

//------ StructureSearch.h ------
#define RUNNINGMODE_WEBSERVER	0
#define RUNNINGMODE_STANDALONE	1

//#define	SEARCH_HMMFILTER		0
//#define	SEARCH_SUBSTRUCTUREFILTER	1
//#define	SEARCH_WHOLESTRUCTURE		2
//#define	SEARCH_STREAMLINE		3

//------ DynamicP.h ------
//#define HMMFILTER_STRUCTURE_SEARCH	3

//adding the empty-stem model. zbhuang 20090211
#define NEGATIVE_ONE	-1
#define ZERO			0
#define ONE				1
#define	EMPTY_STEM_LABEL	"^"

#define MAXNUMPT		2
#define MAXNUM_CHILD	6

//limiting the number of empty-stem. zbhuang 20090213
#define MAXNUMEMPTYSTEM		10	//2

//#define MAXTREE		2

#define NONCONNECTED	0
#define NONDIRECTED		1
#define DIRECTED		2
#define WEIGHTZERO		3
#define BOTH			4

#define SMALLEST -1.0E37//-1.0E10//

#define MYTHRESHOLD		-10000.0
#define INITSCORE		0.0
#define ACTUAL_CANDIDATE_NUM 	1000000

#define MAX_STRUCTURE_ALIGN_LENGTH	2000

//#define STEPSIZE_IN_CANDIDATE		2
//#define CANDIDATE_SCORE_THRESHOLD	15
#define CANDIDATE_POSITION_THRESHOLD 3

//20081019 initialization
#define NEGATIVE_ONE	-1

#define	TWO_PASTALINE	2
#define WIDTH 10

#define TRUE  1
#define FALSE 0

#define ZERO 0

#define	SCORE_STEM			0
#define	SCORE_STEM_POSLOOP	1
#define	SCORE_STEM_LOOP		2

//add extend_dir zbhuang
#define EXTLEFT		'L'
#define EXTRIGHT	'R'
#define EXTNA		'N'
#define EXTNO		' '

#endif
