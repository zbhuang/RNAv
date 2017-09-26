#ifndef _SEARCHOUTPUT_H_
#define _SEARCHOUTPUT_H_

//#include "constdef.h"
//#include "timecalc.h"

#include "GenomeData.h"
#include "ConstantsInSearch.h"

#include <iomanip>

#define GENOME_TAG_HIT_EXTENSION		"Extension of the hit"
#define GENOME_TAG_HIT_EXTENSION_TAIL	"--------"
#define GENOME_TAG_TOTAL_TIME_USED		"Total time used"
#define GENOME_TAG_TOTAL_NUM_HIT		"Total no of hits"

//#include <string.h>	//gcc3.4
#include <string>
using std::string;

char * strrpl(char *s, const char s1, const char s2);
//char * strdel(char *s, const char s1);
char * strappend(char *s, const char s1);
char * strextract(char *s, int start_pos, int end_pos);

//zbhuang 20100201
int countNumChar(char *str, char ch);

string extractfilename(string s);
string printYN(int iFlag);
string printDirection(int iFlag);
void printGenomeSegments(char * c_seq, int num);
void printGenomeSegments(char * c_seq, int num, ofstream &outfile);

//
int isHmmFilterFile(char * filename);
/*
void printValuesForParameters(double	pseudocount,
							  int		top_k,
							  double	threshold,
							  int		num_nts_overlap_between_stem,
							  int		iFlagMergeCandInPreprocess,
							  int		iFlagCandwithShortestLength,
							  int		iShiftNumMergeCand,
							  int		iAllowedNullLoopInsNum,
							  double	pcoeff,
							  int		iJumpStrategy,
							  int		iStepSize,
							  double	dScoreThresholdInJump,
							  int		iFlagPrintStrAlignInfo,
							  int		iFlagSearchReverse,
							  int		iFlagPrintDebugInfo,
							  int		iFlagPrintDebugInputInfo,
							  int		iFlagPrintDebugTreeInfo,
							  int		iFlagPrintDebugDPWinInfo,
							  int		iFlagPrintDebugTDInfo,
							  int		iFlagPrintDebugDPSearchInfo);
*/
void printValuesForParameters(double	pseudocount,
							  int		top_k,
							  double	threshold,
							  int		num_nts_overlap_between_stem,
							  int		iFlagMergeCandInPreprocess,
							  int		iFlagCandwithShortestLength,
							  int		iShiftNumMergeCand,
							  int		iAllowedNullLoopInsNum,
							  double	pcoeff,
							  int		iJumpStrategy,
							  int		iStepSize,
							  int		iFlagSearchReverse,
							  int		max_empty_stem_num);

void printCurrentTime();
void printOneStarLine();
void printAdditionalInfo();
void printFilterSearchResultHeader(string	trainfile,
								   int		profilelength,
								   int		filterauto, 
								   int		filtertype, 
								   int		filterbegin, int filterend, 
								   string	genomefile,	//char * genomefile,
								   int		genomenum,
								   long		total_genome_length);
void printFilterHitIndex(int index);
void printHitPos(int start, int end);
void printHitPos(int start, int end, ofstream &outfile);
void printHMMHitScore(double score);
void printHMMHitScore(double score, ofstream &outfile);
void printWholeStructureHitScore(double score);
void printWholeStructureHitScore(double score, ofstream &outfile);
void printHitScore(double score);
void printHitScore(double score, ofstream &outfile);
void printFilterHitAlignment(char *seq, char *align);
void printFilterHitAlignment(char *uSeq, char *align, ofstream &outfile);
void printFilterHitExtenPos(int start, int end);
void printFilterHitExtenPos(int start, int end, ofstream &outfile);
void printFilterHitExtenNtsHeader();
void printFilterHitExtenNtsHeader(ofstream &outfile);
//void printFilterHitExtenNtsTail();
void printTotalHitNum(int hit_num);
void printTotalTimeInfo(double dTime);
void printSearchEndingTimeInfo();
void printWholeSearchResultHeader(string	trainfile,
								  int		profilelength,
								  string	genomefile,
								  int		genomenum,
								  long		total_genome_length,
								  double	pseudocount,
								  int		top_k,
								  double	threshold,
								  int		num_nts_overlap_between_stem,
								  int		iFlagMergeCandInPreprocess,
								  int		iFlagCandwithShortestLength,
								  int		iShiftNumMergeCand,
								  int		iAllowedNullLoopInsNum,
								  double	pcoeff,
								  int		iJumpStrategy,
								  int		iStepSize,
								  int		iFlagSearchReverse,
								  //int		streamline,
								  int		searchtype,
								  int		filterbegin,
								  int		filterend,
								  int		max_empty_stem_num
								  );
int	getProfileLength(char * trainfile);
void printWholeSearchHitIndex(int index);
void printWholeSearchHitIndex(int index, ofstream &outfile);
void printEValue(double evalue);
void formatStringOutput(string illustration, string line_one, string line_two, string line_three, int start, int end, int flag, ofstream &outfile);
void formatOneOutput(string line_one_unit, string line_two_unit, string line_thr_unit, int index, int start, int flag, ofstream &outfile);
void formatStringOutput(string illustration, string line_one, string line_two, string line_three, int start, int end, int flag);
void formatOneOutput(string line_one_unit, string line_two_unit, string line_thr_unit, int index, int start, int flag);

#endif
