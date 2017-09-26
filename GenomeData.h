// GenomeData.h: interface for the GenomeData class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GENOMEDATA_H__E09901D8_4A6A_4AF3_ACBC_AEF08B96E963__INCLUDED_)
#define AFX_GENOMEDATA_H__E09901D8_4A6A_4AF3_ACBC_AEF08B96E963__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "common.h"
#include "SearchOutput.h"

#define NUM_NTS_PER_LINE		60
#define MAX_NUM_NTS_PER_LINE	200

#define ORIGINAL_GENOME	0
#define GENOME_SEGMENT	1
#define GENOME_NONE		2

//#define GENOME_SEGMENT_TAG "#Frequency of ACGT"
#define GENOME_TAG_COMMENT "#"
#define GENOME_TAG_FILTERSEARCH_DIRECTION_PLUS "Plus"
#define GENOME_TAG_FILTERSEARCH_DIRECTION_MINUS "Minus"

//#define GENOME_SEARCH_DIRECTION_PLUS		-1
#define GENOME_SEARCH_DIRECTION_PLUS		0
#define GENOME_SEARCH_DIRECTION_MINUS		1

#define GENOME_SEARCH_DIRECTION_ONE			2
#define GENOME_SEARCH_DIRECTION_BOTH		3

#define GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS		"Plus search result"	//"#Plus search result"
#define GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS	"Minus search result"	//"#Minus search result"

#define GENOME_TAG_SEARCH_DIRECTION "search result"

#define GENOME_TAG_FREQUENCY "#Frequency of ACGT"
#define GENOME_TAG_NORESULT "No results have been found"

#define GENOME_TAG_FILTERTIME "#Time for total filter search "
#define GENOME_TAG_FILTERTIMEEND " hours"

#include <string.h>	//gcc3.4

/*
GenomeData supports fasta-format genome file
The first line is genome description including genome id
From the second line is genome data, ACGT, GenomeData filters out non-ACGT
GenomeData supports multiple genomes file, that is, multiple genomes in one fasta file.
Zhibin Huang
*/
class GenomeData  
{
private:
	//single genome file
	char *	SGenomeName;
	char *	SGenomeDataC;
	int	*	SGenomeDataI;

	char *	SGenomeReverseDataC;	//reverse genome character
	int	*	SGenomeReverseDataI;	//reverse genome integer

	int		SGenomeLength;
	double  SGenomeBaseFreq[4];

	//multiple genome files
	int		MGenomeNum;
	char **	MGenomeName;
	char **	MGenomeDataC;
	char **	MGenomeReverseDataC;	//reverse genome character
	int	*	MGenomeLength;
	double **MGenomeBaseFreq;

	//file type: filtered genome segment (GENOME_SEGMENT) or original genome (ORIGINAL_GENOME)
	int	*	MGenomeType;

	//filter search direction (Plus | Minus)
	int *	MSearchDirect;

	int *	MGenomeExtLeft;
	int *	MGenomeExtRight;

	float	TotalFilterSearchTime;
public:
	GenomeData();
	virtual ~GenomeData();

	//single genome file
	void	loadSingleGenomeFile(char * genomefileName);
	void	loadSingleGenomeFile_improved(char * pseudoGenomefileName);
	void	loadSingleGenomeFile_improved2(char * pseudoGenomefileName);
	char *	getSingleGenomeName();
	char *	getSingleGenomeDataC();
	int	*	getSingleGenomeDataI();

	char *	getSingleGenomeReverseDataC();	//reverse genome character
	int	*	getSingleGenomeReverseDataI();	//reverse genome integer
	void	generateSingleGenomeReverseData();
	void	freeSingleGenomeReverseData();

	int		getSingleGenomeLength();
	double *getSingleGenomeBaseFreq();
	void	printSingleGenomeBaseFreq();
	void	printSingleGenomeData();
	void	freeSingleGenomeName();
	void	freeSingleGenomeData();

//	void	freeSingleGenome();

	//multiple genome file
	
	void	loadMultipleGenomes(char * genomefileName);		//deprecated
	void	loadMultipleGenomes2(char * genomefileName);
	
	int		getMultipleGenomeNum();
	int *	getMultipleGenomeLength();
	int *	getMultipleGenomeType();		//
	int *	getMultipleSearchDirect();
	int *	getMultipleGenomeExtLeft();	//
	int *	getMultipleGenomeExtRight();	//
	char **	getMultipleGenomeName();
	char **	getMultipleGenomeDataC();
	double **getMultipleGenomeBaseFreq();
	float	getTotalFilterSearchTime();
	
	//20080909
	long	getTotalGenomeLength();

	void	freeMultipleGenomeLength();
	void	freeMultipleGenomeType();
	void	freeMultipleSearchDirect();
	void	freeMultipleGenomeExtLeft();	//
	void	freeMultipleGenomeExtRight();	//
	void	freeMultipleGenomeName();
	void	freeMultipleGenomeData();
	void	freeMultipleGenomeBaseFreq();

//	void	freeMultipleGenome();

	void	printMultipleGenomeName();
	void	printMultipleGenomeData();
	void	printMultipleGenomeLength();
	void	printMultipleGenomeType();
	void	printMultipleSearchDirect();
	void	printMultipleGenomeExtLeft();	//
	void	printMultipleGenomeExtRight();	//
	void	printMultipleGenomeBaseFreq();

	//reverse
	void	generateMultipleGenomeReverseData();
	char **	getMultipleGenomeReverseDataC();	//reverse genome character
	void	freeMultipleGenomeReverseData();
	void	printMultipleGenomeReverseData();

};

#endif // !defined(AFX_GENOMEDATA_H__E09901D8_4A6A_4AF3_ACBC_AEF08B96E963__INCLUDED_)
