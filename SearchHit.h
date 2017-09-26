// SearchHit.h: interface for the Hit class.
//
//////////////////////////////////////////////////////////////////////
//solve the problem of 'class' type redefinition	20081121
#pragma once

#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include "common.h"
#include "SearchOutput.h"

#define	TWO_PASTALINE	2

class SearchHit  
{
private:
	string	genome_name;
	double	align_score;
	string	direction;
	int		pos_start;
	int		pos_end;
	string	folded_stru_seq;	//Folded structure	  - seq
	string	folded_stru_align;	//Folded structure	  - alignment
	string	folded_stru_pasta_charid;
	string	stru_seq;			//structure alignment - seq
	string	stru_align;			//structure alignment - alignment
	string	stru_pasta_charid;
//	string	pasta_charid;
	int		pasta_line_num;

	//limiting the number of empty-stem. zbhuang
	int		empty_stem_num;
public:
	SearchHit();
	virtual ~SearchHit();

	//get method
	string	getGenomeName();
	double	getAlignScore();
	string	getDirection();
	int		getPosStart();
	int		getPosEnd();
	string	getFoldedStruSeq();
	string	getFoldedStruAlign();
	string	getStruSeq();
	string	getStruAlign();
	int		getPastaLineNum();
//	string	getPastaCharId();
	string	getFoldedStruPastaCharId();
	string	getStruPastaCharId();
	int		getEmptyStemNum();		//limiting the number of empty-stem. zbhuang

	//set method
	void	setGenomeName(string str);
	void	setAlignScore(double score);
	void	setDirection(string dir);
	void	setPosStart(int pos);
	void	setPosEnd(int pos);
	void	setFoldedStruSeq(string str);
	void	setFoldedStruAlign(string str);
	void	setStruSeq(string str);
	void	setStruAlign(string str);
	void	setPastaLineNum(int num);
//	void	setPastaCharId(string str);
	void	setFoldedStruPastaCharId(string str);
	void	setStruPastaCharId(string str);
	void	setEmptyStemNum(int num);	//limiting the number of empty-stem. zbhuang

	//print out for debuging
	void	print();
	void	formattedPrint();

	//sorting by align_score
	bool operator< (const SearchHit& rhs) const
	{ 
		return align_score > rhs.align_score; 
	}
};
