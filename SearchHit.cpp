// SearchHit.cpp: implementation of the SearchHit class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning(disable: 4786)

#include "SearchHit.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SearchHit::SearchHit()
{
	this->setPosStart(NEGATIVE_ONE);
	this->setPosEnd(NEGATIVE_ONE);
}

SearchHit::~SearchHit()
{
}

string SearchHit::getGenomeName()
{
	return this->genome_name;
}
double SearchHit::getAlignScore()
{
	return this->align_score;
}

string SearchHit::getDirection()
{
	return this->direction;
}
int	SearchHit::getPosStart()
{
	return this->pos_start;
}
int	SearchHit::getPosEnd()
{
	return this->pos_end;
}
string SearchHit::getFoldedStruSeq()
{
	return this->folded_stru_seq;
}
string SearchHit::getFoldedStruAlign()
{
	return this->folded_stru_align;
}
string SearchHit::getStruSeq()
{
	return this->stru_seq;
}
string SearchHit::getStruAlign()
{
	return this->stru_align;
}
int	SearchHit::getPastaLineNum()
{
	return this->pasta_line_num;
}
//string Hit::getPastaCharId()
//{
//	return this->pasta_charid;
//}
string SearchHit::getFoldedStruPastaCharId()
{
	return this->folded_stru_pasta_charid;
}
string SearchHit::getStruPastaCharId()
{
	return this->stru_pasta_charid;
}
int	SearchHit::getEmptyStemNum()
{
	return this->empty_stem_num;
}

void SearchHit::setGenomeName(string str)
{
	this->genome_name = str;
}
void SearchHit::setAlignScore(double score)
{
	this->align_score = score;
}
void SearchHit::setDirection(string dir)
{
	this->direction = dir;
}
void SearchHit::setPosStart(int pos)
{
	this->pos_start = pos;
}
void SearchHit::setPosEnd(int pos)
{
	this->pos_end = pos;
}
void SearchHit::setFoldedStruSeq(string str)
{
	this->folded_stru_seq = str;
}
void SearchHit::setFoldedStruAlign(string str)
{
	this->folded_stru_align = str;
}
void SearchHit::setStruSeq(string str)
{
	this->stru_seq = str;
}
void SearchHit::setStruAlign(string str)
{
	this->stru_align = str;
}
void SearchHit::setPastaLineNum(int num)
{
	this->pasta_line_num = num;
}
//void SearchHit::setPastaCharId(string str)
//{
//	this->pasta_charid = str;
//}
void SearchHit::setFoldedStruPastaCharId(string str)
{
	this->folded_stru_pasta_charid = str;
}
void SearchHit::setStruPastaCharId(string str)
{
	this->stru_pasta_charid = str;
}
void SearchHit::setEmptyStemNum(int num)
{
	this->empty_stem_num = num;
}

void SearchHit::print()
{
	cout<<"Alignment score   "<<align_score<<endl
		<<"empty stem num    "<<empty_stem_num<<endl	//limiting the number of empty-stem. zbhuang
		<<"direction         "<<direction<<endl
		<<"pos_start         "<<pos_start<<endl
		<<"pos_end           "<<pos_end<<endl
		<<"folded_stru_seq   "<<folded_stru_seq<<endl
		<<"folded_stru_align "<<folded_stru_align<<endl
		<<"stru_seq          "<<stru_seq<<endl
		<<"stru_align        "<<stru_align<<endl;
}

void SearchHit::formattedPrint()
{
	cout<<this->genome_name<<endl
		<<this->direction<<endl;
	printHitPos(this->pos_start, this->pos_end);
	printWholeStructureHitScore(this->align_score);
	cout<<"Num of empty stem   = "<<this->getEmptyStemNum()<<endl;	//limiting the number of empty-stem. zbhuang
	if(this->getPastaLineNum() == TWO_PASTALINE)
	{
		formatStringOutput(	string("Folded structure"), 
							this->getFoldedStruAlign(),
							this->getFoldedStruPastaCharId(),		
							this->getFoldedStruSeq(),
							this->getPosStart(),
							this->getPosEnd(),
							1);

		formatStringOutput(	string("Structure alignment"), 
							this->getStruAlign(),
							this->getStruPastaCharId(),
							this->getStruSeq(),
							this->getPosStart(),
							this->getPosEnd(),
							1);
	} else {
		formatStringOutput(	string("Folded structure"), 
							this->getFoldedStruAlign(),
							this->getFoldedStruPastaCharId(),
							this->getFoldedStruSeq(),
							this->getPosStart(),
							this->getPosEnd(),
							0);

		formatStringOutput(	string("Structure alignment"), 
							this->getStruAlign(),
							this->getStruPastaCharId(),
							this->getStruSeq(),
							this->getPosStart(),
							this->getPosEnd(),
							0);
	}
}
