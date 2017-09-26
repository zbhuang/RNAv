#pragma warning( disable: 4786 )

#include "Score.h"
#include <math.h>
//#include <stdio.h>

Score::Score(const vector <string> &oridata, string thestandard, 
			 int NSize, int WinSize, double GapThreshold)
{
	int i, j, k, seqnum, length, category, reallength, letter_acc[6];
	vector <int> CurComposition;
	string curchar;
	copydata = oridata;
	standardstr=thestandard;//"_ACGU";
	NL = standardstr.length();
	NeighbourSize = NSize;//7;
	SmoothingWinSize = WinSize;//3;
	lambda = 0.5;//////////////////////
	length = oridata[0].length();
	seqnum = oridata.size();
	//composition.resize(length);
    for (i = 0; i < length; i++)
		CurComposition.resize(6);
	for (i = 0; i < length; i++) 
	{
		for (k = 0; k < 6; k++)
			CurComposition[k]=0;
			//composition[i][k]=0; //initializaion
		for (j = 1; j < seqnum; j++) //skip the first line which is algnment information like "AAA...aaa"
		{
			curchar = oridata[j].substr(i, 1);
			category = standardstr.find(curchar);
			if (category != string::npos)
				CurComposition[category]++;
			//else
			//	j=j;
				//composition [i][category]++;
		}
		for (k = 0; k < 5; k++)
			CurComposition[5] += CurComposition[k];
		    //composition [i][5] += composition [i][k];
		if (((double)CurComposition[5]/(seqnum-1))>=GapThreshold)//!!!!!
		{
			composition.push_back(CurComposition);	
			recordcol.push_back(i);
		} 
	}
	reallength = composition.size();
	for (i=0; i<6; i++)
		letter_acc[i]=0;
	for (i=0; i<reallength; i++)
		for (j=0; j<6; j++)
			letter_acc[j] += composition[i][j];

	for (i=0; i<5; i++)
		q[i] = (double)letter_acc[i] / letter_acc[5]; 
}

Score::~Score(void)
{
}

double Score::Freq(int ColId, string letter)
{
	double frequency=0;
	int category;
	category = standardstr.find (letter);
	if (category != string::npos)
		frequency = (double)composition[ColId][category]/composition[ColId][5];
	return frequency;
}

double Score::Freq(int ColId, int LetterId)
{
	double frequency=0;
	frequency = (double)composition[ColId][LetterId]/composition[ColId][5];
	return frequency;
}

int Score::GetRecordCol(int index)
{
	return recordcol[index];
}

double Score::RelativeEntropy (double *p, double *q)
{
	double value=0;
	int i;
	for (i=0; i<5; i++)
		if ((p[i]>0)&&(q[i]>0))
			value += p[i]*log(p[i]/q[i]);
	return value;
}

double Score::JSD (int ColId)
{
	double pc[5], r[5], thejsd;
	int i;
	for (i=0; i<5; i++)
	{
		pc[i] = (double)composition[ColId][i]/composition[ColId][5];
		r[i]  = lambda * pc[i] + (1-lambda) * q[i];
	}
	thejsd = lambda*RelativeEntropy(pc,r)+lambda*RelativeEntropy(q,r);
	return thejsd;
}

int Score::Min(int firstvalue, int secondvalue)
{
	int value;
	if (firstvalue < secondvalue)
		value = firstvalue;
	else 
		value = secondvalue;
	return value;
}

int Score::Max(int firstvalue, int secondvalue)
{
	int value;
	if (firstvalue > secondvalue)
		value = firstvalue;
	else 
		value = secondvalue;
	return value;
}

vector <double> Score::Smooth_Score ()
{
	int i, j, length, size, begincol, endcol;
	double Value;
	vector <double> S_Jsd; // smooth arcs
	vector <double> JSDScore;
	length = composition.size();
	JSDScore.resize(length);
	S_Jsd.resize(length);
	if (SmoothingWinSize%2 == 0)
		size = SmoothingWinSize - 1;
	else
		size = SmoothingWinSize;
	for (i=0; i<length; i++)
		JSDScore[i] = JSD(i);
	for (i=0; i<length; i++)
	{
		Value=0;
		begincol = Max (0, i-(size-1)/2);
		endcol = Min (length-1, i+(size-1)/2);
		for (j=begincol; j<=endcol; j++)
			Value += JSDScore[j];
		Value = Value / (endcol - begincol + 1);
		S_Jsd[i] = Value;
	}
	return S_Jsd;
}
/*

double Score::ColValue (int ColId)
{
	double value=0, curfreq;
	int i, catenum;
	catenum = standardstr.length();
	for (i=0; i<catenum; i++)
	{
		curfreq = Freq(ColId, standardstr.substr(i, 1));
		if (curfreq > 0)
			value += (-1)*curfreq*log(curfreq)/log((double)2);
	}
	return value;
}



double Score::Logos (int ColId)
{
	double logosvalue, h_max, h;
	h_max = log((double)Min(NL, composition.size()))/log((double)2);
	h=ColValue(ColId);
	logosvalue = h_max - h;
	return logosvalue;
}

int Score::RelatedCount(int Col_i, int Col_j, char char_p, char char_q)
{
	int line, Count=0;
	for (line=1; line<(int)copydata.size(); line++)
		if((copydata[line][recordcol[Col_i]]==char_p)&&(copydata[line][recordcol[Col_j]]==char_q))
			Count++;
	return Count;
}

double Score::InfoDepen(int Col_i, int Col_j)
{
	double h_ij=0;
	int CharNum, p_id, q_id, relatedcount;
	CharNum = standardstr.length();
	for (p_id=0; p_id<CharNum; p_id++)
		for (q_id=0; q_id<CharNum; q_id++)
		{
			relatedcount=RelatedCount(Col_i, Col_j, 
				standardstr[p_id], standardstr[q_id]);
			if (relatedcount>0)
				h_ij += (double)composition[Col_i][p_id]/(double)composition[Col_i][5]
						*((double)((double)relatedcount/composition[Col_i][p_id]))*
						(log((double)(composition[Col_i][p_id]/(double)relatedcount))
						/log((double)2));
		}
	return h_ij;
}

double Score::FunctionDepen(int Col_i, int Col_j)
{
	double Fd;
	Fd = 1- InfoDepen(Col_i, Col_j)/(log((double)(copydata.size()-1))/log((double) 2));
	return Fd;
}

double Score::ARCS(int Col_i)
{
	double value =0;
	int neighbour, begincol, endcol, Col_j;
	if (NeighbourSize%2 == 0)
		neighbour = NeighbourSize - 1;
	else
		neighbour = NeighbourSize;
	begincol = Max(0, Col_i-(neighbour-1)/2);
	endcol   = Min(composition.size()-1, Col_i+(neighbour-1)/2);
	for (Col_j=begincol; Col_j<=endcol; Col_j++)
		value +=  FunctionDepen(Col_j, Col_i) * Logos(Col_j);
	return value;
}

vector <double> Score::Smooth_ARCS()
{
	int i, j, length, size, begincol, endcol;
	double Value;
	vector <double> S_ARCS; // smooth arcs
	vector <double> ArcsScore;
	length = composition.size();
	ArcsScore.resize(length);
	S_ARCS.resize(length);
	if (SmoothingWinSize%2 == 0)
		size = SmoothingWinSize - 1;
	else
		size = SmoothingWinSize;
	for (i=0; i<length; i++)
		ArcsScore[i] = ARCS(i);
	for (i=0; i<length; i++)
	{
		Value=0;
		begincol = Max (0, i-(size-1)/2);
		endcol = Min (length-1, i+(size-1)/2);
		for (j=begincol; j<=endcol; j++)
			Value += ArcsScore[j];
		Value = Value / (endcol - begincol + 1);
		S_ARCS[i] = Value;
	}
	return S_ARCS;
}
*/


