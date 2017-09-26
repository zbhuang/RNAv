#include "PeakFinding.h"

PeakFinding::PeakFinding(const vector <double> &oriscore, 
						 double theTmh, int theTew)
{
	T_mh = theTmh;
	T_ew = theTew;
	Filter_Beg=-1;
	Filter_End=-1;
	SmoothScore=Smooth(oriscore);
	//
}

PeakFinding::~PeakFinding(void)
{
}

vector <double> PeakFinding::Smooth(const vector <double> &oriscore)
{
	int i, j, length, begincol, endcol;
	double Value;
	vector <double> smoothresult;

	length = oriscore.size();
	smoothresult.resize(length);
	for (i=0; i<length; i++)
	{
		Value=0;
		begincol = Max (0, i-T_ew);
		endcol = Min (length-1, i+T_ew);
		for (j=begincol; j<=endcol; j++)
			Value += oriscore[j];
		Value = Value / (endcol - begincol + 1);
		smoothresult[i] = Value;
	}

	return smoothresult;
}

int PeakFinding::Min(int firstvalue, int secondvalue)
{
	int value;
	if (firstvalue < secondvalue)
		value = firstvalue;
	else 
		value = secondvalue;
	return value;
}

int PeakFinding::Max(int firstvalue, int secondvalue)
{
	int value;
	if (firstvalue > secondvalue)
		value = firstvalue;
	else 
		value = secondvalue;
	return value;
}

void PeakFinding::getresult()
{
	//int i_test;
	//LocalMinMax();
	RegionFind();
	//cout<<Filter_Beg<<endl;
	//cout<<Filter_End<<endl;
	//for (i_test=0; i_test<(int)SmoothScore.size(); i_test++)
	//	cout<<LocalValue[i_test]<<endl;
}

int PeakFinding::LocalMinMax()
{
	int i;
	double HighestValue=-100;
	int HighestPosi;
	vector <int> record(2);// position, statur(min(-1) or max(+1))
	//LocalValue.resize(SmoothScore.size());
	for (i=0; i<(int)SmoothScore.size(); i++)
	{
		if (((i==0)&&(SmoothScore[i]<=SmoothScore[i+1] ))||
			((i==((int)SmoothScore.size()-1))&&
			(SmoothScore[i-1]>=SmoothScore[i])))
		{
			record[0] = i;
			record[1] = -1;
			LocalValue.push_back(record);
			continue;
		}
		if (((i==0)&&(SmoothScore[i]>=SmoothScore[i+1] ))||
			((i==((int)SmoothScore.size()-1))&&
			(SmoothScore[i-1]<=SmoothScore[i])))
		{
			record[0] = i;
			record[1] = 1;
			LocalValue.push_back(record);
			if (SmoothScore[i]>HighestValue)
			{
				HighestValue = SmoothScore[i];
				HighestPosi  = i;
			}
			continue;
		}
		if ((SmoothScore[i]<=SmoothScore[i-1])&&      //minimum
			(SmoothScore[i]<=SmoothScore[i+1]))
		{
			if ((SmoothScore[i]==SmoothScore[i-1])&&
				(SmoothScore[i]==SmoothScore[i+1]))
				continue;
			else
			{
				record[0] = i;
				record[1] = -1;
				LocalValue.push_back(record);
			}
		} else if ((SmoothScore[i]>=SmoothScore[i-1])&&    //maximum
					(SmoothScore[i]>=SmoothScore[i+1]))
		{
			if ((SmoothScore[i]==SmoothScore[i-1])&&
				(SmoothScore[i]==SmoothScore[i+1]))
				continue;
			else
			{
				record[0] = i;
				record[1] = 1;
				LocalValue.push_back(record);
				if (SmoothScore[i]>HighestValue)
				{
					HighestValue = SmoothScore[i];
					HighestPosi  = i;
				}
			}
		} 
	}
	return HighestPosi;
}

void PeakFinding::RegionFind()
{
	double PeakValue;
	int i, length, PeakPosi, RecordPosi, left, right;
	PeakPosi=LocalMinMax(); //Get minimum and maximum
	length = LocalValue.size();
	PeakValue = SmoothScore [PeakPosi];
	//cout<<PeakValue<<endl;
	for (i=0; i<length; i++)
		if (LocalValue[i][0] == PeakPosi)
		{
			RecordPosi = i;
			break;
		}
	if (RecordPosi == 0)
		Filter_Beg = LocalValue[0][0];
	else
	{
		for (left=(RecordPosi-1); left>=0; left--)
			if ((LocalValue[left][1]==-1)&&
				((PeakValue-SmoothScore[LocalValue[left][0]])>=T_mh))
			{
				Filter_Beg = LocalValue[left][0];
				break;
			}
		if (left == -1)
			Filter_Beg = 0;
	}
	if (RecordPosi == (length-1))
		Filter_End = LocalValue[length-1][0];
	else
	{
		for (right=(RecordPosi+1); right<length; right++)
			if ((LocalValue[right][1]==-1)&&
				((PeakValue-SmoothScore[LocalValue[right][0]])>=T_mh))
			{
				Filter_End = LocalValue[right][0];
				break;
			}
		if (right == length)
			Filter_End = SmoothScore.size()-1;
	}
}

int PeakFinding::GetBeg()
{
	return Filter_Beg;
}

int PeakFinding::GetEnd()
{
	return Filter_End;
}
