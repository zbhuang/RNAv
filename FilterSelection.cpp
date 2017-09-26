#pragma warning( disable: 4786 )

#include "PeakFinding.h"
#include "Score.h"
#include "FilterSelection.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>

using namespace std;
FilterSelection::FilterSelection(string str, double theTmh, string thestandard,
								  int NSize, int WinSize,  int theTew,
								 double GapThreshold)
{
	thefilename    = str;
	mystandard     = thestandard;
	myNSize        = NSize;
    myWinSize      = WinSize;
	myTmh          = theTmh;
	myTew          = theTew;
	myGapThreshold = GapThreshold;
}

FilterSelection::~FilterSelection(void)
{
}

vector <string> FilterSelection::getFilter()
{
	int i, length, count=1;
	string str_firstline="";
	readfile(); //read the input file;
	vector <double> myscore;
	vector <string> outputvector;
	length = ori_alignment.size();
	Score score(ori_alignment, mystandard, myNSize, 
		myWinSize, myGapThreshold);
	myscore = score.Smooth_Score();
	PeakFinding pf(myscore, myTmh, myTew);
	pf.getresult();
///////////////////////////////////////////////////////
//	debug	//int index_zhibin=-1000, value=20000;
/*
	ofstream f_info;
	f_info.open("data_3.txt", ios::out);
	for (i=0; i<myscore.size(); i++)
		f_info<<score.JSD(i)<<endl;
	f_info<<"**********************"<<endl;
	for (i=0; i<myscore.size(); i++)
		f_info<<myscore[i]<<endl;
	f_info<<"**********************"<<endl;
	for (i=0; i<myscore.size(); i++)
		f_info<<pf.SmoothScore[i]<<endl;
	f_info<<"my region"<<endl;
	f_info<<pf.GetBeg()<<" to "<<pf.GetEnd()<<endl;
	f_info<<"zhibin region"<<endl;
	for (i=0; i<myscore.size(); i++)
		if (score.GetRecordCol(i)==777)
		{
			f_info<<i<<" to ";
			break;
		}
	for (i=0; i<myscore.size(); i++)
		if (score.GetRecordCol(i)==833)
		{
			f_info<<i<<endl;
			break;
		}
	f_info.close();
*/
/////////////////////
	begposi=score.GetRecordCol(pf.GetBeg());
	endposi=score.GetRecordCol(pf.GetEnd());
	if ((begposi>=0)&&(endposi>=0)&&((endposi-begposi)<(int)ori_alignment[0].length()))
	{
		filterflag = 1;
		outputvector.resize(length*2-1);
		for (i=begposi; i<=endposi; i++) // for HMM
			str_firstline += ".";
		outputvector[0]=str_firstline;
		for (i=1; i<length; i++)
		{ 
			if (i<=(int)commentline.size())
			{
				outputvector[count]= commentline[i-1];
				count++;
			}
			outputvector[count]=ori_alignment[i].substr(begposi,endposi-begposi+1);
			count++;
		}
	}else
		filterflag = 0;
	return outputvector;
}

void FilterSelection::readfile ()
{
	int length, index;
	string str_line;
	ifstream f_in;
	f_in.open(thefilename.c_str(), ios::in);
	if (!f_in)
	{
		cerr<<"can't find input file"<<endl;
		exit(1);
	}
	getline (f_in, str_line);
	length = str_line.length();
	ori_alignment.push_back(str_line);
	while ((!f_in.eof()))
	{
		getline (f_in, str_line);
		if (str_line == "")
			continue;
		if (str_line[0]=='>')
		{
			commentline.push_back(str_line);
			continue;
		}
		if ((int)str_line.length()<length)
			break;	
		if (str_line == "!")
			break;
		index = str_line.find("1");
		if (index != string::npos)
			continue;
		ori_alignment.push_back(str_line);

		str_line="!";
	}
	f_in.close();
}

int FilterSelection::getbegposi()
{
	return begposi;
}

int FilterSelection::getendposi()
{
	return endposi;
}
	
int FilterSelection::getfilterflag()
{
	return filterflag;
}
