#pragma warning( disable: 4786 )

//#pragma once
#include <vector>
#include <string>

using namespace std;

class FilterSelection
{
private:
	//vector <string> aa;//for test;
	vector <string>  ori_alignment, commentline;
	string thefilename, mystandard;
	int myNSize, myWinSize, myTew, begposi, endposi, filterflag;//filterflag =0 means no filter;
	double myTmh, myGapThreshold;

public:
	FilterSelection(string str, double theTmh=0.05, string thestandard="TACGU",
					 int NSize=7, int WinSize=3,  int theTew=10, 
					double GapThreshold=0.5);
	~FilterSelection(void);
	vector <string>  getFilter();
    void readfile();
	int getbegposi();
	int getendposi();
	int getfilterflag();
};
