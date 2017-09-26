//#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class PeakFinding
{
private:
	double T_mh; //T_mh is the minimum height of the local maximum from the local minimum
	int    T_ew; //T_ew is the half window size 
	int    Filter_Beg, Filter_End;
	vector <double> SmoothScore; //Score after the second smooth;
	vector <vector<int> > LocalValue; // whether it is local minimum (-1), local maximum (1), or nothing (0); 
public:
	PeakFinding(const vector <double> &oriscore, double theTmh, int theTew);
	~PeakFinding(void);
	vector <double> Smooth(const vector <double> &);
	void getresult();
	int Min(int, int);
	int Max(int, int);
	int LocalMinMax();
	void RegionFind();
	int GetBeg();
	int GetEnd();
//	vector <double> SmoothScore; //Score after the second smooth , for debug;
};
