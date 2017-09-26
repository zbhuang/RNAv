#pragma warning( disable: 4786 )

//#pragma once
#include <vector>
#include <string>

using namespace std;

class Score
{
private:
	vector <string> copydata;
	vector <vector<int> > composition;
	vector <int> recordcol;
	string standardstr;
	int NL, NeighbourSize, SmoothingWinSize;
	double lambda, q[5]; //q[] stores backgroud distribution
public:
	Score(const vector <string> &oridata, string thestandard, 
		int NSize, int WinSize, double GapThreshold);
	~Score(void);
	double Freq(int ColId, string letter);
	double Freq(int ColId, int LetterId);
	int GetRecordCol(int);
	double RelativeEntropy (double *p, double *q); // probability and backgroud probability
    double JSD (int); //return score
	vector <double> Smooth_Score ();
	int Min(int, int);
	int Max(int, int);
	/*
	double ColValue (int);
	double Logos (int);
	int RelatedCount(int Col_i, int Col_j, char str_p, char str_q);
	double InfoDepen(int Col_i, int Col_j);
	double FunctionDepen(int Col_i, int Col_j);
	double ARCS(int);
	vector <double> Smooth_ARCS ();
	*/
};

