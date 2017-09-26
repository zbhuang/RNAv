#include "SearchOutput.h"

char * strrpl(char *s, const char s1, const char s2)
{
	int m;
	int iLength = strlen(s);
	char *newChar = new char[iLength+1];
	for(m=0; m<iLength; m++)
	{
		if(s[m] == s1)
			newChar[m] = s2;
		else
			newChar[m] = s[m];
	}
	newChar[m] = '\0';
	return newChar;
}
/*
char * strdel(char *s, const char s1)
{
	int m;
	int iLength = strlen(s);
	int count = 0;
	char *newChar = new char[iLength+1];
	for(m=0; m<iLength; m++)
	{
		if(s[m] != s1)
			newChar[count++] = s[m];
	}
	newChar[count] = '\0';
	return newChar;
}
*/
char * strappend(char *s, const char s1)
{
    int m;
	int length = strlen(s);
	char * new_sequence = new char[length + 2];
	for(m=0; m<length; m++)
	{
		new_sequence[m] = s[m];
	}
	new_sequence[m++] = s1;
	new_sequence[m] = '\0';
	return new_sequence;
}
/*
 * From start_pos to end_pos including char on the position of end_pos
 */
char * strextract(char *s, int start_pos, int end_pos)
{
    int m;
	char * new_sequence = new char[end_pos - start_pos + 2];
	for(m=start_pos; m<=end_pos; m++) {
		new_sequence[m-start_pos] = s[m];
	}
	new_sequence[m-start_pos] = '\0';
	return new_sequence;
}

//zbhuang 20100201
int countNumChar(char *str, char ch)
{
	int length = strlen(str);
	int count = 0;
	for(int i=0; i<length; i++) {
		if(str[i] == ch)
			count++;
	}
	return count;
}

string extractfilename(string s)
{
	int b_onlyfilename = 1;
	int ipos = s.find_first_of("/");
	if(ipos != -1)
		b_onlyfilename = 0;
	if(b_onlyfilename == 1)
		return s;
	else {
		ipos = s.find_last_of("/");
		return s.substr(ipos+1, s.length());
	}
}
string printYN(int iFlag)
{
	if(iFlag)
		return string("Yes");
	else
		return string("No");
}
string printDirection(int iFlag)
{
	if(iFlag == GENOME_SEARCH_DIRECTION_BOTH)
		return string("Yes");
	else
		return string("No");
}
/*
void printGenomeSegments(char * c_seq, int num)
{
	string s_seq(c_seq);
	int length = s_seq.length();
	for(int i=0; i*num<length; i++)
	{
		cout<<s_seq.substr(i*num, num)<<endl;
	}
}
*/
void printGenomeSegments(char * uSeq, int num)
{
	char * tSeq = strrpl(uSeq, 'U', 'T');
	string s_seq(tSeq);
	int length = s_seq.length();
	for(int i=0; i*num<length; i++) {
		cout<<s_seq.substr(i*num, num)<<endl;
	}
	if(tSeq != NULL) {
		delete [] tSeq;
		tSeq = NULL;
	}
}
void printGenomeSegments(char * uSeq, int num, ofstream &outfile)
{
	char * tSeq = strrpl(uSeq, 'U', 'T');
	string s_seq(tSeq);
	int length = s_seq.length();
	for(int i=0; i*num<length; i++) {
		outfile<<s_seq.substr(i*num, num)<<endl;
	}
	if(tSeq != NULL) {
		delete [] tSeq;
		tSeq = NULL;
	}
}
//------------------------------------------------------------------
int isHmmFilterFile(char * filename)
{
	int iResult = 1;
	FILE *  pFile = fopen(filename, "r");
    if( pFile == NULL ) {
        printf(" Fail to open %s\n", filename);
        exit( 0 );
    }
	fseek(pFile, 0, SEEK_SET);

	const int line_max_num = 200;
	char c_line[line_max_num-1];
	fgets(c_line, line_max_num, pFile);
//	string s_line(c_line);
//	for(int i=0; i<s_line.length() && (iResult==0); i++)
	for(unsigned int i=0; i<strlen(c_line) && (iResult==1); i++)
	{
		if(c_line[i] != '.' && c_line[i] != '\n')
			iResult = 0;
	}
	fclose(pFile);

	return iResult;
}
//deprecated
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
							  int		iFlagPrintDebugDPSearchInfo)
{
	cout<<"-------------Parameter settings-------------"<<endl;
	cout<<"pseudocount="<<pseudocount<<endl;
	cout<<"k="<<top_k<<endl;
	cout<<"Threshold="<<threshold<<endl;
	cout<<"Number of overlap between stems="<<num_nts_overlap_between_stem<<endl;
	cout<<"Take the merge-candidate strategy in preprocessing="<<printYN(iFlagMergeCandInPreprocess)<<endl;
	cout<<"About the merge strategy, take the candidate with the (s)hortest length="<<printYN(iFlagCandwithShortestLength)<<endl;
	cout<<"iShiftNumMergeCand="<<printYN(iShiftNumMergeCand)<<endl;
	cout<<"iAllowedNullLoopInsNum="<<iAllowedNullLoopInsNum<<endl;
	cout<<"pcoeff="<<pcoeff<<endl;
	cout<<"iJumpStrategy="<<printYN(iJumpStrategy)<<endl;
	cout<<"iStepSize="<<iStepSize<<endl;
	cout<<"iScoreThresholdInJump="<<dScoreThresholdInJump<<endl;
	cout<<"Print structure alignment option="<<printYN(iFlagPrintStrAlignInfo)<<endl;
	cout<<"reverse complement search="<<printYN(iFlagSearchReverse)<<endl;
//	cout<<"Debug option="<<printYN(iFlagPrintDebugInfo)<<endl;
//	cout<<"Print the debug info of input-checking in dp search="<<printYN(iFlagPrintDebugInputInfo)<<endl;
//	cout<<"Print the debug info of tree in dp search="<<printYN(iFlagPrintDebugTreeInfo)<<endl;
//	cout<<"Print the debug info of window size in dp search="<<printYN(iFlagPrintDebugDPWinInfo)<<endl;
//	cout<<"Print the debug info of td in dp search="<<printYN(iFlagPrintDebugTDInfo)<<endl;
//	cout<<"Print the debug info of the dp search="<<printYN(iFlagPrintDebugDPSearchInfo)<<endl;
	cout<<"--------------------------------------------"<<endl;
}
*/
void printCurrentTime()
{
	time_t     now;
	struct tm  *ts;
	char       buf[80];
 
    // Get the current time
	now = time(NULL);
 
    // Format and print the time, "ddd yyyy-mm-dd hh:mm:ss zzz"
	ts = localtime(&now);
//	strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", ts);
	strftime(buf, sizeof(buf), "%H:%M:%S %Z %Y-%m-%d", ts);
    cout<<left<<setw(40)<<buf;	//<<endl;
}

void printSearchEndingTimeInfo()
{
	cout<<"* Time: ";
	printCurrentTime();
	cout<<right<<setw(NUM_NTS_PER_LINE-48)<<"*"<<endl;
}

void printOneStarLine()
{
	cout<<"************************************************************"<<endl;
}
void printAdditionalInfo()
{
	printOneStarLine();
//	cout<<"* Searched done : with RNATOPS V1.2                        *"<<endl
	cout<<"* Searched done : with RNATOPS V2.0                        *"<<endl
		<<"* By RNA-Informatics @ UGA                                 *"<<endl;
}

void printFilterSearchResultHeader(string	trainfile,
								   int		profilelength,
								   int		filterauto, 
								   int		filtertype, 
								   int		filterbegin, int filterend, 
								   string	genomefile,	//char * genomefile,
								   int		genomenum,
								   long		total_genome_length)
{
	printOneStarLine();
	cout<<"* Filtering Result                                         *"<<endl;
	cout<<"*                                                          *"<<endl;
	cout<<"* Profile file : "<<left<<setw(NUM_NTS_PER_LINE - 18)<<trainfile<<"*"<<endl;
	cout<<"* Profile length : "<<left<<setw(NUM_NTS_PER_LINE - 20)<<profilelength<<"*"<<endl;
	cout<<"*                                                          *"<<endl;
	if(filterauto)
		cout<<"* Filter generation: automatic seleted                     *"<<endl;
	else
		cout<<"* Filter generation: manual seleted                        *"<<endl;
	if(filtertype == SEARCH_HMMFILTER)
		cout<<"* Filter type : HMM                                        *"<<endl;
	else
		cout<<"* Filter type : substructure                               *"<<endl;
	cout<<"* Filter info : positions from "<<left<<setw(6)<<filterbegin<<" to "<<left<<setw(NUM_NTS_PER_LINE-42)<<filterend<<"*"<<endl;
//	cout<<"* : regions from "<<left<<setw(6)<<filterbegin<<" to "<<left<<setw(NUM_NTS_PER_LINE-28)<<filterend<<"*"<<endl;
	cout<<"* Genome file : "<<left<<setw(NUM_NTS_PER_LINE-17)<<genomefile<<"*"<<endl;
	cout<<"* Number of sequences : "<<left<<setw(NUM_NTS_PER_LINE-25)<<genomenum<<"*"<<endl;
	cout<<"* Total length of sequences : "<<left<<setw(NUM_NTS_PER_LINE-31)<<total_genome_length<<"*"<<endl;
	printOneStarLine();
}

void printFilterHitIndex(int index)
{
	cout<<endl<<"Filtering hit "<<index<<endl
		<<"---------------"<<endl;
}
void printHitPos(int start, int end)
{
	cout<<"Hit Positions: "<<start<<"-"<<end<<endl;
}
void printHitPos(int start, int end, ofstream &outfile)
{
	outfile<<"Hit Positions: "<<start<<"-"<<end<<endl;
}
void printFilterHitAlignment(char *uSeq, char *align)
{
	char * tSeq = strrpl(uSeq, 'U', 'T');
	cout<<"Alignment to the filter"<<endl
		<<tSeq<<endl
		<<align<<endl;
	if(tSeq != NULL) {
		delete [] tSeq;
		tSeq = NULL;
	}
}
void printFilterHitAlignment(char *uSeq, char *align, ofstream &outfile)
{
	char * tSeq = strrpl(uSeq, 'U', 'T');
	outfile<<"Alignment to the filter"<<endl
		<<tSeq<<endl
		<<align<<endl;
	if(tSeq != NULL) {
		delete [] tSeq;
		tSeq = NULL;
	}
}
void printHMMHitScore(double score)
{
	cout<<"HMM ";
	printHitScore(score);
}
void printHMMHitScore(double score, ofstream &outfile)
{
	outfile<<"HMM ";
	printHitScore(score, outfile);
}

void printWholeStructureHitScore(double score)
{
	cout<<"HIT ";
	printHitScore(score);
}
void printWholeStructureHitScore(double score, ofstream &outfile)
{
	outfile<<"HIT ";
	printHitScore(score, outfile);
}

void printHitScore(double score)
{
	cout<<"Alignment score = "<<score<<endl;
}
void printHitScore(double score, ofstream &outfile)
{
	outfile<<"Alignment score = "<<score<<endl;
}

void printFilterHitExtenPos(int start, int end)
{
	cout<<"Extension positions: "<<start<<"-"<<end<<endl;
}
void printFilterHitExtenPos(int start, int end, ofstream &outfile)
{
	outfile<<"Extension positions: "<<start<<"-"<<end<<endl;
}
void printFilterHitExtenNtsHeader()
{
	cout<<GENOME_TAG_HIT_EXTENSION<<endl;
}
void printFilterHitExtenNtsHeader(ofstream &outfile)
{
	outfile<<GENOME_TAG_HIT_EXTENSION<<endl;
}

//void printFilterHitExtenNtsTail()
//{
//	cout<<GENOME_TAG_HIT_EXTENSION_TAIL<<endl;
//}

void printTotalHitNum(int hit_num)
{
	cout<<"* "<<GENOME_TAG_TOTAL_NUM_HIT<<": "<<left<<setw(15)<<hit_num<<right<<setw(NUM_NTS_PER_LINE-35)<<"*"<<endl;
}

void printTotalTimeInfo(double dTime)
{
	cout<<"* "<<GENOME_TAG_TOTAL_TIME_USED<<" "<<left<<setw(15)<<dTime<<" hours "<<right<<setw(NUM_NTS_PER_LINE-40)<<"*"<<endl;
}

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
								  int		filtertype, 
								  int		filterbegin, 
								  int		filterend,
								  int		max_empty_stem_num
								  )
{
	printOneStarLine();
	cout<<"* Whole Structure Result                                   *"<<endl;
	cout<<"*                                                          *"<<endl;
	cout<<"* Profile file : "<<left<<setw(NUM_NTS_PER_LINE - 18)<<trainfile<<"*"<<endl;
	cout<<"* Profile length : "<<left<<setw(NUM_NTS_PER_LINE - 20)<<profilelength<<"*"<<endl;
	cout<<"*                                                          *"<<endl;
	//if(streamline)	//adding some filter info
	//{
		if(filtertype != SEARCH_WHOLESTRUCTURE)
		{
			if(filtertype == SEARCH_HMMFILTER)
				cout<<"* Filter type : HMM                                        *"<<endl;
			else
				cout<<"* Filter type : SUBSTRUCTURE                               *"<<endl;
			cout<<"* Filter info : positions from "<<left<<setw(6)<<filterbegin<<" to "<<left<<setw(NUM_NTS_PER_LINE-42)<<filterend<<"*"<<endl;
			cout<<"*                                                          *"<<endl;
		} else {
			//whole structure search
		}
	//}
	cout<<"* Genome file : "<<left<<setw(NUM_NTS_PER_LINE - 17)<<genomefile<<"*"<<endl;
	cout<<"* Number of sequences : "<<left<<setw(NUM_NTS_PER_LINE - 25)<<genomenum<<"*"<<endl;
	cout<<"* Total length of sequences : "<<left<<setw(NUM_NTS_PER_LINE - 31)<<total_genome_length<<"*"<<endl;
	cout<<"*                                                          *"<<endl;
	printValuesForParameters(pseudocount,
							top_k,
							threshold,
							num_nts_overlap_between_stem,
							iFlagMergeCandInPreprocess,
							iFlagCandwithShortestLength,
							iShiftNumMergeCand,
							iAllowedNullLoopInsNum,
							pcoeff,
							iJumpStrategy,
							iStepSize,
							iFlagSearchReverse,
							max_empty_stem_num);	//limiting the number of empty-stem. zbhuang
	printOneStarLine();
}
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
							  int		max_empty_stem_num)	//limiting the number of empty-stem. zbhuang
{
	cout<<"* Search parameters setting:                               *"<<endl;
	cout<<"* Pseudocount = "<<left<<setw(NUM_NTS_PER_LINE - 17)<<pseudocount<<"*"<<endl;
	cout<<"* Num of stem candidates = "<<left<<setw(NUM_NTS_PER_LINE - 28)<<top_k<<"*"<<endl;
	cout<<"* Score threshold for hits = "<<left<<setw(NUM_NTS_PER_LINE - 30)<<threshold<<"*"<<endl;
	cout<<"* Num of nt overlap between stems = "<<left<<setw(NUM_NTS_PER_LINE - 37)<<num_nts_overlap_between_stem<<"*"<<endl;
	cout<<"* Candidate representatives only = "<<left<<setw(NUM_NTS_PER_LINE - 36)<<printYN(iFlagMergeCandInPreprocess)<<"*"<<endl;
	cout<<"* Shortest candidate representatives = "<<left<<setw(NUM_NTS_PER_LINE - 40)<<printYN(iFlagCandwithShortestLength)<<"*"<<endl;
	cout<<"* IShiftNumMergeCand = "<<left<<setw(NUM_NTS_PER_LINE - 24)<<printYN(iShiftNumMergeCand)<<"*"<<endl;
	cout<<"* Nts allowed in null loops = "<<left<<setw(NUM_NTS_PER_LINE - 31)<<iAllowedNullLoopInsNum<<"*"<<endl;
	cout<<"* Pcoeff = "<<left<<setw(NUM_NTS_PER_LINE - 12)<<pcoeff<<"*"<<endl;
	cout<<"* Search with jump = "<<left<<setw(NUM_NTS_PER_LINE - 22)<<printYN(iJumpStrategy)<<"*"<<endl;
	cout<<"* Search step size = "<<left<<setw(NUM_NTS_PER_LINE - 22)<<iStepSize<<"*"<<endl;
	cout<<"* Search reversed complement sequence = "<<left<<setw(NUM_NTS_PER_LINE - 41)<<printDirection(iFlagSearchReverse)<<"*"<<endl;
	cout<<"* Max num of empty stem = "<<left<<setw(NUM_NTS_PER_LINE - 27)<<max_empty_stem_num<<"*"<<endl;
}

int getProfileLength(char * trainfile)
{
	int length = 0;
	FILE *  pFile = fopen(trainfile, "r");
    if( pFile == NULL ) {
        printf(" Fail to open %s\n", trainfile);
        exit( 0 );
    }
	fseek(pFile, 0, SEEK_SET);

	const int line_max_num = 5000;
	char c_line[line_max_num-1];
	fgets(c_line, line_max_num, pFile);
	fclose(pFile);

	//Yf Wang report this bug. 20081029
//	length = strlen(c_line);
	string str_line(c_line);
	length = str_line.find_first_of("\n");
	return length;
}

void printWholeSearchHitIndex(int index)
{
	cout<<endl<<"Whole structure search hit "<<index<<endl
		<<"----------------------------"<<endl;
}
void printWholeSearchHitIndex(int index, ofstream &outfile)
{
	outfile<<"Whole structure search hit "<<index<<endl
		<<"----------------------------"<<endl;
}

void printEValue(double evalue)
{
	cout<<"E-value = "<<evalue<<endl;
}

void formatStringOutput(string illustration, string line_one, string line_two, string line_thr, int start, int end, int flag, ofstream &outfile)
{
	//NUM_NTS_PER_LINE
	outfile<<endl<<illustration<<endl;

	string line_one_unit;
	string line_two_unit;
	string line_thr_unit;
	int length = line_one.length();
	for(int i=0; i*NUM_NTS_PER_LINE<length; i++) {
		line_one_unit.assign(line_one.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
		if(flag)
			line_two_unit.assign(line_two.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
		line_thr_unit.assign(line_thr.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
/*
		if((i+1)*NUM_NTS_PER_LINE<length) {
			line_one_unit.assign(line_one.substr(i*NUM_NTS_PER_LINE, (i+1)*NUM_NTS_PER_LINE));
			if(flag)
				line_two_unit.assign(line_two.substr(i*NUM_NTS_PER_LINE, (i+1)*NUM_NTS_PER_LINE));
			line_thr_unit.assign(line_thr.substr(i*NUM_NTS_PER_LINE, (i+1)*NUM_NTS_PER_LINE));
		} else {
			line_one_unit.assign(line_one.substr(i*NUM_NTS_PER_LINE, length));
			if(flag)
				line_two_unit.assign(line_two.substr(i*NUM_NTS_PER_LINE, length));
			line_thr_unit.assign(line_thr.substr(i*NUM_NTS_PER_LINE, length));
		}
*/
		formatOneOutput(line_one_unit, line_two_unit, line_thr_unit, i, start, flag, outfile);
	}
}

void formatOneOutput(string line_one_unit, string line_two_unit, string line_thr_unit, int index, int start, int flag, ofstream &outfile)
{
	int length_one = line_one_unit.length();
	int length_two = 0;
	if(flag)
		length_two = line_two_unit.length();
	int length_thr = line_thr_unit.length();
	if(flag) {
		outfile<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_one_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_one)<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_two_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_two)<<endl
			<<right<<setw(10)<<(start+1+index*NUM_NTS_PER_LINE)<<" "<<line_thr_unit<<" "<<(start+1+index*NUM_NTS_PER_LINE+length_thr)<<endl;
	} else {
		outfile<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_one_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_one)<<endl
			<<right<<setw(10)<<(start+1+index*NUM_NTS_PER_LINE)<<" "<<line_thr_unit<<" "<<(start+1+index*NUM_NTS_PER_LINE+length_thr)<<endl;
	}
}

void formatStringOutput(string illustration, string line_one, string line_two, string line_thr, int start, int end, int flag)
{
	//NUM_NTS_PER_LINE
	cout<<endl<<illustration<<endl;

	string line_one_unit;
	string line_two_unit;
	string line_thr_unit;
	int length = line_one.length();
	for(int i=0; i*NUM_NTS_PER_LINE<length; i++) {
		line_one_unit.assign(line_one.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
		if(flag)
			line_two_unit.assign(line_two.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
		line_thr_unit.assign(line_thr.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
		formatOneOutput(line_one_unit, line_two_unit, line_thr_unit, i, start, flag);
	}
}

void formatOneOutput(string line_one_unit, string line_two_unit, string line_thr_unit, int index, int start, int flag)
{
	int length_one = line_one_unit.length();
	int length_two = 0;
	if(flag)
		length_two = line_two_unit.length();
	int length_thr = line_thr_unit.length();
	if(flag) {
		cout<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_one_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_one)<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_two_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_two)<<endl
			<<right<<setw(10)<<(start+1+index*NUM_NTS_PER_LINE)<<" "<<line_thr_unit<<" "<<(start+1+index*NUM_NTS_PER_LINE+length_thr)<<endl;
	} else {
		cout<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_one_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_one)<<endl
			<<right<<setw(10)<<(start+1+index*NUM_NTS_PER_LINE)<<" "<<line_thr_unit<<" "<<(start+1+index*NUM_NTS_PER_LINE+length_thr)<<endl;
	}
}

