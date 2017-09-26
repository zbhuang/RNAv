// GenomeData.cpp: implementation of the GenomeData class.
//
//////////////////////////////////////////////////////////////////////

#include "GenomeData.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GenomeData::GenomeData()
{
	//single genome file
    SGenomeName			= NULL;
	SGenomeDataC		= NULL;
	SGenomeDataI		= NULL;
	SGenomeReverseDataC = NULL;
	SGenomeReverseDataI = NULL;

	SGenomeLength=0;
    for(int i=0; i<4; i++ )
        SGenomeBaseFreq[i]=-1.0;

	//multiple genome files
	MGenomeName			= NULL;
	MGenomeDataC		= NULL;
	MGenomeReverseDataC = NULL;
	MGenomeLength		= NULL;
	MGenomeType			= NULL;
	MSearchDirect		= NULL;
	MGenomeBaseFreq		= NULL;

	TotalFilterSearchTime = 0.0;
}

GenomeData::~GenomeData()
{
	//single genome file
	this->freeSingleGenomeName();
	this->freeSingleGenomeData();
	this->freeSingleGenomeReverseData();

	//multiple genome files
	this->freeMultipleGenomeName();
	this->freeMultipleGenomeData();
	this->freeMultipleGenomeReverseData();
	this->freeMultipleGenomeBaseFreq();
  	this->freeMultipleGenomeLength();
  	this->freeMultipleGenomeType();		//
	this->freeMultipleSearchDirect();
	this->freeMultipleGenomeExtLeft();	//
	this->freeMultipleGenomeExtRight();	//
}

char * GenomeData::getSingleGenomeName()
{
  return this->SGenomeName;
}

int  GenomeData::getSingleGenomeLength()
{
  return this->SGenomeLength;
}

char * GenomeData::getSingleGenomeDataC()
{
  return this->SGenomeDataC;
}

int  * GenomeData::getSingleGenomeDataI()
{
  return this->SGenomeDataI;
}
//reverse genome data
char * GenomeData::getSingleGenomeReverseDataC()
{
  return this->SGenomeReverseDataC;
}

int  * GenomeData::getSingleGenomeReverseDataI()
{
  return this->SGenomeReverseDataI;
}

double * GenomeData::getSingleGenomeBaseFreq()
{
  return this->SGenomeBaseFreq;
}

//generate single genome's reverse sequence
void GenomeData::generateSingleGenomeReverseData()
{
	int iLength = this->SGenomeLength;
	//reverse genome
	this->SGenomeReverseDataI = new int[iLength];
	this->SGenomeReverseDataC = new char[iLength + 1];

	char  ch;
	int i;
	//reverse genome
	for(i=0; i<iLength; i++)
	{
		ch = this->SGenomeDataC[i];
		if(ch=='a' || ch=='A')
			ch = 'T';
		else if(ch=='c' || ch=='C')
			ch = 'G';
		else if(ch=='g' || ch=='G')
			ch = 'C';
		else
			ch = 'A';
	
		this->SGenomeReverseDataC[iLength-1-i] = ch;
		this->SGenomeReverseDataI[iLength-1-i] = chartonum(ch);
	}
	this->SGenomeReverseDataC[iLength] = '\0';
}

//1. Count the number of amino acid in the genome file
//2. Count the frequency of amino acid in the genome file
//3. Load the genome file data
void GenomeData::loadSingleGenomeFile(char * genomeFileName)
{
	int iLength = strlen(genomeFileName);
	this->SGenomeName = new char[iLength + 1];
	strcpy(this->SGenomeName, genomeFileName);

	FILE *  fp;
	fp = fopen(genomeFileName, "r");
	if(fp == NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}

	char  ch;
	int    numbase[5]={0,0,0,0,0};
	iLength = 0;

	//dealing with the fasta genome file
	ch=fgetc(fp);
	if(ch == '>')
	{
		//skip the first line
		char   array[500];	//assuming the length of the first line will not be 500 nts long.
		fgets(array,100,fp);
	}
	while((ch=fgetc(fp)) != EOF)
	{
		if(	ch=='a' || ch=='A' ||
			ch=='c' || ch=='C' ||
			ch=='g' || ch=='G' ||
			ch=='u' || ch=='U' ||
			ch=='t' || ch=='T') 
		{
			numbase[chartonum(ch)]++;
			iLength++;
		}
	}
	fclose(fp);

	int i;
	int allbasenum = numbase[1]+numbase[2]+numbase[3]+numbase[4];
	for(i=1; i<5; i++) {
		SGenomeBaseFreq[i-1]=(double)numbase[i]/allbasenum;
	}  
  
	this->SGenomeLength = iLength;
	this->SGenomeDataI = new int[iLength];
	this->SGenomeDataC = new char[iLength + 1];

	fp = fopen(genomeFileName, "r");
	if(fp == NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}

	int  index = 0;
	ch=fgetc(fp);
	if(ch == '>')
	{
		//skip the first line
		char   array[500];	//assuming the length of the first line will not be 500 nts long.
		fgets(array,100,fp);
	}
	while((ch=fgetc(fp)) != EOF)
	{
		if(	ch=='a' || ch=='A' ||
			ch=='c' || ch=='C' ||
			ch=='g' || ch=='G' ||
			ch=='u' || ch=='U' ||
			ch=='t' || ch=='T')
		{
			this->SGenomeDataC[index] = ch;
			this->SGenomeDataI[index] = chartonum(ch);
			index++;
		}
	}
	this->SGenomeDataC[index] = '\0';
	fclose(fp);

	this->freeSingleGenomeName();
}

void GenomeData::loadSingleGenomeFile_improved(char * pseudoGenomeFileName)
{
	string str_line;
	string genomename;

	ifstream f_in;
	f_in.open(pseudoGenomeFileName, ios::in);
	if (!f_in)
	{
		cerr<<"can't open genome file "<<pseudoGenomeFileName<<endl;
		exit(1);
	}

	//dealing with the fasta format genome file
	//the first line is genome name
	getline(f_in, str_line);
	if(str_line.find_first_of(">") == 0)
		genomename = str_line.substr(1, str_line.find_first_of(" ")-1);
	else
		cerr<<"This genome file "<<pseudoGenomeFileName<<" is not fasta format, pls check it"<<endl;

	int iLength = strlen(genomename.c_str());
	this->SGenomeName = new char[iLength + 1];
	strcpy(this->SGenomeName, genomename.c_str());

	char  ch;
	int    numbase[5]={0,0,0,0,0};
	iLength = 0;

	unsigned int i=0;
	while ((!f_in.eof())) {
		getline(f_in, str_line);
		for(i=0; i<str_line.length(); i++)
		{
			ch = str_line[i];
			//filter out nts not A/C/G/T/U
			if(	ch=='a' || ch=='A' ||
				ch=='c' || ch=='C' ||
				ch=='g' || ch=='G' ||
				ch=='u' || ch=='U' ||
				ch=='t' || ch=='T') 
			{
				numbase[chartonum(ch)]++;
				iLength++;
			}
		}
	}
	f_in.close();
	f_in.clear();

	int allbasenum = numbase[1]+numbase[2]+numbase[3]+numbase[4];
	for(i=1; i<5; i++) {
		SGenomeBaseFreq[i-1]=(double)numbase[i]/allbasenum;
	}  
  
	this->SGenomeLength = iLength;
	this->SGenomeDataI = new int[iLength];
	this->SGenomeDataC = new char[iLength + 1];

	f_in.open(pseudoGenomeFileName, ios::in);
	if (!f_in)
	{
		cerr<<"can't open genome file "<<pseudoGenomeFileName<<endl;
		exit(1);
	}

	//skip the first fasta line
	getline(f_in, str_line);

	unsigned int index = 0;
	while ((!f_in.eof())) {
		getline(f_in, str_line);
		for(i=0; i<str_line.length(); i++)
		{
			ch = str_line[i];
			//filter out nts not A/C/G/T/U
			if(	ch=='a' || ch=='A' ||
				ch=='c' || ch=='C' ||
				ch=='g' || ch=='G' ||
				ch=='u' || ch=='U' ||
				ch=='t' || ch=='T') 
			{
				this->SGenomeDataC[index] = ch;
				this->SGenomeDataI[index] = chartonum(ch);
				index++;
			}
		}
	}
	this->SGenomeDataC[index] = '\0';
	f_in.close();

	this->freeSingleGenomeName();
}
void GenomeData::loadSingleGenomeFile_improved2(char * pseudoGenomeFileName)
{
	const int line_max_num = 200;
	char	c_line[line_max_num-1];
	FILE *  fp;
	fp = fopen(pseudoGenomeFileName, "r");
	if(fp == NULL) {
		cerr<<"Cannot open file"<<pseudoGenomeFileName<<endl;
		exit(1);
	}

	//dealing with the fasta genome file
	//get the first line
	fgets(c_line,line_max_num,fp);
	//extract the genome name
	string genomename;
	string str_line(c_line);

	if(str_line.find_first_of(">") == 0)
		genomename = str_line.substr(1, str_line.find_first_of(" ")-1);
	else
		cerr<<"This genome file "<<pseudoGenomeFileName<<" is not fasta format, pls check it"<<endl;

	int iLength = strlen(genomename.c_str());
	this->SGenomeName = new char[iLength + 1];
	strcpy(this->SGenomeName, genomename.c_str());

	int    numbase[5]={0,0,0,0,0};
	iLength = 0;

	//dealing with the other genome files
	int i = 0, seql = 0;
	char  ch;
	memset (c_line,'\0',line_max_num);
	fgets(c_line, line_max_num, fp);
	while(!feof(fp))
	{
        seql=strlen(c_line);
        for(i=0; i<seql; i++)
        {
			ch = c_line[i];
			if(	ch=='a' || ch=='A' ||
				ch=='c' || ch=='C' ||
				ch=='g' || ch=='G' ||
				ch=='u' || ch=='U' ||
				ch=='t' || ch=='T') 
			{
				numbase[chartonum(ch)]++;
				iLength++;
			}
        }
		memset (c_line,'\0',line_max_num);
		//fgets(str, seqlength+100, pfdata);
		fgets(c_line, line_max_num, fp);
	}
//	fclose(fp);

	int allbasenum = numbase[1]+numbase[2]+numbase[3]+numbase[4];
	for(i=1; i<5; i++) {
		SGenomeBaseFreq[i-1]=(double)numbase[i]/allbasenum;
	}  
  
	this->SGenomeLength = iLength;
	this->SGenomeDataI = new int[iLength];
	this->SGenomeDataC = new char[iLength + 1];

//	fp = fopen(genomeFileName, "r");
	fseek(fp, 0, SEEK_SET);

	//skip the first line
	fgets(c_line,line_max_num,fp);

	int  index = 0;
	memset (c_line,'\0',line_max_num);
	fgets(c_line, line_max_num, fp);
	while(!feof(fp))
	{
        seql=strlen(c_line);
        for(i=0; i<seql; i++)
        {
			ch = c_line[i];
			if(	ch=='a' || ch=='A' ||
				ch=='c' || ch=='C' ||
				ch=='g' || ch=='G' ||
				ch=='u' || ch=='U' ||
				ch=='t' || ch=='T') 
			{
				this->SGenomeDataC[index] = ch;
				this->SGenomeDataI[index] = chartonum(ch);
				index++;
			}
		}
		memset (c_line,'\0',line_max_num);
		fgets(c_line, line_max_num, fp);
	}
	this->SGenomeDataC[index] = '\0';
	fclose(fp);

	this->freeSingleGenomeName();
}

/*
Scan the pseudoGenomeFile three times.
The first time is to count the num of genome sequences, that is the num of ">";
The second time is to count the num of sequences for every genome;
The third time is to store the genome sequences;
*/
void GenomeData::loadMultipleGenomes(char * pseudoGenomeFileName)
{
	const int line_max_num = MAX_NUM_NTS_PER_LINE;
	char	c_line[line_max_num-1];
	FILE *  fp;
	fp = fopen(pseudoGenomeFileName, "r");
	if(fp == NULL) {
		cerr<<"Cannot open file"<<pseudoGenomeFileName<<endl;
		exit(1);
	}

	//the first time: count the number of genome
	int genome_num = 0;
	fgets(c_line,line_max_num,fp);
	while(!feof(fp))
	{
		if(c_line[0] == '>')
			genome_num++;
		fgets(c_line,line_max_num,fp);
	}
	this->MGenomeNum = genome_num;

	//the second time: count the num of nts, ACGT
	fseek(fp, 0, SEEK_SET);
	this->MGenomeLength		= new int[genome_num];
	this->MGenomeType		= new int[genome_num];	//
	this->MGenomeBaseFreq	= new double*[genome_num];
	this->MSearchDirect		= new int[genome_num];
	int index = -1;
	int i, seql, iLength;
	char ch;
	fgets(c_line,line_max_num,fp);
	while(!feof(fp))
	{
		if(c_line[0] == '>')	//the beginning of one new genome data
		{
			index++;
			this->MGenomeBaseFreq[index] = new double[4];
			iLength = 0;
			fgets(c_line,line_max_num,fp);	//the second line
			string s_line(c_line);
			if((signed int)s_line.find(GENOME_TAG_NORESULT) == -1)
			{
//				if((signed int)s_line.find(GENOME_TAG_FREQUENCY) != -1)	//contain that frequency for the original genome
				if((signed int)s_line.find(GENOME_TAG_COMMENT) != -1)	//contain that frequency for the original genome
				{
					this->MGenomeType[index] = GENOME_SEGMENT;
					//Plus/Minus filter-search result
					if((signed int)s_line.find(GENOME_TAG_FILTERSEARCH_DIRECTION_PLUS) != -1)
						this->MSearchDirect[index]= GENOME_SEARCH_DIRECTION_PLUS;
					else
						this->MSearchDirect[index]= GENOME_SEARCH_DIRECTION_MINUS;

					fgets(c_line,line_max_num,fp);
/*
					int pos_start	= s_line.find_last_of("(")+1;
					int pos_end		= s_line.find_last_of(")");
					string s_freq_info = s_line.substr(pos_start, pos_end - pos_start);
					pos_end = s_freq_info.find_first_of(", ");
					string s_freqA = s_freq_info.substr(0, pos_end);
					float d_freqA = atof(s_freqA.c_str());

					s_freq_info = s_freq_info.substr(pos_end+2, s_freq_info.length());
					pos_end = s_freq_info.find_first_of(", ");
					string s_freqC = s_freq_info.substr(0, pos_end);
					float d_freqC = atof(s_freqC.c_str());

					s_freq_info = s_freq_info.substr(pos_end+2, s_freq_info.length());
					pos_end = s_freq_info.find_first_of(", ");
					string s_freqG = s_freq_info.substr(0, pos_end);
					float d_freqG = atof(s_freqG.c_str());

					s_freq_info = s_freq_info.substr(pos_end+2, s_freq_info.length());
					pos_end = s_freq_info.length();
					string s_freqT = s_freq_info.substr(0, pos_end);
					float d_freqT = atof(s_freqT.c_str());

					//loading the frequency
					this->MGenomeBaseFreq[index][0] = d_freqA;
					this->MGenomeBaseFreq[index][1] = d_freqC;
					this->MGenomeBaseFreq[index][2] = d_freqG;
					this->MGenomeBaseFreq[index][3] = d_freqT;
*/
					fgets(c_line,line_max_num,fp);
				}
				else
				{
					this->MGenomeType[index]	= ORIGINAL_GENOME;
					this->MSearchDirect[index]	= -1;
				}

//				while(c_line[0] != '>' && !feof(fp))
				while(c_line[0] != '>' && c_line[0] != '#' && !feof(fp))
				{
					seql=strlen(c_line);
					for(i=0; i<seql; i++)
					{
						ch = c_line[i];
						if(	ch=='a' || ch=='A' ||
							ch=='c' || ch=='C' ||
							ch=='g' || ch=='G' ||
							ch=='u' || ch=='U' ||
							ch=='t' || ch=='T') 
						{
							iLength++;
						}
					}
					memset (c_line,'\0',line_max_num);
					fgets(c_line, line_max_num, fp);
				}
				this->MGenomeLength[index] = iLength;

				//skip those lines with the beginning of '#'
				while(c_line[0] == '#' && !feof(fp))
					fgets(c_line, line_max_num, fp);
			}
			else
			{
				this->MGenomeType[index] = GENOME_NONE;
				this->MSearchDirect[index]= -1;
				//skip the current line
				fgets(c_line,line_max_num,fp);
			}
		}
		else
		{
			//skip the unexpected blank line
			fgets(c_line,line_max_num,fp);
		}
	}
	
	//the third time
	int j;
	string genomename;
	fseek(fp, 0, SEEK_SET);
	this->MGenomeName		= new char*[genome_num];
	this->MGenomeDataC		= new char*[genome_num];
	this->MGenomeExtLeft	= new int[genome_num];
	this->MGenomeExtRight	= new int[genome_num];

	int ipos = -1;
	char *ptr_filtertime = NULL;
	//get the genome name
	fgets(c_line,line_max_num,fp);
	for(i=0; i<genome_num; i++)
	{
		int numbase[5]={0,0,0,0,0};

		if(this->MGenomeType[i] != GENOME_NONE)
		{
			while(c_line[0] != '>')
				fgets(c_line,line_max_num,fp);	//skip blank line

			string str_line(c_line);
//			ipos = str_line.find_first_of("(");
			ipos = str_line.find_first_of("[");
			if(ipos == -1) {
//				genomename = str_line.substr(1, str_line.find_first_of(" ")-1);
//				genomename = str_line.substr(1, str_line.find_first_of("\n")-1);
				//20080902
				if((signed int)str_line.find_first_of(" ") == -1) {
					genomename = str_line.substr(1, str_line.find_first_of("\n")-1);
				} else {
					genomename = str_line.substr(1, str_line.find_first_of(" ")-1);
				}

				this->MGenomeExtLeft[i]  = 0;
				this->MGenomeExtRight[i] = 0;
			} else {
				//
				ipos = str_line.find_first_of("(");
				genomename = str_line.substr(1, ipos-1);

				//
				ipos = str_line.find_first_of("[");
				string s_pos_info = str_line.substr(ipos+1, str_line.length()-3-ipos);
				ipos = s_pos_info.find_first_of("-");
				string s_shift_left = s_pos_info.substr(0, ipos);
				int i_shift_left = atoi(s_shift_left.c_str());

				s_pos_info = s_pos_info.substr(ipos+1, s_pos_info.length());
				string s_shift_right = s_pos_info;
				int i_shift_right = atoi(s_shift_right.c_str());

				this->MGenomeExtLeft[i]  = i_shift_left;
				this->MGenomeExtRight[i] = i_shift_right;
			}

			this->MGenomeName[i] = new char[genomename.length()+1];
			strcpy(this->MGenomeName[i], genomename.c_str());

			this->MGenomeDataC[i] = new char[this->MGenomeLength[i]+1];
			iLength = 0;
			fgets(c_line,line_max_num,fp);

			//skip the line of frequency if it is the genome segment
			if(this->MGenomeType[i] == GENOME_SEGMENT)
			{
				fgets(c_line,line_max_num,fp);	//skip search direction line
				fgets(c_line,line_max_num,fp);	//skip the line of frequency
			}
			
//			while(c_line[0] != '>' && !feof(fp))
			while(c_line[0] != '>' && c_line[0] != '#' && !feof(fp))
			{
				seql=strlen(c_line);
				for(j=0; j<seql; j++)
				{
					ch = c_line[j];
					if(	ch=='a' || ch=='A' ||
						ch=='c' || ch=='C' ||
						ch=='g' || ch=='G' ||
						ch=='u' || ch=='U' ||
						ch=='t' || ch=='T') 
					{
//						if(this->MGenomeType[i] == ORIGINAL_GENOME)
							numbase[chartonum(ch)]++;
						this->MGenomeDataC[i][iLength++] = ch;
					}
				}
				memset (c_line,'\0',line_max_num);
				fgets(c_line, line_max_num, fp);
			}
			this->MGenomeDataC[i][iLength] = '\0';

//			if(this->MGenomeType[i] == ORIGINAL_GENOME)
//			{
				int allbasenum = numbase[1]+numbase[2]+numbase[3]+numbase[4];
				for(j=1; j<5; j++) {
					this->MGenomeBaseFreq[i][j-1]=(double)numbase[j]/allbasenum;
				}  
//			}

			//skip those lines with the beginning of '#'
			while(c_line[0] == '#' && !feof(fp)) {
				ptr_filtertime = strstr(c_line, GENOME_TAG_FILTERTIME);
				if(ptr_filtertime != NULL) {
					string filtertime_line(c_line);
					ipos = filtertime_line.find(GENOME_TAG_FILTERTIMEEND);
					int length = strlen(GENOME_TAG_FILTERTIME);
					string s_filtertime_info = filtertime_line.substr(length, ipos-length);
					TotalFilterSearchTime = atof(s_filtertime_info.c_str());
				}
				fgets(c_line, line_max_num, fp);
			}
		}
		else
		{
			//
			string str_line(c_line);
			genomename = str_line.substr(1, str_line.length()-1);
			this->MGenomeExtLeft[i]  = 0;
			this->MGenomeExtRight[i] = 0;

			this->MGenomeName[i] = new char[genomename.length()+1];
			strcpy(this->MGenomeName[i], genomename.c_str());

			//
			fgets(c_line,line_max_num,fp);
			//skip the current lines until meeting the new line with ">"
			while(c_line[0] != '>' && !feof(fp))
				fgets(c_line,line_max_num,fp);
		}
	}

	fclose(fp);
}

/*
Scan the pseudoGenomeFile three times.
The first time is to count the num of genome sequences, that is the num of ">";
The second time is to count the num of sequences for every genome;
The third time is to store the genome sequences;
*/
void GenomeData::loadMultipleGenomes2(char * pseudoGenomeFileName)
{
	const int line_max_num = MAX_NUM_NTS_PER_LINE;
	char	c_line[line_max_num-1];
	FILE *  fp;
	fp = fopen(pseudoGenomeFileName, "r");
	if(fp == NULL) {
		cerr<<"Cannot open file"<<pseudoGenomeFileName<<endl;
		exit(1);
	}

	//the first time: count the number of genome
	int genome_num = 0;
	fgets(c_line,line_max_num,fp);
	while(!feof(fp))
	{
		if(c_line[0] == '>')
			genome_num++;
		fgets(c_line,line_max_num,fp);
	}
	this->MGenomeNum = genome_num;

	//the second time: count the num of nts, ACGT
	fseek(fp, 0, SEEK_SET);
	this->MGenomeLength		= new int[genome_num];
	this->MGenomeType		= new int[genome_num];	//
	this->MGenomeBaseFreq	= new double*[genome_num];
	this->MSearchDirect		= new int[genome_num];
	int index = -1;
	int i, seql, iLength;
	char ch;
	string s_line;
	fgets(c_line,line_max_num,fp);
	while(!feof(fp))
	{
		//until reaching the line with the beginning of ">"
		if(c_line[0] == '>')	//the beginning of one new genome data
		{
			index++;
			this->MGenomeBaseFreq[index] = new double[4];
			iLength = 0;
			//Plus/Minus filter-search result or genome nts
			fgets(c_line,line_max_num,fp);	//the second line
			s_line.assign(c_line);
			if((signed int)s_line.find(GENOME_TAG_SEARCH_DIRECTION) != -1)
			{
				this->MGenomeType[index] = GENOME_SEGMENT;
				//Plus/Minus filter-search result
				if((signed int)s_line.find(GENOME_TAG_FILTERSEARCH_DIRECTION_PLUS) != -1)
					this->MSearchDirect[index]= GENOME_SEARCH_DIRECTION_PLUS;
				else
					this->MSearchDirect[index]= GENOME_SEARCH_DIRECTION_MINUS;

				while(s_line.find(GENOME_TAG_HIT_EXTENSION) == -1)
				{
					fgets(c_line,line_max_num,fp);
					s_line.assign(c_line);
				}
				//until reach it
				fgets(c_line,line_max_num,fp);
			}
			else
			{
				this->MGenomeType[index]	= ORIGINAL_GENOME;
				this->MSearchDirect[index]	= -1;
			}

//			while(c_line[0] != '>' && c_line[0] != '-' && !feof(fp))
			while(c_line[0] != '>' && c_line[0] != '\n' && !feof(fp))
			{
				seql=strlen(c_line);
				for(i=0; i<seql; i++)
				{
					ch = c_line[i];
					if(	ch=='a' || ch=='A' ||
						ch=='c' || ch=='C' ||
						ch=='g' || ch=='G' ||
						ch=='u' || ch=='U' ||
						ch=='t' || ch=='T') 
					{
						iLength++;
					}
				}
				memset (c_line,'\0',line_max_num);
				fgets(c_line, line_max_num, fp);
			}
			this->MGenomeLength[index] = iLength;

			while(c_line[0] != '>' && !feof(fp))
				fgets(c_line, line_max_num, fp);
		}
		else
		{
			//skip the unexpected blank line
			fgets(c_line,line_max_num,fp);
		}
	}
	
	//the third time
	int j;
	string genomename;
	fseek(fp, 0, SEEK_SET);
	this->MGenomeName		= new char*[genome_num];
	this->MGenomeDataC		= new char*[genome_num];
	this->MGenomeExtLeft	= new int[genome_num];
	this->MGenomeExtRight	= new int[genome_num];

	int ipos = -1;
	char *ptr_filtertime = NULL;
	//get the genome name
	fgets(c_line,line_max_num,fp);
	for(i=0; i<genome_num; i++)
	{
		int numbase[5]={0,0,0,0,0};

		while(c_line[0] != '>')
			fgets(c_line,line_max_num,fp);	//skip blank line

		//parse the genome name
		s_line.assign(c_line);
		ipos = s_line.find_first_of("[");
		if(ipos == -1) {
			//20080902
			if((signed int)s_line.find_first_of(" ") == -1) {
				genomename = s_line.substr(1, s_line.find_first_of("\n")-1);
			} else {
				genomename = s_line.substr(1, s_line.find_first_of(" ")-1);
			}

			this->MGenomeExtLeft[i]  = 0;
			this->MGenomeExtRight[i] = 0;
		} else {
			//
			ipos = s_line.find_first_of("(");
			genomename = s_line.substr(1, ipos-1);

			//
			ipos = s_line.find_first_of("[");
			string s_pos_info = s_line.substr(ipos+1, s_line.length()-3-ipos);
			ipos = s_pos_info.find_first_of("-");
			string s_shift_left = s_pos_info.substr(0, ipos);
			int i_shift_left = atoi(s_shift_left.c_str());

			s_pos_info = s_pos_info.substr(ipos+1, s_pos_info.length());
			string s_shift_right = s_pos_info;
			int i_shift_right = atoi(s_shift_right.c_str());

			this->MGenomeExtLeft[i]  = i_shift_left;
			this->MGenomeExtRight[i] = i_shift_right;
		}

		this->MGenomeName[i] = new char[genomename.length()+1];
		strcpy(this->MGenomeName[i], genomename.c_str());

		this->MGenomeDataC[i] = new char[this->MGenomeLength[i]+1];
		iLength = 0;

		//locate the hit extension
		fgets(c_line,line_max_num,fp);

		//skip the line of frequency if it is the genome segment
		if(this->MGenomeType[i] == GENOME_SEGMENT)
		{
			while(s_line.find(GENOME_TAG_HIT_EXTENSION) == -1)
			{
				fgets(c_line,line_max_num,fp);
				s_line.assign(c_line);
			}
			//until reach it
			fgets(c_line,line_max_num,fp);
		}
			
//		while(c_line[0] != '>' && c_line[0] != '-' && !feof(fp))
		while(c_line[0] != '>' && c_line[0] != '\n' && !feof(fp))
		{
			seql=strlen(c_line);
			for(j=0; j<seql; j++)
			{
				ch = c_line[j];
				if(	ch=='a' || ch=='A' ||
					ch=='c' || ch=='C' ||
					ch=='g' || ch=='G' ||
					ch=='u' || ch=='U' ||
					ch=='t' || ch=='T') 
				{
					numbase[chartonum(ch)]++;
					this->MGenomeDataC[i][iLength++] = ch;
				}
			}
			memset (c_line,'\0',line_max_num);
			fgets(c_line, line_max_num, fp);
		}
		this->MGenomeDataC[i][iLength] = '\0';

		int allbasenum = numbase[1]+numbase[2]+numbase[3]+numbase[4];
		for(j=1; j<5; j++) {
			this->MGenomeBaseFreq[i][j-1]=(double)numbase[j]/allbasenum;
		}  
	}

	//skip those lines until it reaches the time for filter search
	s_line.assign(c_line);
	while(s_line.find(GENOME_TAG_TOTAL_TIME_USED) == -1 && !feof(fp))
		fgets(c_line, line_max_num, fp);

	if(s_line.find(GENOME_TAG_TOTAL_TIME_USED) != -1 && !feof(fp)) {
		ipos = s_line.find(GENOME_TAG_FILTERTIMEEND);
		int length = s_line.find(GENOME_TAG_TOTAL_TIME_USED) + strlen(GENOME_TAG_FILTERTIME);
		string s_filtertime_info = s_line.substr(length, ipos-length);
		TotalFilterSearchTime = atof(s_filtertime_info.c_str());
	}

	fclose(fp);
}

void GenomeData::generateMultipleGenomeReverseData()
{
	int genome_num = this->getMultipleGenomeNum();
	this->MGenomeReverseDataC = new char*[genome_num];
	int i, j, iLength = 0;
	char  ch;
	for(i=0; i<genome_num; i++)
	{
		iLength = this->MGenomeLength[i];
		this->MGenomeReverseDataC[i] = new char[iLength+1];
		//reverse genome
		for(j=0; j<iLength; j++)
		{
			ch = this->MGenomeDataC[i][j];
			if(ch=='a' || ch=='A')
				ch = 'T';
			else if(ch=='c' || ch=='C')
				ch = 'G';
			else if(ch=='g' || ch=='G')
				ch = 'C';
			else
				ch = 'A';
		
			this->MGenomeReverseDataC[i][iLength-1-j] = ch;
		}
		this->MGenomeReverseDataC[i][iLength] = '\0';
	}
}

void GenomeData::freeSingleGenomeReverseData()
{
  if(SGenomeReverseDataI != NULL)
  {
    delete [ ] SGenomeReverseDataI;
    SGenomeReverseDataI = NULL;
  }

  if(SGenomeReverseDataC != NULL)
  {
    delete [ ] SGenomeReverseDataC;
    SGenomeReverseDataC = NULL;
  }
}

void GenomeData::printSingleGenomeBaseFreq()
{
	cout<<"Percentage in genome sequence"<<endl;
	int i=0;
	for(i=0; i<4; i++) {
		printf("%c:\t%f\n", numtochar(i+1), this->SGenomeBaseFreq[i]);
	} 
}

void GenomeData::printSingleGenomeData()
{
	int i, iLength = this->getSingleGenomeLength();
	for(i=0; i<iLength; i++)
	{
		cout<<this->SGenomeDataC[i];
	}
	cout<<endl;
}

void GenomeData::freeSingleGenomeName()
{
	if(SGenomeName != NULL)
	{
		delete [ ] SGenomeName;
		SGenomeName = NULL;
	}
}

void GenomeData::freeSingleGenomeData()
{
	if(SGenomeDataI != NULL)
	{
		delete [ ] SGenomeDataI;
		SGenomeDataI = NULL;
	}

	if(SGenomeDataC != NULL)
	{
		delete [ ] SGenomeDataC;
		SGenomeDataC = NULL;
	}
}

float GenomeData::getTotalFilterSearchTime()
{
	return this->TotalFilterSearchTime;
}
int	GenomeData::getMultipleGenomeNum()
{
	return this->MGenomeNum;
}
int * GenomeData::getMultipleGenomeLength()
{
	return this->MGenomeLength;
}
int * GenomeData::getMultipleGenomeType()
{
	return this->MGenomeType;
}
int * GenomeData::getMultipleSearchDirect()
{
	return this->MSearchDirect;
}
int * GenomeData::getMultipleGenomeExtLeft()
{
	return this->MGenomeExtLeft;
}
int * GenomeData::getMultipleGenomeExtRight()
{
	return this->MGenomeExtRight;
}
char ** GenomeData::getMultipleGenomeName()
{
	return this->MGenomeName;
}
char ** GenomeData::getMultipleGenomeDataC()
{
	return this->MGenomeDataC;
}
double ** GenomeData::getMultipleGenomeBaseFreq()
{
	return this->MGenomeBaseFreq;
}
char ** GenomeData::getMultipleGenomeReverseDataC()
{
	return this->MGenomeReverseDataC;
}

void GenomeData::freeMultipleGenomeLength()
{
	if(MGenomeLength != NULL)
	{
		delete [] MGenomeLength;
		MGenomeLength = NULL;
	}
}
void GenomeData::freeMultipleGenomeType()
{
	if(MGenomeType != NULL)
	{
		delete [] MGenomeType;
		MGenomeType = NULL;
	}
}
void GenomeData::freeMultipleSearchDirect()
{
	if(MSearchDirect != NULL)
	{
		delete [] MSearchDirect;
		MSearchDirect = NULL;
	}
}
void GenomeData::freeMultipleGenomeExtLeft()
{
	if(MGenomeExtLeft != NULL)
	{
		delete [] MGenomeExtLeft;
		MGenomeExtLeft = NULL;
	}
}
void GenomeData::freeMultipleGenomeExtRight()
{
	if(MGenomeExtRight != NULL)
	{
		delete [] MGenomeExtRight;
		MGenomeExtRight = NULL;
	}
}
void GenomeData::freeMultipleGenomeName()
{
	int i;
	if(MGenomeName != NULL) {
		for(i=0; i<this->getMultipleGenomeNum(); i++) {
			delete [] MGenomeName[i];
			MGenomeName[i] = NULL;
		}
		delete [] MGenomeName;
		MGenomeName = NULL;
	}
}
void GenomeData::freeMultipleGenomeData()
{
	int i;
	if(MGenomeDataC != NULL) {
		for(i=0; i<this->getMultipleGenomeNum(); i++) {
			//if GENOME_NONE then no memeory allocated for MGenomeDataC[i]
			if(MGenomeType[i] != GENOME_NONE)
			{
				delete [] MGenomeDataC[i];
				MGenomeDataC[i] = NULL;
			}
		}
		delete [] MGenomeDataC;
		MGenomeDataC = NULL;
	}
}
void GenomeData::freeMultipleGenomeBaseFreq()
{
	int i;
	if(MGenomeBaseFreq != NULL) {
		for(i=0; i<this->getMultipleGenomeNum(); i++) {
			delete [] MGenomeBaseFreq[i];
			MGenomeBaseFreq[i] = NULL;
		}
		delete [] MGenomeBaseFreq;
		MGenomeBaseFreq = NULL;
	}
}
void GenomeData::freeMultipleGenomeReverseData()
{
	int i;
	if(MGenomeReverseDataC != NULL) {
		for(i=0; i<this->getMultipleGenomeNum(); i++) {
			delete [] MGenomeReverseDataC[i];
			MGenomeReverseDataC[i] = NULL;
		}
		delete [] MGenomeReverseDataC;
		MGenomeReverseDataC = NULL;
	}
}

void GenomeData::printMultipleGenomeName()
{
	int i;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
		cout<<MGenomeName[i]<<endl;
}
void GenomeData::printMultipleGenomeData()
{
	int i;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
		cout<<MGenomeDataC[i]<<endl;
}
void GenomeData::printMultipleGenomeLength()
{
	int i;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
		cout<<MGenomeLength[i]<<endl;
}
void GenomeData::printMultipleGenomeType()
{
	int i;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
		cout<<MGenomeType[i]<<endl;
}
void GenomeData::printMultipleGenomeExtLeft()
{
	int i;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
		cout<<MGenomeExtLeft[i]<<endl;
}
void GenomeData::printMultipleGenomeExtRight()
{
	int i;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
		cout<<MGenomeExtRight[i]<<endl;
}
void GenomeData::printMultipleGenomeBaseFreq()
{
	int i, j;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
	{
		cout<<"Percentage in genome "<<MGenomeName[i]<<endl;
		for(j=0; j<4; j++) {
			cout<<numtochar(j+1)<<":\t"<<this->MGenomeBaseFreq[i][j]<<endl;
		}
		cout<<endl;
	}
}
void GenomeData::printMultipleGenomeReverseData()
{
	int i;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
		cout<<MGenomeReverseDataC[i]<<endl;
}

long GenomeData::getTotalGenomeLength()
{
	long total_genome_length = 0;
	int i;
	for(i=0; i<this->getMultipleGenomeNum(); i++)
		total_genome_length += MGenomeLength[i];
	return total_genome_length;
}
