// StructureSearch.cpp: implementation of the StructureSearch class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning(disable: 4786)

#include "StructureSearch.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

StructureSearch::StructureSearch()
{

}

StructureSearch::~StructureSearch()
{

}

/*
Input parameters:
				train data
				genome data
				hmm model
*/
void StructureSearch::hmm_search(
								Tree_bag * &	wRoot,
								DynamicP &		wDP,
								Stem * &		wPStemModel,
								Loop * &		wPLoopModel,
								int				winsize,
								DataLoader &	trainLoader,
								int		iFlagSearchReverse,
								char *	train_file,
								int		iFlagAutoFilter,
								char *	filter_file,	//
								int		filterbegin,	//
								int		filterend,		//
								int		genome_num,
								int	*	genome_length,
								char **	genome_name,
								char **	genome_sequence,
								char **	rc_genome_sequence,
								double **genome_baseFreq,
								int		runmode,
								int		iFlagPrintFilterInfo,	//
								char *	output_file,
								int		max_empty_stem_num
								)
{
	int		iterate = 0;

	float	total_cpu_time = 0.0;
	float	hmm_search_time = 0.0;
	int		search_count = 0;

	//support saving output into files 20081121
	ofstream outfile;
	outfile.open(output_file, ios::app | ios::out | ios::in);

	for(iterate=0; iterate<genome_num; iterate++)
	{
		//support hit sorting 20081121
		vector<SearchHit> vHit;

		//nc search
		hmm_search_time = 0.0;
		this->hmm_search_oneway(
							wRoot,
							wDP,
							wPStemModel,
							wPLoopModel,
							winsize,
							trainLoader,
							GENOME_SEARCH_DIRECTION_PLUS,
							train_file,
							iFlagAutoFilter,
							filter_file,	//
							filterbegin,	//
							filterend,		//
							genome_num,
							genome_length[iterate],
							genome_name[iterate],
							genome_sequence[iterate],
							genome_baseFreq[iterate][0],	//A
							genome_baseFreq[iterate][1],	//C
							genome_baseFreq[iterate][2],	//G
							genome_baseFreq[iterate][3],	//U
							runmode,
							iFlagPrintFilterInfo,	//
							hmm_search_time,
							search_count,
							vHit,
							outfile,
							max_empty_stem_num
							);
		total_cpu_time += hmm_search_time;

		if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_BOTH)
		{
			hmm_search_time = 0.0;
			this->hmm_search_oneway(
								wRoot,
								wDP,
								wPStemModel,
								wPLoopModel,
								winsize,
								trainLoader,
								GENOME_SEARCH_DIRECTION_MINUS,	//iFlagSearchReverse,
								train_file,
								iFlagAutoFilter,
								filter_file,	//
								filterbegin,	//
								filterend,		//
								genome_num,
								genome_length[iterate],
								genome_name[iterate],
								rc_genome_sequence[iterate],
								genome_baseFreq[iterate][3],	//U
								genome_baseFreq[iterate][2],	//G
								genome_baseFreq[iterate][1],	//C
								genome_baseFreq[iterate][0],	//A
								runmode,
								iFlagPrintFilterInfo,	//
								hmm_search_time,
								search_count,
								vHit,
								outfile,
								max_empty_stem_num
								);
			total_cpu_time += hmm_search_time;
		}
		//support hit sorting 20081121
		sort(vHit.begin(), vHit.end());
		search_count = vHit.size();
		for(int i=0; i<search_count; i++)
		{
			cout<<endl
				<<"Whole structure search hit "<<(i+1)<<endl
				<<"----------------------------"<<endl;
			SearchHit one = vHit[i];
//			one.print();
			one.formattedPrint();
		}

	}

	//support saving output into files 20081121
	outfile.close();

	//if(iFlagShowTimeInfo)
	{
		cout<<endl;
		printAdditionalInfo();
		printTotalHitNum(search_count);	//20080915
		printTotalTimeInfo(total_cpu_time);
		printSearchEndingTimeInfo();
		printOneStarLine();
	}
}

void StructureSearch::hmm_search_oneway(	
					Tree_bag * &	wRoot,
					DynamicP &		wDP,
					Stem * &		wPStemModel,
					Loop * &		wPLoopModel,
					int				wWinsize,
					DataLoader &	trainLoader,
					int		iFlagSearchReverse,
					char *	train_file,
					int		iFlagAutoFilter,
					char *	filter_file,	//
					int		filterbegin,	//
					int		filterend,		//
					int		genome_num,
					int	 	genome_length,
					char *	genome_name,
					char *	genome_sequence,
					double	genome_baseFreq_A,
					double	genome_baseFreq_C,
					double	genome_baseFreq_G,
					double	genome_baseFreq_U,
					int		runmode,
					int		iFlagPrintFilterInfo,	//
					float &	hmm_search_time,
					int	 &	count,
					vector<SearchHit> &vHit,
					ofstream &outfile,
					int		max_empty_stem_num
					)
{
	int		i, j;
	clock_t start = 0, end = 0;
	float total_cpu_time = 0.0;
	float structure_search_time = 0.0;

	//if(iFlagShowTimeInfo)
		start = clock();

	double	pdcntgap= 0.001;
	int		numlp	= 0;
	double	probv = 0.0;

	DataLoader filterLoader;

	if(iFlagAutoFilter == TRUE)
	{
		//calling hmm filter selection
		vector <string> cc;
		string thestandard="TACGU";//"_ACGU";

		FilterSelection filterselection(train_file);
		cc=filterselection.getFilter();

		if (filterselection.getfilterflag()== 0) {
			cout<<"No filter is selected, pls checking the filter selection. Thanks."<<endl;
			exit(1);
		}
		else
		{
			filterbegin = filterselection.getbegposi();
			filterend	= filterselection.getendposi();

			filterLoader.inputData(cc);
		}
	} else {//zbhuang 20100209
		//
		filterLoader.inputData(filter_file);
	}

	int algpos[2], loopNo=0;

	int winsize = 0;
	int original_winsize = filterLoader.getScanWinLength(3, '-');
	winsize = original_winsize*2;

	//should call other functions to get these two values
	int extension_left	= 0;
	int extension_right = 0;
	int extension_pos_start = 0;
	int extension_pos_end	= 0;

	filterLoader.setBaseFreq(genome_baseFreq_A, 
							genome_baseFreq_C, 
							genome_baseFreq_G, 
							genome_baseFreq_U);

	HMMBuilder loopbd;
//	loopbd.setPseudocnt(pdcnt, pdcntgap);
	loopbd.setPseudocnt(PSEUDOCNT, pdcntgap);
	loopbd.buildHMM(filterLoader);
	Loop *plp=loopbd.getAllLoops(numlp);
    
	loopNo	= 0;
	probv	= INVLDPROB;
	for(i=0; i<2; i++) {
		algpos[i]=-1;
	}
   
	char *pStr = new char[winsize+1];
    
	double	optimal_score;
	int		pre_end_pos = -1;
	int		optimal_pos_start = -1, optimal_pos_end = -1;	//for hmm-filter
	optimal_score = SMALLEST;

	FVTBSearch vtber(1);
	vtber.setHitThreshold(0.0);
	char *retAlg=NULL, *retSeq=NULL;
	char *pre_retAlg=NULL, *pre_retSeq=NULL;	//

	//store the subgenome based on hmmfilter result
	vector <string> subgenome_based_on_hmmfilter;

	char * c_subgenome_seq = NULL;
	int tmp_count = 0;	//just for the dfilter

//just for the debug convenience
//int i_start = 0, i_end = 0;
//i_start = 130;
//i_end	= 135;
//	for(i=i_start; i<=i_end; i++)
//	for(i=i_start; i<=i_end-winsize; i++)
	for(i=0; i<=genome_length-winsize; i++)
	{
		//1. set the search sequence
		for(j=0; j<winsize; j++)
			pStr[j] = genome_sequence[i+j];
		pStr[j] = '\0';

		//2. search for all candidates
		FVTBCandidate * fcd = vtber.SearchOneNoMerge(&plp[loopNo], pStr);
        
		if(fcd == NULL)  //no hit
			continue;
		else 
		{             //find hit
			fcd->getPosition(algpos);
			fcd->getBothShortSeqs(&retSeq, &retAlg);
			algpos[0]=algpos[0]+i;
			algpos[1]=algpos[1]+i;
			probv = fcd->getScore();

			//merging the result
			if(optimal_score == SMALLEST)
			{
				//the first time we meet the candidate hit.
				if(optimal_score < probv) {
					optimal_score		= probv;
					optimal_pos_start	= algpos[0];
					optimal_pos_end		= algpos[1];
					pre_end_pos			= algpos[1];	//

					pre_retAlg = new char[strlen(retAlg)+1];
					strcpy(pre_retAlg, retAlg);
					pre_retSeq = new char[strlen(retSeq)+1];
					strcpy(pre_retSeq, retSeq);

					//detect memory leak problem in rnatops zbhuang 20090815
					if(retSeq != NULL) { delete [] retSeq; retSeq = NULL; }
					if(retAlg != NULL) { delete [] retAlg; retAlg = NULL; }
				}
			}
			else
			{
				//not the first time of meeting the candidate hit
				if(algpos[0] < pre_end_pos)	//pointing to the same hit
				{
					//belongs to the same candidate hit
					if(optimal_score < probv) {
						optimal_score		= probv;
						optimal_pos_start	= algpos[0];
						optimal_pos_end		= algpos[1];
						pre_end_pos			= algpos[1];

						if( pre_retSeq != NULL ) {   delete [ ] pre_retSeq;  pre_retSeq = NULL;  }
						if( pre_retAlg != NULL ) {   delete [ ] pre_retAlg;  pre_retAlg = NULL;  }
						pre_retAlg = new char[strlen(retAlg)+1];
						strcpy(pre_retAlg, retAlg);
						pre_retSeq = new char[strlen(retSeq)+1];
						strcpy(pre_retSeq, retSeq);

					}
					//detect memory leak problem in rnatops zbhuang 20090815
					if(retSeq != NULL) { delete [] retSeq; retSeq = NULL; }
					if(retAlg != NULL) { delete [] retAlg; retAlg = NULL; }
				}
				else	//other candidate hit is coming
				{
					trainLoader.getBothExtLength(extension_left, extension_right);	//support no-scanning 20081029
					count++;

					if(optimal_pos_start > extension_left)
					{
						extension_pos_start = optimal_pos_start-extension_left;
						extension_pos_end	= optimal_pos_end+extension_right;

						int length_genome_seq = strlen(genome_sequence);
						if(length_genome_seq <= extension_pos_end)
							extension_pos_end = length_genome_seq-1;
						c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
					}
					else
					{
						extension_pos_start = 0;
						extension_pos_end	= optimal_pos_end+extension_right;
//						c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
						int length_genome_seq = strlen(genome_sequence);
						if(length_genome_seq <= extension_pos_end)
							extension_pos_end = length_genome_seq-1;
						c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
					}
					//do the whole structure search
					//if(iFlagStreamLine) 
					{
						count -= 1;
						if(iFlagPrintFilterInfo) {
//							tmp_count = count;

							outfile<<endl;
							printHMMHitScore(optimal_score, outfile);
							outfile<<"Hmm hit positions: "<<optimal_pos_start<<"-"<<optimal_pos_end<<endl;
							printFilterHitAlignment(pre_retSeq, pre_retAlg, outfile);
							printFilterHitExtenPos(extension_pos_start, extension_pos_end, outfile);
							printFilterHitExtenNtsHeader(outfile);
							printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE, outfile);
						}
						end = clock();
						float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
						cpu_time = cpu_time/3600;
						total_cpu_time += cpu_time;

						structure_search_time = 0.0;
						//support hit sorting 20081121
						SearchHit onehit = this->structure_search_based_on_one_hmmfilterresult(
																			wRoot,
																			wDP,
																			wPStemModel,
																			wPLoopModel,
																			winsize,
																			iFlagSearchReverse,	//PLUS | MINUS
																			train_file,
																			iFlagAutoFilter,
																			filter_file,	//
																			filterbegin,	//support no-scanning 20081106
																			filterend,		//support no-scanning 20081106
																			genome_name,
																			c_subgenome_seq,
																			(extension_pos_end-extension_pos_start+1),
																			optimal_pos_start,
																			extension_left,
																			optimal_pos_end,
																			extension_right,
																			runmode,
																			count,	//pre_hit_num,			//
																			structure_search_time,
																			outfile,
																			max_empty_stem_num
																			);	//
						total_cpu_time += structure_search_time;
						start = clock();
						if(onehit.getPosStart() != -1)
							vHit.push_back(onehit);
					}

					if(c_subgenome_seq != NULL) {
						delete [] c_subgenome_seq;
						c_subgenome_seq = NULL;
					}

					if( retSeq != NULL ) {   delete [ ] retSeq;  retSeq = NULL;  }
					if( retAlg != NULL ) {   delete [ ] retAlg;  retAlg = NULL;  }

					optimal_score		= probv;
					optimal_pos_start	= algpos[0];
					optimal_pos_end		= algpos[1];
					pre_end_pos			= algpos[1];
				}
			}
			delete [] fcd;
		}
	}
	if(pStr != NULL) { delete [] pStr; pStr = NULL; }

	//20080709 segmentation fault
	if(optimal_pos_start != -1 && optimal_pos_end != -1)
	{
		trainLoader.getBothExtLength(extension_left, extension_right);	//support no-scanning 20081029
		count++;

		if(optimal_pos_start > extension_left)
		{
			extension_pos_start = optimal_pos_start-extension_left;
			extension_pos_end	= optimal_pos_end+extension_right;
//			cout<<extension_pos_start<<endl<<extension_pos_end<<endl;
//			cout<<genome_sequence<<endl<<strlen(genome_sequence)<<endl;
			
			//here one problem is the genome_sequence is not long enough for the extension_pos_start~extension_pos_end
			//zbhuang 20100125
			int length_genome_seq = strlen(genome_sequence);
			if(length_genome_seq <= extension_pos_end)
				extension_pos_end = length_genome_seq-1;
			c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
//			cout<<c_subgenome_seq<<endl<<strlen(c_subgenome_seq)<<endl;
		}
		else
		{
			extension_pos_start = 0;
			extension_pos_end	= optimal_pos_end+extension_right;
//			c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
			int length_genome_seq = strlen(genome_sequence);
			if(length_genome_seq <= extension_pos_end)
				extension_pos_end = length_genome_seq-1;
			c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
		}
		//do the whole structure search
		//if(iFlagStreamLine) 
		{
			count -= 1;
			if(iFlagPrintFilterInfo) {
//				tmp_count = count;

				outfile<<endl;
				printHMMHitScore(optimal_score, outfile);
				outfile<<"Hmm hit positions: "<<optimal_pos_start<<"-"<<optimal_pos_end<<endl;
				printFilterHitAlignment(pre_retSeq, pre_retAlg, outfile);
				printFilterHitExtenPos(extension_pos_start, extension_pos_end, outfile);
				printFilterHitExtenNtsHeader(outfile);
				printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE, outfile);
			}
			end = clock();
			float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
			cpu_time = cpu_time/3600;
			total_cpu_time += cpu_time;

			structure_search_time = 0.0;
			//support hit sorting 20081121
			SearchHit onehit = this->structure_search_based_on_one_hmmfilterresult(
																wRoot,
																wDP,
																wPStemModel,
																wPLoopModel,
																winsize,
																iFlagSearchReverse,	//PLUS | MINUS
																train_file,
																iFlagAutoFilter,
																filter_file,	//
																filterbegin,	//support no-scanning 20081106
																filterend,		//support no-scanning 20081106
																genome_name,
																c_subgenome_seq,
																(extension_pos_end-extension_pos_start+1),//c_subgenome_seq_length,
																optimal_pos_start,
																extension_left,
																optimal_pos_end,
																extension_right,
																runmode,
																count,	//pre_hit_num,			//
																structure_search_time,
																outfile,
																max_empty_stem_num
																);	//
			total_cpu_time += structure_search_time;
			start = clock();
			if(onehit.getPosStart() != -1)
				vHit.push_back(onehit);
		}

		if(c_subgenome_seq != NULL) {
			delete [] c_subgenome_seq;
			c_subgenome_seq = NULL;
		}
	}

	if( retSeq != NULL ) {   delete [ ] retSeq;  retSeq = NULL;  }
	if( retAlg != NULL ) {   delete [ ] retAlg;  retAlg = NULL;  }

	if( pre_retSeq != NULL ) {   delete [ ] pre_retSeq;  pre_retSeq = NULL;  }
	if( pre_retAlg != NULL ) {   delete [ ] pre_retAlg;  pre_retAlg = NULL;  }

	loopbd.freeAllLoops();

	//if(iFlagShowTimeInfo)
	{
		end = clock();
		float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
		cpu_time = cpu_time/3600;
		total_cpu_time += cpu_time;
		hmm_search_time += total_cpu_time;	//20080906
	}
}

SearchHit StructureSearch::structure_search_based_on_one_hmmfilterresult(
																	Tree_bag * &	wRoot,
																	DynamicP &		wDP,
																	Stem * &		wPStemModel,
																	Loop * &		wPLoopModel,
																	int		winsize,
																	int		iFlagSearchReverse,
																	char *	train_file,
																	int		iFlagAutoFilter,
																	char *	filter_file,	//
																	int		filterbegin,	//
																	int		filterend,		//
																	char *	genome_name,
																	char *	genome_segment,
																	int		genome_segment_length,
																	int		optimal_pos_start,
																	int		extension_left,
																	int		optimal_pos_end,
																	int		extension_right,
																	int		runmode,
																	int	 &	pre_hit_num,
																	float &	structure_hmm_search_time,
																	ofstream &outfile,
																	int		max_empty_stem_num
																	)
{
	clock_t start_all, end_all;
	float total_time = 0.0;
	start_all = clock();

	int j, genome_length;
	//prepare for the whole structure search
	int		inside_genome_num = 1;
				
	int	*	inside_genome_type;
	inside_genome_type = new int[1];
	inside_genome_type[0] = GENOME_SEGMENT;

	//support no-scanning 20081106
	int *	filterhitbeginpos;
	filterhitbeginpos = new int[1];
	filterhitbeginpos[0] = extension_left;

	int *	inside_genome_extleft;
	inside_genome_extleft = new int[1];

//	int *	inside_genome_extright;
//	inside_genome_extright = new int[1];

	if(optimal_pos_start > extension_left) {
		inside_genome_extleft[0] = optimal_pos_start-extension_left;
//		inside_genome_extright[0] = optimal_pos_end+extension_right;
	} else {
		inside_genome_extleft[0] = 0;
//		inside_genome_extright[0] = optimal_pos_end+extension_right;
	}

//	genome_length = inside_genome_extright[0] - inside_genome_extleft[0] + 1;
	genome_length = genome_segment_length;

	int	*	inside_genome_length;
	inside_genome_length = new int[1];
	inside_genome_length[0] = genome_length;

	char **	inside_genome_name;
	inside_genome_name = new char*[1];
	inside_genome_name[0] = new char[strlen(genome_name)+1];
	strcpy(inside_genome_name[0], genome_name);

	char **	inside_genome_sequence;
	inside_genome_sequence = new char*[1];
	inside_genome_sequence[0] = new char[strlen(genome_segment)+1];
	strcpy(inside_genome_sequence[0], genome_segment);

	char ch;

	double ** inside_genome_baseFreq;
	inside_genome_baseFreq = new double*[1];
	inside_genome_baseFreq[0] = new double[4];
	int numbase[5]={0,0,0,0,0};
	for(j=0; j<genome_length; j++)
	{
		ch = genome_segment[j];
		if(	ch=='a' || ch=='A' ||
			ch=='c' || ch=='C' ||
			ch=='g' || ch=='G' ||
			ch=='u' || ch=='U' ||
			ch=='t' || ch=='T') 
		{
			numbase[chartonum(ch)]++;
		}
	}
	int allbasenum = numbase[1]+numbase[2]+numbase[3]+numbase[4];
	for(j=1; j<5; j++) {
		inside_genome_baseFreq[0][j-1]=(double)numbase[j]/allbasenum;
	}  

	end_all = clock();
	total_time = ((float) (end_all - start_all)) / CLOCKS_PER_SEC;
	total_time = total_time/3600;

	float structure_search_time = 0.0;
	//here the whole structure search should be executed. 20080803
	SearchHit onehit = this->SearchStructureHitInGenomeSegment(
							wRoot,
							wDP,
							wPStemModel,
							wPLoopModel,
							winsize,
							iFlagSearchReverse,
							train_file,
							iFlagAutoFilter,
							filter_file,	//
							filterbegin,	//
							filterend,		//
							inside_genome_num,
							inside_genome_type,
							NULL,
							inside_genome_length,
							inside_genome_name,
							inside_genome_sequence,
							NULL,	//inside_rc_genome_sequence,
							inside_genome_baseFreq,
							filterhitbeginpos,
							inside_genome_extleft,
							0.0,	//filtertime
							runmode,	//runmode,
							SEARCH_WHOLESTRUCTURE,	//SEARCH_STREAMLINE);//
							1,				//streamline_flag
							pre_hit_num,
							structure_search_time,	//
							SEARCH_FILTERHIT_BASED_WHOLESTRUCTURE,	//20081029
							outfile,
							max_empty_stem_num
							);
	
	total_time += structure_search_time;
	structure_hmm_search_time += total_time;

	if(inside_genome_type != NULL) {
		delete [] inside_genome_type;
		inside_genome_type = NULL;
	}

	if(inside_genome_length != NULL) {
		delete [] inside_genome_length;
		inside_genome_length = NULL;
	}

	if(inside_genome_name != NULL) {
		delete [] inside_genome_name[0];
		inside_genome_name[0] = NULL;
		delete [] inside_genome_name;
		inside_genome_name = NULL;
	}

	if(inside_genome_sequence != NULL) {
		delete [] inside_genome_sequence[0];
		inside_genome_sequence[0] = NULL;
		delete [] inside_genome_sequence;
		inside_genome_sequence = NULL;
	}

	if(inside_genome_baseFreq != NULL) {
		delete [] inside_genome_baseFreq[0];
		inside_genome_baseFreq[0] = NULL;
		delete [] inside_genome_baseFreq;
		inside_genome_baseFreq = NULL;
	}

	if(filterhitbeginpos != NULL) {
		delete [] filterhitbeginpos;
		filterhitbeginpos = NULL;
	}

	if(inside_genome_extleft != NULL) {
		delete [] inside_genome_extleft;
		inside_genome_extleft = NULL;
	}

//	if(inside_genome_extright != NULL) {
//		delete [] inside_genome_extright;
//		inside_genome_extright = NULL;
//	}

	return onehit;
}

/*
Compared to hmm_search, here trainLoader is totally different from the previous trainLoader
Here filterLoader is the same as the trainLoader in hmm_search.
*/
SearchHit StructureSearch::SearchStructureHitInGenomeSegment(	
								Tree_bag * &	root,
								DynamicP &		dp,
								Stem * &		pStemModel,
								Loop * &		pLoopModel,
								int		winsize,
								int		iFlagSearchReverse,
								char *	train_file,
								int		iFlagAutoFilter,
								char *	filter_file,	//
								int		filterbegin,	//
								int		filterend,		//
								int		genome_num,
								int	*	genome_type,
								int	*	genome_direction,
								int	*	genome_length, 
								char **	genome_name,
								char **	genome_sequence, 
								char **	rc_genome_sequence, 
								double ** genome_baseFreq,
								int *	filterhitbeginpos,
								int *	genome_extleft,
								float	filtertime,
								int		runmode,
								int		searchtype,
								int		streamline_flag,
								int	&	pre_hit_num,		//
								float & struct_search_time,	//
								int		flag_filterhit_based_wholestruct_search,	//20081029
								ofstream &outfile,
								int		max_empty_stem_num
								)
{

	SearchHit onehit;

	dp.setPreHitNum(pre_hit_num);	//

	int minsize = 0;
	float total_searching_time = 0.0;
	for(int iterate=0; iterate<genome_num; iterate++)
	{
		if(genome_type[iterate] != GENOME_NONE)
		{
			//
			if(dp.getFlagPrintDebugDPWinInfo()) {
				cout<<genome_baseFreq[iterate][0]<<" | "
					<<genome_baseFreq[iterate][1]<<" | "
					<<genome_baseFreq[iterate][2]<<" | "
					<<genome_baseFreq[iterate][3]<<endl;
			}
/*
			//support no-scanning 20081029
			if(flag_filterhit_based_wholestruct_search == SEARCH_FILTERHIT_BASED_WHOLESTRUCTURE)
			{
				winsize = genome_length[iterate];
				dp.setFlagFilterSearchBased(TRUE);
				dp.setFilterHitBeginPos(filterhitbeginpos[iterate]);
			} else {
				winsize = genome_length[iterate];
				dp.setFlagFilterSearchBased(TRUE);
				dp.setFilterHitBeginPos(filterhitbeginpos[iterate]);
			}
*/
			//do local structure alignment zbhuang 20090723
			if(flag_filterhit_based_wholestruct_search != SEARCH_LOCALSTRUCTUREALIGNMENT) {
				winsize = genome_length[iterate];
				dp.setFlagFilterSearchBased(TRUE);
				dp.setFilterHitBeginPos(filterhitbeginpos[iterate]);
				dp.setFlagLocalStructureAlign(FALSE);
			} else {
				winsize = genome_length[iterate];
				dp.setFlagFilterSearchBased(FALSE);
				dp.setFlagLocalStructureAlign(TRUE);
				cout<<"Local Structure Alignment"<<endl;
			}

			//load the random/genome sequence file
			dp.setGenomeSequence(genome_length[iterate], genome_sequence[iterate]);

			//load the reverse genome sequence.
			if(genome_type[iterate] == ORIGINAL_GENOME)
			{
				if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_BOTH) {
					dp.setFlagSearchReverseGenome(iFlagSearchReverse);
					dp.setRGenomeSequence(genome_length[iterate], rc_genome_sequence[iterate]);
				}
			} else {
				//GENOME_SEGMENT. filter search result
				dp.setFlagSearchReverseGenome(iFlagSearchReverse);
			}

			//20100121
			dp.initHitObject();

			dp.searchPK(winsize, 
						minsize, 
						pStemModel, 
						pLoopModel, 
						root, 
						genome_extleft[iterate], 
						genome_name[iterate], 
						searchtype, 
						genome_baseFreq[iterate],
						outfile);

			onehit = dp.getHitObject();

			total_searching_time += dp.getCpuTimeAll();
		}//if(genome_type[iterate] != GENOME_NONE)
	}//for iterate

	total_searching_time += filtertime;
	struct_search_time = total_searching_time;	//
	
	pre_hit_num += dp.getNumHit();

	if(!streamline_flag)
	{
		if(searchtype == SEARCH_WHOLESTRUCTURE) {	//20080813
			cout<<endl;
			printAdditionalInfo();
			printTotalHitNum(pre_hit_num);	//20080915
			printTotalTimeInfo(total_searching_time);

			printSearchEndingTimeInfo();
			printOneStarLine();
		}
	}

	return onehit;
}

void StructureSearch::preprocessTDDP(DataLoader & trainLoader, Tree_bag * & overall_root, DynamicP & dp)
{
	int numOfTrainSeq = trainLoader.getNumOfTrainSeq();
	int length = trainLoader.getSeqLength();
	int ** trainingSeqs = trainLoader.getPointerOfSeqs();

	dp.setTrainingSeqs(numOfTrainSeq, length, trainingSeqs);
	if(dp.getFlagPrintDebugInputInfo())
		dp.printTrainingSeqs();
					
	dp.setStems(trainLoader.getNumstem(), trainLoader.getStemarray());
	if(dp.getFlagPrintDebugInputInfo())
		dp.printStems();

	dp.setLoops(trainLoader.getNumloop(), trainLoader.getLooparray());
	if(dp.getFlagPrintDebugInputInfo())
		dp.printLoops();

	//Tree Decomposition Generator
	dp.build_tree();
	if(dp.getFlagPrintDebugInputInfo())
		dp.print_tree();

	dp.identify_loopneighbor();

	//20100120 debug. zbhuang
	//dp.printLoopNeighbor();
	
	dp.build_node_mapping();

	//add the pParent in every tree node.
	overall_root = dp.getMyTDRoot();
	overall_root->pParent = overall_root;
	dp.setTreeBagParent(overall_root);

	dp.buildStemIdxArray();
	dp.buildStemIdxFullArray();
	if(dp.getFlagPrintDebugInputInfo())
		dp.printStemIdxArray();

//	cout<<"----------------"<<endl;
//	dp.printTreebags(overall_root);

	dp.arrangeNodelist(overall_root);	
//	cout<<"----------------"<<endl;
//	dp.printTreebags(overall_root);

	if(dp.getFlagPrintDebugTreeInfo()) {
		printf("<tree decomposition>\n");
		dp.printTreebags(overall_root);
		printf("</tree decomposition>\n");
		exit(0);
	}

	if(dp.getFlagPrintDebugTDInfo())
		dp.printStemIdxArray();
	dp.allocateStemCandIdxArray();

	dp.preprocess_tree(overall_root);

	//dp.print_stem_image(overall_root);	//just for debug 20090211

	dp.buildTreeNodePath(overall_root);
}

void StructureSearch::search_structure(
								Tree_bag * &	wRoot,
								DynamicP &		wDP,
								Stem * &		wPStemModel,
								Loop * &		wPLoopModel,
								int				wWinsize,
								Tree_bag * &	sRoot,
								DynamicP &		sDP,
								Stem * &		sPStemModel,
								Loop * &		sPLoopModel,
								int				sWinsize,
								DataLoader &	trainLoader,
								int		iFlagSearchReverse,
								char *	train_file,
								int		iFlagAutoFilter,
								char *	filter_file,	//
								int		filterbegin,	//
								int		filterend,		//
								int		genome_num,
								int	*	genome_length,
								char **	genome_name,
								char **	genome_sequence,
								char **	rc_genome_sequence,
								double **genome_baseFreq,
								int		runmode,
								int		searchtype,	//supporting substructure-filter search 20090630
								int		iFlagPrintFilterInfo,	//
								char *	output_file,
								int		max_empty_stem_num
								)
{
	int		iterate = 0;

	float	total_cpu_time = 0.0;
	float	structure_search_time = 0.0;
	int		search_count = 0;

	//support saving output into files 20081121
	ofstream outfile;
	outfile.open(output_file, ios::app | ios::out | ios::in);

	for(iterate=0; iterate<genome_num; iterate++)
	{
		//support hit sorting 20081121
		vector<SearchHit> vHit;
		//count the hit with plus direction and minus direction to merge the hit. zbhuang 20090630
		int CountPlus	= 0;
		int CountMinus	= 0;

		//nc search
		structure_search_time = 0.0;
		this->search_structure_oneway(
							wRoot,
							wDP,
							wPStemModel,
							wPLoopModel,
							wWinsize,
							sRoot,
							sDP,
							sPStemModel,
							sPLoopModel,
							sWinsize,
							trainLoader,
							GENOME_SEARCH_DIRECTION_PLUS,
							train_file,
							iFlagAutoFilter,
							filter_file,	//
							filterbegin,	//
							filterend,		//
							genome_num,
							genome_length[iterate],
							genome_name[iterate],
							genome_sequence[iterate],
							genome_baseFreq[iterate][0],	//A
							genome_baseFreq[iterate][1],	//C
							genome_baseFreq[iterate][2],	//G
							genome_baseFreq[iterate][3],	//U
							runmode,
							searchtype,	//supporting substructure-filter search 20090630
							iFlagPrintFilterInfo,	//
							structure_search_time,
							search_count,
							CountPlus,
							CountMinus,
							vHit,
							outfile,
							max_empty_stem_num
							);
		total_cpu_time += structure_search_time;

		if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_BOTH)
		{
			structure_search_time = 0.0;
			this->search_structure_oneway(
								wRoot,
								wDP,
								wPStemModel,
								wPLoopModel,
								wWinsize,
								sRoot,
								sDP,
								sPStemModel,
								sPLoopModel,
								sWinsize,
								trainLoader,
								GENOME_SEARCH_DIRECTION_MINUS,	//iFlagSearchReverse,
								train_file,
								iFlagAutoFilter,
								filter_file,	//
								filterbegin,	//
								filterend,		//
								genome_num,
								genome_length[iterate],
								genome_name[iterate],
								rc_genome_sequence[iterate],
								genome_baseFreq[iterate][3],	//U
								genome_baseFreq[iterate][2],	//G
								genome_baseFreq[iterate][1],	//C
								genome_baseFreq[iterate][0],	//A
								runmode,
								searchtype,	//supporting substructure-filter search 20090630
								iFlagPrintFilterInfo,	//
								structure_search_time,
								search_count,
								CountPlus,
								CountMinus,
								vHit,
								outfile,
								max_empty_stem_num
								);
			total_cpu_time += structure_search_time;
		}
		//support hit sorting 20081121
		sort(vHit.begin(), vHit.end());
		search_count = vHit.size();
		for(int i=0; i<search_count; i++)
		{
			cout<<endl
				<<"Whole structure search hit "<<(i+1)<<endl
				<<"----------------------------"<<endl;
			SearchHit one = vHit[i];
			one.formattedPrint();
		}
	}

	//support saving output into files 20081121
	outfile.close();

	//if(iFlagShowTimeInfo)
	{
		cout<<endl;
		printAdditionalInfo();
		printTotalHitNum(search_count);	//20080915
		printTotalTimeInfo(total_cpu_time);
		printSearchEndingTimeInfo();
		printOneStarLine();
	}
}

void StructureSearch::search_structure_oneway(	
					Tree_bag * &	wRoot,
					DynamicP &		wDP,
					Stem * &		wPStemModel,
					Loop * &		wPLoopModel,
					int				wWinsize,
					Tree_bag * &	sRoot,
					DynamicP &		sDP,
					Stem * &		sPStemModel,
					Loop * &		sPLoopModel,
					int				sWinsize,
					DataLoader &	trainLoader,
					int		iFlagSearchReverse,
					char *	train_file,
					int		iFlagAutoFilter,
					char *	filter_file,	//
					int		filterbegin,	//
					int		filterend,		//
					int		genome_num,
					int	 	genome_length,
					char *	genome_name,
					char *	genome_sequence,
					double	genome_baseFreq_A,
					double	genome_baseFreq_C,
					double	genome_baseFreq_G,
					double	genome_baseFreq_U,
					int		runmode,
					int		searchtype,	//supporting substructure-filter search 20090630
					int		iFlagPrintFilterInfo,	//
					float &	structure_search_time,
					int	 &	pre_hit_num,
					int  &	CountPlus,
					int  &	CountMinus,
					vector<SearchHit> &vHit,
					ofstream &outfile,
					int		max_empty_stem_num
					)
{
	int		i, j;
	clock_t start = 0, end = 0;
	float total_cpu_time = 0.0;
	float one_structure_search_time = 0.0;

	start = clock();


	int winsize = 0;	
	int iStepSize = 0;
	if(searchtype == SEARCH_WHOLESTRUCTURE) {
		winsize = wWinsize;
		iStepSize = wDP.getStepSize();
	} else {
		winsize = sWinsize;
		iStepSize = sDP.getStepSize();
	}

	cout<<"winsize="<<winsize<<endl;

	int	*	outside_genome_type;
	outside_genome_type = new int[1];
	outside_genome_type[0] = GENOME_SEGMENT;

	int	*	outside_genome_length;
//	outside_genome_length = new int[1];
//	outside_genome_length[0] = winsize+1;

	char **	outside_genome_name;
	outside_genome_name = new char*[1];
	outside_genome_name[0] = new char[strlen(genome_name)+1];
	strcpy(outside_genome_name[0], genome_name);

	double ** outside_genome_baseFreq;
	outside_genome_baseFreq = new double*[1];
	outside_genome_baseFreq[0] = new double[4];
	outside_genome_baseFreq[0][0] = genome_baseFreq_A;
	outside_genome_baseFreq[0][1] = genome_baseFreq_C;
	outside_genome_baseFreq[0][2] = genome_baseFreq_G;
	outside_genome_baseFreq[0][3] = genome_baseFreq_U;

	int * outside_filterhitbeginpos;
	outside_filterhitbeginpos = new int[1];
	outside_filterhitbeginpos[0] = ZERO;

	if(winsize > genome_length)
	{
		//genome_length is shorter than sliding window size, then
		//do local structure alignment zbhuang 20090723
		wDP.setFlagFilterSearchBased(FALSE);

		end = clock();
		float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
		cpu_time = cpu_time/3600;
		total_cpu_time += cpu_time;

		char **	outside_genome_sequence;
		outside_genome_sequence = new char*[1];
		outside_genome_sequence[0] = new char[genome_length+1];
		strcpy(outside_genome_sequence[0], genome_sequence);

		int * outside_genome_extleft;
		outside_genome_extleft = new int[1];
		outside_genome_extleft[0] = ZERO;

		outside_genome_length = new int[1];
		outside_genome_length[0] = genome_length+1;

		one_structure_search_time = 0.0;
		SearchHit oneHit_outside;
		//here the whole structure search should be executed. 20080803
		if(searchtype == SEARCH_WHOLESTRUCTURE) {
			oneHit_outside = this->SearchStructureHitInGenomeSegment(
								wRoot,
								wDP,
								wPStemModel,
								wPLoopModel,
								winsize,
								iFlagSearchReverse,
								train_file,
								iFlagAutoFilter,
								filter_file,	//
								filterbegin,	//
								filterend,		//
								ONE,	//genome_num
								outside_genome_type,
								NULL,
								outside_genome_length,
								outside_genome_name,
								outside_genome_sequence,
								NULL,
								outside_genome_baseFreq,
								outside_filterhitbeginpos,	//filterhitbeginpos[]
								outside_genome_extleft,	//ext_left[]
								0.0,	//filtertime
								runmode,	//runmode,
								searchtype,
								TRUE,				//streamline_flag
								pre_hit_num,
								one_structure_search_time,	//
								SEARCH_LOCALSTRUCTUREALIGNMENT,	//do local structure alignment zbhuang 20090723
								outfile,
								max_empty_stem_num
								);
			if(oneHit_outside.getPosStart() != NEGATIVE_ONE)
				this->mergeHits(oneHit_outside, vHit, CountPlus, CountMinus);
		}

		total_cpu_time += one_structure_search_time;
		start = clock();

		//release the memory
		if(outside_genome_extleft != NULL) {
			delete [] outside_genome_extleft;
			outside_genome_extleft = NULL;
		}

		if(outside_genome_sequence != NULL) {
			delete [] outside_genome_sequence[0];
			outside_genome_sequence[0] = NULL;
			delete [] outside_genome_sequence;
			outside_genome_sequence = NULL;
		}
	}
	else	//genome_length is longer than sliding window size, then do the scanning structure search
	{
		outside_genome_length = new int[1];
		outside_genome_length[0] = winsize+1;

		char *pStr = new char[winsize+1];
		//pStr = new char[winsize+1];

//just for the debug convenience
//int i_start = 0, i_end = 0;
//i_start = 133;
//i_end	= 133;
//for(i=i_start; i<=i_end; i++)
		for(i=0; i<=genome_length-winsize; i+=iStepSize)
		{
			//easy for debugging
			outfile<<endl<<"------------------- pos="<<i<<" -------------------"<<endl;

			//1. set the search sequence
			for(j=0; j<winsize; j++)
				pStr[j] = genome_sequence[i+j];
			pStr[j] = '\0';

			end = clock();
			float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
			cpu_time = cpu_time/3600;
			total_cpu_time += cpu_time;

			char **	outside_genome_sequence;
			outside_genome_sequence = new char*[1];
			outside_genome_sequence[0] = new char[winsize+1];
			strcpy(outside_genome_sequence[0], pStr);

			int * outside_genome_extleft;
			outside_genome_extleft = new int[1];
			outside_genome_extleft[0] = i;

			one_structure_search_time = 0.0;
			SearchHit oneHit_outside;
			//here the whole structure search should be executed. 20080803
			if(searchtype == SEARCH_WHOLESTRUCTURE)
				oneHit_outside = this->SearchStructureHitInGenomeSegment(
								wRoot,
								wDP,
								wPStemModel,
								wPLoopModel,
								winsize,
								iFlagSearchReverse,
								train_file,
								iFlagAutoFilter,
								filter_file,	//
								filterbegin,	//
								filterend,		//
								ONE,	//genome_num
								outside_genome_type,
								NULL,
								outside_genome_length,
								outside_genome_name,
								outside_genome_sequence,
								NULL,
								outside_genome_baseFreq,
								outside_filterhitbeginpos,	//filterhitbeginpos[]
								outside_genome_extleft,	//ext_left[]
								//ext_right[],
								0.0,	//filtertime
								runmode,	//runmode,
								searchtype,
								TRUE,				//streamline_flag
								pre_hit_num,
								one_structure_search_time,	//
								SEARCH_WHOLESTRUCTURE,	//20081029
								outfile,
								max_empty_stem_num
								);
			else
				oneHit_outside = this->SearchStructureHitInGenomeSegment(
								sRoot,
								sDP,
								sPStemModel,
								sPLoopModel,
								winsize,
								iFlagSearchReverse,
								train_file,
								iFlagAutoFilter,
								filter_file,	//
								filterbegin,	//
								filterend,		//
								ONE,	//genome_num
								outside_genome_type,
								NULL,
								outside_genome_length,
								outside_genome_name,
								outside_genome_sequence,
								NULL,
								outside_genome_baseFreq,
								outside_filterhitbeginpos,	//filterhitbeginpos[]
								outside_genome_extleft,	//ext_left[]
								//ext_right[],
								0.0,	//filtertime
								runmode,	//runmode,
								searchtype,
								TRUE,				//streamline_flag
								pre_hit_num,
								one_structure_search_time,	//
								SEARCH_WHOLESTRUCTURE,	//20081029
								outfile,
								max_empty_stem_num
								);

			total_cpu_time += one_structure_search_time;
			start = clock();

			//release the memory
			if(outside_genome_extleft != NULL) {
				delete [] outside_genome_extleft;
				outside_genome_extleft = NULL;
			}

			if(outside_genome_sequence != NULL) {
				delete [] outside_genome_sequence[0];
				outside_genome_sequence[0] = NULL;
				delete [] outside_genome_sequence;
				outside_genome_sequence = NULL;
			}

			if(oneHit_outside.getPosStart() != NEGATIVE_ONE)
			{
				if(searchtype == SEARCH_WHOLESTRUCTURE)
					this->mergeHits(oneHit_outside, vHit, CountPlus, CountMinus);
	/**/
				else if(searchtype == SEARCH_SUBSTRUCTUREFILTER)
				{
					//prepare for the whole structure search
					int extension_left	= 0;
					int extension_right = 0;
					trainLoader.getBothExtLength(extension_left, extension_right);
					//
					int filter_hit_length		= oneHit_outside.getPosEnd()-oneHit_outside.getPosStart()+1;
					int while_search_pos_start	= oneHit_outside.getPosStart() - extension_left;// + 1;
					int while_search_length		= 0;

					int * inside_filterhitbeginpos;
					inside_filterhitbeginpos = new int[1];

					if(while_search_pos_start >= 0) {
						while_search_length		= filter_hit_length+extension_left+extension_right;
						inside_filterhitbeginpos[0] = extension_left;
					} else {
						//while_search_pos_start < 0
						while_search_pos_start = 0;
						while_search_length		= filter_hit_length+oneHit_outside.getPosStart()+extension_right;
						inside_filterhitbeginpos[0] = oneHit_outside.getPosStart();
					}

					int	*	inside_genome_length;
					inside_genome_length = new int[1];
					inside_genome_length[0] = while_search_length;

					char **	inside_genome_sequence;
					inside_genome_sequence = new char*[1];
					inside_genome_sequence[0] = new char[while_search_length+1];
					for(j=0; j<while_search_length; j++)
						inside_genome_sequence[0][j] = genome_sequence[while_search_pos_start+j];
					inside_genome_sequence[0][j] = '\0';

					int * inside_genome_extleft;
					inside_genome_extleft = new int[1];
					inside_genome_extleft[0] = while_search_pos_start;

					end = clock();
					float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
					cpu_time = cpu_time/3600;
					total_cpu_time += cpu_time;

					one_structure_search_time = 0.0;
					SearchHit oneHit_inside = this->SearchStructureHitInGenomeSegment(
								wRoot,
								wDP,
								wPStemModel,
								wPLoopModel,
								wWinsize,
								iFlagSearchReverse,
								train_file,
								iFlagAutoFilter,
								filter_file,	//
								filterbegin,	//
								filterend,		//
								ONE,	//genome_num
								outside_genome_type,
								NULL,
								inside_genome_length,
								outside_genome_name,
								inside_genome_sequence,
								NULL,
								outside_genome_baseFreq,
								inside_filterhitbeginpos,	//filterhitbeginpos[]
								inside_genome_extleft,	//ext_left[]
								//ext_right[],
								0.0,	//filtertime
								runmode,	//runmode,
								SEARCH_WHOLESTRUCTURE,	//SEARCH_SUBSTRUCTUREFILTER,
								TRUE,
								pre_hit_num,
								one_structure_search_time,	//
								SEARCH_FILTERHIT_BASED_WHOLESTRUCTURE,	//zbhuang 20090715
								outfile,
								max_empty_stem_num
								);

					total_cpu_time += one_structure_search_time;
					start = clock();

					if(inside_genome_extleft != NULL) {
						delete [] inside_genome_extleft;
						inside_genome_extleft = NULL;
					}

					if(inside_genome_sequence != NULL) {
						delete [] inside_genome_sequence[0];
						inside_genome_sequence[0] = NULL;
						delete [] inside_genome_sequence;
						inside_genome_sequence = NULL;
					}

					if(inside_filterhitbeginpos != NULL) {
						delete [] inside_filterhitbeginpos;
						inside_filterhitbeginpos = NULL;
					}

					if(inside_genome_length != NULL) {
						delete [] inside_genome_length;
						inside_genome_length = NULL;
					}
					if(oneHit_inside.getPosStart() != NEGATIVE_ONE)
						this->mergeHits(oneHit_inside, vHit, CountPlus, CountMinus);

				}//else if(searchtype == SEARCH_SUBSTRUCTUREFILTER)
			}
		}
		if(pStr != NULL) {
			delete [] pStr;
			pStr = NULL;
		}
	}

	if(outside_filterhitbeginpos != NULL) {
		delete [] outside_filterhitbeginpos;
		outside_filterhitbeginpos = NULL;
	}

	if(outside_genome_baseFreq != NULL) {
		delete [] outside_genome_baseFreq[0];
		outside_genome_baseFreq[0] = NULL;
		delete [] outside_genome_baseFreq;
		outside_genome_baseFreq = NULL;
	}

	if(outside_genome_name != NULL) {
		delete [] outside_genome_name[0];
		outside_genome_name[0] = NULL;
		delete [] outside_genome_name;
		outside_genome_name = NULL;
	}

	if(outside_genome_length != NULL) {
		delete [] outside_genome_length;
		outside_genome_length = NULL;
	}

	if(outside_genome_type != NULL) {
		delete [] outside_genome_type;
		outside_genome_type = NULL;
	}

	//if(iFlagShowTimeInfo)
	{
		end = clock();
		float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
		cpu_time = cpu_time/3600;
		total_cpu_time += cpu_time;
		structure_search_time += total_cpu_time;
	}
	
}

int StructureSearch::mergeHits(SearchHit oneHit_outside, vector<SearchHit> & vHit, int & cntPlus, int & cntMinus)
{
	//vector<SearchHit>::iterator iter;
	if(oneHit_outside.getPosStart() != NEGATIVE_ONE)
	{
		if(cntPlus == 0 && oneHit_outside.getDirection() == GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS)
		{
			cntPlus = 1;
			vHit.push_back(oneHit_outside);
		} else if(cntMinus == 0 && oneHit_outside.getDirection() == GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS) {
			cntMinus = 1;
			vHit.push_back(oneHit_outside);
		} else {
			bool bReplace = false;
			bool bAdd = true;
		
			vector<SearchHit>::iterator iter;
			for(iter=vHit.begin(); iter<vHit.end() && !bReplace; iter++)
			{
				SearchHit tmpHit = *iter;
				if( tmpHit.getDirection() == oneHit_outside.getDirection() )
				{
					if((tmpHit.getPosStart() <= oneHit_outside.getPosStart() && oneHit_outside.getPosEnd() <= tmpHit.getPosEnd())
					|| (oneHit_outside.getPosStart() <= tmpHit.getPosStart() && tmpHit.getPosEnd() <= oneHit_outside.getPosEnd())
					|| (tmpHit.getPosStart() <= oneHit_outside.getPosStart() && oneHit_outside.getPosStart() <= tmpHit.getPosEnd() && tmpHit.getPosEnd() <= oneHit_outside.getPosEnd())
					|| (oneHit_outside.getPosStart() <= tmpHit.getPosStart() && tmpHit.getPosStart() <= oneHit_outside.getPosEnd() && oneHit_outside.getPosEnd() <= tmpHit.getPosEnd())
					)
					{
						if(tmpHit.getAlignScore() < oneHit_outside.getAlignScore())
						{
							bReplace = true;
						} else {
							bAdd = false;
						}
					}
				}
			}
			if(bReplace)
			{
				if(vHit.size() == ONE) {
					vHit.pop_back();
				} else {
					vHit.erase(--iter);
				}

				vHit.push_back(oneHit_outside);
			} else if(bAdd) {
				vHit.push_back(oneHit_outside);
			}
		}
	}
	return 1;
}

void StructureSearch::search(int	top_k,
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
							int		iFlagSearchReverse,
							int		iFlagPrintStrAlignInfo,
							int		iFlagPrintScoreInfo,
							int		iFlagPrintDebugInfo,
							int		iFlagPrintDebugFilterInfo,
							int		iFlagPrintDebugInputInfo,
							int		iFlagPrintDebugTreeInfo,
							int		iFlagPrintDebugDPWinInfo,
							int		iFlagPrintDebugTDInfo,
							int		iFlagPrintDebugDPSearchInfo,
							int		iSplitThreshold,
							char *	train_file,
							int		iFlagAutoFilter,
							char *	filter_file,	//
							int		filterbegin,	//
							int		filterend,		//
							double	pdcnt,			//pseudocount
							double	priorh,
							char *	prior_file,
							int		genome_num,
							int	*	genome_type,
							int	*	genome_direction,
							int	*	genome_length, 
							char **	genome_name,
							char **	genome_sequence, 
							char **	rc_genome_sequence, 
							double ** genome_baseFreq,
							int *	genome_extleft,
							int *	genome_extright,
							float	filtertime,
							int		runmode,
							int		searchtype,
							char *	output_file,
							int		max_empty_stem_num,
							int		s_option,		//score option zbhuang 20090304
							double	weight_stem,	//weight factor for stem zbhuang 20090304
							double	weight_loop		//weight factor for loop zbhuang 20090304
							)
{
	//call tree decomposition once. zbhuang 20090720
	DataLoader trainLoader;
	trainLoader.setSplitThreshold(iSplitThreshold);

	//if(searchtype == SEARCH_HMMFILTER)
	if(searchtype != SEARCH_WHOLESTRUCTURE)
		trainLoader.setFilterEndPosition(filterbegin, filterend, -1);

	trainLoader.inputData(train_file);

//	cout<<genome_baseFreq[0][0]<<" | "<<genome_baseFreq[0][1]<<" | "
//		<<genome_baseFreq[0][2]<<" | "<<genome_baseFreq[0][3]<<endl;

	trainLoader.setBaseFreq(genome_baseFreq[0][0],
							genome_baseFreq[0][1],
							genome_baseFreq[0][2],
							genome_baseFreq[0][3]);

	Tree_bag *	wRoot;
	DynamicP	wDP;
	wDP.setSearchParams(top_k,					//[k value]
						threshold,				//[threshold]
						num_nts_overlap_between_stem,	//[number of nts in stem when overlap is allowed]
						iFlagMergeCandInPreprocess,		//whether taking the merge-candidate strategy in preprocessing
						iFlagCandwithShortestLength,
						iShiftNumMergeCand,		//[num of shift when merge candidate in preprocessing]
						iAllowedNullLoopInsNum,	//[number of insertion allowed in null loop]
						pcoeff,					//pcoeff
						iJumpStrategy,			//[skip-and-jump strategy]
						iStepSize,				//[stepsize in skip-and-jump strategy]
						dScoreThresholdInJump,	//[score_filtering threshold in jump strategy]
						iFlagPrintStrAlignInfo,	//[structure alignment info -s|-n]
						iFlagPrintScoreInfo,	//
						iFlagPrintDebugInfo,	//[debug info -n|-d]
						iFlagPrintDebugInputInfo,
						iFlagPrintDebugTreeInfo,
						iFlagPrintDebugDPWinInfo,
						iFlagPrintDebugTDInfo,
						iFlagPrintDebugDPSearchInfo,
						trainLoader.getPastaLines(),
						s_option,					
						weight_stem,				
						weight_loop,				
						max_empty_stem_num
						);
	wDP.setPreHitNum(ZERO);
	this->preprocessTDDP(trainLoader, wRoot, wDP);

	//sliding window size
	int wWinsize = trainLoader.getScanWinLength(3, '-') * 2;

	ScfgBuilder	wSCFGBuilder;
	HMMBuilder	wLOOPBuilder;

	if(searchtype == SEARCH_HMMFILTER) {
		int numsm=0,  numlp=0;

		//wSCFGBuilder
		wSCFGBuilder.setPseudocnt(pdcnt, 0.001);
		wSCFGBuilder.setPriorCoef(priorh);			//coefficient of prior matrix
		wSCFGBuilder.setPriorFilePath(prior_file);	
		wSCFGBuilder.buildSCFG(trainLoader);
		Stem *	wPStemModel = wSCFGBuilder.getAllStems(numsm);

		//wLOOPBuilder
		wLOOPBuilder.setPseudocnt(pdcnt, 0.001);
		wLOOPBuilder.buildHMM(trainLoader);
		Loop *	wPLoopModel = wLOOPBuilder.getAllLoops(numlp);

		this->hmm_search(
						wRoot, 
						wDP,
						wPStemModel,
						wPLoopModel,
						wWinsize,
						trainLoader,
						iFlagSearchReverse,
						train_file,
						iFlagAutoFilter,
						filter_file,	//
						filterbegin,	//
						filterend,		//
						genome_num,
						genome_length,
						genome_name,
						genome_sequence,
						rc_genome_sequence,
						genome_baseFreq,
						runmode,
						iFlagPrintDebugFilterInfo,	//
						output_file,
						max_empty_stem_num
						);

	} else {
		int pre_hit_num = 0;
		float time = 0.0;

		//support saving output into files 20081121
		ofstream outfile;
		outfile.open(output_file, ios::app | ios::out | ios::in);

		int numsm=0,  numlp=0;
		
		//model building for the whole structure search
		//wSCFGBuilder;
		wSCFGBuilder.setPseudocnt(pdcnt, 0.001);
		wSCFGBuilder.setPriorCoef(priorh);			//coefficient of prior matrix
		wSCFGBuilder.setPriorFilePath(prior_file);	
		wSCFGBuilder.buildSCFG(trainLoader);
		Stem *	wPStemModel = wSCFGBuilder.getAllStems(numsm);

		//wLOOPBuilder;
		wLOOPBuilder.setPseudocnt(pdcnt, 0.001);
		wLOOPBuilder.buildHMM(trainLoader);
		Loop *	wPLoopModel = wLOOPBuilder.getAllLoops(numlp);

		//model building for the sub structure search
		DataLoader filterLoader;
		ScfgBuilder	sSCFGBuilder;
		HMMBuilder	sLOOPBuilder;
		Tree_bag *	sRoot	= NULL;
		DynamicP	sDP;
		int sWinsize;
		Stem *	sPStemModel = NULL;
		Loop *	sPLoopModel = NULL;

		if(searchtype == SEARCH_SUBSTRUCTUREFILTER) {
			filterLoader.setSplitThreshold(iSplitThreshold);
			filterLoader.inputData(filter_file);
			filterLoader.setBaseFreq(genome_baseFreq[0][0],
									genome_baseFreq[0][1],
									genome_baseFreq[0][2],
									genome_baseFreq[0][3]);

			sDP.setSearchParams(top_k,					//[k value]
								threshold,				//[threshold]
								num_nts_overlap_between_stem,	//[number of nts in stem when overlap is allowed]
								iFlagMergeCandInPreprocess,		//whether taking the merge-candidate strategy in preprocessing
								iFlagCandwithShortestLength,
								iShiftNumMergeCand,		//[num of shift when merge candidate in preprocessing]
								iAllowedNullLoopInsNum,	//[number of insertion allowed in null loop]
								pcoeff,					//pcoeff
								iJumpStrategy,			//[skip-and-jump strategy]
								iStepSize,				//[stepsize in skip-and-jump strategy]
								dScoreThresholdInJump,	//[score_filtering threshold in jump strategy]
								iFlagPrintStrAlignInfo,	//[structure alignment info -s|-n]
								iFlagPrintScoreInfo,	//
								iFlagPrintDebugInfo,	//[debug info -n|-d]
								iFlagPrintDebugInputInfo,
								iFlagPrintDebugTreeInfo,
								iFlagPrintDebugDPWinInfo,
								iFlagPrintDebugTDInfo,
								iFlagPrintDebugDPSearchInfo,
								filterLoader.getPastaLines(),
								s_option,					
								weight_stem,				
								weight_loop,				
								ZERO	//max_empty_stem_num
								);
			sDP.setPreHitNum(ZERO);
			this->preprocessTDDP(filterLoader, sRoot, sDP);

			//sliding window size
			sWinsize = filterLoader.getScanWinLength(3, '-') * 2;

			int snumsm=0,  snumlp=0;
		
			//model building for the sub-structure search
			//sSCFGBuilder;
			sSCFGBuilder.setPseudocnt(pdcnt, 0.001);
			sSCFGBuilder.setPriorCoef(priorh);			//coefficient of prior matrix
			sSCFGBuilder.setPriorFilePath(prior_file);	
			sSCFGBuilder.buildSCFG(filterLoader);
			sPStemModel = sSCFGBuilder.getAllStems(snumsm);

			//wLOOPBuilder;
			sLOOPBuilder.setPseudocnt(pdcnt, 0.001);
			sLOOPBuilder.buildHMM(filterLoader);
			sPLoopModel = sLOOPBuilder.getAllLoops(snumlp);
		}

		this->search_structure(
						wRoot, 
						wDP,
						wPStemModel,
						wPLoopModel,
						wWinsize,
						sRoot, 
						sDP,
						sPStemModel,
						sPLoopModel,
						sWinsize,
						trainLoader,
						iFlagSearchReverse,
						train_file,
						iFlagAutoFilter,
						filter_file,	//
						filterbegin,	//
						filterend,		//
						genome_num,
						genome_length,
						genome_name,
						genome_sequence,
						rc_genome_sequence,
						genome_baseFreq,
						runmode,
						searchtype,	//supporting substructure-filter search 20090630
						iFlagPrintDebugFilterInfo,	//
						output_file,
						max_empty_stem_num
						);
		if(searchtype == SEARCH_SUBSTRUCTUREFILTER) {
			sLOOPBuilder.freeAllLoops();
			sSCFGBuilder.freeAllStems();

			//detect memory leak problem in rnatops zbhuang 20090815
			//sDP.postprocess_tree(sRoot);
			sDP.postprocessTDDP(sRoot);
		}
	}

	wLOOPBuilder.freeAllLoops();
	wSCFGBuilder.freeAllStems();

	//wDP.postprocess_tree(wRoot);
	wDP.postprocessTDDP(wRoot);
}


