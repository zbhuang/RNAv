// StructureSearch.h: interface for the StructureSearch class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_STRUCTURESEARCH_H__94292247_F0A1_479D_83C9_93D1EADD5AFA__INCLUDED_)
#define AFX_STRUCTURESEARCH_H__94292247_F0A1_479D_83C9_93D1EADD5AFA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "ScfgBuilder.h"
#include "HMMBuilder.h"

#include "SearchOutput.h"

#include "FViterbi.h"
#include "DynamicP.h"
#include "GenomeData.h"
#include "FilterSelection.h"

#include "SearchHit.h"
#include "ConstantsInSearch.h"

#include <algorithm>
#include <vector>

class StructureSearch  
{
private:
	void hmm_search(
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
					int	 *	genome_length,
					char **	genome_name,
					char **	genome_sequence,
					char **	rc_genome_sequence,
					double **genome_baseFreq,
					int		runmode,
					int		iFlagPrintDebugFilterInfo,	//
					char *	output_file,
					int		max_empty_stem_num
					);

	void hmm_search_oneway(	
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
						int	 	genome_length,
						char *	genome_name,
						char *	genome_sequence,
						double	genome_baseFreq_A,
						double	genome_baseFreq_C,
						double	genome_baseFreq_G,
						double	genome_baseFreq_U,
						int		runmode,
						int		iFlagPrintDebugFilterInfo,	//
						float & hmm_search_time,
						int	 &	count,
						vector<SearchHit> &vHit,
						ofstream &outfile,
						int		max_empty_stem_num
						);

	void search_structure(
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
					int	 *	genome_length,
					char **	genome_name,
					char **	genome_sequence,
					char **	rc_genome_sequence,
					double **genome_baseFreq,
					int		runmode,
					int		searchtype,	//supporting substructure-filter search 20090630
					int		iFlagPrintDebugFilterInfo,	//
					char *	output_file,
					int		max_empty_stem_num
					);

	void search_structure_oneway(	
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
						int		iFlagPrintDebugFilterInfo,	//
						float & hmm_search_time,
						int	 &	count,
						int  &	CountPlus,
						int  &	CountMinus,
						vector<SearchHit> &vHit,
						ofstream &outfile,
						int		max_empty_stem_num
						);



	SearchHit structure_search_based_on_one_hmmfilterresult(
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
														int	&	total_hit,
														float & time,
														ofstream &outfile,
														int		max_empty_stem_num
														);

	SearchHit SearchStructureHitInGenomeSegment(
							Tree_bag * &	wRoot,
							DynamicP &		wDP,
							Stem * &		wPStemModel,
							Loop * &		wPLoopModel,
							int		wWinsize,
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
							int	&	pre_hit_num,			//
							float & struct_search_time,		//
							int		flag_filterhit_based_wholestruct_search,	//20081029
							ofstream &outfile,
							int		max_empty_stem_num
							);

	void preprocessTDDP(DataLoader & dataloader, Tree_bag * & root, DynamicP & dp);

	int mergeHits(SearchHit oneHit_outside, vector<SearchHit> &vHit, int & cntPlus, int & cntMinus);
public:
	StructureSearch();
	virtual ~StructureSearch();

	void search(	int		top_k,
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
					double	pdcnt,
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
					int		searchmode,
					int		searchtype,
					char *	output_file,
					int		max_empty_stem_num,
					int		s_option,		//score option zbhuang 20090304
					double	weight_stem,	//weight factor for stem zbhuang 20090304
					double	weight_loop		//weight factor for loop zbhuang 20090304
					);

};

#endif // !defined(AFX_STRUCTURESEARCH_H__94292247_F0A1_479D_83C9_93D1EADD5AFA__INCLUDED_)
