// DynamicP.h: interface for the DynamicP class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning(disable: 4786)

#if !defined(AFX_DYNAMICP_H__EDB6E100_2AE9_4858_9E31_52417CBCFC6D__INCLUDED_)
#define AFX_DYNAMICP_H__EDB6E100_2AE9_4858_9E31_52417CBCFC6D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vector>
#include <map>
#include "TRnaGraphTd.h"

#include "DataLoader.h"
#include "common.h"
#include "ScfgBuilder.h"
#include "HMMBuilder.h"
#include "CandScanner.h"
#include <iostream>
using std::cout;
using std::endl;

#include <string.h>	//gcc3.4

//support SearchHit sorting 20081121
#include "SearchHit.h"

#include "ConstantsInSearch.h"

//A macro used to suggest that a pair of vertices are connected with both directed and nondirected edge.    
typedef struct Stemloc {
	int		a;		//the starting point of the left part of the stem.  
	int		b;		//the ending point of the left part of the stem.
	int		c;		//the starting point of the right part of the stem.
	int		d;		//the ending point of the right part of the stem.
	double	p_score;//the profiling score associated with the stem. 

	//add extend_dir zbhuang
	char	a_ext_dir;
	char	b_ext_dir;
	char	c_ext_dir;
	char	d_ext_dir;
} Stemloc;

typedef struct StructureAlign {
	char *	sequence;
	char *	alignment;
	char *	pastaline_idx;	//support stemid index in two pastalines 20081015
} StructureAlign;

typedef struct CandidateStem {
	int				num;
	Stemloc *		s_locinfo;
	StructureAlign *	leftArm;	//
	StructureAlign *	rightArm;	//
} CandidateStem;

typedef struct StemInfo {
	char	charid;			//character used to label the stem(structure)
	char	charid_idx;		//index of charid, i.e. A1, A2, a1, a2.	20081015
	int		start1;			//The starting location of the left stem in the alignment. 
	int		end1;			//The ending location of the left stem in the alignment.
	int		start2;			//The starting location of the right stem in the alignment.
	int		end2;			//The ending location of the right stem in the alignment.
	int		lg_id;			//The graph node number of the left half of the stem.
	int		rg_id;			//The graph node number of the right half of the stem.
} StemInfo;

typedef struct LoopInfo {
	int		start;			//The starting location of the loop in the concensus structure.
	int		end;			//The ending location of the loop in the concensus structure.
	int		left_stemid;	//The identity of the stem where the left boundary point belongs -1 if it is the source.      
	int		right_stemid;	//The identity of the stem where the right boundary point belongs -1 if it is the sink. 
	int		lg_id;			//The left graph identity for the left end of the loop.     
	int		rg_id;			//The right graph identity for the right end of the loop.
} LoopInfo;

//This part is responsible for generating the graphs needed to be tree decomposed for searching. 
typedef struct C_graphnode {
	int		node_id;  
	int		node_type;  //node_type=0, it is a node of type sink or source. 
						//node_type=1, it is an intermediate node.
	int		mapping_id; //the mapped component id of each node in the graph.      
} C_graphnode; 

typedef struct Strunit {
	int		type;       //type=0, a loop; type=1, a stem.
						//type=1, stem; type=0, half-stem.
	int		g_id1;      //the graph node of the left end of loop(stem).  
	int		g_id2;      //the graph node of the right end of the loop(stem).
	int		loop_id;    //the id of the loop represented by the structure unit.    
	int		stem_id;    //the id of the stem represented by the structure unit.
	int		start1;		//The leftmost location for the left part of the stem.      
	int		end1;		//The rightmost location for the left part of the stem.  
	int		start2;		//The leftmost location for the right part of the stem.
	int		end2;		//The rightmost location for the right part of the stem.
	
	int		left_nid;	//The left node id of the left part of the stem. 
						//left_nid means "index" of this node in this node-list. zbhuang

	int		right_nid;	//The right node id of the right part of the stem.
						//right_nid also means "index" of this node in this node-list. zbhuang

	struct	Strunit * lstem;	//The pointer that points to stem information for the left side of the loop.
	struct	Strunit * rstem;	//The pointer that points to stem information for the right side of the loop.

	struct	Strunit * next;		//the next pointer that points to the next strunit in the linked list. 
} Strunit;

typedef struct Node_list {
	int		g_node;
	int		visited;		//A flag used to construct the stem_list.
    struct	Node_list *next;
} Node_list;

typedef struct New_node_combination {
	double	score;
	int		optimal;					//-1 or 1. 1 means optimal result.
	int		childIndex[MAXNUM_CHILD];	//remove the block memory function in dp part. zbhuang 20090208
	int		overall_empty_stem_num;		//limiting the number of empty-stem. zbhuang 20090211
/**/
	//expend the memory for the dp table when supporting empty-stem-model 20090219
	int	*	lpos;
	int *	rpos;

/**/
	//add extend_dir zbhuang
	//L - left  or  R - right  or  N(initialization value) - not determined
	char *	l_ext_dir;
	char *	r_ext_dir;

} New_node_combination;

//I add lpnum to count the number of su_loop, stnum for su_stem.
typedef struct Tree_bag {
	int		nodenum;		//number of node in this tree_bag.
	int		childnum;		//number of child
	int		stnum;			//number of stems (half stems) in the list.
	int		lpnum;			//number of loops in the list.

	int	 *	stem_image[MAXNUMPT];	//The stem image needed for computing probabilities.
	int	 *	loop_image[MAXNUMPT];	//The loop image needed for computing probabilities.
	int	 *	enum_stem;	//enumerate all possible combinations of structural units(stem) in the tree-node.
	int	 *	enum_node;	//enumerate all possible combinations of elements in the tree-node.

	int	 *	set_up;		//The array really used when the profiling score is extracted.
						//It is based on the state_array while remove some useless nodes.
						//i.e. source/sink node and take one from two nodes forming one stem.

	double	prob_score[MAXNUMPT];	//The probabilities that are shared among the parent node and each of the children.

	Strunit *	su_stem;	//The list of structural units which include both stems and half stems.
	Strunit *	su_loop;	//The list of structural units which include all loops.

	Node_list *	nhead;	//head of tree node
	New_node_combination *	t_head;	//The list of pointers used to store the dynamic programming table results.

	struct	Tree_bag *	pParent;
	struct	Tree_bag *	pChild[MAXNUM_CHILD];
} Tree_bag;

//Declare a structure variable to accomplish the conversion between a graph node and a stem.  
typedef struct Sid_info{
	int		s_id;           //stem id. 
	int		left_or_right;  //whether the node is the left part or the right part.
							//left: 0, right: 1. 
} Sid_info;

class DynamicP  
{
private:
	int		iLeading;
	//training sequences
	int		seqlength;
	int		numtrainingseqs;
	int	**	pTrainingSeqs;

	//stem info
	int			numstems;
	StemInfo *	pStem;

	//loop info
	int			numloops;
	LoopInfo *	pLoop;

	//The graph used to represent the concensus structure.
	int **		pCongraph;

	//every node in the graph belongs to which region. used when calling split_stemcomponents()
	C_graphnode *	node_mapping;

	int		numnode;
	int		nodestart;

	int *		stemstate;

	Sid_info *	gidtosid;

	CandidateStem *	candidateStems;

	//genome sequence
	int		genomeSequenceLength;
	char *	cGenomeSequence;
	//reverse genome sequence
	char *	cRGenomeSequence;

	int		scanLength;
	char *	scanSequence;
	char *	cWindowSequence;
	double *r_scores;	//array of storing the alignment scores.
	int *	final_result_from_array;
	int *	final_result_to_array;

	int		numhit;//numpoint;
	double	threshold;	//
	void	sort_print_FinalResultVariables();

	void	get_stemcomponents(int num, int& stop);	//reference
	int		searchlist(Node_list * n_head, int nodeid);

	//free some memory before preprocess the tree
	void	free_pTrainingSeqs();
	void	free_pStem();
	void	free_pLoop();
	void	free_pCongraph();

	//preprocess_tree
	void	get_conversion_gtos();
	void	print_gidtosid();
	void	allocate_dptable(Tree_bag * root);
	void	construct_stemlist(Tree_bag * root);
	void	tree_stemlist(Tree_bag * root);
	void	construct_looplist(Tree_bag * root);
	void	tree_looplist(Tree_bag * root);
	void	generate_sharedstemlist(Tree_bag * parent, Tree_bag * child, int);  
	void	generate_sharedlooplist(Tree_bag * parent, Tree_bag * child, int);
	void	tree_sharedstemloop(Tree_bag * root);
	void	construct_mappingarray(Tree_bag * parent, Tree_bag * child);   
	void	tree_mappingarray(Tree_bag * root);	

	//postprocess_tree(free memory of sub-unit in the tree node)
	void	free_gtos();
	void	free_dptable(Tree_bag * root);
	Strunit * find_stemtail(Strunit *cur, int index);
	void	free_stem(Tree_bag * root);
	void	free_stemlist(Tree_bag * root);
	Strunit * find_looptail(Strunit *cur, int index);
	void	free_loop(Tree_bag * root);
	void	free_looplist(Tree_bag * root);
	void	free_sharedstemlist(Tree_bag * parent, Tree_bag * child, int c_id);
	void	free_sharedlooplist(Tree_bag * parent, Tree_bag * child, int c_id);
	void	free_sharedstemloop(Tree_bag * root);
	void	deconstruct_mappingarray(Tree_bag * parent, Tree_bag * child);   
	void	free_mappingarray(Tree_bag * root);

	//free memory opposite to the construct_treedecomp(free node->nhead and tree node)
	Node_list * find_treebag_nhead_tail(Node_list *cur, int index);
	void		free_treebag_nhead(Node_list *nodehead, int nodenum);

	//calculate the time.
	clock_t start_all, end_all;
	double	cpu_time_hours_all;
	clock_t start_preprocessing, end_preprocessing;
	double	cpu_time_hours_preprocessing;
	clock_t start_dp, end_dp;
	double	cpu_time_hours_dp;

	void	initDPTable(Tree_bag *root);

    //allow overlap in the adjacent different part of the stems
    int     iStemOverlapNum;

	int		iStemStartPosInStructureLine, iStemEndPosInStructureLine;

	//flag of searching reverse genome sequence
	int		iFlagSearchReverseGenome;

	//flag to print out debug info
	int		flagPrintDebugInfo;
	int		flagPrintDebugInputInfo;
	int		flagPrintDebugTreeInfo;
	int		flagPrintDebugDPWinInfo;
	int		flagPrintDebugTDInfo;
	int		flagPrintDebugDPSearchInfo;

	//adding the empty-stem model. zbhuang 20090211
	int		max_empty_stem_num;
	int		flagSupportEmptyStem;

	//limiting the number of empty-stem. zbhuang 20090213
//	int		empty_stem_num_in_current_treenode;
//	int		emptystemid[MAXNUMEMPTYSTEM];

public:
	DynamicP();
	virtual ~DynamicP();

	void	setLeadingNum(int num);
	int		getLeadingNum();

	//training sequences
	int		getSeqlength();
	int		getNumtrainingseqs();
	void	setSeqlength(int length);
	void	setNumtrainingseqs(int num);
	void	setTrainingSeqs(int length, int num, int ** trainingseqs);
	void	printTrainingSeqs();

	//genome sequence
	void	setGenomeSequence(int ilength, char * cGenomeSequence);
	int		getGenomeSequenceLength();
	char *	getGenomeSequence();
	void	printGenomeSequence();
	void	freeGenomeSequence();

	//reverse genome sequence
	void	setRGenomeSequence(int ilength, char * cRGenomeSequence);
	int		getRGenomeSequenceLength();
	char *	getRGenomeSequence();
	void	printRGenomeSequence();
	void	freeRGenomeSequence();

	//stem and loop info
	int		getNumstems();
	int		getNumloops();
	void	setNumstems(int numstems);
	void	setNumloops(int numloops);
	StemInfo *	getPStem();
	LoopInfo *	getPLoop();
	void	setStems(int numstems, StemLocation *stemarray);
	void	setLoops(int numloops, LoopLocation *looparray);
	void	printStems();
	void	printLoops();

	int		getNumnode();
	int		getNodestart();
	void	build_tree();
	void	print_tree();

	void	printTreenode(Tree_bag *node);
	void	printTreebags(Tree_bag *node);

	void	identify_loopneighbor();
	void	printLoopNeighbor();

	void	build_node_mapping();
	void	free_node_mapping();

	void	allocateStemState(int num);
	void	deallocateStemState();

	void	print_outstructuregraph();

	void	preprocess_tree(Tree_bag * root);
	void	postprocess_tree(Tree_bag * root);

	void	searchPK(int iWinSize, 
					int iMinSize,
					Stem * pStemModel, 
					Loop * pLoopModel, 
					Tree_bag * root, 
					int shift_left, 
					char * genome_name, 
					int searchtype, 
					double * genome_baseFreq,
					ofstream &outfile);

	void	scan_genome(int		iWinSize, 
						int		iMinSize, 
						Stem *	pStemModel, 
						Loop *	pLoopModel, 
						int		type,	//0: random sequence; 1: real genome sequence
						Tree_bag * root, 
						int		genomelength, 
						char *	genome_sequence, 
						int		shift_left, 
						char *	genome_name, 
						int		searchtype, 
						double * genome_baseFreq,
						int		search_direction,
						ofstream &outfile);

	void	compute_profilescore(int i, Tree_bag * root, char * cSequence, Loop *pLoopModel);

	double	getMaxProfilescore(Tree_bag * root);
	double	getMaxProfilescore2(Tree_bag * root);
	double	extract_profilescore(int childIndex, Tree_bag * child, Tree_bag * parent, 
								int * parent_empty_stemid_array, 
								int * parent_emptystemnum,
								int	emptystemnum	//do full combination zbhuang 20090914
								);

	void	tb_indexingfetch(Tree_bag * tr_node, int * index_array, int num_index, int child_empty_stem_num, double * score);//, int * empty_stem_num);
	void	tb_indexingstore(Tree_bag * tr_node, int * index_array, int num_index, double score, int empty_stem_num);
	void	setThreshold(double dThreshold);
	double	getThreshold();


	//just for testing the recursively visiting all treebag algorithm
	void	print_treenode(Tree_bag * node);
	void	print_tree(Tree_bag * root);
	void	printSuStemLoopInfo(Tree_bag * root);
	void	print_treenodelist(Node_list * node_list);
	
	void	free_treenode(Tree_bag * node);
	void	free_tree(Tree_bag * root);

	void	free_scan_genome();

	int		getFlagPrintDebugInfo();
	void	setFlagPrintDebugInfo(int flag);
	int		getFlagPrintDebugInputInfo();
	void	setFlagPrintDebugInputInfo(int flag);
	int		getFlagPrintDebugTreeInfo();
	void	setFlagPrintDebugTreeInfo(int flag);
	int		getFlagPrintDebugDPWinInfo();
	void	setFlagPrintDebugDPWinInfo(int flag);
	int		getFlagPrintDebugTDInfo();
	void	setFlagPrintDebugTDInfo(int flag);
	int		getFlagPrintDebugDPSearchInfo();
	void	setFlagPrintDebugDPSearchInfo(int flag);

	int		getFlagSkipScanRandomSeq();
	void	setFlagSkipScanRandomSeq(int flag);
	int		getFlagCompareTargetPosInfo();
	void	setFlagCompareTargetPosInfo(int flag);

	void	setTreeBagParent(Tree_bag *root);
	void	printTreeBagParent(Tree_bag *root);
	void	removeTreeBagParent(Tree_bag *root);

	void	print_enum_stem(Tree_bag * tr_node);
	void	print_enum_node(int *index_array, int num_index);
	void	print_stem_image(Tree_bag * tr_node);
	void	print_loop_image(Tree_bag * tr_node);

	void	print_treebag_state_array(Tree_bag *root);
	void	print_treebag_set_up(Tree_bag *root);

	void	printTreeBagAllInfo(Tree_bag *root);
	void	printChildNum(Tree_bag *root);

	//20090209
	void	printTraceBackTreeInfo(Tree_bag *root);
	void	printTreeNodeInfo(Tree_bag *root);

	//print the pos for the current candidate
	void	printTreeNodeInfoIncludingPos(Tree_bag *root);

	void	printRootNodeInfo(Tree_bag *root);

	void	printTreeBagChildIdInfo(Tree_bag *root);

	int		findOptimalinRoot(Tree_bag *root);
	void	getOptimalPath(Tree_bag *root);
//	void	printOptimalPath(Tree_bag *root, int currentIndex);

	int		iNumOftheWholeTreeNode;
	int *	pTheWholeTreeNode;
	int		iTheWholeTreeNodeIndex;
	int		countNumOftheWholeTreeNode(Tree_bag *root);
	void	getTheWholeTreeNode(Tree_bag *root);
	int *	pOptimalTheWholeTreeNode;
	int		iOptimalTheWholeTreeNodeIndex;

	void	getOptimalTheWholeTreeNode(Tree_bag *root, int currentIndex);
//	void	getOptimalTheWholeTreeNode(Tree_bag *root, int currentIndex, int empty_stem_num);
	//void	getOptimalTheWholeTreeNode(Tree_bag *tr_node, int * pParentOptimalNode, int optimal_node_num);

	//20090208
	void	printOptimalTheWholeTreeNode();

	//
	int	*	stemPosArray;
	int *	stemIdxArray;	//not including source and sink
	void	buildStemIdxArray();
	void	printStemIdxArray();
	void	freeStemIdxArray();

	//Sort by the nodeid in treebag 20090219
	int *	stemIdxFullArray;			//including source and sink
	void	buildStemIdxFullArray();	//including source and sink
	void	freeStemIdxFullArray();

	//store index of stem candidate in this array
	int *	stemCandIdxArray;
	void	allocateStemCandIdxArray();
	void	buildStemCandIdxArray();
	void	printStemCandIdxArray();
	void	freeStemCandIdxArray();

	//limiting the number of empty-stem. zbhuang
	int		countEmptyStemNumInStemCandIdxArray();

	//
	int		iStemStartPosIdx;
	int		iStemEndPosIdx;
	void	locateStemStartEndPosIdx();

	//deprecated
	int		getIdxInNodeList(Node_list *n_head, int nodeid);

    //allow overlap in the adjacent different part of the stems
    int     iNumOfNtsOverlapBetweenStem;
    void	setNumOfNtsOverlapBetweenStem(int num);
    int		getNumOfNtsOverlapBetweenStem();

	//whether taking the merge-candidate strategy in preprocessing
	int		flagMergeCandInPreprocess;
	int		getFlagMergeCandInPreprocess();
	void	setFlagMergeCandInPreprocess(int pflag);

	//take the sequence with the max score or one with the shortest length 
	//when merge candidate in the preprocessing
	int		flagCandwithShortestLength;
	int		getFlagCandwithShortestLength();
	void	setFlagCandwithShortestLength(int pflag);

	//num of shift when merge candidate in the preprocessing
	int		iShiftNumMergeCand;
	int		getShiftNumMergeCand();
	void	setShiftNumMergeCand(int num);

	//number of insertion allowed in null loop
	int		iAllowedNullLoopInsNum;
	int		getAllowedNullLoopInsNum();
	void	setAllowedNullLoopInsNum(int num);

	TRnaGraphTreeBag<int, Tree_bag>* pTDRoot;
	//top-down
	void	IterateTD_topdown(TRnaGraphTreeBag<int, Tree_bag>* pRoot);
	void	IterateTD_Payload_topdown(TRnaGraphTreeBag<int, Tree_bag>* pRoot);
	//bottom-up
	void	IterateTD_bottomup(TRnaGraphTreeBag<int, Tree_bag>* pRoot);
	void	IterateTD_Payload_bottomup(TRnaGraphTreeBag<int, Tree_bag>* pRoot);
	
	Tree_bag *	mypTDRoot;
	Tree_bag *	copy_treedecomp(TRnaGraphTreeBag<int, Tree_bag>* pRoot);
	Tree_bag *	getMyTDRoot();
	void		freeMyTDRoot();

	//pcoeff
	double	pcoeff;
	double	getPcoeff();
	void	setPcoeff(double _pcoeff);

	//structure alignment info
	int		flagPrintScoreInfo;
	int		getFlagPrintScoreInfo();
	void	setFlagPrintScoreInfo(int pflag);

	//structure alignment info
	int		flagPrintStrAlignInfo;
	int		getFlagPrintStrAlignInfo();
	void	setFlagPrintStrAlignInfo(int pflag);

	//structure alignment info
	char	structureSequence[MAX_STRUCTURE_ALIGN_LENGTH];
	char	strucutreAlignment[MAX_STRUCTURE_ALIGN_LENGTH];
	//support stemid index in two pastalines 20081015
	char	strucutreAlignment_charid_idx[MAX_STRUCTURE_ALIGN_LENGTH];

	//folded structure
	char	foldedstrucutre[MAX_STRUCTURE_ALIGN_LENGTH];
	//support stemid index in two pastalines 20081015
	char	foldedstrucutre_charid_idx[MAX_STRUCTURE_ALIGN_LENGTH];

	//adding the empty-stem model. zbhuang 20090211
	//void	buildStructureAlignInfo(Loop *pLoopModel);
	void	buildStructureAlignInfo(Loop *pLoopModel, int, int);
//	void	printStructureAlignInfo(bool bDebug);
	void	printStructureAlignInfo(int real_start, int real_end, ofstream &outfile);
	void	printStructureAlignInfo2(int real_start, int real_end, ofstream &outfile);

	//for the purpose of debug 20081024
	void	printCandidateStemInfo(ofstream &outfile);
	//re-arrange the debug info module 20081111
	void	printCandidateStemInfo_charid(ofstream &outfile);
	void	printCandidateStemInfo_pos(int index, ofstream &outfile);
	void	printCandidateStemInfo_score(int index, ofstream &outfile);
	void	printHitStemInfo_pos(ofstream &outfile);
	void	printHitStemInfo_score(ofstream &outfile);

	void	printCandidateStemInfo_FailedCase(ofstream &outfile);

	//zbhuang 20100215
	void	printCandidateStemInfo_preprocessingInfo(ofstream &outfile);
	int		checkValidation();

	void	buildFoldedStructureInfo();
//	void	printFoldedStructureInfo(int start, int end, int real_start, int real_end);
	void	printFoldedStructureInfo(int start, int end, int real_start, int real_end, ofstream &outfile);

	char *	buildstemfoldedstructure(int stemidx, int rank, int arm);

	void	printHeader(int startPos, int endPos, bool bDebug);

	int		identifyLoopId(int left_gid, int right_gid);

	void	printScoreForEveryComponent(Loop *pLoopModel, bool bDebug, ofstream &outfile);

	//jump strategy
	int		iJumpStrategy;
	int		getJumpStrategy();
	void	setJumpStrategy(int strategy);
	//stepsize in jump strategy
	int		iStepSize;
	int		getStepSize();
	void	setStepSize(int stepsize);

	double	searchCurrentPos(int iWinSize, Stem *pStemModel, 
							 Loop *pLoopModel, int type, Tree_bag * root,
							 int curPos, 
							 int genomelength,
							 ofstream &outfile);
	void	freeMemForCandidateStems();
	double	getCpuTimeAll();
	double	getCpuTimePreprocessing();
	double	getCpuTimeDP();

	double	dScoreThresInJump;
	double	getScoreThresInJump();
	void	setScoreThresInJump(double dScoreThres);

	int		bNoCandidateInCurrentWindow;

	int		getFlagSearchReverseGenome();
	void	setFlagSearchReverseGenome(int flag);

	void	buildTreeNodePath(Tree_bag * root);
	void	printTreeNodePath();	//20090209
	void	freeTreeNodePath();

	void	setSearchParams(int		top_k,
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
							int		iFlagPrintScoreInfo,
							int		iFlagPrintDebugInfo,
							int		iFlagPrintDebugInputInfo,
							int		iFlagPrintDebugTreeInfo,
							int		iFlagPrintDebugDPWinInfo,
							int		iFlagPrintDebugTDInfo,
							int		iFlagPrintDebugDPSearchInfo,
							int		PastaLineNum,	//20090628
							int		score_option,	//20090628
							double	weight_stem,	//20090628
							double	weight_loop,	//20090628
							int		max_empty_stem_num
							);

	int		getNumHit();

	int		pre_hit_num;
	int		getPreHitNum();
	void	setPreHitNum(int num);

	//support stemid index in two pastalines 20081015
	int		num_pastaline;
	int		getNumPastaLine();
	void	setNumPastaLine(int num);

	//flag of whole structure search based on filter search
	int		flag_filtersearch_based;
	int		getFlagFilterSearchBased();
	void	setFlagFilterSearchBased(int bValue);
	int		filterHitBeginPos;
	int		getFilterHitBeginPos();
	void	setFilterHitBeginPos(int begin_pos);

	//score option zbhuang 20090304
	int		score_option;
	int		getScoreOption();
	void	setScoreOption(int option);

	//support weight factor for stem and loop	20081109
	double	stem_weight;
	double	getStemWeight();
	void	setStemWeight(double weight);
	double	loop_weight;
	double	getLoopWeight();
	void	setLoopWeight(double weight);

	//support stub test 20081116
	int		flagStubTest;
	int		getFlagStubTest();
	void	setFlagStubTest(int flag);

	//mark the invalid combination 20081116
	//output only one invalid combination 20081116
	int		flagInvalidCombination;
	int		getFlagInvalidCombination();
	void	setFlagInvalidCombination(int flag);

	//support hit sorting 20081121
	SearchHit	oneHit;
	void		initHitObject();
	SearchHit	getHitObject();

	//adding the empty-stem model. zbhuang 20090211
	void	setMaxEmptyStemNum(int num);
	int		getMaxEmptyStemNum();

	int		getFlagSupportEmptyStem();
	void	setFlagSupportEmptyStem(int flag);

	int		getStemIdbyLeftNodeId(int nodeid);

    //check and remove the position confliction 20081209
	//also need to minus the score of removed stem
	void	printStemCandIdxPos();
	void	checkStemCandIdxConfliction(double & score);
	int		getIdxInCandIdxArray(int gid);

	//Sort by the nodeid in treebag 20090219
	void	arrangeNodelist(Tree_bag * root);

	void	restoreEnumNodeCombination(Tree_bag * tr_node);
	void	transferEmptyStemPosToParent(Tree_bag * parent, int * parent_empty_stemid_array, 
										Tree_bag * child, int optimallChildIndex);
	int		getNodeIndexInTreenode(Node_list *n_head, int nodeid);
	void	printNumInBase(int total, int base, int num, int overall_empty_stem_num, bool bRoot);

	void	checkSingleEmptyStemPos(Tree_bag * tr_node, double * score);
	void	computeNeighborEmptyStemPos(Tree_bag * tr_node, double * score);
	void	checkNeighborEmptyStemPos(Tree_bag * tr_node, double * score);
	//20090402 for debug
	void	printNeighborStemPos(Tree_bag * tr_node, double score);
	//
	int		_parentIndex[MAXNUM_CHILD];
	int		_childOptimalIndex[MAXNUM_CHILD];

	int		bConflict;

	void	postprocessTDDP(Tree_bag * root);

	//do local structure alignment zbhuang 20090723
	int		flagLocalStructureAlign;
	int		getFlagLocalStructureAlign();
	void	setFlagLocalStructureAlign(int flag);
};

#endif // !defined(AFX_DYNAMICP_H__EDB6E100_2AE9_4858_9E31_52417CBCFC6D__INCLUDED_)
