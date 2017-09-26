// DynamicP.cpp: implementation of the DynamicP class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning( disable: 4786 )

#include "DynamicP.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DynamicP::DynamicP()
{
	pStem			= NULL;
	pLoop			= NULL;
	pCongraph		= NULL;
	pTrainingSeqs	= NULL;
	node_mapping	= NULL;
	//
	stemstate		= NULL;

	seqlength		= 0;
	numtrainingseqs = 0;

	numnode			= 0;
	nodestart		= 0;

	genomeSequenceLength = 0;
	cGenomeSequence = NULL;
	cRGenomeSequence = NULL;

	scanLength = 0;
	scanSequence = NULL;

	candidateStems = NULL;

	r_scores = NULL;

	this->numhit = 0;
	this->pre_hit_num = 0;	//

	iTheWholeTreeNodeIndex			= 0;
	iOptimalTheWholeTreeNodeIndex	= 0;
	iFlagSearchReverseGenome		= 0;

	flagPrintDebugInfo			= 0;
	flagPrintDebugInputInfo		= 0;
	flagPrintDebugTreeInfo		= 0;
	flagPrintDebugDPWinInfo		= 0;
	flagPrintDebugTDInfo		= 0;
	flagPrintDebugDPSearchInfo	= 0;

	//support no-scanning 20081106
	flag_filtersearch_based = 0;
	filterHitBeginPos = 0;

	//support stub test 20081116
	flagStubTest = FALSE;
}

DynamicP::~DynamicP()
{
	this->free_pTrainingSeqs();

	this->free_pStem();

	this->free_pLoop();

	this->free_pCongraph();

	this->free_node_mapping();

	//free memeory for genome sequence
	this->freeGenomeSequence();
	this->freeRGenomeSequence();

}

void DynamicP::setLeadingNum(int num)
{
	this->iLeading = num;
}

int DynamicP::getLeadingNum()
{
	return this->iLeading;
}

void DynamicP::free_pTrainingSeqs()
{
	int i;
	if(pTrainingSeqs != NULL) {
		for(i=0; i<this->getNumtrainingseqs(); i++) {
			delete [] pTrainingSeqs[i];
			pTrainingSeqs[i] = NULL;
		}
		delete [] pTrainingSeqs;
		pTrainingSeqs = NULL;
	}
}

void DynamicP::free_pStem()
{
	if(pStem != NULL)
		delete [] pStem;
	pStem = NULL;
}
void DynamicP::free_pLoop()
{
	if(pLoop != NULL)
		delete [] pLoop;
	pLoop = NULL;
}
void DynamicP::free_pCongraph()
{
	int i;
	if(pCongraph != NULL) {
		for(i=0; i<this->getNumnode(); i++) {
			delete [] pCongraph[i];
			pCongraph[i] = NULL;
		}
		delete [] pCongraph;
		pCongraph = NULL;
	}
}
void DynamicP::free_node_mapping()
{
	if(node_mapping != NULL)
		delete [] node_mapping;
	node_mapping = NULL;
}

void DynamicP::freeGenomeSequence()
{
	if(cGenomeSequence != NULL)
		delete [] cGenomeSequence;
	cGenomeSequence = NULL;
}

void DynamicP::freeRGenomeSequence()
{
	if(cRGenomeSequence != NULL)
		delete [] cRGenomeSequence;
	cRGenomeSequence = NULL;
}

int DynamicP::getSeqlength()
{
	return this->seqlength;
}

int DynamicP::getNumtrainingseqs()
{
	return this->numtrainingseqs;
}

void DynamicP::setTrainingSeqs(int numtrainseqs, int trainingseqlength, int ** trainingseqs)
{
	int i, j;

	this->seqlength = trainingseqlength;
	this->numtrainingseqs = numtrainseqs;

	this->pTrainingSeqs = new int* [numtrainseqs];
	if(this->pTrainingSeqs == 0) {
		printf("Fail to allocate memory!\n");
		exit(0);
    }

	for(i=0; i<numtrainseqs; i++)
	{
		this->pTrainingSeqs[i] = new int[trainingseqlength];
		for(j=0; j<trainingseqlength; j++)
		{
			this->pTrainingSeqs[i][j] = trainingseqs[i][j];
		}
	}
}

void DynamicP::printTrainingSeqs()
{
	int i, j;
	int num = this->getNumtrainingseqs();
	int length = this->getSeqlength();

	for(i=0; i<num; i++)
	{
		for(j=0; j<length; j++)
			cout<<numtochar(this->pTrainingSeqs[i][j]);
		cout<<endl;
	}
}

double DynamicP::getThreshold()
{
	return this->threshold;
}
void DynamicP::setThreshold(double dThreshold)
{
	this->threshold = dThreshold;
}
//genome sequence
void DynamicP::setGenomeSequence(int ilength, char * cGenomeSequence)
{
	this->genomeSequenceLength = ilength;
	this->cGenomeSequence = new char[ilength+1];
	strcpy(this->cGenomeSequence, cGenomeSequence);
}
int	DynamicP::getGenomeSequenceLength()
{
	return this->genomeSequenceLength;
}
char * DynamicP::getGenomeSequence()
{
	return this->cGenomeSequence;
}
void DynamicP::printGenomeSequence()
{
	cout<<"printGenomeSequence()"<<endl;
	cout<<this->cGenomeSequence<<endl;
}

//reverse genome sequence
void DynamicP::setRGenomeSequence(int ilength, char * cRGenomeSequence)
{
	this->cRGenomeSequence = new char[ilength+1];
	strcpy(this->cRGenomeSequence, cRGenomeSequence);
}
char * DynamicP::getRGenomeSequence()
{
	return this->cRGenomeSequence;
}
void DynamicP::printRGenomeSequence()
{
	cout<<"printRGenomeSequence()"<<endl;
	cout<<this->cRGenomeSequence<<endl;
}

void DynamicP::setNumstems(int numstems)
{
	this->numstems = numstems;
}

void DynamicP::setNumloops(int numloops)
{
	this->numloops = numloops;
}
int DynamicP::getNumstems()
{
	return this->numstems;
}

int DynamicP::getNumloops()
{
	return this->numloops;
}

StemInfo * DynamicP::getPStem()
{
	return this->pStem;
}

LoopInfo * DynamicP::getPLoop()
{
	return this->pLoop;
}

void DynamicP::setStems(int numstems, StemLocation *stemarray)
{
	int i;

	this->setNumstems(numstems);

	this->pStem = new StemInfo[numstems];

	//20100319
	this->iStemEndPosInStructureLine = -1;

    for(i=0; i<numstems; i++) {
		this->pStem[i].charid	= stemarray[i].charid;	//
		if(this->getNumPastaLine() == TWO_PASTALINE)
		{
			//support stemid index in two pastalines i.e. A1, A2, a1, a2. 20081015
			string stemstrid	= stemarray[i].getStemStrID();
			if(this->pStem[i].charid == stemstrid[0])
				this->pStem[i].charid_idx	= stemstrid[1];
			else
				cout<<"Something wrong in parsing stemid and stemstrid"<<endl;
		} else {
			this->pStem[i].charid_idx	= NULL;
		}
		this->pStem[i].start1	= stemarray[i].Pos[0];
		this->pStem[i].end1		= stemarray[i].Pos[1];
		this->pStem[i].start2	= stemarray[i].Pos[2];
		this->pStem[i].end2		= stemarray[i].Pos[3];
		if(i == 0) {
			this->iStemStartPosInStructureLine = stemarray[i].Pos[0];
		}
//		if(i == (numstems-1)) {
//			this->iStemEndPosInStructureLine = stemarray[i].Pos[3];
//		}
		if(this->iStemEndPosInStructureLine < stemarray[i].Pos[3]) {
			this->iStemEndPosInStructureLine = stemarray[i].Pos[3];
		}
	}
}

void DynamicP::printStems()
{
	cout<<"----Stem info----"<<endl;
	int i;
	int num = this->getNumstems();
	cout<<"- num of stem: "<<num<<" -"<<endl;
	StemInfo *	_pStem = this->getPStem();
    for(i=0; i<num; i++) {
		cout<<_pStem[i].start1<<"\t"
			<<_pStem[i].end1<<"\t"
			<<_pStem[i].start2<<"\t"
			<<_pStem[i].end2<<endl;
    }
}

void DynamicP::setLoops(int numloops, LoopLocation *looparray)
{
	int i;
	this->setNumloops(numloops);
	this->pLoop = new LoopInfo[numloops];

	for(i=0; i<numloops; i++) {
		this->pLoop[i].start= looparray[i].Pos[0];
		this->pLoop[i].end	= looparray[i].Pos[1];
	}
}

void DynamicP::printLoops()
{
	cout<<"----Loop info----"<<endl;
	int i;
	int num = this->getNumloops();
	cout<<"- num of loop: "<<num<<" -"<<endl;
	LoopInfo *	_pLoop = this->getPLoop();
    for(i=0; i<num; i++) {
		cout<<_pLoop[i].start<<"\t"
			<<_pLoop[i].end<<endl;
    }
}

int	DynamicP::getNumnode()
{
	return this->numnode;
}

int	DynamicP::getNodestart()
{
	return this->nodestart;
}

//construct the graph to represent the consensus structure.
void DynamicP::build_tree()
{
	std::vector< TRnaGraphNode<int> > v;
	std::map< int, std::pair<int, TRnaGraphNode<int> > > m;
	int i, j;
	int numstem = this->getNumstems();
	int numloop = this->getNumloops();
	int numnode = 2*numstem + 2;
	this->numnode = numnode;

	//allocate memory for pCongraph
	this->pCongraph = new int* [numnode];
    for(i=0; i<numnode; i++) {
		this->pCongraph[i] = new int[numnode];
    }

    for(i=0; i<numnode; i++){
		for(j=0; j<numnode; j++){
			this->pCongraph[i][j] = NONCONNECTED;
		}
    }

    this->nodestart = 1;
    for(i=0; i<numstem; i++)
	{
		//Assign the nodes i+1 and i+2 to the stem i.     
		this->pStem[i].lg_id = nodestart;  
		m.insert(std::pair<int,std::pair<int, TRnaGraphNode<int> > >(pStem[i].start1, std::pair<int, TRnaGraphNode<int> >(pStem[i].start2,TRnaGraphNode<int>(-1,-1,-1,pStem[i].lg_id))));
		this->pStem[i].rg_id = nodestart+1;
		m.insert(std::pair<int,std::pair<int, TRnaGraphNode<int> > >(pStem[i].start2, std::pair<int, TRnaGraphNode<int> >(pStem[i].start1,TRnaGraphNode<int>(-1,-1,-1,pStem[i].rg_id))));
		this->pCongraph[nodestart][nodestart+1] = NONDIRECTED;
		this->pCongraph[nodestart+1][nodestart] = NONDIRECTED;
		this->nodestart += 2;
    }

	i = 0;
	std::map<int,std::pair<int, TRnaGraphNode<int> > >::iterator it = m.begin();
	while( it != m.end() )
	{
		it->second.second.ID(i);
		it->second.second.Next(i+1);
		m.find(it->second.first)->second.second.Pair(i);
		it++;
		i++;
	}

	it = m.begin();
	while( it != m.end() )
	{
		v.push_back( it->second.second );
		it++;
	}
/*
	cout << "<--------------" << endl;
	for( i = 0; i < v.size(); i++ )
	{
		cout << v[i].ID() << ' ' << v[i].Next() << ' ' << v[i].Pair() << ' ' << v[i].Payload() << endl;
	}
	cout << "-------------->" << endl;
*/
	TRnaGraphTd< int, Tree_bag > rgtd( v );
	rgtd.Construct();
//	cout << "###############" << endl;
//	rgtd.Print();
	//get the root of the tree-bag
	this->pTDRoot = rgtd.Decomposition();

//	this->IterateTD(this->pTDRoot);
//	this->IterateTD_Payload_topdown(this->pTDRoot);
//	this->IterateTD_Payload_bottomup(this->pTDRoot);
	mypTDRoot = this->copy_treedecomp(this->pTDRoot);
/*
	//20080630
	if(this->pTDRoot != NULL) {
		delete this->pTDRoot;
		this->pTDRoot = NULL;
	}
*/
	//20080702
	rgtd.Release();

//	this->printTreebags(mypTDRoot);

    //Construct the edges in the graph, go through all the
    //loops in the consensus structure and make connections between nodes.
	int source = -1, sink = -1;
	int loop_start, loop_end, start_found, end_found;

    for(i=0; i<numloop; i++)
	{
		loop_start	= this->pLoop[i].start;
		loop_end	= this->pLoop[i].end;
		start_found	= 0;
		end_found	= 0;

		for(j=0; j<numstem; j++)
		{
			//Go through the list of stems to check for the corresponding stem; 
			if(loop_start == this->pStem[j].end1+1) {
				source		= this->pStem[j].lg_id;
				start_found = 1;
			} else if(loop_start == this->pStem[j].end2+1) {
				source		= this->pStem[j].rg_id;
				start_found = 1;
			}
			if(loop_end == this->pStem[j].start1-1) {
				sink		= pStem[j].lg_id;
				end_found	= 1;
			} else if(loop_end == this->pStem[j].start2-1) {
				sink		= pStem[j].rg_id;
				end_found	= 1;
			}
		}//for j.
		if(start_found == 0) {
			source = 0;
		}

		if(end_found == 0) {
			sink = this->getNodestart();
		}

		this->pCongraph[source][sink] = DIRECTED;
	}//for i.
/*
	//Now start checking if there are contiguous stems without a loop in between.
	int s1_start1, s1_end1, s1_start2, s1_end2;
	int s2_start1, s2_end1, s2_start2, s2_end2;

	for(i=0; i<numstem; i++)
	{
		s1_start1 = this->pStem[i].start1;
		s1_end1 = this->pStem[i].end1;
		s1_start2 = this->pStem[i].start2;
		s1_end2 = this->pStem[i].end2;

		for(j=0; j<numstem; j++)
		{
			s2_start1	= this->pStem[j].start1;
			s2_end1		= this->pStem[j].end1;
			s2_start2	= this->pStem[j].start2;
			s2_end2		= this->pStem[j].end2;

			if(s1_end1+1 == s2_start1)
			{
				//If the end of the stem1 is the starting of the start2.      
				source	= this->pStem[i].lg_id;
				sink	= this->pStem[j].lg_id;
				this->pCongraph[source][sink] = WEIGHTZERO;  //Set the weight of the edge to be zero.
			}
			if(s2_end2+1 == s1_start2)
			{
				//If the end of the right stem2 is the starting of the right half of stem1.
				source	= this->pStem[j].rg_id;
				sink	= this->pStem[i].rg_id;
				this->pCongraph[source][sink] = WEIGHTZERO;  //Set the weight of the edge to be zero.
			}
			if(s2_end1+1 == s1_start2)
			{
				//If the end of the left stem2 is the starting of the right half of stem1. 
				source	= this->pStem[j].lg_id;
				sink	= this->pStem[i].rg_id;
				this->pCongraph[source][sink] = WEIGHTZERO; //Set the weight of the edge to be zero.
			}
			if(s1_end2+1 == s2_start1)
			{
				//If the end of the right stem1 is the starting of the left half of stem2. 
				source	= this->pStem[i].rg_id;
				sink	= this->pStem[j].lg_id;
				this->pCongraph[source][sink] = WEIGHTZERO; //Set the weight of the edge to be zero.
			}//if
		}//for j.
	}//for i.
*/
	//Finally, the source and sink are checked to see if they should be connected to the stems.     
	for(i=0; i<numstem; i++)
	{
//		if(this->pStem[i].start1 == 0)
		if(this->pStem[i].start1 == this->iStemStartPosInStructureLine)
		{
			//If there is a stem starting with the source.   
			this->pCongraph[0][this->pStem[i].lg_id] = WEIGHTZERO;
		}

//		if(this->pStem[i].end2 == this->getSeqlength()-1) {
		if(this->pStem[i].end2 == this->iStemEndPosInStructureLine) {
			//If there is a stem ending with the sink.
			this->pCongraph[this->pStem[i].rg_id][nodestart] = WEIGHTZERO;
		}
	}
}

void DynamicP::print_tree()
{
	cout<<"----Graph Info----"<<endl;
	int i, j;
	int numnode = this->getNumnode();
	for(i=0; i<numnode; i++)
	{
		for(j=0; j<numnode; j++)
			cout<<this->pCongraph[i][j]<<" ";
		cout<<endl;
	}
	cout<<"numnode="<<numnode<<endl;
}

//Find the left stem id and the right stem id (allloops.left_id/right_id/lg_id/rg_id)
//int left_id;	- The identity of the stem where the left boundary point belongs -1 if it is the source.      
//int right_id; - The identity of the stem where the right boundary point belongs -1 if it is the sink. 
//int lg_id;	- The left graph identity for the left end of the loop.     
//int rg_id;	- The right graph identity for the right end of the loop.
void DynamicP::identify_loopneighbor()
{
	int numloop		= this->getNumloops();
	int numstem		= this->getNumstems();
	int nodestart	= this->getNodestart();

	int i, j;
	int l_start, l_end;
	int l_sid, r_sid;
	int l_gid = -1, r_gid = -1;

	for(i=0; i<numloop; i++)
	{
		//Go through every loop in the structure.
		l_start = this->pLoop[i].start;
		l_end	= this->pLoop[i].end;

		l_sid = -1;
		r_sid = -1;

		for(j=0; j<numstem; j++)
		{
			//Go through every stem in the secondary structure to check for the ending points of loops. 
            if(this->pStem[j].end1+1 == l_start || this->pStem[j].end2+1 == l_start)
			{
				l_sid = j;
				if(this->pStem[j].end1+1 == l_start)
					l_gid = this->pStem[j].lg_id;
				else
					l_gid = this->pStem[j].rg_id;
			}
            if(l_end+1 == this->pStem[j].start1 || l_end+1 == this->pStem[j].start2)
			{
				r_sid = j;
				if(this->pStem[j].start1 == l_end+1)
					r_gid = this->pStem[j].lg_id;
				else
					r_gid = this->pStem[j].rg_id;
            }
		}//for j.
		this->pLoop[i].left_stemid	 = l_sid;
		this->pLoop[i].right_stemid = r_sid;
		this->pLoop[i].lg_id	 = l_gid;
		this->pLoop[i].rg_id	 = r_gid;

		if(l_sid == -1) {
			//If the leftmost end of the loop is the source.
			this->pLoop[i].lg_id = 0;
		}
		if(r_sid == -1) {
			//If the rightmost end of the loop is the sink. 
			this->pLoop[i].rg_id = nodestart;
		}
	}//for i.
}

void DynamicP::printLoopNeighbor()
{
	cout<<"----Loop Neighbor Info----"<<endl;
	int numloop = this->getNumloops();

	int i;
	for(i=0; i<numloop; i++) {
        cout<<"this->pLoop["<<i<<"].left_stemid is: "<<this->pLoop[i].left_stemid<<endl;
        cout<<"this->pLoop["<<i<<"].right_stemid is: "<<this->pLoop[i].right_stemid<<endl;
	}

	for(i=0; i<numloop; i++) {
        cout<<"this->pLoop["<<i<<"].lg_id is: "<<this->pLoop[i].lg_id<<endl;
        cout<<"this->pLoop["<<i<<"].rg_id is: "<<this->pLoop[i].rg_id<<endl;
	}
}

//identify the loop id given the left stem id and right stem id
//it is used in the printing out the whole structure alignment info
//including stem and loop
int DynamicP::identifyLoopId(int left_gid, int right_gid)
{
	int		iLoopId = -1;
	int		i;
	int		numloop = this->getNumloops();
	bool	bFound = false;

	for(i=0; i<numloop && !bFound; i++) {
		if((this->pLoop[i].lg_id == left_gid) 
		&& (this->pLoop[i].rg_id == right_gid))
		{
			iLoopId = i;
			bFound	= true;
		}
	}
	return iLoopId;
}
void DynamicP::build_node_mapping()
{
	int nodestart = this->getNodestart();
	node_mapping = new C_graphnode[nodestart+1];

    //Initialize the type and mapping_id of the source and the sink.
    node_mapping[0].node_type = 0;
    node_mapping[0].mapping_id = 0;

	//do the initialization
	if(nodestart>1) {
		for(int i=0; i<nodestart; i++) {
			node_mapping[i].node_type = 1;
		}
	}

    node_mapping[nodestart].node_type = 0;
    node_mapping[nodestart].mapping_id = 0;
}

void DynamicP::allocateStemState(int num)
{
	int i;
	stemstate = new int[num];
	for(i=0; i<num; i++) {
		stemstate[i] = 0;
	}
}

void DynamicP::deallocateStemState()
{
	delete [] stemstate;
	stemstate = NULL;
}

//Compute the number of stems.    
void DynamicP::get_stemcomponents(int num, int& stop) 
{
	int numstem = this->getNumstems();

    int i, j;
    int s1_start1, s1_end1, s1_start2, s1_end2;
    int s2_start1, s2_end1, s2_start2, s2_end2;

    int maxnum, maxstem = 0;

	//Store the number of stems that cross with each given stem in the consensus structure.
    int * numcross = new int[numstem];  
    for(i=0; i<numstem; i++) {
		numcross[i]=0;
    }
 
    for(i=0; i<numstem; i++)
	{
		//Go through all possible stems. 
		if(stemstate[i] == num)
		{
			s1_start1	= this->pStem[i].start1;
			s1_end1		= this->pStem[i].end1;
			s1_start2	= this->pStem[i].start2;
			s1_end2		= this->pStem[i].end2;

			for(j=0; j<numstem; j++)
			{
				//Go through the stems with an integer index greater than i. 
				if(stemstate[j] == num)
				{
					s2_start1	= this->pStem[j].start1; 
					s2_end1		= this->pStem[j].end1;
					s2_start2	= this->pStem[j].start2;
					s2_end2		= this->pStem[j].end2;
					if(s1_start1 < s2_start1 && s1_end2 < s2_end2 && s2_end1 < s1_start2 ) {
						//Increase the number of crossing stems of stem i by 1.
						numcross[i]++;
					} else if(s1_start1 > s2_start1 && s1_end2 > s2_end2 && s1_end1<s2_start2) {
						//Increase the number of crossing stems of stem i by 1.
						numcross[i]++;
					}
				}//if stemstate[j]
			}//for j.
		}//if stemstate[i]
	}//for i.

    //Find the stem with the maximum number of crossing stems. 
    maxnum = 0;

    for(i=0; i<numstem; i++) {
		if(maxnum < numcross[i]) {
			maxnum	= numcross[i];
			maxstem	= i;
        }
	}

	if(maxnum == 0) {
		//Stop the procedure if there is no crossing stems in the array.   
		stop = 1;
	} else {
		stemstate[maxstem] = num+1;
	}

	delete [] numcross;
	numcross = NULL;
}
/*
//Separate the crossing stems into different components those will be considered 
//when the tree decomposition of the graph needs to be constructed.        
int DynamicP::split_stemcomponents()
{
	int numstem = this->getNumstems();

	int		i;
	int		s_id, mappedstem;
	int		stop = 0;

	//this->allocateStemState(numstem);

	s_id = 0;
	mappedstem = 0;

	while(mappedstem < numstem)
	{
		while(stop == 0) {
			this->get_stemcomponents(s_id, stop);
        }

        for(i=0; i<numstem; i++)
		{
			if(stemstate[i] == s_id)
			{
				mappedstem++;
				node_mapping[this->pStem[i].lg_id].node_type	= 1;	//set the node type to be 1.
				node_mapping[this->pStem[i].lg_id].mapping_id	= s_id;
				node_mapping[this->pStem[i].rg_id].node_type	= 1;	//set the node type to be 1.
				node_mapping[this->pStem[i].rg_id].mapping_id	= s_id;
			}//if
		}//for
		s_id++;
		stop=0; ///get stop back to 0
	}//while

	//this->deallocateStemState();
	//s_id is the number of trees we have obtained for the mapping.
	return s_id;
}
*/
//Print out the nodes in one tree bag. 
void DynamicP::printTreenode(Tree_bag *node)
{
    //int i;
    Node_list *nodehead;

	if(node == NULL)
		return;

    nodehead = node->nhead;
    printf("{");
    while(nodehead) {
		//Go through the list of nodes, print out the head information.
		printf("%d ", nodehead->g_node);
		nodehead = nodehead->next;
    }
    printf("}");
}

//Print out the tree_bag, given the root of the tree node. 
void DynamicP::printTreebags(Tree_bag *node)
{
	int i;
	//node_list *nodehead;
	if(node == NULL)
		return;
	printTreenode(node);

	printf(":=");
	for(i=0; i<node->childnum; i++)
		printTreenode(node->pChild[i]);

	printf("\n");

	for(i=0; i<node->childnum; i++) {
		//printf("node->ptnum is: %d\n", node->ptnum);
		printTreebags(node->pChild[i]);
	}
}

//Search for a given node in a linked list. 
//It returns 1 if it is found and returns 0 if it is not.
int DynamicP::searchlist(Node_list *n_head, int nodeid) 
{
	Node_list *nhead;
	nhead = n_head;
	while(nhead) {
		if(nhead->g_node == nodeid) {
			return 1;
		}
		nhead = nhead->next;
	}
	return 0; 
}

void DynamicP::setTreeBagParent(Tree_bag *root) {
	int i;
	if(root == NULL)
		return;

	for(i=0; i<root->childnum; i++) {
		(root->pChild[i])->pParent = root;
	}

	for(i=0; i<root->childnum; i++) {
		setTreeBagParent(root->pChild[i]);
	}
}
void DynamicP::printTreeBagParent(Tree_bag *root) {
	int i;
	if(root == NULL)
		return;

	for(i=0; i<root->childnum; i++) {
		printTreenode(root->pChild[i]);
		cout<<"'s parent=";
		printTreenode((root->pChild[i])->pParent);
		cout<<endl;
	}

	for(i=0; i<root->childnum; i++) {
		printTreeBagParent(root->pChild[i]);
	}
}

void DynamicP::removeTreeBagParent(Tree_bag *root) {
	int i;
	if(root == NULL)
		return;

	for(i=0; i<root->childnum; i++) {
		removeTreeBagParent(root->pChild[i]);
	}

	root->pParent = NULL;

//	for(i=0; i<root->childnum; i++) {
//		(root->pChild[i])->pParent = NULL;
//	}

}
//Printing out the conformational graph (only the connectivity so far).  
void DynamicP::print_outstructuregraph()
{
	int		numnode = this->getNumnode();

    int		i, j;

    printf("<graph>\n"); 
    for(i=0; i<numnode; i++)
	{
        for(j=0; j<numnode; j++) {
			printf("%d ", this->pCongraph[i][j]);
		}
		printf("\n");
	}
	printf("</graph>\n"); 
}

//Convert between all the vertices in the graph and the stems in the secondary structure.
void DynamicP::get_conversion_gtos()
{
	int numstem = this->getNumstems();
	int i;
	int lg_id, rg_id;

	//Allocate enough enough memory to store the information .
	gidtosid = new Sid_info[2*numstem + 2];

	gidtosid[0].s_id			= -1;	//The starting node (source).
	gidtosid[0].left_or_right	= -1;	//

	gidtosid[nodestart].s_id			= -2;	//The ending node (sink).
	gidtosid[nodestart].left_or_right	= -1;	//

	for(i=0; i<numstem; i++)
	{
		//Now, go through all the stems in the secondary structure and store the conversion relationship in the table.
        lg_id = this->pStem[i].lg_id;
        rg_id = this->pStem[i].rg_id;

        gidtosid[lg_id].left_or_right = 0;
        gidtosid[lg_id].s_id = i;

        gidtosid[rg_id].left_or_right = 1;
        gidtosid[rg_id].s_id = i;
    }
}

void DynamicP::print_gidtosid()
{
	printf("--- print out gidtosid array ---\n");
	int numstem = this->getNumstems();
	int i;

	for(i=0; i<2*numstem+2; i++) {
		cout<<"gidtosid["<<i<<"].s_id="<<gidtosid[i].s_id<<", left_or_right="<<gidtosid[i].left_or_right<<endl;
	}
}

void DynamicP::allocate_dptable(Tree_bag *root)
{
	if(root == NULL)
		return;

//	this->printTreenode(root);	cout<<endl;
	int		i, nodenum, stnum, iTotal;
	nodenum = root->nodenum;
	stnum	= root->stnum;

	//Now start allocating memory needed for storing intermediate computational results.
	root->enum_stem = new int[stnum];
	root->enum_node = new int[nodenum];

	if(root->pParent == root)	//the decision of being the root
	{
		iTotal = 1;
	} else {
		iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), 
									static_cast<double>(nodenum)));
	}

	//support max empty stem num 
	//scale the dp table by MaxEmptyStemNum
	if(getFlagSupportEmptyStem() == TRUE)
		iTotal *= (this->getMaxEmptyStemNum() + 1);
	
	root->t_head = new New_node_combination[iTotal];
/**/
	//expend the memory for the dp table when supporting empty-stem-model
	if(getFlagSupportEmptyStem() == TRUE) {
		for(int j=0; j<iTotal; j++) {
			root->t_head[j].lpos = new int[nodenum];
			root->t_head[j].rpos = new int[nodenum];

			//add extend_dir 
			root->t_head[j].l_ext_dir = new char[nodenum];
			root->t_head[j].r_ext_dir = new char[nodenum];
		}
	}

	for(i=0; i<root->childnum; i++)
		this->allocate_dptable(root->pChild[i]); 
}

void DynamicP::initDPTable(Tree_bag *root)
{
	int		i, nodenum, iTotal;
	nodenum = root->nodenum;

	if(root->pParent == root)	//the decision of being the root
	{
		iTotal = 1;
	} else {
		iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), 
									static_cast<double>(nodenum)));
	}

	//support max empty stem num 
	if(getFlagSupportEmptyStem() == TRUE)
		iTotal *= (this->getMaxEmptyStemNum()+1);

	for(i=0; i<iTotal; i++) {
		root->t_head[i].score	= SMALLEST;
		root->t_head[i].optimal	= -1;

		//limiting the number of empty-stem.  20090211
		root->t_head[i].overall_empty_stem_num	= 0;
/**/
		//expend the memory for the dp table when supporting empty-stem-model
		if(getFlagSupportEmptyStem() == TRUE)
		{
			for(int j=0; j<nodenum; j++)
			{
				root->t_head[i].lpos[j] = NEGATIVE_ONE;
				root->t_head[i].rpos[j] = NEGATIVE_ONE;

				//add extend_dir 
				root->t_head[i].l_ext_dir[j] = EXTNA;
				root->t_head[i].r_ext_dir[j] = EXTNA;
			}
		}

		//remove the block memory function in dp part.  20090208
		for(int j=0; j<root->childnum; j++)
			root->t_head[i].childIndex[j] = 0;	
	}

	for(i=0; i<root->childnum; i++)
		this->initDPTable(root->pChild[i]); 
}

//Construct the list of stems or half stems for all the nodes in the tree.
void DynamicP::construct_stemlist(Tree_bag *root) 
{
	int		node1, node2 = -1;
	int		node_cnt1, node_cnt2; 
	int		found_pair;
	int		stem_num;
	Strunit *	head = NULL, *tail = NULL;
	Node_list *	n_head1 = NULL, *n_head2 = NULL;
     
	n_head1 = root->nhead;
	head	= NULL;

	while(n_head1) {
		n_head1->visited = 0;	//Initialize the visited bit to be zero.
		n_head1 = n_head1->next;
	}

	node_cnt1	= 0;
	stem_num	= 0;
	n_head1		= root->nhead;

	while(n_head1)
	{
		if(n_head1->visited == 0)
		{
			node1 = n_head1->g_node;
			//just for trace. 
//			printf("node1=%d\n", node1);
			found_pair	= 0;
			n_head2		= n_head1->next; 
			node_cnt2	= node_cnt1 + 1;
             
			while(n_head2 != NULL && n_head2->visited != 0) {
				n_head2 = n_head2->next;
				node_cnt2++;
			}

			///if(n_head2!=NULL && n_head2->visited==0){
			if((n_head2 != NULL && n_head2->visited == 0) 
				|| (n_head2 == NULL && n_head1->visited == 0))
			{
                //If node2 has not been visited yet.
                while(n_head2)
				{
					node2 = n_head2->g_node;
					//just for trace. 
//					printf("\tnode2=%d\n", node2);
					if(gidtosid[node1].s_id == gidtosid[node2].s_id)
					{
						//If node1 and node2 belongs to the same stem.
						found_pair = 1;
						n_head2->visited = 1;
						break;
					}
					node_cnt2++;
					n_head2 = n_head2->next;
				}

                if(head == NULL && gidtosid[node1].s_id >= 0)
				{
					head = new Strunit[1];	//allocate memory
					if(found_pair == 1)
					{
						head->type = 1; 
						if(gidtosid[node1].left_or_right == 0)
						{ 
							head->g_id1 = node1;
							head->g_id2 = node2;
							head->left_nid = node_cnt1;
							head->right_nid = node_cnt2;
							head->stem_id = gidtosid[node1].s_id;
						} else {
							head->g_id1 = node2;
							head->g_id2 = node1;
							head->left_nid = node_cnt2;
							head->right_nid = node_cnt1;
							head->stem_id = gidtosid[node1].s_id;
						}
						stem_num++;
					} else {
						head->type = 2;
						head->g_id1 = node1;
						head->g_id2 = NEGATIVE_ONE;		//20081019
						head->left_nid = node_cnt1;
						head->right_nid= NEGATIVE_ONE;	//20081019
						head->stem_id = gidtosid[node1].s_id;
						stem_num++;
					}
					head->next = NULL;
					tail = head;
				}//if
                else if(gidtosid[node1].s_id>=0)
				{
					tail->next = new Strunit[1];	//allocate memory
					tail = tail->next;

					if(found_pair == 1)
					{
						tail->type = 1;
						if(gidtosid[node1].left_or_right == 0)
						{
							tail->g_id1 = node1;
							tail->g_id2 = node2;
							tail->left_nid = node_cnt1;
							tail->right_nid = node_cnt2;
							tail->stem_id = gidtosid[node1].s_id;
						} else {
							tail->g_id1 = node2;
							tail->g_id2 = node1;
							tail->left_nid = node_cnt2;
							tail->right_nid = node_cnt1;
							tail->stem_id = gidtosid[node1].s_id;
						}
						stem_num++;
					} else {
						tail->type = 2;
						tail->g_id1 = node1;
						tail->g_id2 = NEGATIVE_ONE;		//20081019
						tail->left_nid = node_cnt1;
						tail->right_nid= NEGATIVE_ONE;	//20081019
						tail->stem_id = gidtosid[node1].s_id;
						stem_num++;
					}
					tail->next = NULL;
				}//else if
			}//if(n_head2!=NULL && n_head2->visited==0)
		}//if
		node_cnt1++;
		n_head1 = n_head1->next;
	}//while
	root->stnum = stem_num; 
	root->su_stem = head;
}

void DynamicP::tree_stemlist(Tree_bag *root)
{
	int		i;
	this->construct_stemlist(root);

	for(i=0; i<root->childnum; i++)
		this->tree_stemlist(root->pChild[i]);
}

//Construct the list of loops for a given node in the tree.
//At the end, also incorporate stem information into the data structure.
void DynamicP::construct_looplist(Tree_bag *root)
{
	int numloop = this->getNumloops();
	int		i;
	int		node1, node2;
	int		loop_num;
	Node_list * n_head1 = NULL, * n_head2 = NULL;
	Strunit * head = NULL, * tail = NULL;
	Strunit * s_head = NULL;

	n_head1 = root->nhead;
	loop_num = 0;
	head = NULL;

	while(n_head1)
	{
		node1 = n_head1->g_node;
		//just for trace. 
//		printf("node1=%d\n", node1);
		n_head2 = n_head1->next;
        while(n_head2)
		{
			node2 = n_head2->g_node;
			//just for trace. 
//			printf("\tnode2=%d\n", node2);
			if(this->pCongraph[node1][node2] == DIRECTED || this->pCongraph[node2][node1] == DIRECTED)
			{
				for(i=0; i<numloop; i++)
				{
					if((node1 == this->pLoop[i].lg_id && node2 == this->pLoop[i].rg_id)
						|| (node1 == this->pLoop[i].rg_id && node2 == this->pLoop[i].lg_id))
					{
						if(head == NULL)
						{
							head = new Strunit[1];	//allocate memory
							head->type = 0;
							if(this->pCongraph[node1][node2] == DIRECTED) {
								head->g_id1 = node1;
								head->g_id2 = node2;
							} else {
								head->g_id1 = node2;
								head->g_id2 = node1;
							}
							head->loop_id = i;
							head->next = NULL;
							tail = head;
						} else {
							tail->next = new Strunit[1];	//allocate memory
							tail = tail->next;
							tail->type = 0;

							if(this->pCongraph[node1][node2] == DIRECTED){
								tail->g_id1 = node1;
								tail->g_id2 = node2;
							} else {
								tail->g_id1 = node2;
								tail->g_id2 = node1;
							}
							tail->loop_id = i;
							tail->next = NULL;
						}
						loop_num++;
					}//if
				}//for i.
			}//if
			n_head2 = n_head2->next;
		}//while
		n_head1 = n_head1->next;
	}//while
	root->lpnum = loop_num;
	root->su_loop = head;

	//Now, go through the linked list to incorporate more information into the data structure.
	while(head) {
		s_head = root->su_stem;
		//Initialize the values of lstem and rstem to be NULL.   
		head->lstem = NULL;
		head->rstem = NULL;

		while(s_head)
		{
			if(s_head->g_id1 == head->g_id1 || s_head->g_id2 == head->g_id1) {
				head->lstem = s_head;   
			}

			if(s_head->g_id2 == head->g_id2 || s_head->g_id1 == head->g_id2) {
				head->rstem = s_head;
			}

			s_head = s_head->next;
		}
		head = head->next;
	}
}

//Generate the loop list for all the nodes in the tree.  
void DynamicP::tree_looplist(Tree_bag *root)
{
    int		i;
/*
	cout<<"----tree_looplist----"<<endl;
	this->printTreenode(root);
	cout<<endl;
*/
	this->construct_looplist(root);

    for(i=0; i<root->childnum; i++)
			this->tree_looplist(root->pChild[i]);
}

//Generate the list of stems that is shared by both a parent node and a child node. 
//Also the array stem_image in the root node indicates whether a
//a particular stem is also contained in a child node or not.  
void DynamicP::generate_sharedstemlist(Tree_bag * parent, Tree_bag * child, int c_id)
{
//	this->printTreenode(parent);	cout<<"\t";
//	this->printTreenode(child);		cout<<endl;

	int		i;
	int		prog_id;
	int		node1, node2;
	int		n_found1, n_found2;
	Strunit		*s_head;
	Node_list	*n_head;
    
	//Check the list of stems.
	s_head = parent->su_stem;

	//Now allocate sufficient memory for the array. 
	parent->stem_image[c_id] = new int[parent->nodenum];	//allocate memory

	for(i=0; i<parent->nodenum; i++) {
		parent->stem_image[c_id][i] = -1; 
	}

	//Now, go through the list of stems in the stem list.   
	prog_id = 0;

	while(s_head)
	{
		n_found1 = 0;   ///
		n_found2 = 0;
		node1 = -1;
		node2 = -1;

		if(s_head->type == 1) {
           //If the s_head is a whole stem.
           node1 = s_head->g_id1;
           node2 = s_head->g_id2;
           n_found1 = 0;
           n_found2 = 0;
        } else {
           //Otherwise, s_head contains a half stem
			//and only the graph node for that particular half stem is searched in the following procedure.
           node1 = s_head->g_id1;
           n_found1 = 0;
        }

        n_head = child->nhead;

		while(n_head) {
			if(n_head->g_node == node1) {
				n_found1 = 1;
			}
			if(n_head->g_node == node2) {
				n_found2 = 1;
			}
			n_head = n_head->next;
		}

		if(s_head->type == 1 && n_found1 == 1 && n_found2 == 1) {
			//Set the corresponding flag value to be one.
			parent->stem_image[c_id][prog_id] = 1;
		}

		if(s_head->type == 2 && n_found1 == 1) {
			parent->stem_image[c_id][prog_id] = 1;
		}
		prog_id++;
		s_head = s_head->next;
    }//while
}

//Similarly, Generate the list of loops shared by a child and a root. 
//The loop_image array in the node root will takes the value of either 1 or -1, 
//when the value is 1, it means that the child also contains the loop.   
void DynamicP::generate_sharedlooplist(Tree_bag * parent, Tree_bag * child, int c_id)
{
     int	i;
     int	prog_id;
     int	node1, node2;
     int	n_found1, n_found2;
     Strunit *		l_head;
     Node_list *	n_head;

     //Now, allocate memory to the array.
     l_head = parent->su_loop;

     parent->loop_image[c_id] = new int[parent->nodenum];

	//Initialize all the array elements to be -1.
	for(i=0; i<parent->nodenum; i++) {
		parent->loop_image[c_id][i] = -1;
	}

	//Now go through the list of loops.
	prog_id = 0;

	while(l_head)
	{
        node1 = l_head->g_id1;
        node2 = l_head->g_id2;

        n_found1 = 0;
        n_found2 = 0;

        n_head = child->nhead;

		while(n_head)
		{
			if(n_head->g_node == node1) {
				n_found1 = 1;
			}
			if(n_head->g_node == node2) {
				n_found2 = 1;
			}
			n_head = n_head->next;
		}

		if(n_found1 == 1 && n_found2 == 1) {
			//Set up the loop id to be 1.
			parent->loop_image[c_id][prog_id] = 1;
        }

        prog_id++;
        l_head = l_head->next;
	}//while
}

//Generate the shared stem and loop list for all the tree nodes in the tree. 
void DynamicP::tree_sharedstemloop(Tree_bag * root)
{
	int		i;

	//Initialize both arrays to be NULL.
	for(i=0; i<MAXNUMPT; i++) {
		root->stem_image[i] = NULL;
		root->loop_image[i] = NULL;
	}

	for(i=0; i<root->childnum; i++) {
		this->generate_sharedstemlist(root, root->pChild[i], i);
		this->generate_sharedlooplist(root, root->pChild[i], i);
	}

	//Call the function recursively to generate the list for all the tree nodes. 
	for(i=0; i<root->childnum; i++) {
		this->tree_sharedstemloop(root->pChild[i]);
	}
}

//Generate the mapping array between a node and the parent for all the stem nodes in the tree bag node.      
void DynamicP::construct_mappingarray(Tree_bag * parent, Tree_bag * child)
{
	int		i;
	int		c_node, p_node;
	int		c_loc, p_loc;
	int		f_flag, s_count;
	Node_list	*c_head;
	Node_list	*p_head;
	Strunit		*s_head;

	//Allocate sufficient memory space for the array state_array.
	//make state_array as a local variable
	int * state_array = new int[child->nodenum];
	child->set_up = new int[child->stnum];

	//Initialize the state_array to be all -1.
	for(i=0; i<child->nodenum; i++) {
		state_array[i] = -1;
	}

	//Now start looking at every node in the child and
	//check if it is also included in the tree bag node of parent. 
	c_head = child->nhead;
	c_loc = 0;

	while(c_head)
	{
		c_node = c_head->g_node;
        p_loc = 0;
        p_head = parent->nhead;
        f_flag = 0;

        while(p_head)
		{
			p_node = p_head->g_node;
			if(p_node == c_node) {
				f_flag = 1;
				break;
			}
			p_loc++;
			p_head = p_head->next;
		}

		if(f_flag == 1) {
			//If the found flag is found to be one.
			state_array[c_loc] = p_loc;
		}

		c_loc++;
		c_head = c_head->next;
	}//while
 
	//Now, go through the list of stems in the child and generate the array of set_up in the child node 
    //such that a mapping between the stems of the child node and the nodes in the parent node is set up.  
    s_head = child->su_stem;
    s_count = 0;
    while(s_head)
	{
		if(s_head->type == 1)
		{
			if(state_array[s_head->left_nid] != -1) {
				child->set_up[s_count] = state_array[s_head->left_nid];
			} else if(state_array[s_head->right_nid] != -1) {
				child->set_up[s_count] = state_array[s_head->right_nid];
			} else {
				child->set_up[s_count] = -1;
			}
		} else if(s_head->type == 2) {
			if(state_array[s_head->left_nid] != -1){
				child->set_up[s_count] = state_array[s_head->left_nid];
			} else {
				child->set_up[s_count] = -1;
			}
		}
		s_count++;
		s_head = s_head->next;
	}//while

	if(state_array != NULL) {
		delete [] state_array;	
		state_array = NULL;
	}
}

//Constructing the mapping array for all the nodes in the tree recursively. 
void DynamicP::tree_mappingarray(Tree_bag * root) 
{
    int		i;
    int		childnum;

    //If root is found to be NULL.
    if(root == NULL)
		return;

    childnum = root->childnum;
    
	for(i=0; i<childnum; i++)
		this->construct_mappingarray(root, root->pChild[i]);

    for(i=0; i<childnum; i++) {
        this->tree_mappingarray(root->pChild[i]);
    }
}

void DynamicP::print_treebag_set_up(Tree_bag *root)
{
	int i, j;
	if(root == NULL)
		return;

	this->printTreenode(root);
	cout<<endl;
	//print out its child's set_up
    for(i=0; i<root->childnum; i++) {
			cout<<"\t";
			this->printTreenode(root->pChild[i]);
//			cout<<endl;
		for(j=0; j<root->pChild[i]->stnum; j++) {
			cout<<"\t"<<root->pChild[i]->set_up[j]<<" ";
		}
		cout<<endl;
	}

	//recursively call itself
    for(i=0; i<root->childnum; i++) {
        this->print_treebag_set_up(root->pChild[i]);
    }
}

void DynamicP::preprocess_tree(Tree_bag *root)
{
	//1.Perform a conversation between the graph node id and the stem id represented in the secondary structure.  
	this->get_conversion_gtos();
//	this->print_gidtosid();

//	this->free_pStem();

	//3.Construct the list of stems and loops for all nodes in the tree.
//	cout<<"---- this->tree_stemlist(root) ----"<<endl;
	this->tree_stemlist(root);

//	cout<<"---- this->tree_looplist(root) ----"<<endl;
	this->tree_looplist(root);

//	this->printSuStemLoopInfo(root);
//	exit(0);

	this->free_pCongraph();
//	this->free_pLoop();

	//2.Allocate enough memory to perform dynamic programming and store intermediate results.
	this->allocate_dptable(root);

	//4.Construct the list of shared stems and loops for all nodes in the tree.
	this->tree_sharedstemloop(root);

	//5.Construct the arrays used for mapping purposes. 
	this->tree_mappingarray(root);

//	cout<<"print_treebag_set_up"<<endl;
//	this->print_treebag_set_up(root);

}

void DynamicP::free_gtos()
{
	if(gidtosid != NULL)
		delete [] gidtosid;
	gidtosid = NULL;
}
void DynamicP::free_dptable(Tree_bag *root)
{
	int i;
	int nodenum = root->nodenum;
	//20080630
    for(i=0; i<root->childnum; i++)
		this->free_dptable(root->pChild[i]);

	//just for debug. 20080701
//	this->printTreenode(root);
//	cout<<endl;

	int iTotal;
	if(root->pParent == root)	//the decision of being the root
	{
		iTotal = 1;
	} else {
		iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), 
									static_cast<double>(nodenum)));
	}

	//support max empty stem num 
	if(getFlagSupportEmptyStem() == TRUE)
		iTotal *= (this->getMaxEmptyStemNum()+1);

	if(root->t_head != NULL) 
	{
/**/
		if(getFlagSupportEmptyStem() == TRUE)
		{
			for(int j=0; j<iTotal; j++)
			{
				if(root->t_head[j].lpos != NULL) {
					delete [] root->t_head[j].lpos;
					root->t_head[j].lpos = NULL;
				}
				if(root->t_head[j].rpos != NULL) {
					delete [] root->t_head[j].rpos;
					root->t_head[j].rpos = NULL;
				}

				//add extend_dir 
				if(root->t_head[j].l_ext_dir != NULL) {
					delete [] root->t_head[j].l_ext_dir;
					root->t_head[j].l_ext_dir = NULL;
				}
				if(root->t_head[j].r_ext_dir != NULL) {
					delete [] root->t_head[j].r_ext_dir;
					root->t_head[j].r_ext_dir = NULL;
				}
			}
		}

		delete [] root->t_head;			
		root->t_head = NULL;
	}
	if(root->enum_stem != NULL) {
		delete [] root->enum_stem;	
		root->enum_stem = NULL;
	}
	if(root->enum_node != NULL) {
		delete [] root->enum_node;	
		root->enum_node = NULL;
	}

}

Strunit * DynamicP::find_stemtail(Strunit *cur, int index)
{
	int i;
	Strunit * tail = cur;
	if(tail != NULL) {
		i = 0;
		while(i < index) {
			tail = tail->next;
			i++;
		}
	}
	return tail;
}

void DynamicP::free_stem(Tree_bag *root) 
{
	int i;
	//detect memory leak problem in rnatops (linux-valgrind)  20090815
	//Strunit * head, * tail;
	Strunit * head = NULL, * tail = NULL;
	head = root->su_stem;

	for(i=root->stnum; i>=0; i--)
	{
		tail = this->find_stemtail(head, i);
		if(tail != NULL)
		{
//			cout<<"delete "<<tail->g_id1<<","<<tail->g_id2<<endl;
//			if(tail->lstem != NULL) {
				tail->lstem = NULL;
//			}
//			if(tail->rstem != NULL) {
				tail->rstem = NULL;
//			}
			//detect memory leak problem in rnatops (linux-valgrind)  20090815
			//delete tail;
			delete [] tail;
			tail = NULL;
		}
	}
}

void DynamicP::free_stemlist(Tree_bag *root)
{
	int		i;

	//starts from the child node
	for(i=0; i<root->childnum; i++)
		this->free_stemlist(root->pChild[i]);

//	this->printTreenode(root);
//	cout<<endl;

	this->free_stem(root);
}

Strunit * DynamicP::find_looptail(Strunit *cur, int index)
{
	int i;
	Strunit * tail = cur;
	if(tail != NULL) {
		i = 0;
		while(i < index) {
			tail = tail->next;
			i++;
		}
	}
	return tail;
}

void DynamicP::free_loop(Tree_bag *root) 
{
	int i;
	Strunit * head, * tail;
	head = root->su_loop;

	for(i=root->lpnum; i>=0; i--)
	{
		tail = this->find_looptail(head, i);
		if(tail != NULL)
		{
//			cout<<"delete "<<tail->g_id1<<","<<tail->g_id2<<endl;
			delete [] tail;
			tail = NULL;
		}
	}

	//detect memory leak problem in rnatops (linux-valgrind)  20090815
//	if(root->su_loop != NULL) {
//		root->su_loop = NULL;
//	}
}

void DynamicP::free_looplist(Tree_bag * root)
{
    int		i;

    for(i=0; i<root->childnum; i++)
		this->free_looplist(root->pChild[i]);

//	this->printTreenode(root);
//	cout<<endl;

	this->free_loop(root);
}

void DynamicP::free_sharedstemlist(Tree_bag * root, Tree_bag * child, int c_id)
{
	if(root->stem_image[c_id] != NULL) {
		delete [] root->stem_image[c_id];
		root->stem_image[c_id] = NULL;
	}
}

//Similarly, Generate the list of loops shared by a child and a root. 
//The loop_image array in the node root will takes the value of either 1 or -1, 
//when the value is 1, it means that the child also contains the loop.   
void DynamicP::free_sharedlooplist(Tree_bag * root, Tree_bag * child, int c_id)
{
	if(root->loop_image[c_id] != NULL) {
		delete [] root->loop_image[c_id];  
		root->loop_image[c_id] = NULL;
	}
}

void DynamicP::free_sharedstemloop(Tree_bag * root)
{
	int		i;

	//Call the function recursively to generate the list for all the tree nodes. 
	for(i=0; i<root->childnum; i++)
		this->free_sharedstemloop(root->pChild[i]);

//	this->printTreenode(root);
//	cout<<endl;

	for(i=0; i<root->childnum; i++) {
			this->free_sharedstemlist(root, root->pChild[i], i);
			this->free_sharedlooplist(root, root->pChild[i], i);
	}
}

//Generate the mapping array between a node and the parent for all the stem nodes in the tree bag node.      
void DynamicP::deconstruct_mappingarray(Tree_bag * parent, Tree_bag * child)
{
	if(child->set_up != NULL) {
		delete [] child->set_up;
		child->set_up = NULL;
	}
}

//Constructing the mapping array for all the nodes in the tree recursively. 
void DynamicP::free_mappingarray(Tree_bag * root) 
{
//	cout<<"---- free_mappingarray ----"<<endl;
    int		i;

    //If root is found to be NULL.
    if(root == NULL)
		return;

    for(i=0; i<root->childnum; i++) {
        this->free_mappingarray(root->pChild[i]);
    }

//	this->printTreenode(root);
//	cout<<endl;

	for(i=0; i<root->childnum; i++)
		this->deconstruct_mappingarray(root, root->pChild[i]);
}

//Do not forget to free some pointers in Tree_bag, e.g. * nhead, * t_head, * t_pointers[MAXNUMPT].
void DynamicP::postprocess_tree(Tree_bag *root)
{
	this->free_gtos();

	this->free_dptable(root);

	this->free_stemlist(root);	//Do not forget to delete su_stem/su_loop->lstem/rstem. 

	this->free_looplist(root);

	this->free_sharedstemloop(root);

	this->free_mappingarray(root);
}
//------------------------------------------------------
void DynamicP::print_treenode(Tree_bag * node)
{
    Node_list *nodehead;

    nodehead = node->nhead;
    cout<<"{";
    while(nodehead) {
		//Go through the list of nodes, print out the head information.
		cout<<nodehead->g_node<<" ";
		nodehead = nodehead->next;
    }
    cout<<"}"<<endl;
}

void DynamicP::print_treenodelist(Node_list * node_list)
{
    Node_list *nodehead;

    nodehead = node_list;
    cout<<"{";
    while(nodehead) {
		//Go through the list of nodes, print out the head information.
		cout<<nodehead->g_node<<" ";
		nodehead = nodehead->next;
    }
    cout<<"}"<<endl;
}

void DynamicP::print_tree(Tree_bag * root)
{
	cout<<"---- print_tree ----"<<endl;
    int		i;

    //If root is found to be NULL.
    if(root == NULL)
		return;

    for(i=0; i<root->childnum; i++) {
        this->print_tree(root->pChild[i]);
    }

	this->print_treenode(root);
}

void DynamicP::printSuStemLoopInfo(Tree_bag * root)
{
	int		i;
	Strunit *	s_head;
	Strunit *	l_head;

	for(i=0; i<root->childnum; i++) {
		this->printSuStemLoopInfo(root->pChild[i]);
	}

	cout<<endl<<"compnode:\t";
	this->print_treenode(root);

	//1. printing out the info of su_stem
	cout<<"Info of su_stem"<<endl;
	s_head	= root->su_stem;
	while(s_head)
	{
		cout<<"\tg_id1:"<<s_head->g_id1<<"\t\tg_id2:"<<s_head->g_id2<<endl;
		cout<<"\tleft_nid:"<<s_head->left_nid<<"\tright_nid:"<<s_head->right_nid<<endl;
		s_head = s_head->next;
	}

	//2. printing out the info of su_loop
	cout<<"Info of su_loop"<<endl;
	l_head = root->su_loop;
	while(l_head)
	{
		cout<<"\tg_id1:"<<l_head->g_id1<<"\t\tg_id2:"<<l_head->g_id2<<endl;
//		cout<<"\tleft_nid:"<<l_head->left_nid<<"\tright_nid:"<<l_head->right_nid<<endl;
		l_head = l_head->next;
	}
}

Node_list * DynamicP::find_treebag_nhead_tail(Node_list *cur, int index)
{
	int i;
	Node_list * tail = cur;
	if(tail != NULL) {
		i = 0;
		while(i < index) {
			tail = tail->next;
			i++;
		}
	}
	return tail;
}

void DynamicP::free_treebag_nhead(Node_list *nodehead, int nodenum) 
{
	int i;
	Node_list * tail;

//	cout<<"delete:\t{";
	for(i=nodenum-1; i>=0; i--)
	{
		tail = this->find_treebag_nhead_tail(nodehead, i);
		if(tail != NULL)
		{
//			cout<<tail->g_node<<" ";

			//detect memory leak problem in rnatops  20090815
			//delete tail;
			//tail = NULL;
			if(tail->next != NULL) {
				delete [] tail->next;
				tail->next = NULL;
			}
		}
	}
//	cout<<"}"<<endl;
	//detect memory leak problem in rnatops  20090815
	//nodehead = NULL;
	if(nodehead != NULL) {
		delete [] nodehead;
		nodehead = NULL;
	}
}

void DynamicP::free_treenode(Tree_bag * node)
{
	this->free_treebag_nhead(node->nhead, node->nodenum);

	//the treenode
	if(node != NULL) {
		delete [] node;
		node = NULL;
	}
}

void DynamicP::free_tree(Tree_bag * root)
{
    int		i;
    
	for(i=0; i<root->childnum; i++) {
        this->free_tree(root->pChild[i]);
    }

//	this->printTreenode(root);
//	cout<<endl;

	this->free_treenode(root);
}

void DynamicP::buildTreeNodePath(Tree_bag * root)
{
	//get the whole path of all tree node
	iNumOftheWholeTreeNode = this->countNumOftheWholeTreeNode(root);
//	if(this->getFlagPrintDebugTreeInfo())
//		cout<<"iNumOftheWholeTreeNode="<<iNumOftheWholeTreeNode<<endl;

	pTheWholeTreeNode = new int[iNumOftheWholeTreeNode];	//
	memset(pTheWholeTreeNode, 0, iNumOftheWholeTreeNode);
	this->getTheWholeTreeNode(root);

//	if(this->getFlagPrintDebugTreeInfo()) {
//		for(int j=0; j<iTheWholeTreeNodeIndex; j++) {
//			cout<<pTheWholeTreeNode[j]<<"-";
//		}
//	}

	pOptimalTheWholeTreeNode = new int[iNumOftheWholeTreeNode];	//
}

void DynamicP::printTreeNodePath()
{
	cout<<endl<<"================ TreePathInfo ================"<<endl;
	cout<<"iNumOftheWholeTreeNode="<<iNumOftheWholeTreeNode<<endl;
	for(int j=0; j<iTheWholeTreeNodeIndex; j++) {
		cout<<pTheWholeTreeNode[j]<<"-";
	}
	cout<<endl;
}

void DynamicP::freeTreeNodePath()
{
	this->iOptimalTheWholeTreeNodeIndex = 0;
	if(pOptimalTheWholeTreeNode != NULL) {
		delete [] pOptimalTheWholeTreeNode;
		pOptimalTheWholeTreeNode = NULL;
	}

	this->iTheWholeTreeNodeIndex = 0;
	if(pTheWholeTreeNode != NULL) {
		delete [] pTheWholeTreeNode;
		pTheWholeTreeNode = NULL;
	}
}

void DynamicP::searchPK(int iWinSize, 
						int iMinSize, 
						Stem * pStemModel, 
						Loop * pLoopModel, 
						Tree_bag * root, 
						int shift_left, 
						char * genome_name, 
						int searchtype, 
						double * genome_baseFreq,
						ofstream &outfile)
{
	this->r_scores = new double[ACTUAL_CANDIDATE_NUM];
	this->cWindowSequence = new char[iWinSize+1];

	this->scan_genome(	iWinSize, 
						iMinSize, 
						pStemModel, 
						pLoopModel, 
						1,	//0: random sequence; 1: real genome sequence
						root, 
						this->getGenomeSequenceLength(), 
						this->getGenomeSequence(), 
						shift_left, 
						genome_name, 
						searchtype, 
						genome_baseFreq,
						GENOME_SEARCH_DIRECTION_PLUS,
						outfile);
	this->freeGenomeSequence();

	if(cWindowSequence != NULL) {
		delete [] this->cWindowSequence;
		this->cWindowSequence = NULL;
	}
	if(r_scores != NULL) {
		delete [] this->r_scores;
		this->r_scores = NULL;
	}

	if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_BOTH)
	{
		this->r_scores = new double[ACTUAL_CANDIDATE_NUM];
		this->cWindowSequence = new char[iWinSize+1];

		this->scan_genome(	iWinSize, 
							iMinSize, 
							pStemModel, 
							pLoopModel, 
							1,	//0: random sequence; 1: real genome sequence
							root, 
							this->getGenomeSequenceLength(), 
							this->getRGenomeSequence(), 
							shift_left, 
							genome_name, 
							searchtype, 
							genome_baseFreq,
							GENOME_SEARCH_DIRECTION_MINUS,
							outfile);
		this->freeRGenomeSequence();

		if(cWindowSequence != NULL) {
			delete [] this->cWindowSequence;
			this->cWindowSequence = NULL;
		}
		if(r_scores != NULL) {
			delete [] this->r_scores;
			this->r_scores = NULL;
		}
	}
}

double DynamicP::searchCurrentPos(	int iWinSize, 
									Stem *pStemModel, 
									Loop *pLoopModel, 
									int type, 
									Tree_bag * root, 
									int curPos, 
									int genomelength,
									ofstream &outfile)
{
	int j, k;
	int	iSize = iWinSize;

	if((curPos+iSize) < genomelength)
	{
		for(j=curPos; j<curPos+iSize; j++)
			cWindowSequence[j-curPos] = scanSequence[j];
		cWindowSequence[iSize] = '\0';
	}
	else
	{
		//
		for(j=curPos; j<genomelength; j++)
			cWindowSequence[j-curPos] = scanSequence[j];
		cWindowSequence[genomelength-curPos] = '\0';
	}
//	cout<<cWindowSequence<<endl;

	int		numstem = this->getNumstems();
	int		numloop = this->getNumloops();
	double	profilescore;

	//--------------------------------
	//start computing the time for preprocessing
	start_preprocessing = clock();

	CandCollect	col(numstem, numloop);

	//do local structure alignment  20090723
	if(this->getFlagLocalStructureAlign() == TRUE) {
		col.setTwoRegionScan(false);
		col.enableOffset(false);
	} else {
		col.setTwoRegionScan(true);
		col.enableOffset(true);
	}

//	bool bTwoReg = true;
//	col.setTwoRegionScan(bTwoReg);
//	bool offsetvld = true;
//	col.enableOffset( offsetvld );  //enable or disable offset

	double droprate = -1;
	col.setCYKMergeParamter(this->getFlagMergeCandInPreprocess(),
							this->getFlagCandwithShortestLength(),
							this->getShiftNumMergeCand(),
							droprate);

	double timesSD = 3.0;
	col.setCoeffStandardDev(timesSD);	//set the coefficient of standarad deviation

    double plowerbound=1.0;				//start to add penalty out of this bound
	col.setCYKLengthPenalty(plowerbound, this->getPcoeff());	//length penalty parameters

	//---------------------------------------------------------------------------
	if(this->getFlagPrintDebugDPWinInfo())
		cout<<"searchCandids function begins ..."<<endl;

	//support no-scanning 20081106
	if(this->getFlagFilterSearchBased() == TRUE) {
		//adding the empty-stem model.  20090211
		if(this->getFlagSupportEmptyStem() == TRUE) {
			col.searchCandids(pStemModel, numstem, (this->iLeading - 1), INVLDPROB, cWindowSequence, this->getFilterHitBeginPos());
		} else {
			col.searchCandids(pStemModel, numstem, this->iLeading, INVLDPROB, cWindowSequence, this->getFilterHitBeginPos());
		}
	} else {
		//adding the empty-stem model.  20090211
		if(this->getFlagSupportEmptyStem() == TRUE) {
			col.searchCandids(pStemModel, numstem, (this->iLeading - 1), INVLDPROB, cWindowSequence);
		} else {
			col.searchCandids(pStemModel, numstem, this->iLeading, INVLDPROB, cWindowSequence);
		}
	}

	if(this->getFlagPrintDebugDPWinInfo())
		cout<<"searchCandids function ends ..."<<endl;

	//support stub test 20081116
	this->setFlagStubTest(FALSE);
/*
	if(this->getFlagPrintDebugInfo()) {
		this->setFlagStubTest(TRUE);
	} else {
		this->setFlagStubTest(FALSE);
	}
*/

	int iCandidateStemNum;

	//detect memory leak problem in rnatops  20090815
	//this->candidateStems = new CandidateStem[numstem];

	//make sure every stem has candidates, if one of them has none, then move to the next pos
	bNoCandidateInCurrentWindow = 1;
	for(j=0; j<numstem; j++)
	{
		iCandidateStemNum = col.stemgrp[j].num;
		if(iCandidateStemNum > 0)
			bNoCandidateInCurrentWindow *= 1;
		else
			bNoCandidateInCurrentWindow *= 0;
	}

	//check the candidate alignment info.  20100201
	int num_charid_left = 0, num_charid_right = 0;
	bool bSkip = false;

	//20090103
	if(bNoCandidateInCurrentWindow == TRUE)
	{
		//detect memory leak problem in rnatops  20090815
		this->candidateStems = new CandidateStem[numstem];

		for(j=0; j<numstem; j++)
		{
/*
			if(this->getFlagPrintDebugInfo()) {
				cout<<"---------- stem["<<j<<"] ----------"<<endl;
			}
*/
			iCandidateStemNum = col.stemgrp[j].num;

			//add candidate number to be k
			candidateStems[j].num = this->getLeadingNum();
			candidateStems[j].s_locinfo = new Stemloc[this->getLeadingNum()];	//

			if(this->getFlagSupportEmptyStem() == TRUE)
			{
				if(this->getFlagPrintStrAlignInfo()) {
					candidateStems[j].leftArm	= new StructureAlign[this->getLeadingNum()-1];
					candidateStems[j].rightArm	= new StructureAlign[this->getLeadingNum()-1];
				}
			} else {
				if(this->getFlagPrintStrAlignInfo()) {
					candidateStems[j].leftArm	= new StructureAlign[this->getLeadingNum()];
					candidateStems[j].rightArm	= new StructureAlign[this->getLeadingNum()];
				}
			}

			for(k=0; k<iCandidateStemNum; k++)
			{
//				//support stub test 20081116
				if(this->getFlagStubTest() == TRUE)	//stub test
				{
				} 
				else 
				{
					//if not stub test
					candidateStems[j].s_locinfo[k].a = col.stemgrp[j].pMember[k].stempos[0][0];
					candidateStems[j].s_locinfo[k].b = col.stemgrp[j].pMember[k].stempos[0][1];
					candidateStems[j].s_locinfo[k].c = col.stemgrp[j].pMember[k].stempos[1][0];
					candidateStems[j].s_locinfo[k].d = col.stemgrp[j].pMember[k].stempos[1][1];
					candidateStems[j].s_locinfo[k].p_score = this->getStemWeight() * col.stemgrp[j].pMember[k].prob;
				}
/**/
				//add extend_dir 
				candidateStems[j].s_locinfo[k].a_ext_dir = EXTNO;
				candidateStems[j].s_locinfo[k].b_ext_dir = EXTNO;
				candidateStems[j].s_locinfo[k].c_ext_dir = EXTNO;
				candidateStems[j].s_locinfo[k].d_ext_dir = EXTNO;

				if(this->getFlagPrintStrAlignInfo()) {
					//
					char	*cLeftSequence = NULL, *cRightSequence = NULL, 
							*cLeftStructureAlign = NULL, *cRightStructureAlign = NULL;
					char	*leftAlign = NULL, *rightAlign = NULL;
					//support stemid index in two pastalines 20081015
					char	*leftAlign_charid_idx = NULL, *rightAlign_charid_idx = NULL;
					int		iLeftLength=0, iRightLength=0;

					col.getAlignmentStr(j, k, 
										cLeftSequence, cLeftStructureAlign, iLeftLength, 
										cRightSequence, cRightStructureAlign, iRightLength);

					candidateStems[j].leftArm[k].sequence	= new char[iLeftLength + 1];
					candidateStems[j].leftArm[k].alignment	= new char[iLeftLength + 1];

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
						candidateStems[j].leftArm[k].pastaline_idx	= new char[iLeftLength + 1];

					candidateStems[j].rightArm[k].sequence	= new char[iRightLength + 1];
					candidateStems[j].rightArm[k].alignment	= new char[iRightLength + 1];

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
						candidateStems[j].rightArm[k].pastaline_idx	= new char[iRightLength + 1];

					//check the candidate alignment info.  20100201
					num_charid_left = countNumChar(cLeftStructureAlign, '[');
					num_charid_right = countNumChar(cRightStructureAlign, ']');
					if(num_charid_left != num_charid_right) {
						bSkip = true;
					}

					//replace cLeftStructureAlign/cRightStructureAlign "[" and "]" with this->pStem[j].charid
					//toupper/tolower	 20100315
					leftAlign	= strrpl(cLeftStructureAlign, '[', toupper(this->pStem[j].charid));	//left
					rightAlign	= strrpl(cRightStructureAlign, ']', tolower(this->pStem[j].charid));	//right
					strcpy(candidateStems[j].leftArm[k].sequence,	cLeftSequence);
					strcpy(candidateStems[j].leftArm[k].alignment,	leftAlign);
					strcpy(candidateStems[j].rightArm[k].sequence,	cRightSequence);
					strcpy(candidateStems[j].rightArm[k].alignment, rightAlign);

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
					{
						leftAlign_charid_idx	= strrpl(cLeftStructureAlign, '[', this->pStem[j].charid_idx);
						rightAlign_charid_idx	= strrpl(cRightStructureAlign, ']', this->pStem[j].charid_idx);
						strcpy(candidateStems[j].leftArm[k].pastaline_idx,	leftAlign_charid_idx);
						strcpy(candidateStems[j].rightArm[k].pastaline_idx, rightAlign_charid_idx);
					}

					if(cLeftSequence != NULL) {
						delete [] cLeftSequence;
						cLeftSequence = NULL;
					}
					if(cRightSequence != NULL) {
						delete [] cRightSequence;
						cRightSequence = NULL;
					}
					if(cLeftStructureAlign != NULL) {
						delete [] cLeftStructureAlign;
						cLeftStructureAlign = NULL;
					}
					if(cRightStructureAlign != NULL) {
						delete [] cRightStructureAlign;
						cRightStructureAlign = NULL;
					}
					if(leftAlign != NULL) {
						delete [] leftAlign;
						leftAlign = NULL;
					}
					if(rightAlign != NULL) {
						delete [] rightAlign;
						rightAlign = NULL;
					}
					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
					{
						if(leftAlign_charid_idx != NULL) {
							delete [] leftAlign_charid_idx;
							leftAlign_charid_idx = NULL;
						}
						if(rightAlign_charid_idx != NULL) {
							delete [] rightAlign_charid_idx;
							rightAlign_charid_idx = NULL;
						}
					}
				}//if(this->getFlagPrintStrAlignInfo())
			}//for k

			//add candidate number to be k
			if(this->getFlagSupportEmptyStem() == TRUE)
			{
				// 20090614
				//empty stem model. the last one is reserved for empty stem candidate
				for(k=iCandidateStemNum; k<this->getLeadingNum()-1; k++)
				{
					candidateStems[j].s_locinfo[k].a = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].b = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].c = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].d = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].p_score = MYTHRESHOLD;

					candidateStems[j].leftArm[k].sequence	= NULL;
					candidateStems[j].leftArm[k].alignment	= NULL;
					candidateStems[j].rightArm[k].sequence	= NULL;
					candidateStems[j].rightArm[k].alignment	= NULL;
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						candidateStems[j].leftArm[k].pastaline_idx	= NULL;
						candidateStems[j].rightArm[k].pastaline_idx	= NULL;
					}

					//add extend_dir 
					candidateStems[j].s_locinfo[k].a_ext_dir = EXTNO;
					candidateStems[j].s_locinfo[k].b_ext_dir = EXTNO;
					candidateStems[j].s_locinfo[k].c_ext_dir = EXTNO;
					candidateStems[j].s_locinfo[k].d_ext_dir = EXTNO;
				}

				if(k == this->getLeadingNum()-1)
				{
					candidateStems[j].s_locinfo[k].a = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].b = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].c = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].d = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].p_score = ZERO;
				}

				//add extend_dir 
				candidateStems[j].s_locinfo[k].a_ext_dir = EXTNA;
				candidateStems[j].s_locinfo[k].b_ext_dir = EXTNA;
				candidateStems[j].s_locinfo[k].c_ext_dir = EXTNA;
				candidateStems[j].s_locinfo[k].d_ext_dir = EXTNA;

			} else {
				// 20090614
				//non-empty stem model
				for(k=iCandidateStemNum; k<this->getLeadingNum(); k++)
				{
					candidateStems[j].s_locinfo[k].a = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].b = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].c = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].d = NEGATIVE_ONE;
					candidateStems[j].s_locinfo[k].p_score = MYTHRESHOLD;

					candidateStems[j].leftArm[k].sequence	= NULL;
					candidateStems[j].leftArm[k].alignment	= NULL;
					candidateStems[j].rightArm[k].sequence	= NULL;
					candidateStems[j].rightArm[k].alignment	= NULL;
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						candidateStems[j].leftArm[k].pastaline_idx	= NULL;
						candidateStems[j].rightArm[k].pastaline_idx	= NULL;
					}

					//add extend_dir 
					candidateStems[j].s_locinfo[k].a_ext_dir = EXTNO;
					candidateStems[j].s_locinfo[k].b_ext_dir = EXTNO;
					candidateStems[j].s_locinfo[k].c_ext_dir = EXTNO;
					candidateStems[j].s_locinfo[k].d_ext_dir = EXTNO;
				}
			}
		}//for j

		//print out all the candidates info
		if(this->getFlagPrintDebugInfo()) 
		{
			for(j=0; j<numstem; j++)
			{
				cout<<"---------- stem["<<j<<"] ----------"<<endl;

				if(this->getFlagSupportEmptyStem() == TRUE)
				{
					//for(k=0; k<iCandidateStemNum+1; k++)
					//add candidate number to be k
					for(k=0; k<this->getLeadingNum(); k++)
					{
						cout<<"("<<candidateStems[j].s_locinfo[k].a<<"~"<<candidateStems[j].s_locinfo[k].b<<")"
							<<"-("<<candidateStems[j].s_locinfo[k].c<<"~"<<candidateStems[j].s_locinfo[k].d<<")"
							<<" \t"<<candidateStems[j].s_locinfo[k].p_score<<endl;
					}
				} else {
					//for(k=0; k<iCandidateStemNum; k++)
					//add candidate number to be k
					for(k=0; k<this->getLeadingNum(); k++)
					{
						cout<<"("<<candidateStems[j].s_locinfo[k].a<<"~"<<candidateStems[j].s_locinfo[k].b<<")"
							<<"-("<<candidateStems[j].s_locinfo[k].c<<"~"<<candidateStems[j].s_locinfo[k].d<<")"
							<<" \t"<<candidateStems[j].s_locinfo[k].p_score<<endl;
					}
				}
			}
		}//if(this->getFlagPrintDebugInfo()) 
	}
	//-------------------
	col.freeAllCandids();

	if(bSkip) {
		profilescore = SMALLEST;
	} else {
		if(bNoCandidateInCurrentWindow == 1)
		{
			// 20100215
			//printCandidateStemInfo_preprocessingInfo(outfile);

			end_preprocessing = clock();
			double	cpu_time_preprocessing = ((double) (end_preprocessing - start_preprocessing)) / CLOCKS_PER_SEC;
			cpu_time_hours_preprocessing += cpu_time_preprocessing/3600;

			start_dp = clock();
			this->compute_profilescore(curPos, root, cWindowSequence, pLoopModel);
	//		profilescore = this->getMaxProfilescore(root);
			profilescore = this->getMaxProfilescore2(root);
			end_dp = clock();
			double	cpu_time_dp = ((double) (end_dp - start_dp)) / CLOCKS_PER_SEC;
			cpu_time_hours_dp += cpu_time_dp/3600;
		} else {
			profilescore = SMALLEST;
		}
	}
	return profilescore;
}
void DynamicP::freeMemForCandidateStems()
{
	int		j, k;
	int		numstem = this->getNumstems();
	
	//adding the empty-stem model.  20090211
	int		iCandidateStemNum = 0;

	//start: free the memory for candidateStems
	for(j=0; j<numstem; j++)
	{
		//
		if(this->getFlagPrintStrAlignInfo()) 
		{
			if((this->candidateStems[j].leftArm	!= NULL)
				&& (this->candidateStems[j].rightArm != NULL))
			{

				//adding the empty-stem model.  20090211
				iCandidateStemNum = candidateStems[j].num;

				if(this->getFlagSupportEmptyStem() == TRUE)
					iCandidateStemNum -= 1;

				//adding the empty-stem model.  20090211
				//for(k=0; k<candidateStems[j].num; k++) {
				for(k=0; k<iCandidateStemNum; k++) {
					//add candidate number to be k
					if(candidateStems[j].leftArm[k].sequence != NULL) {
						delete [] candidateStems[j].leftArm[k].sequence;
						candidateStems[j].leftArm[k].sequence	= NULL;
					}
					if(candidateStems[j].leftArm[k].alignment != NULL) {
						delete [] candidateStems[j].leftArm[k].alignment;
						candidateStems[j].leftArm[k].alignment	= NULL;
					}
					if(candidateStems[j].rightArm[k].sequence != NULL) {
						delete [] candidateStems[j].rightArm[k].sequence;
						candidateStems[j].rightArm[k].sequence	= NULL;
					}
					if(candidateStems[j].rightArm[k].alignment != NULL) {
						delete [] candidateStems[j].rightArm[k].alignment;
						candidateStems[j].rightArm[k].alignment	= NULL;
					}
					// 20100124
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						if(candidateStems[j].leftArm[k].pastaline_idx != NULL) {
							delete [] candidateStems[j].leftArm[k].pastaline_idx;
							candidateStems[j].leftArm[k].pastaline_idx	= NULL;
						}
						if(candidateStems[j].rightArm[k].pastaline_idx != NULL) {
							delete [] candidateStems[j].rightArm[k].pastaline_idx;
							candidateStems[j].rightArm[k].pastaline_idx	= NULL;
						}
					}

				}
				delete [] this->candidateStems[j].leftArm;
				delete [] this->candidateStems[j].rightArm;
				this->candidateStems[j].leftArm	 = NULL;
				this->candidateStems[j].rightArm = NULL;
			}
		}
		if(this->candidateStems[j].s_locinfo != NULL)
		{
			delete [] this->candidateStems[j].s_locinfo;
			this->candidateStems[j].s_locinfo = NULL;
		}
	}
	delete [] this->candidateStems;
	this->candidateStems = NULL;
	//end:  free the memory for candidateStems
}
void DynamicP::scan_genome(	int		iWinSize, 
							int		iMinSize, 
							Stem *	pStemModel, 
							Loop *	pLoopModel, 
							int		type,	//0: random sequence; 1: real genome sequence
							Tree_bag * root, 
							int		genomelength, 
							char *	genome_sequence, 
							int		shift_left, 
							char *	genome_name, 
							int		search_type, 
							double * genome_baseFreq,
							int		search_direction,
							ofstream &outfile)
{
	int		iSize = iWinSize;
	int		numstem = this->getNumstems();
	double	profilescore;
	int		segstartpos, segendpos;	//the absolute value of the ending location of a sequence segment

	//1. First of all, decide to take the random sequence or the genome sequence.
	if(type == 0)
	{
		//Firstly, deal with random sequence then computing the threshold
	} else {
		//scan the genome sequence.
		this->scanLength = genomelength;
		this->scanSequence = new char[scanLength+1];
		strcpy(scanSequence, genome_sequence);
	}

	if(this->getFlagPrintDebugDPWinInfo()) {
		cout<<"scanLength="<<scanLength<<endl;
		cout<<"window size="<<iSize<<endl;
	}

	memset(cWindowSequence,0,iSize+1);

	//skip-and-jump strategy. 
//	int		stepsize = this->getStepSize();	//incremented stepsize

	int		start_pos, cur_start_pos, next_start_pos;
	int		pre_end_pos, cur_end_pos;
	double	optimal_score;
	int		optimal_pos = ZERO;//, optimal_pos_start, optimal_pos_end;

	int		curPos = 0;	//start position
	int		endPos = 0;
/*
	if(search_type == SEARCH_SUBSTRUCTUREFILTER)
		endPos = scanLength - iSize + 1;
	else	//if(search_type == SEARCH_WHOLESTRUCTURE)
		endPos = scanLength - iMinSize + 1;
*/
	//support no-scanning 20081029
	if(iMinSize == -1)
		endPos = 0;
	
	if(this->getFlagPrintDebugDPWinInfo())
		cout<<"endPos="<<endPos<<endl;

	//compute the starting time
	double	cpu_time_all = 0.0;

	start_all = clock();
	cpu_time_hours_all = 0.0;
	cpu_time_hours_preprocessing = 0.0;
	cpu_time_hours_dp = 0.0;

	optimal_score = SMALLEST;
	pre_end_pos = cur_end_pos = -1;
	start_pos = cur_start_pos = next_start_pos = -1;
	int lastStemId = (stemIdxArray[2*numstem-1]-1)/2;

//	int extension_left	= 0;
//	int extension_right = 0;
//	int	hitlength		= 0;
//	char * c_subgenome_seq = NULL;

	int hit_pos_start	= 0;
	int hit_pos_end		= 0;

	do
	{
		//output only one invalid combination 20081116
		//for debugging
		//this->setFlagInvalidCombination(FALSE);
		this->initDPTable(root);

		//profilescore = searchCurrentPos(iWinSize, pStemModel, pLoopModel, type, root, curPos, genomelength);
		profilescore = searchCurrentPos(iWinSize, pStemModel, pLoopModel, type, root, curPos, genomelength, outfile);

		if(this->getFlagPrintDebugDPWinInfo()) {
			cout<<curPos<<"\t";
			cout<<profilescore<<endl;
		}

		if(type == ZERO) {
			//do nothing
		} else {
			if(profilescore > this->getThreshold()) 
			{
				memset(pOptimalTheWholeTreeNode, 0, iNumOftheWholeTreeNode);
				iOptimalTheWholeTreeNodeIndex = 0;

				//if(this->getFlagSupportEmptyStem() == TRUE)
					if(this->getFlagPrintDebugDPSearchInfo())
						this->printTraceBackTreeInfo(root);	//just for debug 20090209
				
				//this->printTreeNodePath();	//just for debug 20090209

				this->getOptimalPath(root);
				//this->printOptimalTheWholeTreeNode();	//just for debug 20090209

				this->buildStemCandIdxArray();
//				this->printStemCandIdxArray();	//just for debug 20090209

				int empty_stem_num = 0;
				
				// 20100215
				if(this->checkValidation() == TRUE) {
					if(this->getFlagSupportEmptyStem() == TRUE)
						empty_stem_num = this->countEmptyStemNumInStemCandIdxArray();
					
					this->locateStemStartEndPosIdx();
					segstartpos = candidateStems[0].s_locinfo[iStemStartPosIdx].a;
					segendpos	= candidateStems[lastStemId].s_locinfo[iStemEndPosIdx].d;

					//adding the empty-stem model.  20090211
					if(segstartpos == NEGATIVE_ONE) {
						segstartpos		= 0;
					}
					if(segendpos == NEGATIVE_ONE) {
						//loop until it reach one candidate with its .d > 0
						int k = 0;
						int index = 2*numstem-1-k;
						int cur_stem_index	= (stemIdxArray[index]-1)/2;
						int	candidate_id	= stemCandIdxArray[index];//[cur_stem_index];
						while(candidateStems[cur_stem_index].s_locinfo[candidate_id].d == NEGATIVE_ONE)
						{
							k++;
							index = 2*numstem-1-k;
							cur_stem_index	= (stemIdxArray[index]-1)/2;
							candidate_id	= stemCandIdxArray[index];//[cur_stem_index];
						}
						segendpos	= candidateStems[cur_stem_index].s_locinfo[candidate_id].d;
					}

					hit_pos_start	= optimal_pos+segstartpos+shift_left;
					hit_pos_end		= optimal_pos+segendpos+shift_left;

					r_scores[numhit++] = profilescore;

					printWholeSearchHitIndex(numhit+this->getPreHitNum(), outfile);
					outfile<<genome_name<<endl;

					//support hit sorting 20081121
					oneHit.setAlignScore(profilescore);
					oneHit.setGenomeName(genome_name);
					oneHit.setEmptyStemNum(empty_stem_num);	//limiting the number of empty-stem. 

					if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_PLUS)
					{
						outfile<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
						oneHit.setDirection(GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS);
					}
					else if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_MINUS)
					{
						outfile<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
						oneHit.setDirection(GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS);
					}
					else if(search_direction == GENOME_SEARCH_DIRECTION_PLUS)
					{
						outfile<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
						oneHit.setDirection(GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS);
					}
					else if(search_direction == GENOME_SEARCH_DIRECTION_MINUS)
					{
						outfile<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
						oneHit.setDirection(GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS);
					}

	//				printHitPos(hit_pos_start, hit_pos_end);
	//				printWholeStructureHitScore(profilescore);
					printHitPos(hit_pos_start, hit_pos_end, outfile);
					printWholeStructureHitScore(profilescore, outfile);

					//support hit sorting 20081121
					oneHit.setPosStart(hit_pos_start);
					oneHit.setPosEnd(hit_pos_end);
					oneHit.setPastaLineNum(this->getNumPastaLine());

					this->buildFoldedStructureInfo();
					this->printFoldedStructureInfo(segstartpos, segendpos, hit_pos_start, hit_pos_end, outfile);

					this->buildStructureAlignInfo(pLoopModel, segstartpos, segendpos);
					this->printStructureAlignInfo2(hit_pos_start, hit_pos_end, outfile);

					if(this->getFlagPrintScoreInfo()) {
						this->printScoreForEveryComponent(pLoopModel, false, outfile);
						outfile<<"-------------------"<<endl;
						this->printCandidateStemInfo(outfile);
					}
				} else {
					if(bNoCandidateInCurrentWindow == TRUE) {
						//just for the debug 20081111
						if(this->getFlagPrintScoreInfo()) {
							outfile<<"No hit with score = "<<profilescore<<endl;
							this->printCandidateStemInfo_FailedCase(outfile);
						}
					}
				}

			}
			else
			{
//				//20090103
				if(bNoCandidateInCurrentWindow == TRUE) {
					//just for the debug 20081111
					if(this->getFlagPrintScoreInfo()) {
						outfile<<"No hit with score = "<<profilescore<<endl;
						this->printCandidateStemInfo_FailedCase(outfile);
					}
				}
			}

//			//20090103
			if(bNoCandidateInCurrentWindow == TRUE)
				freeMemForCandidateStems();
		}
	} while(curPos < endPos);

	end_all = clock();
	cpu_time_all = ((double) (end_all - start_all)) / CLOCKS_PER_SEC;
	cpu_time_hours_all += cpu_time_all/3600;

	if(this->scanSequence != NULL) {
		delete [] this->scanSequence;
		this->scanSequence = NULL;
	}
}

int	DynamicP::getStemIdbyLeftNodeId(int leftnodeid)
{
	int stemid = (leftnodeid - 1) / 2;
	return stemid;
}

void DynamicP::compute_profilescore(int j, Tree_bag * tr_node, char * cSequence, Loop *pLoopModel)
{
	int		i;
	int		sum;
	int		iValue, iIndex;
	double	score, child_score;
	double	l_score;

	int		l_id;
	int		l_end;//, l_start;
	int		r_start;//, r_end;
	int		loop_end;
	int		loop_len;

	int		iTmp;
	iTmp =	strlen(cSequence);

	Strunit *	s_head = NULL;
	Strunit *	l_head = NULL;
	Strunit *	tail_head = NULL;

	for(i=0; i<tr_node->childnum; i++) {
		this->compute_profilescore(j, tr_node->pChild[i], cSequence, pLoopModel);
	}

	if(this->getFlagPrintDebugDPSearchInfo()) {
		cout<<endl
			<<"compnode:\t";
		this->print_treenode(tr_node);

		//debug  20100325
		if(tr_node->nhead->g_node == 5) {
			if(tr_node->nhead->next->g_node == 16) {
				cout<<"hello"<<endl;
			}
		}
	}

	//enum_stem:	all possible combinations of structural units(stem).
	for(i=0; i<tr_node->stnum; i++) {
		tr_node->enum_stem[i] = 0;
	}

	//enum_node:	all possible combinations of nodes.
	for(i=0; i<tr_node->nodenum; i++) {
		tr_node->enum_node[i] = 0;
	}

	VTBSearch vs(1);
	vs.setAllowedInsNum(this->getAllowedNullLoopInsNum());

	//mark the invalid combination 20081116
	int	invalid_combination_mark = 0;

	//limiting the number of empty-stem.  20090211
	int current_empty_stem_num = 0;
	int	current_empty_stemid_array[MAXNUMEMPTYSTEM];

	//Compute the maximum value of max_end that can be reached for this particular location on the genome.
	while(1)
	{
		//limiting the number of empty-stem.  20090211
		if(this->getFlagSupportEmptyStem() == TRUE) {
			current_empty_stem_num = 0;
		
			//limiting the number of empty-stem.  20090213
			for(i=0; i<MAXNUMEMPTYSTEM; i++) {
				current_empty_stemid_array[i] = NEGATIVE_ONE;
			}
		}

		for(i=0; i<tr_node->childnum; i++) {
			tr_node->prob_score[i] = 0.0;
		}

		//1. dealing with stems.
		//The list of structural units which include both stems and half stems.
		s_head	= tr_node->su_stem;
		iIndex	= 0;
		score	= INITSCORE;//0.0;

		while(s_head)
		{
			if(s_head->type == 1)	//type=1, a stem.
			{
				//If the stem unit contains a complete stem.
				iValue			= tr_node->enum_stem[iIndex];
				score			+= this->candidateStems[s_head->stem_id].s_locinfo[iValue].p_score;

				//Compute the real locations of the stems on the genome sequence.
				s_head->start1	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].a;
				s_head->end1	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].b;
				s_head->start2	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].c;
				s_head->end2	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].d;

				//consider the score of parts that is shared with children.
				//save them and substract them when possible.
				for(i=0; i<tr_node->childnum; i++)
				{
					if(tr_node->pChild[i] != NULL && tr_node->stem_image[i][iIndex] == 1)
					{
						//If it is included in the child node.  
						tr_node->prob_score[i] += this->candidateStems[s_head->stem_id].s_locinfo[iValue].p_score;               
					}
				}

				tr_node->enum_node[s_head->left_nid]	= iValue;
				tr_node->enum_node[s_head->right_nid]	= iValue;

				//limiting the number of empty-stem.  20090211
				if(this->getFlagSupportEmptyStem() == TRUE) {
					if(iValue == (this->getLeadingNum()-1)) {
						//get this emptystem's stemid
						current_empty_stemid_array[current_empty_stem_num] = this->getStemIdbyLeftNodeId(s_head->g_id1);

						current_empty_stem_num++;
					}
				}
			}
			else
			{
				//Otherwise, only part of the stem is included.
				iValue	= tr_node->enum_stem[iIndex];

				if(gidtosid[s_head->g_id1].left_or_right == 0)
				{
					//if it is the left part of the stem.
					s_head->start1	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].a;
					s_head->end1	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].b;

				} else {
					//Otherwise it is the right part of the stem.
					s_head->start2	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].c;
					s_head->end2	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].d;
				}

				tr_node->enum_node[s_head->left_nid] = iValue;
			}
			iIndex++; 
			s_head = s_head->next;
		}
/*
		cout<<"Score after considering stems = "<<score<<endl;
		this->print_enum_node(tr_node->enum_node, tr_node->nodenum);
		cout<<endl;
*/
		//2. dealing with loops.
		l_head	= tr_node->su_loop;
		l_id	= 0;
		int		l_source, l_sink;

		//compute the pos for the empty-stem candidate
		int preNodeId = NEGATIVE_ONE, curNodeId = NEGATIVE_ONE;
		int preNodeIndex = NEGATIVE_ONE, curNodeIndex = NEGATIVE_ONE;
		int preCandId = NEGATIVE_ONE, curCandId = NEGATIVE_ONE;
		bool bGo = true;
		Node_list * nodeList = tr_node->nhead;

//		curNodeIndex = ZERO;
		while(l_head)	// && (this->getScoreOption() != SCORE_STEM)
		{
			bGo = true;

			//avoid computing the loop score for the part with either part of the stem is empty
			if(this->getFlagSupportEmptyStem() == TRUE) {
				preNodeId = l_head->g_id1;
				curNodeId = l_head->g_id2;

				preNodeIndex = this->getNodeIndexInTreenode(nodeList, preNodeId);
				curNodeIndex = this->getNodeIndexInTreenode(nodeList, curNodeId);

				if(preNodeIndex != NEGATIVE_ONE)
					preCandId = tr_node->enum_node[preNodeIndex];
				else {
					cout<<"Something wrong when computing index for nodeid"<<endl;
					exit(0);
				}

				if(curNodeIndex != NEGATIVE_ONE)
					curCandId = tr_node->enum_node[curNodeIndex];
				else {
					cout<<"Something wrong when computing index for nodeid"<<endl;
					exit(0);
				}

				if(preCandId == (this->getLeadingNum()-1) || curCandId == (this->getLeadingNum()-1))
					bGo = false;
			}

			if(bGo)
			{
				l_source = l_sink = 0;

				if(node_mapping[l_head->g_id1].node_type != 0)
				{
					//If the loop doesn't start with the source. 
					if(gidtosid[l_head->g_id1].left_or_right == 0) {
						l_end	= (l_head->lstem)->end1 + 1;
					} else {
						l_end	= (l_head->lstem)->end2 + 1;
					}

					if(node_mapping[l_head->g_id2].node_type != 0) {
						//If the loop doesn't end with the sink. 
						if(gidtosid[l_head->g_id2].left_or_right == 0) {
							r_start = (l_head->rstem)->start1 - 1;
						} else {
							r_start = (l_head->rstem)->start2 - 1;
						}
					} 
					else 
					{
						l_sink = 1;
					}
					if(l_sink == 0)
					{
						//modify this part because overlap in parts of stems should be allowed
						//l_end and r_start all one nts shift, so shift = this->getNumOfNtsOverlapBetweenStem() - 1 + 2
						if(l_end <= r_start + 1 + this->getNumOfNtsOverlapBetweenStem())
	//					if(l_end <= r_start - 1 + this->getNumOfNtsOverlapBetweenStem())
						{
							if(this->getScoreOption() != SCORE_STEM)
							{
								//If there is no overlap occurs, simply take the profile score of the loop.
								loop_end = r_start;
								loop_len = (r_start - l_end) + 1;

								l_score = vs.VTBSearchMax(&pLoopModel[l_head->loop_id], cSequence, l_end, loop_end, NULL, NULL);
								//support weight factor for stem and loop	20081109
								if(l_score != INVLDPROB) {
									//we need to check this condition and 
									//can not JUST simply make l_score = this->getLoopFactor()*l_score;
									//20081109
									l_score = this->getLoopWeight()*l_score;
								} else {
									//if l_score == INVLDPROB, what does that mean ? 20081118
								}
								loop_end += l_end;	//make it the absolute coordinate

								//
								//score += l_score;
								//only consider the non-negative loop score
								if(this->getScoreOption() == SCORE_STEM_POSLOOP)
								{
									if(l_score >= 0.0) {
										score += l_score;
									}
								} else {
									score += l_score;
								}

								//consider the score of parts that is shared with children.
								//save them and substract them when possible.
								for(i=0; i<tr_node->childnum; i++) 
								{
									if(tr_node->pChild[i] != NULL && tr_node->loop_image[i][l_id] == 1) 
									{
										//only consider the non-negative loop score
										if(this->getScoreOption() == SCORE_STEM_POSLOOP)
										{
											if(l_score >= 0.0) {
												tr_node->prob_score[i] += l_score;
											}
										} else {
											tr_node->prob_score[i] += l_score;
										}
									}
								}
							} else {
								score += ZERO;
							}
						}//end if(l_end <= r_start + 1)
						else
						{
							score = SMALLEST;
						}
					}//end if(l_sink == 0)
					else
					{
						//reaching the sink part. --->End
						//do nothing, just skip it
					}
				}
				else {
					//reaching the source part. Start(0)--->
					l_source = 1;
				}
			} else {
				score += ZERO;
			}
			l_id++;
			l_head = l_head->next;
		}
//		cout<<"Score after considering loops = "<<score<<endl;

		//3. Now, start compute the probability score of the node;
		if(score > MYTHRESHOLD) {
			//do full combination  20090914
			int		iEmptyStemNum = 0;
			double	score_backup = score;
			int		current_empty_stem_num_backup = current_empty_stem_num;
			if(tr_node->childnum == 0) {
				score					= score_backup;
				current_empty_stem_num	= current_empty_stem_num_backup;
				bConflict = FALSE;
				for(i=0; i<tr_node->childnum && !bConflict; i++) {
					this->_parentIndex[i] = NEGATIVE_ONE;
					this->_childOptimalIndex[i] = NEGATIVE_ONE;
					child_score = this->extract_profilescore(i, tr_node->pChild[i], tr_node, 
															current_empty_stemid_array, 
															&current_empty_stem_num,
															iEmptyStemNum
															);
					score += child_score - tr_node->prob_score[i];
				}

				if(bConflict == TRUE)
					score = SMALLEST;

/**/
				//change the flowchart  20100215
				//limiting the number of empty-stem. 
				if(this->getFlagSupportEmptyStem() == TRUE)
					if(current_empty_stem_num > this->getMaxEmptyStemNum())
						score = SMALLEST;

				if(getFlagSupportEmptyStem() == TRUE)
				{
					if(score > MYTHRESHOLD)
					{
						//20090402 for debug
						if(this->getFlagPrintDebugDPSearchInfo())
							this->printNeighborStemPos(tr_node, score);
							
						this->checkSingleEmptyStemPos(tr_node, &score);
						this->computeNeighborEmptyStemPos(tr_node, &score);
						this->checkSingleEmptyStemPos(tr_node, &score);

						//20090402 for debug
						if(this->getFlagPrintDebugDPSearchInfo())
							this->printNeighborStemPos(tr_node, score);
							
						this->checkNeighborEmptyStemPos(tr_node, &score);
					}
				}

				if(this->getFlagPrintDebugDPSearchInfo()) {
					if(this->getFlagSupportEmptyStem() == TRUE) {
						cout<<"score for this compnode with empty="<<current_empty_stem_num
							<<" is "
							<<score<<endl;
					} else {
						cout<<"score for this compnode is "
							<<score<<endl;
					}
				}
				if(score > MYTHRESHOLD)
					this->tb_indexingstore(tr_node, tr_node->enum_node, tr_node->nodenum, score, current_empty_stem_num);

				if(this->getFlagSupportEmptyStem() == TRUE) {
					//restore the original enum_node
					this->restoreEnumNodeCombination(tr_node);
				}
			}
			else
			for(iEmptyStemNum=0; iEmptyStemNum<=this->getMaxEmptyStemNum(); iEmptyStemNum++)
			{
				if(tr_node->childnum == 1) {
					score					= score_backup;
					current_empty_stem_num	= current_empty_stem_num_backup;
					bConflict = FALSE;
					for(i=0; i<tr_node->childnum && !bConflict; i++) {
						this->_parentIndex[i] = NEGATIVE_ONE;
						this->_childOptimalIndex[i] = NEGATIVE_ONE;
						child_score = this->extract_profilescore(i, tr_node->pChild[i], tr_node, 
																current_empty_stemid_array, 
																&current_empty_stem_num,
																iEmptyStemNum
																);
						score += child_score - tr_node->prob_score[i];
					}

					if(bConflict == TRUE)
						score = SMALLEST;

/**/
					//change the flowchart  20100215
					//limiting the number of empty-stem. 
					if(this->getFlagSupportEmptyStem() == TRUE)
						if(current_empty_stem_num > this->getMaxEmptyStemNum())
							score = SMALLEST;

					if(getFlagSupportEmptyStem() == TRUE)
					{
						if(score > MYTHRESHOLD)
						{
							//20090402 for debug
							if(this->getFlagPrintDebugDPSearchInfo())
								this->printNeighborStemPos(tr_node, score);
							
							this->checkSingleEmptyStemPos(tr_node, &score);
							this->computeNeighborEmptyStemPos(tr_node, &score);
							this->checkSingleEmptyStemPos(tr_node, &score);

							//20090402 for debug
							if(this->getFlagPrintDebugDPSearchInfo())
								this->printNeighborStemPos(tr_node, score);
							
							this->checkNeighborEmptyStemPos(tr_node, &score);
						}
					}

					if(this->getFlagPrintDebugDPSearchInfo()) {
						if(this->getFlagSupportEmptyStem() == TRUE) {
							cout<<"score for this compnode with empty="<<current_empty_stem_num
								<<" is "
								<<score<<endl;
						} else {
							cout<<"score for this compnode is "
								<<score<<endl;
						}
					}
					if(score > MYTHRESHOLD)
						this->tb_indexingstore(tr_node, tr_node->enum_node, tr_node->nodenum, score, current_empty_stem_num);

					if(this->getFlagSupportEmptyStem() == TRUE) {
						//restore the original enum_node
						this->restoreEnumNodeCombination(tr_node);
					}
				}//if(tr_node->childnum == 1)
				else if(tr_node->childnum == 2) 
				{
					//for(int childEmptyStemNum=0; childEmptyStemNum<=iEmptyStemNum; childEmptyStemNum++)
					for(int childEmptyStemNum=0; childEmptyStemNum<=iEmptyStemNum; childEmptyStemNum++)
					{
						score					= score_backup;
						current_empty_stem_num	= current_empty_stem_num_backup;
						bConflict = FALSE;
						for(i=0; i<tr_node->childnum && !bConflict; i++) {
							this->_parentIndex[i] = NEGATIVE_ONE;
							this->_childOptimalIndex[i] = NEGATIVE_ONE;
							if(i == 0) {
								child_score = this->extract_profilescore(i, tr_node->pChild[i], tr_node, 
																		current_empty_stemid_array, 
																		&current_empty_stem_num,
																		childEmptyStemNum
																		);
							} else {
								child_score = this->extract_profilescore(i, tr_node->pChild[i], tr_node, 
																		current_empty_stemid_array, 
																		&current_empty_stem_num,
																		iEmptyStemNum	//(iEmptyStemNum-childEmptyStemNum)
																		);
							}
							score += child_score - tr_node->prob_score[i];
						}

						if(bConflict == TRUE)
							score = SMALLEST;

/**/
						//change the flowchart  20100215
						//limiting the number of empty-stem. 
						if(this->getFlagSupportEmptyStem() == TRUE)
							if(current_empty_stem_num > this->getMaxEmptyStemNum())
								score = SMALLEST;

						if(getFlagSupportEmptyStem() == TRUE)
						{
							if(score > MYTHRESHOLD)
							{
								//20090402 for debug
								if(this->getFlagPrintDebugDPSearchInfo())
									this->printNeighborStemPos(tr_node, score);
								
								this->checkSingleEmptyStemPos(tr_node, &score);
								this->computeNeighborEmptyStemPos(tr_node, &score);
								this->checkSingleEmptyStemPos(tr_node, &score);

								//20090402 for debug
								if(this->getFlagPrintDebugDPSearchInfo())
									this->printNeighborStemPos(tr_node, score);
								
								this->checkNeighborEmptyStemPos(tr_node, &score);
							}
						}

						if(this->getFlagPrintDebugDPSearchInfo()) {
							if(this->getFlagSupportEmptyStem() == TRUE) {
								cout<<"score for this compnode with empty="<<current_empty_stem_num
									<<" is "
									<<score<<endl;
							} else {
								cout<<"score for this compnode is "
									<<score<<endl;
							}
						}
						if(score > MYTHRESHOLD) {
							this->tb_indexingstore(tr_node, tr_node->enum_node, tr_node->nodenum, score, current_empty_stem_num);
						}

						if(this->getFlagSupportEmptyStem() == TRUE) {
							//restore the original enum_node
							this->restoreEnumNodeCombination(tr_node);
						}
					}//for childEmptyStemNum
				}
			}//for iEmptyStemNum
		}//if
/**/
		if(score == SMALLEST) {
			//mark the invalid combination 20081116
			invalid_combination_mark++;
		}

		int real_candidate_num = 0;
		int k;
		for(i=tr_node->stnum-1; i>=0; i--) {
			tr_node->enum_stem[i]++;
			tail_head = tr_node->su_stem;
			k = i;
			while(k>0) {
				tail_head = tail_head->next;
				k--;
			}
			if(tail_head != NULL)
				real_candidate_num = this->candidateStems[tail_head->stem_id].num;
			if(tr_node->enum_stem[i] < real_candidate_num) {
				break;
			} else {
				tr_node->enum_stem[i] = 0;
			}
		}

		if(this->getFlagPrintDebugDPSearchInfo()) {
			printf("Check tr_node->enum_stem:\t");
			for(i=0; i<tr_node->stnum; i++) {
				printf("%d ", tr_node->enum_stem[i]);
/**/
				if(tr_node->stnum == 4) {
					if(tr_node->enum_stem[0] == 1
						&& tr_node->enum_stem[1] == 10
						&& tr_node->enum_stem[2] == 9) {
						cout<<"hello"<<endl;
					}
				}

			}
			printf("\n");
		}

		sum = 0;

		for(i=0; i<tr_node->stnum; i++) {
			sum += tr_node->enum_stem[i];
		}

		if(sum == 0) {
/*
			//mark the invalid combination 20081111
			//---------------------------------------
			if(invalid_combination_mark == static_cast<int>(pow(static_cast<double>(this->iLeading), 
																static_cast<double>(tr_node->stnum))))
			{
				if(this->getFlagInvalidCombination() == FALSE)
				{
					cout<<"Invalid combination in :";
					l_head = tr_node->su_loop;
					while(l_head)
					{
						if(gidtosid[l_head->g_id1].left_or_right == 0) {
							cout<<this->pStem[(l_head->lstem)->stem_id].charid<<"(L)->";
						} else {
							cout<<this->pStem[(l_head->lstem)->stem_id].charid<<"(R)->";
						}

						if(gidtosid[l_head->g_id2].left_or_right == 0) {
							cout<<this->pStem[(l_head->rstem)->stem_id].charid<<"(L)\t";
						} else {
							cout<<this->pStem[(l_head->rstem)->stem_id].charid<<"(R)\t";
						}
						l_head = l_head->next;
					}
					cout<<endl;
					this->setFlagInvalidCombination(TRUE);	//only output one invalid combination 20081116
				}
			}
*/
			//---------------------------------------
			break;
		}
	}
}

//Extract the information needed to go bottom up from the child node to a parent. 
double DynamicP::extract_profilescore(int childIndex, Tree_bag *child, Tree_bag *parent, 
									  int *parent_empty_stemid_array, 
									  int * parent_empty_stem_num,
									  int	empty_stem_num_bound	//do full combination  20090914
									  )	
{
	int		i;
	int		sum, s_count, iValue;
	double	maxscore, score;
	Strunit		*s_head = NULL;
	Node_list	*n_list = NULL;
	Strunit *	tail_head = NULL;
	maxscore = SMALLEST;
	int optimallChildIndex = 0;

	if(this->getFlagPrintDebugDPSearchInfo())
		cout<<"\t-------- extract_profilescore -------"<<endl;

	//Now we go through the set up list to enumerate all possible combinations.
	for(i=0; i<child->stnum; i++)
    {
		if(child->set_up[i] == -1) {
			//If it is not preset by the parent, set it to be zero. 
			child->enum_stem[i] = 0;
		} else {
			child->enum_stem[i] = parent->enum_node[child->set_up[i]];
		}
	}
//	//for the purpose of debugging
//	cout<<"\t\t";
//	this->print_enum_stem(child);

	//Initialize the aray of enum_node to be zero. 
	for(i=0; i<child->nodenum; i++) {
		child->enum_node[i] = 0;
	}

	int current_empty_stem_num = 0;
	int	current_empty_stemid_array[MAXNUMEMPTYSTEM];

	int	child_empty_stem_num_total = 0;

	//support max empty stem num 
	int	childEmptyStemNum = 0;
	int parent_overall_empty_stem_num = 0;
	int parent_overall_empty_stem_num_opt = 0;

	bool bFound = false;
	//Now start doing the enumeration.
	while(1)
	{
		//limiting the number of empty-stem. 
		if(this->getFlagSupportEmptyStem() == TRUE) {
			current_empty_stem_num = 0;
			for(i=0; i<MAXNUMEMPTYSTEM; i++) {
				current_empty_stemid_array[i] = NEGATIVE_ONE;
			}

			child_empty_stem_num_total = 0;
		}

		//Now based on the generated array, start doing the look up and conversion; 
        s_head = child->su_stem;
        s_count = 0;

        while(s_head) {
			iValue = child->enum_stem[s_count];
			if(s_head->type == 1) {
				//If it is a complete stem.
				child->enum_node[s_head->left_nid]	= iValue;
				child->enum_node[s_head->right_nid]	= iValue;

				//limiting the number of empty-stem.  20090213
				if(this->getFlagSupportEmptyStem() == TRUE) {
					if(iValue == (this->getLeadingNum()-1)) {
						//get this emptystem's stemid
						current_empty_stemid_array[current_empty_stem_num] = this->getStemIdbyLeftNodeId(s_head->g_id1);

						current_empty_stem_num++;
					}
				}
			} else {
				child->enum_node[s_head->left_nid]	= iValue;
			}
			s_count++;
			s_head = s_head->next;
		}

		n_list = child->nhead;

		if(this->getFlagPrintDebugDPSearchInfo()) {
			printf("\tChild node : ");
			this->printTreenode(child);
			cout<<"\t";
			this->print_enum_node(child->enum_node, child->nodenum);
			cout<<endl;
		}

		//support max empty stem num 
		if(getFlagSupportEmptyStem() == TRUE)
		{
			for(childEmptyStemNum=0; childEmptyStemNum <= this->getMaxEmptyStemNum(); childEmptyStemNum++)
			{
				parent_overall_empty_stem_num = * parent_empty_stem_num;
				//compute the overal_empty_stem_num for the parent tree-bag 20090224
				int overlap_empty_stem_num = 0;
				int parent_stemid = NEGATIVE_ONE;
				for(int m=0; m<*parent_empty_stem_num; m++)
				{
					parent_stemid = parent_empty_stemid_array[m];
					for(int n=0; n<current_empty_stem_num; n++)
					{
						if(parent_stemid == current_empty_stemid_array[n])
							overlap_empty_stem_num++;
					}
				}
				parent_overall_empty_stem_num -= overlap_empty_stem_num;
				parent_overall_empty_stem_num += childEmptyStemNum;

				//do full combination  20090914
				if(parent_overall_empty_stem_num == empty_stem_num_bound)
				{
					this->tb_indexingfetch(child, child->enum_node, child->nodenum, childEmptyStemNum, &score);//, &child_empty_stem_num_total);

					if(maxscore < score) {
						bFound = true;

						//Here is the greedy strategy
						maxscore = score;

						//remove the block memory function in dp part.  20090208
						//compute the current child index that cause the maxscore
						int nodenum = child->nodenum;
						optimallChildIndex = 0;
						for(i=0; i<nodenum; i++) {
							int iTmp = child->enum_node[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), 
																				static_cast<double>(nodenum-1-i)));
							optimallChildIndex += iTmp;
						}
						//support max empty stem num 
						optimallChildIndex += childEmptyStemNum * static_cast<int>(pow(static_cast<double>(this->iLeading), 
																						static_cast<double>(nodenum)));

						parent_overall_empty_stem_num_opt = parent_overall_empty_stem_num;
					}
				}
			}//for(childEmptyStemNum=0
		} else {
			this->tb_indexingfetch(child, child->enum_node, child->nodenum, childEmptyStemNum, &score);

			if(this->getFlagPrintDebugDPSearchInfo()) {
				cout<<"\t\tScore of this child node is "<<score<<endl;
			}

			if(maxscore < score) {
				maxscore = score;

				//20090309
				int nodenum = child->nodenum;
				optimallChildIndex = 0;
				for(i=0; i<nodenum; i++) {
					int iTmp = child->enum_node[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), 
																		static_cast<double>(nodenum-1-i)));
					optimallChildIndex += iTmp;
				}
			}
		}
		int real_candidate_num = 0;
		int k=0;
		for(i=child->stnum-1; i>=0; i--) {
			if(child->set_up[i] == -1) {
				child->enum_stem[i]++;
				tail_head = child->su_stem;
				k = i;
				while(k>0) {
					tail_head = tail_head->next;
					k--;
				}
				if(tail_head != NULL)
					real_candidate_num = this->candidateStems[tail_head->stem_id].num;
				if(child->enum_stem[i] < real_candidate_num) {
					break;
				} else {
					child->enum_stem[i] = 0;
				}
			}
		}

		sum = 0;
		for(i=0; i<child->stnum; i++) {
			if(child->set_up[i] == -1) {
				sum += child->enum_stem[i];
			}
		}
		if(sum == 0) break;
	}//while

	if(this->getFlagPrintDebugDPSearchInfo()) {
		cout<<"\tmaxscore = "<<maxscore<<endl;
	}

	//remove the block memory function in dp part.  20090208
	//store the optimalChildIndex to the p_node_enumerate_array
	if(this->getFlagSupportEmptyStem() && bFound) {
		int index = 0;
		if(parent->pParent != parent)	//root node
		{
			int array_num = parent->nodenum;
			for(i=0; i<array_num; i++)
			{
				index += parent->enum_node[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), 
																	static_cast<double>(array_num-1-i)));
			}
		} else {
			index += parent_overall_empty_stem_num_opt;
		}

		//global variable
		this->_parentIndex[childIndex]		= index;
		this->_childOptimalIndex[childIndex]= optimallChildIndex;

		* parent_empty_stem_num = parent_overall_empty_stem_num_opt;

		//access the pos for the empty-stem candidate from its child treebag
		//here we know 
		//currentParentIndex and optimallChildIndex
		//then we need to access the pos for the empty-stem candidate from this child
		this->transferEmptyStemPosToParent(parent, parent_empty_stemid_array, child, optimallChildIndex);
	} else {
		//20090309
		int index = 0;
		if(parent->pParent != parent)	//root node
		{
			int array_num = parent->nodenum;
			for(i=0; i<array_num; i++)
			{
				index += parent->enum_node[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), 
																	static_cast<double>(array_num-1-i)));
			}
		} else {
			//do nothing
		}

		//global variable
		this->_parentIndex[childIndex]		= index;
		this->_childOptimalIndex[childIndex]= optimallChildIndex;
	}
	return maxscore;
}

/**/
//access the pos for the empty-stem candidate from its child treebag
//from parent, we can compute the currentParentIndex
void DynamicP::transferEmptyStemPosToParent(Tree_bag * parent, 
											int * parent_empty_stemid_array, 
											Tree_bag * child, 
											int optimallChildIndex)
{
	int srcNodeId = 0;
	int sinkNodeId = 2*this->getNumstems()+1;

	//calculate the index for the current parent tree-bag combination
	int pNodeNum = parent->nodenum;
	int curParentIndex = 0;
	for(int i=0; i<pNodeNum; i++) {
		curParentIndex += parent->enum_node[i]
						*static_cast<int>(pow(static_cast<double>(this->iLeading), 
												static_cast<double>(pNodeNum-1-i)));
	}

	Node_list * pNodeList = parent->nhead;
	int pNodeId = NEGATIVE_ONE, pStemId = NEGATIVE_ONE, 
		pCandId = NEGATIVE_ONE, pNodeIndex = NEGATIVE_ONE;
	int cNodeId = NEGATIVE_ONE, cStemId = NEGATIVE_ONE, 
		cCandId = NEGATIVE_ONE, cNodeIndex = NEGATIVE_ONE;

	int lpos, rpos;
	while(pNodeList)
	{
		pNodeIndex++;
		pNodeId = pNodeList->g_node;
		if(pNodeId != NEGATIVE_ONE
		&& pNodeId != srcNodeId 
		&& pNodeId != sinkNodeId)
		{
			pStemId = gidtosid[pNodeId].s_id;
			pCandId = parent->enum_node[pNodeIndex];
			if(pCandId == (this->getLeadingNum()-1))
			{
				//check whether we can access this empty-stem pos from current child
				cNodeId = pNodeId;
				cNodeIndex = getNodeIndexInTreenode(child->nhead, cNodeId);
				if(cNodeIndex != NEGATIVE_ONE)
				{
					cStemId = gidtosid[cNodeId].s_id;
					cCandId = child->enum_node[cNodeIndex];
					if(cCandId == (this->getLeadingNum()-1))
					{
						lpos = child->t_head[optimallChildIndex].lpos[cNodeIndex];
						rpos = child->t_head[optimallChildIndex].rpos[cNodeIndex];
						
						//parent->t_head[curParentIndex].lpos[pNodeIndex] = lpos;
						//parent->t_head[curParentIndex].rpos[pNodeIndex] = rpos;
						//add extend_dir 
						//it seems this info has not been used ???
/*
						cout<<parent->t_head[curParentIndex].l_ext_dir[pNodeIndex]
							<<parent->t_head[curParentIndex].r_ext_dir[pNodeIndex]
							<<endl;
*/
						//Here there is a potential error here, if something confliction here
						//We need to check whether there is pos confliction before replace it.  20090227
						if(gidtosid[pNodeId].left_or_right == 0)	//left part
						{
							if(this->candidateStems[pStemId].s_locinfo[pCandId].a == NEGATIVE_ONE) {
								this->candidateStems[pStemId].s_locinfo[pCandId].a = lpos;

								//add extend_dir 
								this->candidateStems[pStemId].s_locinfo[pCandId].a_ext_dir
									= child->t_head[optimallChildIndex].l_ext_dir[cNodeIndex];
							}
							else	//20090305	//20090403	//20090411
							{
								//if(this->candidateStems[pStemId].s_locinfo[pCandId].a > lpos)
								//add extend_dir 
								//debug  20090914
								if(lpos != NEGATIVE_ONE)
								{
									if((this->candidateStems[pStemId].s_locinfo[pCandId].a > lpos
										&& this->candidateStems[pStemId].s_locinfo[pCandId].a_ext_dir == EXTRIGHT)
									|| (this->candidateStems[pStemId].s_locinfo[pCandId].a < lpos
										&& this->candidateStems[pStemId].s_locinfo[pCandId].a_ext_dir == EXTLEFT))
										bConflict = TRUE;
								}
							}

							if(this->candidateStems[pStemId].s_locinfo[pCandId].b == NEGATIVE_ONE) {
								this->candidateStems[pStemId].s_locinfo[pCandId].b = rpos;

								//add extend_dir 
								this->candidateStems[pStemId].s_locinfo[pCandId].b_ext_dir
									= child->t_head[optimallChildIndex].r_ext_dir[cNodeIndex];
							}
							else	//20090305	//20090411
								//if(this->candidateStems[pStemId].s_locinfo[pCandId].b > rpos)
								//add extend_dir 
								//debug  20090914
								if(rpos != NEGATIVE_ONE)
								{
									if((this->candidateStems[pStemId].s_locinfo[pCandId].b > rpos
										&& this->candidateStems[pStemId].s_locinfo[pCandId].b_ext_dir == EXTRIGHT)
									||	(this->candidateStems[pStemId].s_locinfo[pCandId].b < rpos
										&& this->candidateStems[pStemId].s_locinfo[pCandId].b_ext_dir == EXTLEFT))
										bConflict = TRUE;
								}

						} else {	//right part
							if(this->candidateStems[pStemId].s_locinfo[pCandId].c == NEGATIVE_ONE) {
								this->candidateStems[pStemId].s_locinfo[pCandId].c = lpos;

								//add extend_dir 
								//parent->t_head[curParentIndex].l_ext_dir[pNodeIndex]
								this->candidateStems[pStemId].s_locinfo[pCandId].c_ext_dir
									= child->t_head[optimallChildIndex].l_ext_dir[cNodeIndex];
							}
							else	//20090305	//20090403	//20090411
							{
								//if(this->candidateStems[pStemId].s_locinfo[pCandId].c > lpos)
								//add extend_dir 
								//debug  20090914
								if(lpos != NEGATIVE_ONE)
								{
									if((this->candidateStems[pStemId].s_locinfo[pCandId].c > lpos
										&& this->candidateStems[pStemId].s_locinfo[pCandId].c_ext_dir == EXTRIGHT)
									|| (this->candidateStems[pStemId].s_locinfo[pCandId].c < lpos
										&& this->candidateStems[pStemId].s_locinfo[pCandId].c_ext_dir == EXTLEFT))
										bConflict = TRUE;
								}
							}

							if(this->candidateStems[pStemId].s_locinfo[pCandId].d == NEGATIVE_ONE) {
								this->candidateStems[pStemId].s_locinfo[pCandId].d = rpos;

								//add extend_dir 
								this->candidateStems[pStemId].s_locinfo[pCandId].d_ext_dir
									= child->t_head[optimallChildIndex].r_ext_dir[cNodeIndex];
							}
							else	//20090305	//20090411
							{
								//if(this->candidateStems[pStemId].s_locinfo[pCandId].d > rpos)
								//add extend_dir 
								//debug  20090914
								if(rpos != NEGATIVE_ONE)
								{
									if((this->candidateStems[pStemId].s_locinfo[pCandId].d > rpos
										&& this->candidateStems[pStemId].s_locinfo[pCandId].d_ext_dir == EXTRIGHT)
									|| (this->candidateStems[pStemId].s_locinfo[pCandId].d < rpos
										&& this->candidateStems[pStemId].s_locinfo[pCandId].d_ext_dir == EXTLEFT))
										bConflict = TRUE;
								}
							}
						}
					} else {
						cout<<"Something wrong in transferEmptyStemPosToParent!"<<endl;
						exit(0);
					}
				}//if(cNodeIndex != NEGATIVE_ONE)
			}//if(pCandId == (this->getLeadingNum()-1))
		}
		pNodeList = pNodeList->next;
	}//while(pNodeList)
}

//modification by using direction
// 20100312
void DynamicP::checkSingleEmptyStemPos(Tree_bag * tr_node, double * score)
{
	int curNodeId = NEGATIVE_ONE;
	int curNodePosL = NEGATIVE_ONE, curNodePosR = NEGATIVE_ONE;
	char curNodeExtDirL, curNodeExtDirR;	// 20100312
	int curStemId = NEGATIVE_ONE;
	int curNodeIndex = NEGATIVE_ONE;
	int curCandId = NEGATIVE_ONE;
	int srcNodeId = 0;
	int sinkNodeId = 2*this->getNumstems()+1;

	Node_list * node_list;
	//check the position confliction for the two arms of the stem in the treebag
	//note: here "nodeid" is not equal to "nodeindex"
	//"nodeid" means id label for the node
	//"nodeindex" means the index of the nodeid
	if((this->getFlagSupportEmptyStem() == TRUE) && (*score > MYTHRESHOLD)) {
		node_list = tr_node->nhead;
		while(node_list)
		{
			curNodeIndex++;
			curNodeId = node_list->g_node;
			if(curNodeId != NEGATIVE_ONE
				&& curNodeId != srcNodeId 
				&& curNodeId != sinkNodeId)
			{
				curStemId = gidtosid[curNodeId].s_id;
				curCandId = tr_node->enum_node[curNodeIndex];

				//compute curNode's pos
				if(gidtosid[curNodeId].left_or_right == 0)	//left part
				{
					curNodePosL		= this->candidateStems[curStemId].s_locinfo[curCandId].a;
					curNodeExtDirL	= this->candidateStems[curStemId].s_locinfo[curCandId].a_ext_dir;
					curNodePosR		= this->candidateStems[curStemId].s_locinfo[curCandId].b;
					curNodeExtDirR	= this->candidateStems[curStemId].s_locinfo[curCandId].b_ext_dir;
				} else {
					curNodePosL		= this->candidateStems[curStemId].s_locinfo[curCandId].c;
					curNodeExtDirL	= this->candidateStems[curStemId].s_locinfo[curCandId].c_ext_dir;
					curNodePosR		= this->candidateStems[curStemId].s_locinfo[curCandId].d;
					curNodeExtDirR	= this->candidateStems[curStemId].s_locinfo[curCandId].d_ext_dir;
				}

				if((curNodePosL != NEGATIVE_ONE) && (curNodePosR != NEGATIVE_ONE)) {
					//check whether there is position confliction here
//					if(curNodePosL <= curNodePosR) {
					if(curNodePosL <= curNodePosR) {
						if((curNodeExtDirL == EXTLEFT) && (curNodeExtDirR == EXTRIGHT)) {
							*score = SMALLEST;
							break;
						} else if((curNodeExtDirL == EXTRIGHT) && (curNodeExtDirR == EXTRIGHT)) {
							// 20100323
							if(gidtosid[curNodeId].left_or_right == 0)	//left part
							{
								this->candidateStems[curStemId].s_locinfo[curCandId].a			= curNodePosR;
								this->candidateStems[curStemId].s_locinfo[curCandId].a_ext_dir	= EXTRIGHT;
							} else {
								this->candidateStems[curStemId].s_locinfo[curCandId].c			= curNodePosR;
								this->candidateStems[curStemId].s_locinfo[curCandId].c_ext_dir	= EXTRIGHT;
							}
						} else if((curNodeExtDirL == EXTLEFT) && (curNodeExtDirR == EXTLEFT)) {
							// 20100323
							if(gidtosid[curNodeId].left_or_right == 0)	//left part
							{
								this->candidateStems[curStemId].s_locinfo[curCandId].b			= curNodePosL;
								this->candidateStems[curStemId].s_locinfo[curCandId].b_ext_dir	= EXTLEFT;
							} else {
								this->candidateStems[curStemId].s_locinfo[curCandId].d			= curNodePosL;
								this->candidateStems[curStemId].s_locinfo[curCandId].d_ext_dir	= EXTLEFT;
							}
						} else {
							//do nothing
						}
					} else {
						//position confliction
						*score = SMALLEST;
						break;
					}
				}
			}
			node_list = node_list->next;
		}
	}
}

void DynamicP::computeNeighborEmptyStemPos(Tree_bag * tr_node, double * score)
{
	int preNodeId = NEGATIVE_ONE, curNodeId = NEGATIVE_ONE;
	int preNodePos = NEGATIVE_ONE, curNodePos = NEGATIVE_ONE;
	int preStemId = NEGATIVE_ONE, curStemId = NEGATIVE_ONE;
	int preCandId = NEGATIVE_ONE, curCandId = NEGATIVE_ONE;
	int	preNodeIndex = NEGATIVE_ONE, curNodeIndex = NEGATIVE_ONE;

	int srcNodeId = 0;
	int sinkNodeId = 2*this->getNumstems()+1;
	Node_list * node_list;
	int lpos, rpos;

	//compute the pos for the empty-stem candidate based on pos of other nodes in the same treebag
	if((this->getFlagSupportEmptyStem() == TRUE) && (*score > MYTHRESHOLD)) 
	{
		node_list = tr_node->nhead;
		while(node_list)
		{
			curNodeIndex++;
			curNodeId = node_list->g_node;
			if(curNodeId != NEGATIVE_ONE
				&& curNodeId != srcNodeId 
				&& curNodeId != sinkNodeId)
			{
				curStemId = gidtosid[curNodeId].s_id;
				curCandId = tr_node->enum_node[curNodeIndex];
				if(preNodeId != NEGATIVE_ONE
					&& preNodeId != srcNodeId 
					&& preNodeId != sinkNodeId)
				{
					if((preCandId == (this->getLeadingNum()-1) && curCandId != (this->getLeadingNum()-1))
					|| (preCandId != (this->getLeadingNum()-1) && curCandId == (this->getLeadingNum()-1)))
					{
						if(preCandId == (this->getLeadingNum()-1) && curCandId != (this->getLeadingNum()-1))
						{
							//compute curNode's pos
							if(gidtosid[curNodeId].left_or_right == 0)	//left part
								curNodePos = this->candidateStems[curStemId].s_locinfo[curCandId].a;
							else
								curNodePos = this->candidateStems[curStemId].s_locinfo[curCandId].c;
							//
							if(gidtosid[preNodeId].left_or_right == 0)	//left part
							{
								rpos = this->candidateStems[preStemId].s_locinfo[preCandId].b;
								if(rpos == NEGATIVE_ONE)
								{
									//rpos = curNodePos-1+this->getNumOfNtsOverlapBetweenStem();
									rpos = curNodePos+this->getNumOfNtsOverlapBetweenStem();
									this->candidateStems[preStemId].s_locinfo[preCandId].b = rpos;

									//add extend_dir 
									this->candidateStems[preStemId].s_locinfo[preCandId].b_ext_dir = EXTLEFT;
/*
									lpos = this->candidateStems[preStemId].s_locinfo[preCandId].a;
									if(lpos != NEGATIVE_ONE)
										//if(lpos > rpos)	//check position confliction
										if(lpos >= rpos)
 											* score = SMALLEST;
*/
								} else {
									//handle the case of node inserted into the neighbor of empty-stem 
									if((curNodePos < rpos)
									&& (this->candidateStems[preStemId].s_locinfo[preCandId].b_ext_dir == EXTLEFT))
										this->candidateStems[preStemId].s_locinfo[preCandId].b = curNodePos;
								}
							} else {
								rpos = this->candidateStems[preStemId].s_locinfo[preCandId].d;
								if(rpos == NEGATIVE_ONE)
								{
									//rpos = curNodePos-1+this->getNumOfNtsOverlapBetweenStem();
									rpos = curNodePos+this->getNumOfNtsOverlapBetweenStem();
									this->candidateStems[preStemId].s_locinfo[preCandId].d = rpos;

									//add extend_dir 
									this->candidateStems[preStemId].s_locinfo[preCandId].d_ext_dir = EXTLEFT;
/*
									lpos = this->candidateStems[preStemId].s_locinfo[preCandId].c;
									if(lpos != NEGATIVE_ONE)
										//if(lpos > rpos)	//check position confliction
										if(lpos >= rpos)
 											* score = SMALLEST;
*/
								} else {
									//handle the case of node inserted into the neighbor of empty-stem 
									if((curNodePos < rpos)
									&& (this->candidateStems[preStemId].s_locinfo[preCandId].d_ext_dir == EXTLEFT))
										this->candidateStems[preStemId].s_locinfo[preCandId].d = curNodePos;
								}
							}
						}//if(preCandId == (this->getLeadingNum()-1) && curCandId != (this->getLeadingNum()-1))

						if(preCandId != (this->getLeadingNum()-1) && curCandId == (this->getLeadingNum()-1))
						{
							//compute preNode's pos
							if(gidtosid[preNodeId].left_or_right == 0)	//left part
								preNodePos = this->candidateStems[preStemId].s_locinfo[preCandId].b;
							else
								preNodePos = this->candidateStems[preStemId].s_locinfo[preCandId].d;
							//
							if(gidtosid[curNodeId].left_or_right == 0)	//left part
							{
								lpos = this->candidateStems[curStemId].s_locinfo[curCandId].a;
								if(lpos == NEGATIVE_ONE)
								{
									//lpos = preNodePos+1-this->getNumOfNtsOverlapBetweenStem();
									lpos = preNodePos-this->getNumOfNtsOverlapBetweenStem();
									this->candidateStems[curStemId].s_locinfo[curCandId].a = lpos;

									//add extend_dir 
									this->candidateStems[curStemId].s_locinfo[curCandId].a_ext_dir = EXTRIGHT;
/*
									rpos = this->candidateStems[curStemId].s_locinfo[curCandId].b;
									if(rpos != NEGATIVE_ONE)
										if(lpos >= rpos)
 											* score = SMALLEST;
*/ 
								} else {
									//handle the case of node inserted into the neighbor of empty-stem 
									if((lpos < preNodePos)
									&& (this->candidateStems[curStemId].s_locinfo[curCandId].a_ext_dir == EXTRIGHT))
										this->candidateStems[curStemId].s_locinfo[curCandId].a = preNodePos;
								}
							} else {
								lpos = this->candidateStems[curStemId].s_locinfo[curCandId].c;
								if(lpos == NEGATIVE_ONE)
								{
									//lpos = preNodePos+1-this->getNumOfNtsOverlapBetweenStem();
									lpos = preNodePos-this->getNumOfNtsOverlapBetweenStem();
									this->candidateStems[curStemId].s_locinfo[curCandId].c = lpos;

									//add extend_dir 
									this->candidateStems[curStemId].s_locinfo[curCandId].c_ext_dir = EXTRIGHT;
/*
									rpos = this->candidateStems[curStemId].s_locinfo[curCandId].d;
									if(rpos != NEGATIVE_ONE)
										if(lpos >= rpos)
 											* score = SMALLEST;
*/ 
								} else {
									//handle the case of node inserted into the neighbor of empty-stem 
									if((lpos < preNodePos)
									&& (this->candidateStems[curStemId].s_locinfo[curCandId].c_ext_dir == EXTRIGHT))
										this->candidateStems[curStemId].s_locinfo[curCandId].c = preNodePos;
								}
							}
						}//if(preCandId != (this->getLeadingNum()-1) && curCandId == (this->getLeadingNum()-1))
					}//one of them is empty stem candidate

					else if(preCandId == (this->getLeadingNum()-1) && curCandId == (this->getLeadingNum()-1))
					{
						//case of two continuous empty stems
						int preNodePosL = NEGATIVE_ONE, preNodePosR = NEGATIVE_ONE;
						int curNodePosL = NEGATIVE_ONE, curNodePosR = NEGATIVE_ONE;

						//handle the pos confliction in the continuous two empty-stems 
						char preNodeLExtDir, preNodeRExtDir;
						char curNodeLExtDir, curNodeRExtDir;

						//compute preNode's pos
						if(gidtosid[preNodeId].left_or_right == 0)	//left part
						{
							preNodePosL = this->candidateStems[preStemId].s_locinfo[preCandId].a;
							preNodePosR = this->candidateStems[preStemId].s_locinfo[preCandId].b;

							//handle the pos confliction in the continuous two empty-stems 
							preNodeLExtDir	= this->candidateStems[preStemId].s_locinfo[preCandId].a_ext_dir;
							preNodeRExtDir	= this->candidateStems[preStemId].s_locinfo[preCandId].b_ext_dir;
						} else {
							preNodePosL = this->candidateStems[preStemId].s_locinfo[preCandId].c;
							preNodePosR = this->candidateStems[preStemId].s_locinfo[preCandId].d;

							//handle the pos confliction in the continuous two empty-stems 
							preNodeLExtDir = this->candidateStems[preStemId].s_locinfo[preCandId].c_ext_dir;
							preNodeRExtDir = this->candidateStems[preStemId].s_locinfo[preCandId].d_ext_dir;
						}

						//compute curNode's pos
						if(gidtosid[curNodeId].left_or_right == 0)	//left part
						{
							curNodePosL = this->candidateStems[curStemId].s_locinfo[curCandId].a;
							curNodePosR = this->candidateStems[curStemId].s_locinfo[curCandId].b;

							//handle the pos confliction in the continuous two empty-stems 
							curNodeLExtDir = this->candidateStems[curStemId].s_locinfo[curCandId].a_ext_dir;
							curNodeRExtDir = this->candidateStems[curStemId].s_locinfo[curCandId].b_ext_dir;
						} else {
							curNodePosL = this->candidateStems[curStemId].s_locinfo[curCandId].c;
							curNodePosR = this->candidateStems[curStemId].s_locinfo[curCandId].d;

							//handle the pos confliction in the continuous two empty-stems 
							curNodeLExtDir = this->candidateStems[curStemId].s_locinfo[curCandId].c_ext_dir;
							curNodeRExtDir = this->candidateStems[curStemId].s_locinfo[curCandId].d_ext_dir;
						}

						if(preNodePosL != NEGATIVE_ONE && preNodePosR != NEGATIVE_ONE
						&& curNodePosL != NEGATIVE_ONE && curNodePosR != NEGATIVE_ONE)
						{
							if((preNodePosL >= preNodePosR) || (curNodePosL >= curNodePosR))
								*score = SMALLEST;
							else {
								if(preNodePosR >= curNodePosL)
									*score = SMALLEST;
							}
						}
						else 
						{	//				 <-
							//(-1, -1) - (-1, *)  ->  (*, -1) - (-1, *)
							if(preNodePosL == NEGATIVE_ONE && preNodePosR == NEGATIVE_ONE
							&& curNodePosL == NEGATIVE_ONE && curNodePosR != NEGATIVE_ONE)
							{
								// 20100215
								if(curNodeRExtDir == EXTLEFT) {
									//preNode
									if(gidtosid[preNodeId].left_or_right == 0)	//left part
									{
										this->candidateStems[preStemId].s_locinfo[preCandId].a = curNodePosR;

										//add extend_dir 
										this->candidateStems[preStemId].s_locinfo[preCandId].a_ext_dir = EXTLEFT;
									} else {
										this->candidateStems[preStemId].s_locinfo[preCandId].c = curNodePosR;

										//add extend_dir 
										this->candidateStems[preStemId].s_locinfo[preCandId].c_ext_dir = EXTLEFT;
									}
								}
							}
							//			 <-
							//(-1, -1) - (*, -1)  ->  (*, -1) - (*, -1)
							else if(preNodePosL == NEGATIVE_ONE && preNodePosR == NEGATIVE_ONE
							&& curNodePosL != NEGATIVE_ONE && curNodePosR == NEGATIVE_ONE)
							{
								// 20100215
								if(curNodeLExtDir == EXTLEFT) {
									//preNode
									if(gidtosid[preNodeId].left_or_right == 0)	//left part
									{
										this->candidateStems[preStemId].s_locinfo[preCandId].a = curNodePosL;

										//add extend_dir 
										this->candidateStems[preStemId].s_locinfo[preCandId].a_ext_dir = EXTLEFT;
									} else {
										this->candidateStems[preStemId].s_locinfo[preCandId].c = curNodePosL;

										//add extend_dir 
										this->candidateStems[preStemId].s_locinfo[preCandId].c_ext_dir = EXTLEFT;
									}
								}
							}
							//	  ->
							//(-1, *) - (-1, -1)  ->  (-1, *) - (-1, *)
							else if(preNodePosL == NEGATIVE_ONE && preNodePosR != NEGATIVE_ONE
							&& curNodePosL == NEGATIVE_ONE && curNodePosR == NEGATIVE_ONE)
							{
								if(preNodeRExtDir == EXTRIGHT) {
									//curNode
									if(gidtosid[curNodeId].left_or_right == 0)	//left part
									{
										this->candidateStems[curStemId].s_locinfo[curCandId].b = preNodePosR;

										//add extend_dir 
										this->candidateStems[curStemId].s_locinfo[curCandId].b_ext_dir = EXTRIGHT;
									} else {
										this->candidateStems[curStemId].s_locinfo[curCandId].d = preNodePosR;

										//add extend_dir 
										this->candidateStems[curStemId].s_locinfo[curCandId].d_ext_dir = EXTRIGHT;
									}
								}
							}
							//->
							//(*, -1) - (-1, -1)  ->  (*, -1) - (-1, *)
							else if(preNodePosL != NEGATIVE_ONE && preNodePosR == NEGATIVE_ONE
							&& curNodePosL == NEGATIVE_ONE && curNodePosR == NEGATIVE_ONE)
							{
								if(preNodeLExtDir == EXTRIGHT) {
									//curNode
									if(gidtosid[curNodeId].left_or_right == 0)	//left part
									{
										this->candidateStems[curStemId].s_locinfo[curCandId].b = preNodePosL;

										//add extend_dir 
										this->candidateStems[curStemId].s_locinfo[curCandId].b_ext_dir = EXTRIGHT;
									} else {
										this->candidateStems[curStemId].s_locinfo[curCandId].d = preNodePosL;

										//add extend_dir 
										this->candidateStems[curStemId].s_locinfo[curCandId].d_ext_dir = EXTRIGHT;
									}
								}
							}

							//handle the pos confliction in the continuous two empty-stems 
							//(*, *) - (-1, -1)  -> 
							else if(preNodePosL != NEGATIVE_ONE && preNodePosR != NEGATIVE_ONE
								&&  curNodePosL == NEGATIVE_ONE && curNodePosR == NEGATIVE_ONE)
							{
								if(preNodeLExtDir == EXTRIGHT && preNodeRExtDir == EXTLEFT)
								{
									//									 <-
									//(*R, *L) - (-1, -1) -> (*R, *L) - (*L, -1)
									if(curNodePosL == NEGATIVE_ONE)
									{
										//curNode
										if(gidtosid[curNodeId].left_or_right == 0)	//left part
										{
											this->candidateStems[curStemId].s_locinfo[curCandId].a = preNodePosR;
											this->candidateStems[curStemId].s_locinfo[curCandId].a_ext_dir = EXTLEFT;
										} else {
											this->candidateStems[curStemId].s_locinfo[curCandId].c = preNodePosR;
											this->candidateStems[curStemId].s_locinfo[curCandId].c_ext_dir = EXTLEFT;
										}
									}//if(preNodePosR == NEGATIVE_ONE)
								}
							}
							//handle the pos confliction in the continuous two empty-stems 
							//(-1, -1) - (*, *)  -> 
							else if(preNodePosL == NEGATIVE_ONE && preNodePosR == NEGATIVE_ONE
								&&	curNodePosL != NEGATIVE_ONE && curNodePosR != NEGATIVE_ONE)
							{
								if(curNodeLExtDir == EXTRIGHT && curNodeRExtDir == EXTLEFT)
								{
									//							  <-
									//(-1, -1) - (*R, *L) -> (-1, *R) - (*R, *L)
									if(preNodePosR == NEGATIVE_ONE)
									{
										//preNode
										if(gidtosid[preNodeId].left_or_right == 0)	//left part
										{
											this->candidateStems[preStemId].s_locinfo[preCandId].b = curNodePosL;
											this->candidateStems[preStemId].s_locinfo[preCandId].b_ext_dir = EXTRIGHT;
										} else {
											this->candidateStems[preStemId].s_locinfo[preCandId].d = curNodePosL;
											this->candidateStems[preStemId].s_locinfo[preCandId].d_ext_dir = EXTRIGHT;
										}
									}//if(preNodePosR == NEGATIVE_ONE)
								}
								else if(curNodeLExtDir == EXTLEFT && curNodeRExtDir == EXTLEFT)	// 20100312
								{
									//							  <-
									//(-1, -1) - (*L, *L) -> (-1, *L) - (*L, *L)
									if(preNodePosR == NEGATIVE_ONE)
									{
										//preNode
										if(gidtosid[preNodeId].left_or_right == 0)	//left part
										{
											this->candidateStems[preStemId].s_locinfo[preCandId].b = curNodePosL;
											this->candidateStems[preStemId].s_locinfo[preCandId].b_ext_dir = EXTLEFT;
										} else {
											this->candidateStems[preStemId].s_locinfo[preCandId].d = curNodePosL;
											this->candidateStems[preStemId].s_locinfo[preCandId].d_ext_dir = EXTLEFT;
										}
									}//if(preNodePosR == NEGATIVE_ONE)
								}
							}
							//handle the pos confliction in the continuous two empty-stems  20100325
							//(*, -1) - (*, -1)  -> 
							else if(preNodePosL != NEGATIVE_ONE && preNodePosR == NEGATIVE_ONE
								&&	curNodePosL != NEGATIVE_ONE && curNodePosR == NEGATIVE_ONE)
							{
								if(preNodeLExtDir == EXTRIGHT && curNodeLExtDir == EXTLEFT)
								{
									//							  <-
									//(*R, -1) - (*L, -1) -> (*R, *L) - (*L, -1)
									if(preNodePosR == NEGATIVE_ONE)
									{
										//preNode
										if(gidtosid[preNodeId].left_or_right == 0)	//left part
										{
											this->candidateStems[preStemId].s_locinfo[preCandId].b = curNodePosL;
											this->candidateStems[preStemId].s_locinfo[preCandId].b_ext_dir = EXTLEFT;
										} else {
											this->candidateStems[preStemId].s_locinfo[preCandId].d = curNodePosL;
											this->candidateStems[preStemId].s_locinfo[preCandId].d_ext_dir = EXTLEFT;
										}
									}//if(preNodePosR == NEGATIVE_ONE)
								}
							}
							// 20100209
							//(*, -1) - (*, *)  bring the <- from the right side if possible
							else if(preNodePosL != NEGATIVE_ONE && preNodePosR == NEGATIVE_ONE
								&&	curNodePosL != NEGATIVE_ONE && curNodePosR != NEGATIVE_ONE)
							{
								if(curNodeRExtDir == EXTLEFT)// && curNodeRExtDir == EXTLEFT)
								{
									//						    <-
									//(*, -1) - (*R, *L) -> (*, *L) - (*R, *L)
									//preNode
									if(gidtosid[preNodeId].left_or_right == 0)	//left part
									{
										this->candidateStems[preStemId].s_locinfo[preCandId].b = curNodePosR;
										this->candidateStems[preStemId].s_locinfo[preCandId].b_ext_dir = EXTLEFT;
									} else {
										this->candidateStems[preStemId].s_locinfo[preCandId].d = curNodePosR;
										this->candidateStems[preStemId].s_locinfo[preCandId].d_ext_dir = EXTLEFT;
									}
								}
							}
						}////two continuous empty stem candidate case
					}

				}//if(preNodeId != NEGATIVE_ONE && preNodeId != srcNodeId && preNodeId != sinkNodeId)
			}
			preNodeId = curNodeId;
			preStemId = curStemId;
			preCandId = curCandId;
			preNodeIndex = curNodeIndex;

			node_list = node_list->next;
		}
	}
}

//20090402 for debug
void DynamicP::printNeighborStemPos(Tree_bag * tr_node, double score)
{
	int curNodeId = NEGATIVE_ONE;
	int curNodePosL = NEGATIVE_ONE, curNodePosR = NEGATIVE_ONE;
	int curStemId = NEGATIVE_ONE;
	int curNodeIndex = NEGATIVE_ONE;
	int curCandId = NEGATIVE_ONE;
	int srcNodeId = 0;
	int sinkNodeId = 2*this->getNumstems()+1;

	Node_list * node_list;
	//check the position confliction for the two arms of the stem in the treebag
	//note: here "nodeid" is not equal to "nodeindex"
	//"nodeid" means id label for the node
	//"nodeindex" means the index of the nodeid
	if((this->getFlagSupportEmptyStem() == TRUE) && (score > MYTHRESHOLD)) {
		node_list = tr_node->nhead;
		while(node_list)
		{
			curNodeIndex++;
			curNodeId = node_list->g_node;
			if(curNodeId != NEGATIVE_ONE
				&& curNodeId != srcNodeId 
				&& curNodeId != sinkNodeId)
			{
				curStemId = gidtosid[curNodeId].s_id;
				curCandId = tr_node->enum_node[curNodeIndex];

				//compute curNode's pos
				if(gidtosid[curNodeId].left_or_right == 0)	//left part
				{
					curNodePosL = this->candidateStems[curStemId].s_locinfo[curCandId].a;
					curNodePosR = this->candidateStems[curStemId].s_locinfo[curCandId].b;
				} else {
					curNodePosL = this->candidateStems[curStemId].s_locinfo[curCandId].c;
					curNodePosR = this->candidateStems[curStemId].s_locinfo[curCandId].d;
				}

				cout<<"("<<curNodePosL<<","<<curNodePosR<<") ";

			}
			node_list = node_list->next;
		}
	}
	cout<<endl;
}

void DynamicP::checkNeighborEmptyStemPos(Tree_bag * tr_node, double * score)
{
	int preNodeId = NEGATIVE_ONE, curNodeId = NEGATIVE_ONE;
	int preNodePos = NEGATIVE_ONE, curNodePos = NEGATIVE_ONE;
	int preStemId = NEGATIVE_ONE, curStemId = NEGATIVE_ONE;
	int preNodeIndex = NEGATIVE_ONE, curNodeIndex = NEGATIVE_ONE;
	int preCandId = NEGATIVE_ONE, curCandId = NEGATIVE_ONE;
	int srcNodeId = 0;
	int sinkNodeId = 2*this->getNumstems()+1;

	Node_list * node_list;
	//check the position confliction in the same treebag even they are not neighbors
	//note: here "nodeid" is not equal to "nodeindex"
	//"nodeid" means id label for the node
	//"nodeindex" means the index of the nodeid
	if((this->getFlagSupportEmptyStem() == TRUE) && (*score > MYTHRESHOLD)) {
		node_list = tr_node->nhead;
		while(node_list)
		{
			curNodeIndex++;
			curNodeId = node_list->g_node;
			if(curNodeId != NEGATIVE_ONE
				&& curNodeId != srcNodeId 
				&& curNodeId != sinkNodeId)
			{
				curStemId = gidtosid[curNodeId].s_id;
				curCandId = tr_node->enum_node[curNodeIndex];
				if(preNodeId != NEGATIVE_ONE)
				{
					if(preNodeId != srcNodeId && preNodeId != sinkNodeId)
					{
						{
							//compute curNode's pos
							if(gidtosid[curNodeId].left_or_right == 0)	//left part
							{
								curNodePos = this->candidateStems[curStemId].s_locinfo[curCandId].a;
								if(curNodePos == NEGATIVE_ONE)
									curNodePos = this->candidateStems[curStemId].s_locinfo[curCandId].b;
							} else {
								curNodePos = this->candidateStems[curStemId].s_locinfo[curCandId].c;
								if(curNodePos == NEGATIVE_ONE)
									curNodePos = this->candidateStems[curStemId].s_locinfo[curCandId].d;
							}

							//compute preNode's pos
							if(gidtosid[preNodeId].left_or_right == 0)	//left part
							{
								preNodePos = this->candidateStems[preStemId].s_locinfo[preCandId].b;
								if(preNodePos == NEGATIVE_ONE)
									preNodePos = this->candidateStems[preStemId].s_locinfo[preCandId].a;
							} else {
								preNodePos = this->candidateStems[preStemId].s_locinfo[preCandId].d;
								if(preNodePos == NEGATIVE_ONE)
									preNodePos = this->candidateStems[preStemId].s_locinfo[preCandId].c;
							}

							//check whether there is position confliction here
							if(preNodePos <= curNodePos + this->getNumOfNtsOverlapBetweenStem()) {
								//currently do nothing
							} else {
								//position confliction
								*score = SMALLEST;
								break;
							}
						}
					}
				}
			}
			preNodeId = curNodeId;
			preStemId = curStemId;
			preCandId = curCandId;
			preNodeIndex = curNodeIndex;

			node_list = node_list->next;
		}
	}
}

int	DynamicP::getNodeIndexInTreenode(Node_list *n_head, int nodeid)
{
	Node_list *nhead;
	nhead = n_head;
	int index = NEGATIVE_ONE;
	bool bFound = false;
	while(nhead && !bFound) {
		index++;
		if(nhead->g_node == nodeid) {
			bFound = true;
		}
		nhead = nhead->next;
	}
	if(!bFound)
		return NEGATIVE_ONE;
	else
		return index; 
}

double DynamicP::getMaxProfilescore(Tree_bag *root)
{
	double	dOptimalScore = SMALLEST;
	dOptimalScore = root->t_head[0].score;
	return dOptimalScore;
}

double DynamicP::getMaxProfilescore2(Tree_bag *root)
{
	double	dOptimalScore = SMALLEST;
	if(getFlagSupportEmptyStem() == TRUE)
	{
		for(int childEmptyStemNum=0; childEmptyStemNum <= this->getMaxEmptyStemNum(); childEmptyStemNum++)
		{
			if(root->t_head[childEmptyStemNum].optimal == 1)
			{
				double dTmp = root->t_head[childEmptyStemNum].score;
				if(dTmp > dOptimalScore)
					dOptimalScore = dTmp;
			}
		}
	} else {
		dOptimalScore = root->t_head[0].score;
	}
	return dOptimalScore;
}

void DynamicP::tb_indexingfetch(Tree_bag * tr_node, int * index_array, int num_index, int child_empty_stem_num, double * score)//, int * empty_stem_num)
{
	int	i, index;
	int nodenum = tr_node->nodenum;

	//calculate the index
	index = 0;
	for(i=0; i<nodenum; i++) {
		index += tr_node->enum_node[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), 
															static_cast<double>(nodenum-1-i)));
	}

	//support max empty stem num 
	if(getFlagSupportEmptyStem() == TRUE) {
		index += child_empty_stem_num * static_cast<int>(pow(static_cast<double>(this->iLeading), 
															static_cast<double>(nodenum)));
	}

	if(tr_node->t_head[index].optimal == 1) {
		*score	= tr_node->t_head[index].score;
	} else {
		*score	= SMALLEST;
	}
}

void DynamicP::tb_indexingstore(Tree_bag *tr_node, int *index_array, int num_index, double score, int empty_stem_num)
{
	int		i, index;
	int		nodenum = tr_node->nodenum;

	if(tr_node->pParent == tr_node) {
		//root
		index = 0;
		if(getFlagSupportEmptyStem() == TRUE)
			index += empty_stem_num;
	} else {
		//calculate the index
		index = 0;
		for(i=0; i<nodenum; i++) {
			index += tr_node->enum_node[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), 
																static_cast<double>(nodenum-1-i)));
		}
		//support max empty stem num 
		//recompute the index after adding one property - overall_empty_stem_num
		if(getFlagSupportEmptyStem() == TRUE) {
			index += empty_stem_num * static_cast<int>(pow(static_cast<double>(this->iLeading), 
															static_cast<double>(nodenum)));
		}
	}

	// 20100215
	if(score > tr_node->t_head[index].score) {
		tr_node->t_head[index].score		= score;
		tr_node->t_head[index].optimal		= 1;

		// 20100215
		for(i=0; i<tr_node->childnum; i++) {
			tr_node->t_head[index].childIndex[i] = this->_childOptimalIndex[i];
		}

		//limiting the number of empty-stem. 
		if(this->getFlagSupportEmptyStem() == TRUE) {
			tr_node->t_head[index].overall_empty_stem_num = empty_stem_num;

			//store the pos for the current candidate
			Node_list *node_list = NULL;
			int nodeid, stemid;
			int candidateIdx;
			for(i=0; i<nodenum; i++) {
				node_list = tr_node->nhead;
				for(int j=0; j<i; j++)
					node_list = node_list->next;
				nodeid = node_list->g_node;
				if(nodeid != 0 && nodeid != (2*this->getNumstems()+1))
				{
					stemid = gidtosid[nodeid].s_id;
					candidateIdx = tr_node->enum_node[i];
					if(gidtosid[nodeid].left_or_right == 0)
					{
						//if it is the left part of the stem.
						tr_node->t_head[index].lpos[i]		= this->candidateStems[stemid].s_locinfo[candidateIdx].a;
						tr_node->t_head[index].l_ext_dir[i] = this->candidateStems[stemid].s_locinfo[candidateIdx].a_ext_dir;

						tr_node->t_head[index].rpos[i]		= this->candidateStems[stemid].s_locinfo[candidateIdx].b;
						tr_node->t_head[index].r_ext_dir[i] = this->candidateStems[stemid].s_locinfo[candidateIdx].b_ext_dir;
					} else {
						//Otherwise it is the right part of the stem.
						tr_node->t_head[index].lpos[i]		= this->candidateStems[stemid].s_locinfo[candidateIdx].c;
						tr_node->t_head[index].l_ext_dir[i] = this->candidateStems[stemid].s_locinfo[candidateIdx].c_ext_dir;

						tr_node->t_head[index].rpos[i]		= this->candidateStems[stemid].s_locinfo[candidateIdx].d;
						tr_node->t_head[index].r_ext_dir[i] = this->candidateStems[stemid].s_locinfo[candidateIdx].d_ext_dir;
					}
				}
			}
		}
	}
}

void DynamicP::restoreEnumNodeCombination(Tree_bag *tr_node)
{
	int		i, index;
	int		nodenum = tr_node->nodenum;

	if(tr_node->pParent == tr_node) {
		//root
		index = 0;
	} else {
		//calculate the index
		index = 0;
		for(i=0; i<nodenum; i++) {
			index += tr_node->enum_node[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), 
																static_cast<double>(nodenum-1-i)));
		}
	}

	int srcNodeId = 0;
	int sinkNodeId = 2*this->getNumstems() + 1;

	//store the pos for the current candidate
	Node_list *node_list = NULL;
	int nodeid, stemid;
	int candidateIdx;
	for(i=0; i<nodenum; i++) {
		node_list = tr_node->nhead;
		for(int j=0; j<i; j++)
			node_list = node_list->next;
		nodeid = node_list->g_node;
		if(nodeid != srcNodeId && nodeid != sinkNodeId)
		{
			stemid = gidtosid[nodeid].s_id;
			candidateIdx = tr_node->enum_node[i];
			if( candidateIdx == (this->getLeadingNum() - 1) )
			{
				this->candidateStems[stemid].s_locinfo[candidateIdx].a = NEGATIVE_ONE;
				this->candidateStems[stemid].s_locinfo[candidateIdx].b = NEGATIVE_ONE;
				this->candidateStems[stemid].s_locinfo[candidateIdx].c = NEGATIVE_ONE;
				this->candidateStems[stemid].s_locinfo[candidateIdx].d = NEGATIVE_ONE;

				//add extend_dir 
				this->candidateStems[stemid].s_locinfo[candidateIdx].a_ext_dir = EXTNA;
				this->candidateStems[stemid].s_locinfo[candidateIdx].b_ext_dir = EXTNA;
				this->candidateStems[stemid].s_locinfo[candidateIdx].c_ext_dir = EXTNA;
				this->candidateStems[stemid].s_locinfo[candidateIdx].d_ext_dir = EXTNA;
			}
		}
	}
}

int	DynamicP::getFlagPrintDebugInfo()
{
	return this->flagPrintDebugInfo;
}

void DynamicP::setFlagPrintDebugInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugInfo = 1;
	else
		this->flagPrintDebugInfo = 0;
}

int	DynamicP::getFlagPrintDebugInputInfo()
{
	return this->flagPrintDebugInputInfo;
}

void DynamicP::setFlagPrintDebugInputInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugInputInfo = 1;
	else
		this->flagPrintDebugInputInfo = 0;
}

int	DynamicP::getFlagPrintDebugTreeInfo()
{
	return this->flagPrintDebugTreeInfo;
}

void DynamicP::setFlagPrintDebugTreeInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugTreeInfo = 1;
	else
		this->flagPrintDebugTreeInfo = 0;
}

int	DynamicP::getFlagPrintDebugDPWinInfo()
{
	return this->flagPrintDebugDPWinInfo;
}

void DynamicP::setFlagPrintDebugDPWinInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugDPWinInfo = 1;
	else
		this->flagPrintDebugDPWinInfo = 0;
}

int	DynamicP::getFlagPrintDebugTDInfo()
{
	return this->flagPrintDebugTDInfo;
}

void DynamicP::setFlagPrintDebugTDInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugTDInfo = 1;
	else
		this->flagPrintDebugTDInfo = 0;
}

int	DynamicP::getFlagPrintDebugDPSearchInfo()
{
	return this->flagPrintDebugDPSearchInfo;
}

void DynamicP::setFlagPrintDebugDPSearchInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugDPSearchInfo = 1;
	else
		this->flagPrintDebugDPSearchInfo = 0;
}

void DynamicP::print_enum_stem(Tree_bag * tr_node)
{
	cout<<"enum_stem_info\t";
	int i;
	for(i=0; i<tr_node->stnum; i++) {
		cout<<tr_node->enum_stem[i]<<" ";
	}
	cout<<endl;

}
void DynamicP::print_enum_node(int *index_array, int num_index)
{
	cout<<"enum_node = ";
	int i;
	//enum_node:	all possible combinations of nodes.
	for(i=0; i<num_index; i++) {
		cout<<index_array[i]<<" ";
	}
//	cout<<endl;
}

void DynamicP::print_stem_image(Tree_bag * root)
{
	int i, j;
	if(root == NULL)
		return;

	printTreenode(root);
	cout<<endl;

	for(i=0; i<root->childnum; i++) {
		for(j=0; j<root->nodenum; j++) {
			cout<<root->stem_image[i][j]<<" ";
		}
		cout<<endl;
	}

	for(i=0; i<root->childnum; i++) {
		print_stem_image(root->pChild[i]);
	}
}
void DynamicP::print_loop_image(Tree_bag * tr_node)
{
}
void DynamicP::printTreeBagAllInfo(Tree_bag *root)
{
	int i, iTotal, nodenum;
	nodenum = root->nodenum;
	iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), 
								static_cast<double>(nodenum)));

	if(root == NULL)
		return;
	
	this->printTreenode(root);
	cout<<endl;

	cout<<"index\toptimal\tscore\toptimal_c_id"<<endl;
	//print info of the current treenode
	for(i=0; i<iTotal; i++) {

		cout<<i;
		cout<<"\t"<<root->t_head[i].optimal<<"\t"<<root->t_head[i].score;
		cout<<endl;
	}

	//then recursively calls itself
	for(i=0; i<root->childnum; i++) {
		printTreeBagAllInfo(root->pChild[i]);
	}
}
void DynamicP::printChildNum(Tree_bag *root)
{
	int i;
	if(root == NULL)
		return;
	this->printTreenode(root);
	cout<<"\t child num="<<root->childnum<<endl;
	for(i=0; i<root->childnum; i++) {
		printChildNum(root->pChild[i]);
	}
}

void DynamicP::printTraceBackTreeInfo(Tree_bag *root)
{
	cout<<endl
		<<"================ TraceBack ================"
		<<endl;
	
	//this->printTreeNodeInfo(root);

	//this->printRootNodeInfo(root);

	this->printTreeNodeInfoIncludingPos(root);
}

void DynamicP::printNumInBase(int total, int base, int nodenum, int overall_empty_stem_num, bool bRoot)
{
	int * baseArray = new int[nodenum];
	int left;
	if(!bRoot) {
		left = total - overall_empty_stem_num * 
						static_cast<int>(pow(static_cast<double>(this->iLeading), 
											static_cast<double>(nodenum)));
	} else {
		left = total - overall_empty_stem_num;
	}
	int iTmp, curBase;
	int index = 0;
	while(index < nodenum)
	{
		curBase = static_cast<int>(pow(static_cast<double>(base), 
										static_cast<double>(nodenum-1-index)));
		iTmp = left / curBase;
		left = left % curBase;
		baseArray[index++] = iTmp;
	}
// 	baseArray[index-1] = left;
	for(int i=0; i<nodenum; i++)
		cout<<baseArray[i];
	if(baseArray != NULL)
	{
		delete [] baseArray;
		baseArray = NULL;
	}
}

void DynamicP::printTreeNodeInfo(Tree_bag *root)
{
	int i, nodenum;
	int iTotal = 0;
	nodenum = root->nodenum;
	if(root->pParent == root) {
		iTotal = this->getMaxEmptyStemNum() + 1;
	} else {
		iTotal = (this->getMaxEmptyStemNum() + 2) * 
				static_cast<int>(pow(static_cast<double>(this->iLeading), 
									static_cast<double>(root->nodenum)));
	}

	this->printTreenode(root);
	cout<<endl;

	cout<<"index\toptimal\tscore\temptystemnum\tchildIndex"<<endl;
	//print info of the current treenode
	for(i=0; i<iTotal; i++) 
	{
		if(root->t_head[i].optimal == 1)
		{
			cout<<i;
			//printNumInBase(iTotal, this->iLeading, nodenum);
			cout<<"\t"<<root->t_head[i].optimal
				<<"\t"<<root->t_head[i].score
				//limiting the number of empty-stem.  20090211
				<<"\t"<<root->t_head[i].overall_empty_stem_num
				<<"\t";
				//remove the block memory function in dp part.  20090208
			for(int j=0; j<root->childnum; j++)
				cout<<root->t_head[i].childIndex[j]<<"  ";
			cout<<endl;
		}
	}

	//then recursively calls itself
	for(i=0; i<root->childnum; i++) {
		printTreeNodeInfo(root->pChild[i]);
	}
}

void DynamicP::printTreeNodeInfoIncludingPos(Tree_bag *root)
{
	int i, j, iTotal, nodenum;
	nodenum = root->nodenum;
	bool bRoot = false;
	if(root->pParent == root) {
		iTotal = this->getMaxEmptyStemNum() + 1;
		bRoot = true;
	} else {
		iTotal = (this->getMaxEmptyStemNum() + 1) * 
				static_cast<int>(pow(static_cast<double>(this->iLeading), 
									static_cast<double>(nodenum)));
	}

	this->printTreenode(root);
	cout<<endl;

	cout<<"index\toptimal\tscore\tept_num\tchildIdx"<<endl;
	//print info of the current treenode
	for(i=0; i<iTotal; i++) 
	{
		if(root->t_head[i].optimal == 1)
		{
			cout<<i<<"-";
			printNumInBase(i, this->iLeading, nodenum, root->t_head[i].overall_empty_stem_num, bRoot);
			cout<<"\t"<<root->t_head[i].optimal
				<<"\t"<<root->t_head[i].score
				//limiting the number of empty-stem.  20090211
				<<"\t"<<root->t_head[i].overall_empty_stem_num
				<<"\t";
				//remove the block memory function in dp part.  20090208
			for(j=0; j<root->childnum; j++)
				cout<<root->t_head[i].childIndex[j]<<"  ";
/**/
			if(this->getFlagSupportEmptyStem() == TRUE) {
				//print the pos for the current candidate
				for(j=0; j<nodenum; j++)
				{
					//cout<<"("<<root->t_head[i].lpos[j]<<","<<root->t_head[i].rpos[j]<<") ";

					//add extend_dir 
					cout<<"("
						<<root->t_head[i].lpos[j]<<"-"<<root->t_head[i].l_ext_dir[j]
						<<","
						<<root->t_head[i].rpos[j]<<"-"<<root->t_head[i].r_ext_dir[j]
						<<") ";
				}
			}

			cout<<endl;
		}
	}

	//then recursively calls itself
	for(i=0; i<root->childnum; i++) {
		printTreeNodeInfoIncludingPos(root->pChild[i]);
	}
}

void DynamicP::printRootNodeInfo(Tree_bag *root)
{
	int i, iTotal, nodenum;
	nodenum = root->nodenum;
	if(root->pParent == root)
		iTotal = 1;
	else
		iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), 
									static_cast<double>(root->nodenum)));

	this->printTreenode(root);
	cout<<endl;

	cout<<"index\toptimal\tscore\temptystemnum\tchildIndex"<<endl;
	//print info of the current treenode
	for(i=0; i<iTotal; i++) 
	{
		if(root->t_head[i].optimal == 1)
		{
			cout<<i;
			cout<<"\t"<<root->t_head[i].optimal
				<<"\t"<<root->t_head[i].score
				//limiting the number of empty-stem.  20090211
				<<"\t\t"<<root->t_head[i].overall_empty_stem_num;
				//remove the block memory function in dp part.  20090208
			for(int j=0; j<root->childnum; j++)
				cout<<root->t_head[i].childIndex[j]<<"  ";
			cout<<endl;
		}
	}
}

void DynamicP::printTreeBagChildIdInfo(Tree_bag *root)
{
	int i, iTotal, nodenum;
	nodenum = root->nodenum;
	iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), 
								static_cast<double>(nodenum)));

	if(root == NULL)
		return;

	this->printTreenode(root);
	cout<<endl;
	cout<<"index\toptimal\tscore"<<endl;
	//print info of the current treenode
	for(i=0; i<iTotal; i++) {
		if(root->t_head[i].optimal == 1) {
			cout<<i;
			cout<<"\t"<<root->t_head[i].optimal
				<<"\t"<<root->t_head[i].score;
			cout<<endl;
		}
	}

	//then recursively calls itself
	for(i=0; i<root->childnum; i++) {
		printTreeBagChildIdInfo(root->pChild[i]);
	}
}

int DynamicP::countNumOftheWholeTreeNode(Tree_bag *root)
{
	int i, nodenum;
	nodenum = root->nodenum;
//	this->printTreenode(root);
//	cout<<"\t"<<nodenum<<endl;

	for(i=0; i<root->childnum; i++) {
		nodenum += countNumOftheWholeTreeNode(root->pChild[i]);
	}
	return nodenum;
}

void DynamicP::getTheWholeTreeNode(Tree_bag *root)
{
	//get the current node's tree node info
	int i, nodenum;
	nodenum = root->nodenum;
	int	* pCurrentNode = new int[nodenum];
	memset(pCurrentNode, 0, nodenum);
	Node_list * nodeList = root->nhead;
	i = 0;
	while(nodeList != NULL) {
		pCurrentNode[i++] = nodeList->g_node;
		nodeList = nodeList->next;
	}

	for(i=0; i<nodenum; i++)
		this->pTheWholeTreeNode[this->iTheWholeTreeNodeIndex++] = pCurrentNode[i];

	if(pCurrentNode != NULL) {
		delete [] pCurrentNode;
		pCurrentNode = NULL;
	}

	for(i=0; i<root->childnum; i++) {
		getTheWholeTreeNode(root->pChild[i]);
	}
}

int DynamicP::findOptimalinRoot(Tree_bag *root)
{
	int iTotal, nodenum;
	int index = -1;
	nodenum = root->nodenum;
	iTotal = 1 * this->getMaxEmptyStemNum();

	//do full combination  20090914
	double	dOptimalScore = SMALLEST;

	for(int childEmptyStemNum=0; childEmptyStemNum <= iTotal; childEmptyStemNum++)
	{
		if(root->t_head[childEmptyStemNum].optimal == 1)
		{
			double dTmp = root->t_head[childEmptyStemNum].score;
			if(dTmp > dOptimalScore) {
				dOptimalScore = dTmp;
				index = childEmptyStemNum;
			}
		}
	}

	return index;
}
void DynamicP::printOptimalTheWholeTreeNode()
{
	int i;
	for(i=0; i<this->iOptimalTheWholeTreeNodeIndex; i++)
		cout<<this->pOptimalTheWholeTreeNode[i]<<" ";
	cout<<endl;
}
/**/
void DynamicP::getOptimalPath(Tree_bag *root)
{
	//starts from the root
	int iRoot;
	if(this->getFlagSupportEmptyStem() == TRUE)
		iRoot = this->findOptimalinRoot(root);
	else
		iRoot = 0;	//root->t_head[0].optimalIndex;

	if(iRoot != -1) {
		this->getOptimalTheWholeTreeNode(root, iRoot);
	}
}

void DynamicP::getOptimalTheWholeTreeNode(Tree_bag *root, int currentIndex)
{
	int i, count, nodenum;
	nodenum = root->nodenum;
	int	* pCurrentOptimalNode = new int[nodenum];
	memset(pCurrentOptimalNode, 0, nodenum);

//	this->printTreenode(root);
//	cout<<endl;
	//save current optimal index to the pCurrentOptimalNode
	count = 0;
//	int optimal_index = currentIndex - pow(this->iLeading, nodenum) * empty_stem_num;
	int optimal_index;
	int iShift;

	if(this->getFlagSupportEmptyStem() == TRUE)
	{
		if(root->pParent == root) {
			iShift = root->t_head[currentIndex].overall_empty_stem_num;
		} else {
			iShift = root->t_head[currentIndex].overall_empty_stem_num 
					* static_cast<int>(pow(static_cast<double>(this->iLeading), 
											static_cast<double>(nodenum)));
		}
	} else {
		//empty_stem = 0;
		iShift = 0;
	}
	optimal_index = currentIndex - iShift;
	int result = 0;
	while (optimal_index) 
	{
		result = optimal_index % this->iLeading;
		pCurrentOptimalNode[count++] = result;
		optimal_index = optimal_index / this->iLeading;
	}
//	for(i=0; i<count; i++)
//		cout<<pCurrentOptimalNode[i]<<" ";
//	cout<<endl;

	if(count < nodenum) {
		for(i=0; i<nodenum-count; i++)
			this->pOptimalTheWholeTreeNode[this->iOptimalTheWholeTreeNodeIndex++] = 0;
	}
	for(i=0; i<count; i++)
		this->pOptimalTheWholeTreeNode[this->iOptimalTheWholeTreeNodeIndex++] = pCurrentOptimalNode[count-1 - i];

	if(pCurrentOptimalNode != NULL) {
		delete [] pCurrentOptimalNode;
		pCurrentOptimalNode = NULL;
	}

	for(i=0; i<root->childnum; i++) {
		getOptimalTheWholeTreeNode(	root->pChild[i], 
									root->t_head[currentIndex].childIndex[i]);

	}
}

void DynamicP::buildStemIdxArray()
{
	int i, j;
	int stemnum = this->getNumstems();

	stemPosArray = new int[2*stemnum];
	stemIdxArray = new int[2*stemnum];
	for(i=0; i<2*stemnum; i++) {
		stemPosArray[i] = stemIdxArray[i] = 0;
	}

	for(i=0; i<stemnum; i++) {
		j = 2*i;
		stemIdxArray[j]		= j+1;
		stemIdxArray[j+1]	= j+2;
		stemPosArray[j]		= this->pStem[i].start1;
		stemPosArray[j+1]	= this->pStem[i].start2;
	}

	int iTmp;

	for(i=1; i<2*stemnum; i++)
	{
		for(j=2*stemnum-1; j>=i; j--)
		{
			if(stemPosArray[j-1] > stemPosArray[j])
			{
				//swap the element
				iTmp = stemPosArray[j];
				stemPosArray[j] = stemPosArray[j-1];
				stemPosArray[j-1] = iTmp;

				iTmp = stemIdxArray[j];
				stemIdxArray[j] = stemIdxArray[j-1];
				stemIdxArray[j-1] = iTmp;
			}
		}
	}

	if(stemPosArray != NULL) {
		delete [] stemPosArray;
		stemPosArray = NULL;
	}
}

void DynamicP::buildStemIdxFullArray()
{
	int stemIdxNum = 2 * this->getNumstems();

	//adding source and sink two nodes at the beginning and the end;
	stemIdxFullArray = new int[stemIdxNum + 2];

	stemIdxFullArray[0] = 0;
	
	for(int i=0; i<stemIdxNum; i++)
		stemIdxFullArray[i+1] = stemIdxArray[i];
	
	stemIdxFullArray[stemIdxNum + 1] = stemIdxNum + 1;
}

void DynamicP::printStemIdxArray()
{
	int i, j;
	int stemnum = this->getNumstems();

	for(i=0; i<stemnum; i++) {
		j = 2*i;
		cout<<stemIdxArray[j]<<"-"<<stemIdxArray[j+1]<<"-";
	}
	cout<<endl;
}

void DynamicP::freeStemIdxArray()
{
	if(stemIdxArray != NULL)
		delete [] stemIdxArray;
	stemIdxArray = NULL;
}

void DynamicP::freeStemIdxFullArray()
{
	if(stemIdxFullArray != NULL)
		delete [] stemIdxFullArray;
	stemIdxFullArray = NULL;
}

void DynamicP::allocateStemCandIdxArray()
{
	int stemnum = this->getNumstems();
	stemCandIdxArray = new int[2*stemnum];
}
void DynamicP::buildStemCandIdxArray()
{
	int i, index;
	int iElement;
	int stemnum = this->getNumstems();

	//for every element in stemIdxArray, find its corresponding 
	//optimal index in pOptimalTheWholeTreeNode of pTheWholeTreeNode
	for(i=0; i<2*stemnum; i++) {
		iElement = stemIdxArray[i];
		index = 0;
		while(index<iTheWholeTreeNodeIndex) {
			if(iElement == pTheWholeTreeNode[index])
				break;
			index++;
		}
		stemCandIdxArray[i] = pOptimalTheWholeTreeNode[index];
	}
}

void DynamicP::printStemCandIdxArray()
{
	int i;
	int stemnum = this->getNumstems();

	for(i=0; i<2*stemnum; i++) {
		cout<<stemCandIdxArray[i]<<" ";
	}
	cout<<endl;
}

void DynamicP::freeStemCandIdxArray()
{
	if(stemCandIdxArray != NULL)
		delete [] stemCandIdxArray;
	stemCandIdxArray = NULL;
}

//limiting the number of empty-stem. 
int DynamicP::countEmptyStemNumInStemCandIdxArray()
{
	int i;
	int stemnum = this->getNumstems();
	int num = 0;

	for(i=0; i<2*stemnum; i++) {
		if(stemCandIdxArray[i] == (this->getLeadingNum() - 1))
			num++;
	}

	if(num % 2 != 0) {
		cout<<"Something wrong when computing the num of empty stem"<<endl;
		exit(1);
	} else {
		return num/2;
	}
}

void DynamicP::locateStemStartEndPosIdx()
{
	int stemnum = this->getNumstems();
	iStemStartPosIdx	= stemCandIdxArray[0];
	iStemEndPosIdx		= stemCandIdxArray[2*stemnum-1];
}

int	DynamicP::getFlagMergeCandInPreprocess()
{
	return this->flagMergeCandInPreprocess;
}

void DynamicP::setFlagMergeCandInPreprocess(int pflag)
{
	if(pflag >= 1)
		this->flagMergeCandInPreprocess = 1;
	else
		this->flagMergeCandInPreprocess = 0;
}

int	DynamicP::getFlagCandwithShortestLength()
{
	return this->flagCandwithShortestLength;
}

void DynamicP::setFlagCandwithShortestLength(int pflag)
{
	if(pflag >= 1)
		this->flagCandwithShortestLength = 1;
	else
		this->flagCandwithShortestLength = 0;
}

int	DynamicP::getShiftNumMergeCand()
{
	return this->iShiftNumMergeCand;
}

void DynamicP::setShiftNumMergeCand(int num)
{
	this->iShiftNumMergeCand = num;
}

int	DynamicP::getAllowedNullLoopInsNum()
{
	return this->iAllowedNullLoopInsNum;
}

void DynamicP::setAllowedNullLoopInsNum(int num)
{
	this->iAllowedNullLoopInsNum = num;
}

void DynamicP::setNumOfNtsOverlapBetweenStem(int num)
{
	this->iNumOfNtsOverlapBetweenStem = num;
}

int DynamicP::getNumOfNtsOverlapBetweenStem()
{
	return this->iNumOfNtsOverlapBetweenStem;
}

void DynamicP::IterateTD_topdown(TRnaGraphTreeBag<int, Tree_bag>* pRoot)
{
	unsigned int i = 0;
	cout << pRoot->m_nID<< ':';// << pRoot->m_pPayload << " : ";
	while( i < pRoot->m_vecpNodes.size() )
	{
		cout << pRoot->m_vecpNodes[i]->ID() << ' ';// << pRoot->m_vecpNodes[i]->Payload() << ' ';
		i++;
	}
	cout << endl;
	
	i = 0;
	while( i < pRoot->m_vecChildren.size() )
	{
		IterateTD_topdown( pRoot->m_vecChildren[i] );
		i++;
	}
}

void DynamicP::IterateTD_Payload_topdown(TRnaGraphTreeBag<int, Tree_bag>* pRoot)
{
	unsigned int i = 0;
	cout << pRoot->m_nID<< ':';// << pRoot->m_pPayload << " : ";
	while( i < pRoot->m_vecpNodes.size() )
	{
		cout << pRoot->m_vecpNodes[i]->Payload() << ' ';
		i++;
	}
	cout << endl;
	
	i = 0;
	while( i < pRoot->m_vecChildren.size() )
	{
		IterateTD_Payload_topdown( pRoot->m_vecChildren[i] );
		i++;
	}
}

void DynamicP::IterateTD_Payload_bottomup(TRnaGraphTreeBag<int, Tree_bag>* pRoot)
{
	unsigned int i = 0;
	while( i < pRoot->m_vecChildren.size() )
	{
		IterateTD_Payload_bottomup( pRoot->m_vecChildren[i] );
		i++;
	}

	cout << pRoot->m_nID<< ':';// << pRoot->m_pPayload << " : ";
	i = 0;
	while( i < pRoot->m_vecpNodes.size() )
	{
		cout << pRoot->m_vecpNodes[i]->Payload() << ' ';
		i++;
	}
	cout << endl;
	
}

//create my own tree-bag tree by copying the tree-bag
Tree_bag * DynamicP::copy_treedecomp(TRnaGraphTreeBag<int, Tree_bag>* pRoot)
{
	Tree_bag * cur_node = NULL;
	Node_list * nodelist = NULL;
	unsigned int i = 0;
	//create a tree-bag, and copy the tree-node to its node-list
	cur_node = new Tree_bag[1];
	cur_node->nodenum = pRoot->m_vecpNodes.size();
	cur_node->nhead = new Node_list[1];
	nodelist = cur_node->nhead;

	while(i < (pRoot->m_vecpNodes.size()-1))
	{
		nodelist->g_node = pRoot->m_vecpNodes[i]->Payload();
		nodelist->next = new Node_list[1];
		nodelist = nodelist->next;
		i++;
	}
	nodelist->g_node = pRoot->m_vecpNodes[i]->Payload();
	nodelist->next = NULL;

	//then deal with its child tree-bag.
	cur_node->childnum = pRoot->m_vecChildren.size();
	i = 0;
	while( i < pRoot->m_vecChildren.size() )
	{
		cur_node->pChild[i] = copy_treedecomp(pRoot->m_vecChildren[i]);
		i++;
	}
	return cur_node;
}

void DynamicP::freeMyTDRoot()
{
	this->free_treenode(this->mypTDRoot);
}
Tree_bag * DynamicP::getMyTDRoot()
{
	return this->mypTDRoot;
}

double DynamicP::getPcoeff()
{
	return this->pcoeff;
}

void DynamicP::setPcoeff(double _pcoeff)
{
	this->pcoeff = _pcoeff;
}

int	DynamicP::getFlagPrintStrAlignInfo()
{
	return this->flagPrintStrAlignInfo;
}

void DynamicP::setFlagPrintStrAlignInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintStrAlignInfo = 1;
	else
		this->flagPrintStrAlignInfo = 0;
}

int	DynamicP::getFlagPrintScoreInfo()
{
	return this->flagPrintScoreInfo;
}

void DynamicP::setFlagPrintScoreInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintScoreInfo = 1;
	else
		this->flagPrintScoreInfo = 0;
}

/*
 * The idea is to extract all the information from
 *		pTheWholeTreeNode
 *		stemIdxArray
 *		stemCandidIdxArray
 */
void DynamicP::buildStructureAlignInfo(Loop *pLoopModel, int seq_start, int seq_end)
{
	//adding back the sequence in the sequence-structure alignment 20081220
	string winseq(cWindowSequence);

	int m, n;
	int stemIdx;
	int	stemnum = this->getNumstems();

	int preLeftStemPos = -1, preRightStemPos = -1;
	int curLeftStemPos = -1, curRightStemPos = -1;
	int pre_gid = 0, cur_gid = 0;

	VTBSearch vtber(1);
	vtber.setAllowedInsNum(this->getAllowedNullLoopInsNum());
	int loopid, start, end;
//	double dLoopScore = 0.0;

	//structure alignment
	memset(structureSequence, 0, MAX_STRUCTURE_ALIGN_LENGTH);
	memset(strucutreAlignment, 0, MAX_STRUCTURE_ALIGN_LENGTH);

	//support stemid index in two pastalines 20081015
	memset(strucutreAlignment_charid_idx, 0, MAX_STRUCTURE_ALIGN_LENGTH);

	char *	tmp_sub_charid_idx = NULL;

	//adding the empty-stem model.  20090211
	int		gap_len = 0;
	int		empty_start = 0;

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;

			//adding the empty-stem model.  20090211
			if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
			{
				if(stemIdx == 0) {
					strcat(structureSequence, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
					strcat(strucutreAlignment, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment);

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						strcat(strucutreAlignment_charid_idx, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].pastaline_idx);
					}

				} else {
					//adding the empty-stem model.  20090211
					if(preRightStemPos != NEGATIVE_ONE)
					{
						if(curLeftStemPos > preRightStemPos) 
						{
							if(curLeftStemPos > (preRightStemPos+1))
							{
								//add loop structure alignment info here.

								char * loopAlign	= NULL;
								char * loopSequence	= NULL;

								loopid	= this->identifyLoopId(pre_gid, cur_gid);
								if(loopid != NEGATIVE_ONE)	//20100120
								{
									start	= preRightStemPos + 1;
									end		= curLeftStemPos - 1;

									double dLoopScore = vtber.VTBSearchMax(&pLoopModel[loopid], cWindowSequence, start, end, &loopAlign, &loopSequence);

									//only consider the non-negative loop score
									if(dLoopScore != INVLDPROB)
									{
										//add loop sequence/alignment info
										if((this->getAllowedNullLoopInsNum() >= (end-start+1))
											&& (dLoopScore == 0.0))
										{
											for(int n=start; n<=end; n++)
											{
												char *	cTmp = strappend(structureSequence, cWindowSequence[n]);
												strcpy(structureSequence, cTmp);
												delete [] cTmp; cTmp = NULL;
												cTmp = strappend(strucutreAlignment, '.');
												strcpy(strucutreAlignment, cTmp);

												//support stemid index in two pastalines 20081015
												if(this->getNumPastaLine() == TWO_PASTALINE) {
													strcpy(strucutreAlignment_charid_idx, cTmp);
												}
												delete [] cTmp; cTmp = NULL;
											}
										} else {
				//							cout<<loopSequence<<" <-> "<<loopAlign<<endl;
											//postprocessing the loopSequence/loopAlign
											//if '-' exists in loopSequence and 'd' in loopalign, then delete 'd' from loopalign
											strcat(structureSequence, loopSequence);
											strcat(strucutreAlignment, loopAlign);

											//support stemid index in two pastalines 20081015
											if(this->getNumPastaLine() == TWO_PASTALINE) {
												strcat(strucutreAlignment_charid_idx, loopAlign);
											}
										}
									}
								}

								if(loopAlign != NULL) {
									delete [] loopAlign;
									loopAlign = NULL;
								}
								if(loopSequence != NULL) {
									delete [] loopSequence;
									loopSequence = NULL;
								}
							}

							//add the next stem sequence/alignment info
							strcat(structureSequence, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
							strcat(strucutreAlignment, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment);

							//support stemid index in two pastalines 20081015
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								strcat(strucutreAlignment_charid_idx, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].pastaline_idx);	//20081015
							}
						} else {
							//overlap
							int	iOverlap = preRightStemPos - curLeftStemPos + 1;
							int	iTmpLength = strlen(strucutreAlignment);
							for(n=0; n<iOverlap; n++) {
								strucutreAlignment[iTmpLength-1-n] = '*';

								//support stemid index in two pastalines 20081015
								if(this->getNumPastaLine() == TWO_PASTALINE) {
									strucutreAlignment_charid_idx[iTmpLength-1-n] = '*';
								}
							}


							iTmpLength = strlen(candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
							char * cSeqOther = new char[iTmpLength - iOverlap + 1];
							char * cAliOther = new char[iTmpLength - iOverlap + 1];

							for(n=0; n<(iTmpLength - iOverlap); n++) {
								cSeqOther[n] = candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence[n+iOverlap];
								cAliOther[n] = candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment[n+iOverlap];
							}
							cSeqOther[n] = '\0';
							cAliOther[n] = '\0';
							
							strcat(structureSequence, cSeqOther);
							strcat(strucutreAlignment, cAliOther);

							//support stemid index in two pastalines 20081015
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								tmp_sub_charid_idx = strrpl(cAliOther, toupper(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
								strcat(strucutreAlignment_charid_idx, tmp_sub_charid_idx);
								if(tmp_sub_charid_idx != NULL) {
									delete [] tmp_sub_charid_idx;
									tmp_sub_charid_idx = NULL;
								}
							}

							if(cSeqOther != NULL) {
								delete cSeqOther;
								cSeqOther = NULL;
							}
							if(cAliOther != NULL) {
								delete cAliOther;
								cAliOther = NULL;
							}
						}
					}
					else	//current is not empty stem, but at least previous is empty stem
					{
						//adding the empty-stem model.  20090211
						gap_len = curLeftStemPos - empty_start;//preLeftStemPos;
						if(gap_len >= 0) {
							strcat(structureSequence, winseq.substr(empty_start, gap_len).c_str());
							for(n=0; n<gap_len; n++) {
								strcat(strucutreAlignment, EMPTY_STEM_LABEL);
							}
							//add the next stem sequence/alignment info
							strcat(structureSequence, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
							strcat(strucutreAlignment, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment);

							//making k-empty stem model support two pasta line format 20081221
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								for(n=0; n<gap_len; n++) {
									strcat(strucutreAlignment_charid_idx, EMPTY_STEM_LABEL);
								}

								tmp_sub_charid_idx = strrpl(candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment, toupper(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
								strcat(strucutreAlignment_charid_idx, tmp_sub_charid_idx);
								if(tmp_sub_charid_idx != NULL) {
									delete [] tmp_sub_charid_idx;
									tmp_sub_charid_idx = NULL;
								}
							}
						}
						else
						{	//adding allowing overlap when checking position confliction 20081221
							//overlap case
							gap_len = abs(gap_len);

							//alignment structure info
							for(n=0; n<gap_len; n++)
								strucutreAlignment[strlen(strucutreAlignment)-1-n] = '*';
							string strtmp(candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment);
							strtmp = strtmp.substr(gap_len, strtmp.length() - gap_len);
							strcat(strucutreAlignment, strtmp.c_str());

							//alignment sequence info
							string strseq(candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
							strseq = strseq.substr(gap_len, strseq.length() - gap_len);
							strcat(structureSequence, strseq.c_str());

							//making k-empty stem model support two pasta line format 20081221
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								for(n=0; n<gap_len; n++)
									strucutreAlignment_charid_idx[strlen(strucutreAlignment_charid_idx)-1-n] = '*';

								string strcharidx(candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].pastaline_idx);
								strcharidx = strcharidx.substr(gap_len, strcharidx.length() - gap_len);
								strcat(strucutreAlignment_charid_idx, strcharidx.c_str());
							}
						}
/*
						//support stemid index in two pastalines 20081015
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							strcat(strucutreAlignment_charid_idx, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].pastaline_idx);	//20081015
						}
*/
					}
				}
			} //end of if
		}
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;

			//adding the empty-stem model.  20090211
			if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
			{
				if(preRightStemPos != NEGATIVE_ONE)
				{
					if(curLeftStemPos > preRightStemPos) 
					{
						if(curLeftStemPos > (preRightStemPos+1))
						{
							//add loop structure alignment info here.

							char * loopAlign	= NULL;
							char * loopSequence	= NULL;

							loopid	= this->identifyLoopId(pre_gid, cur_gid);
							if(loopid != NEGATIVE_ONE)	//20100120
							{
								start	= preRightStemPos + 1;
								end		= curLeftStemPos - 1;
								double dLoopScore = vtber.VTBSearchMax(&pLoopModel[loopid], cWindowSequence, start, end, &loopAlign, &loopSequence);

								//only consider the non-negative loop score
								if(dLoopScore != INVLDPROB)
								{
									//add loop sequence/alignment info
									if((this->getAllowedNullLoopInsNum() >= (end-start+1))
										&& (dLoopScore == 0.0))
									{
										for(int n=start; n<=end; n++)
										{
											char *	cTmp = strappend(structureSequence, cWindowSequence[n]);
											strcpy(structureSequence, cTmp);
											if(cTmp != NULL) {
												delete [] cTmp; cTmp = NULL;
											}
											cTmp = strappend(strucutreAlignment, '.');
											strcpy(strucutreAlignment, cTmp);

											//support stemid index in two pastalines 20081015
											if(this->getNumPastaLine() == TWO_PASTALINE) {
												strcpy(strucutreAlignment_charid_idx, cTmp);
											}
											if(cTmp != NULL) {
												delete [] cTmp; cTmp = NULL;
											}
										}
									} else {
				//						cout<<loopSequence<<" <-> "<<loopAlign<<endl;
										strcat(structureSequence, loopSequence);
										strcat(strucutreAlignment, loopAlign);

										//support stemid index in two pastalines 20081015
										if(this->getNumPastaLine() == TWO_PASTALINE) {
											strcat(strucutreAlignment_charid_idx, loopAlign);
										}
									}
								}
							}

							if(loopAlign != NULL) {
								delete [] loopAlign;
								loopAlign = NULL;
							}
							if(loopSequence != NULL) {
								delete [] loopSequence;
								loopSequence = NULL;
							}
						}

						//add the next stem sequence/alignment info
						strcat(structureSequence, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
						strcat(strucutreAlignment, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment);

						//support stemid index in two pastalines 20081015
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							strcat(strucutreAlignment_charid_idx, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].pastaline_idx);	//20081015
						}
					} else {
						//overlap
						int	iOverlap = preRightStemPos - curLeftStemPos + 1;

						int	iTmpLength = strlen(strucutreAlignment);
						for(n=0; n<iOverlap; n++) {
							strucutreAlignment[iTmpLength-1-n] = '*';
							//support stemid index in two pastalines 20081015
							strucutreAlignment_charid_idx[iTmpLength-1-n] = '*';
						}

						iTmpLength = strlen(candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
						char * cSeqOther = new char[iTmpLength - iOverlap + 1];
						char * cAliOther = new char[iTmpLength - iOverlap + 1];

						for(n=0; n<(iTmpLength - iOverlap); n++) {
							cSeqOther[n] = candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence[n+iOverlap];
							cAliOther[n] = candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment[n+iOverlap];
						}
						cSeqOther[n] = '\0';
						cAliOther[n] = '\0';
						strcat(structureSequence, cSeqOther);
						strcat(strucutreAlignment, cAliOther);

						//support stemid index in two pastalines 20081015
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							tmp_sub_charid_idx = strrpl(cAliOther, tolower(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
							strcat(strucutreAlignment_charid_idx, tmp_sub_charid_idx);
							if(tmp_sub_charid_idx != NULL) {
								delete [] tmp_sub_charid_idx;
								tmp_sub_charid_idx = NULL;
							}
						}

						if(cSeqOther != NULL) {
							delete cSeqOther;
							cSeqOther = NULL;
						}
						if(cAliOther != NULL) {
							delete cAliOther;
							cAliOther = NULL;
						}
					}
				}
				else	//current is not empty stem, but at least previous is empty stem
				{
/*
					//adding the empty-stem model.  20090211
					gap_len = curLeftStemPos - empty_start;//preLeftStemPos;
					strcat(structureSequence, winseq.substr(empty_start, gap_len).c_str());
					for(n=0; n<gap_len; n++) {
//						strcat(structureSequence, EMPTY_STEM_LABEL);
						strcat(strucutreAlignment, EMPTY_STEM_LABEL);
					}

					//add the next stem sequence/alignment info
					strcat(structureSequence, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
					strcat(strucutreAlignment, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment);
*/
					//adding the empty-stem model.  20090211
					gap_len = curLeftStemPos - empty_start;//preLeftStemPos;
					if(gap_len >= 0) {
						strcat(structureSequence, winseq.substr(empty_start, gap_len).c_str());
						for(n=0; n<gap_len; n++) {
							strcat(strucutreAlignment, EMPTY_STEM_LABEL);
						}
						//add the next stem sequence/alignment info
						strcat(structureSequence, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
						strcat(strucutreAlignment, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment);

						//making k-empty stem model support two pasta line format 20081221
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							for(n=0; n<gap_len; n++) {
								strcat(strucutreAlignment_charid_idx, EMPTY_STEM_LABEL);
							}

							tmp_sub_charid_idx = strrpl(candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment, tolower(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
							strcat(strucutreAlignment_charid_idx, tmp_sub_charid_idx);
							if(tmp_sub_charid_idx != NULL) {
								delete [] tmp_sub_charid_idx;
								tmp_sub_charid_idx = NULL;
							}
						}
					}
					else
					{	//adding allowing overlap when checking position confliction 20081221
						//overlap case
						gap_len = abs(gap_len);

						//alignment structure info
						for(n=0; n<gap_len; n++)
							strucutreAlignment[strlen(strucutreAlignment)-1-n] = '*';
						string strtmp(candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment);
						strtmp = strtmp.substr(gap_len, strtmp.length() - gap_len);
						strcat(strucutreAlignment, strtmp.c_str());

						//alignment sequence info
						string strseq(candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
						strseq = strseq.substr(gap_len, strseq.length() - gap_len);
						strcat(structureSequence, strseq.c_str());

						//making k-empty stem model support two pasta line format 20081221
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							for(n=0; n<gap_len; n++)
								strucutreAlignment_charid_idx[strlen(strucutreAlignment_charid_idx)-1-n] = '*';

							string strcharidx(candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].pastaline_idx);
							strcharidx = strcharidx.substr(gap_len, strcharidx.length() - gap_len);
							strcat(strucutreAlignment_charid_idx, strcharidx.c_str());
						}
					}
/*
					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						strcat(strucutreAlignment_charid_idx, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].pastaline_idx);	//20081015
					}
*/
				}
			}
		}

		//adding the empty-stem model.  20090211
		if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
		{
			preLeftStemPos	= curLeftStemPos;
			preRightStemPos	= curRightStemPos;
		} else {
			if(preLeftStemPos != NEGATIVE_ONE && preRightStemPos != NEGATIVE_ONE)
			{
				empty_start = preRightStemPos + 1;
				preLeftStemPos	= curLeftStemPos;
				preRightStemPos	= curRightStemPos;
			}
			else
			{
				preLeftStemPos	= empty_start;
				preRightStemPos	= curRightStemPos;
			}
		}
		
		pre_gid	= cur_gid;
	}//for
//	cout<<structureSequence<<endl
//		<<strucutreAlignment<<endl;
}

void DynamicP::printStructureAlignInfo(int real_start, int real_end, ofstream &outfile)
{
	//preprocess the structureAlignment
	char * strAlign_01 = strrpl(strucutreAlignment, 'm', '.');
	char * strAlign_02 = strrpl(strAlign_01, 'd', '.');
	char * strAlign_03 = strrpl(strAlign_02, 'i', '-');
	
	char * strAlign_04 = NULL;
	char * strAlign_05 = NULL;
	char * strAlign_06 = NULL;
	if(this->getNumPastaLine() == TWO_PASTALINE) {
		//preprocess the strucutreAlignment_charid_idx
		strAlign_04 = strrpl(strucutreAlignment_charid_idx, 'm', '.');
		strAlign_05 = strrpl(strAlign_04, 'd', '.');
		strAlign_06 = strrpl(strAlign_05, 'i', '.');

		formatStringOutput(	string("Structure alignment"), 
							string(strAlign_03),
							string(strAlign_06),
							string(structureSequence), 
							real_start,
							real_end,
							1,
							outfile);
		oneHit.setStruPastaCharId(string(strAlign_06));
	} else {
		formatStringOutput(	string("Structure alignment"), 
							string(strAlign_03),
							string(""),
							string(structureSequence), 
							real_start,
							real_end,
							0,
							outfile);
		oneHit.setStruPastaCharId("");
	}
	//support hit sorting 20081121
	oneHit.setStruSeq(structureSequence);
	oneHit.setStruAlign(strAlign_03);

	if(strAlign_01 != NULL) {
		delete [] strAlign_01;
		strAlign_01 = NULL;
	}
	if(strAlign_02 != NULL) {
		delete [] strAlign_02;
		strAlign_02 = NULL;
	}
	if(strAlign_03 != NULL) {
		delete [] strAlign_03;
		strAlign_03 = NULL;
	}
	if(strAlign_04 != NULL) {
		delete [] strAlign_04;
		strAlign_04 = NULL;
	}
	if(strAlign_05 != NULL) {
		delete [] strAlign_05;
		strAlign_05 = NULL;
	}
	if(strAlign_06 != NULL) {
		delete [] strAlign_06;
		strAlign_06 = NULL;
	}
}
void DynamicP::printStructureAlignInfo2(int real_start, int real_end, ofstream &outfile)
{
	int i, length = strlen(strucutreAlignment);
	if(length != strlen(strucutreAlignment_charid_idx)) {
		cout<<"Something wrong here DynamicP::printStructureAlignInfo2()"<<endl;
		exit(0);
	}
	for(i=0; i<length; i++) {
		if((strucutreAlignment[i] == 'm' && strucutreAlignment_charid_idx[i]==109)	//ascii
		||(strucutreAlignment[i] == 'd' && strucutreAlignment_charid_idx[i]==100))
		{
			strucutreAlignment[i] = '.';
			strucutreAlignment_charid_idx[i] = '.';
		} else if(strucutreAlignment[i] == 'i' && strucutreAlignment_charid_idx[i]==105) {
			strucutreAlignment[i] = '-';
			strucutreAlignment_charid_idx[i] = '.';
		}
	}

	if(this->getNumPastaLine() == TWO_PASTALINE) {
		formatStringOutput(	string("Structure alignment"), 
							string(strucutreAlignment),
							string(strucutreAlignment_charid_idx),
							string(structureSequence), 
							real_start,
							real_end,
							1,
							outfile);
		oneHit.setStruPastaCharId(string(strucutreAlignment_charid_idx));
	} else {
		formatStringOutput(	string("Structure alignment"), 
							string(strucutreAlignment),
							string(""),
							string(structureSequence), 
							real_start,
							real_end,
							0,
							outfile);
		oneHit.setStruPastaCharId("");
	}
	//support hit sorting 20081121
	oneHit.setStruSeq(structureSequence);
	oneHit.setStruAlign(strucutreAlignment);
}

void DynamicP::printCandidateStemInfo_charid(ofstream &outfile)
{
	int m;
	int cur_gid;
	int stemIdx;
	int	stemnum = this->getNumstems();
	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			outfile<<"       "<<this->pStem[stemIdx].charid<<left<<setw(WIDTH-8)<<"";
		} 
		else 
		{	
			//even - right_arm
			outfile<<"       "<<this->pStem[stemIdx].charid<<left<<setw(WIDTH-8)<<"";
		}
	}
	outfile<<endl;
}
void DynamicP::printCandidateStemInfo_pos(int index, ofstream &outfile)
{
	int m;
	int cur_gid;
	int stemIdx;
	int	stemnum = this->getNumstems();
	int curLeftStemPos = 0, curRightStemPos = 0;

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[index].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[index].b;
			outfile<<"["<<right<<setw(3)<<curLeftStemPos<<","<<left<<setw(3)<<curRightStemPos<<"] ";
		} 
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[index].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[index].d;
			outfile<<"["<<right<<setw(3)<<curLeftStemPos<<","<<left<<setw(3)<<curRightStemPos<<"] ";
		}
	}
	outfile<<endl;
}
void DynamicP::printCandidateStemInfo_score(int index, ofstream &outfile)
{
	int m;
	int cur_gid;
	int stemIdx;
	int	stemnum = this->getNumstems();

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			if(candidateStems[stemIdx].s_locinfo[index].p_score > 1)
				outfile<<"( "<<setw(5)<<setprecision(4)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > 0.1)
				outfile<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > 0.01)
				outfile<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > 0)
				outfile<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > -0.01)
				outfile<<"( "<<setw(5)<<setprecision(0)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > -0.1)
				outfile<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > -1.0)
				outfile<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else
				outfile<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
		} 
		else 
		{	
			//even - right_arm
			if(candidateStems[stemIdx].s_locinfo[index].p_score > 1)
				outfile<<"( "<<setw(5)<<setprecision(4)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > 0.1)
				outfile<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > 0.01)
				outfile<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > 0)
				outfile<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > -0.01)
				outfile<<"( "<<setw(5)<<setprecision(0)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > -0.1)
				outfile<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[index].p_score > -1.0)
				outfile<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
			else
				outfile<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[index].p_score<<" ) ";
		}
	}
	outfile<<endl;
}
void DynamicP::printHitStemInfo_pos(ofstream &outfile)
{
	int m;
	int cur_gid;
	int stemIdx;
	int	stemnum = this->getNumstems();
	int curLeftStemPos = 0, curRightStemPos = 0;
			
	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;
			outfile<<"["<<right<<setw(3)<<curLeftStemPos<<","<<left<<setw(3)<<curRightStemPos<<"] ";
		} 
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;
			outfile<<"["<<right<<setw(3)<<curLeftStemPos<<","<<left<<setw(3)<<curRightStemPos<<"] ";
		}
	}
	outfile<<endl;
}

void DynamicP::printHitStemInfo_score(ofstream &outfile)
{
	int m;
	int cur_gid;
	int stemIdx;
	int	stemnum = this->getNumstems();

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > 1)
				outfile<<"( "<<setw(5)<<setprecision(4)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > 0.1)
				outfile<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > 0.01)
				outfile<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > 0)
				outfile<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > -0.01)
				outfile<<"( "<<setw(5)<<setprecision(0)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > -0.1)
				outfile<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > -1.0)
				outfile<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else
				outfile<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
		} 
		else 
		{	
			//even - right_arm
			if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > 1)
				outfile<<"( "<<setw(5)<<setprecision(4)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > 0.1)
				outfile<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > 0.01)
				outfile<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > 0)
				outfile<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > -0.01)
				outfile<<"( "<<setw(5)<<setprecision(0)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > -0.1)
				outfile<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else if(candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score > -1.0)
				outfile<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
			else
				outfile<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
		}
	}
	outfile<<endl;
}
/*
 * print the candidate stem info
 * 20081024
 */
void DynamicP::printCandidateStemInfo(ofstream &outfile)
{
	int m, n;
	int	stemnum = this->getNumstems();

	for(n=0; n<this->getLeadingNum(); n++)
	{
		if(n == 0)
		{
			outfile<<"Result from Preprocessing"<<endl;
			//print out the stem charid
			this->printCandidateStemInfo_charid(outfile);
		}
		//print out the position of one stem
		outfile<<left<<setw(2)<<n<<":";
		this->printCandidateStemInfo_pos(n, outfile);

		//print out the score for that stem
		outfile<<"   ";
		this->printCandidateStemInfo_score(n, outfile);

		if(n == (this->getLeadingNum()-1))
		{
			outfile<<endl<<"Result from Dynamic Programming"<<endl;
			//print out the candidate stem index in valid combination
			for(m=0; m<2*stemnum; m++)
			{
				outfile<<"       "<<stemCandIdxArray[m]<<left<<setw(WIDTH-8)<<"";
			}
			outfile<<endl;
			outfile<<"   ";
			this->printHitStemInfo_pos(outfile);

			//print out the score for that stem
			outfile<<"   ";
			this->printHitStemInfo_score(outfile);
		}
	}
}

void DynamicP::printCandidateStemInfo_FailedCase(ofstream &outfile)
{
	int n;
//	int	stemnum = this->getNumstems();

	for(n=0; n<this->getLeadingNum(); n++)
	{
		if(n == 0)
		{
			outfile<<"Result from Preprocessing"<<endl;
			//print out the stem charid
			this->printCandidateStemInfo_charid(outfile);
		}
		//print out the position of one stem
		outfile<<left<<setw(2)<<n<<":";
		this->printCandidateStemInfo_pos(n, outfile);

		//print out the score for that stem
		outfile<<"   ";
		this->printCandidateStemInfo_score(n, outfile);
	}
}

void DynamicP::printCandidateStemInfo_preprocessingInfo(ofstream &outfile) {
	printCandidateStemInfo_FailedCase(outfile);
}

int DynamicP::checkValidation()
{
	int iValid = TRUE;
	int m;
	int stemIdx;
	int	stemnum = this->getNumstems();

	int preLeftStemPos = NEGATIVE_ONE, preRightStemPos = NEGATIVE_ONE;
	int curLeftStemPos = NEGATIVE_ONE, curRightStemPos = NEGATIVE_ONE;
	int pre_gid, cur_gid;
	pre_gid = cur_gid = 0;

	int	empty_start = 0;

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;

//			cout<<cur_gid<<"["<<stemIdx<<"]\tL("<<curLeftStemPos<<","<<curRightStemPos<<")"<<endl;

			///adding the empty-stem model.  20090211
			if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
			{
				if(preLeftStemPos != NEGATIVE_ONE && preRightStemPos != NEGATIVE_ONE)
				{
					if(preRightStemPos >= curLeftStemPos) {
						iValid = FALSE;
						break;
					}
				} else {
					if(empty_start > curLeftStemPos) {
						iValid = FALSE;
						break;
					}
				}
			}//end of if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
		} 
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;

//			cout<<cur_gid<<"["<<stemIdx<<"]\tL("<<curLeftStemPos<<","<<curRightStemPos<<")"<<endl;

			//adding the empty-stem model.  20090211
			if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
			{
				if(preLeftStemPos != NEGATIVE_ONE && preRightStemPos != NEGATIVE_ONE)
				{
					if(preRightStemPos >= curLeftStemPos) {
						iValid = FALSE;
						break;
					}
				} else {
					if(empty_start > curLeftStemPos) {
						iValid = FALSE;
						break;
					}
				}
			}
		}

		if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
		{
			preLeftStemPos	= curLeftStemPos;
			preRightStemPos	= curRightStemPos;
		} else {
			if(preLeftStemPos != NEGATIVE_ONE && preRightStemPos != NEGATIVE_ONE)
			{
				//empty_start = preRightStemPos + 1;
				empty_start = preRightStemPos;
				preLeftStemPos	= curLeftStemPos;
				preRightStemPos	= curRightStemPos;
			}
			else
			{
				preLeftStemPos	= empty_start;
				preRightStemPos	= curRightStemPos;
			}
		}
		pre_gid	= cur_gid;
	}//for
	return iValid;
}

void DynamicP::buildFoldedStructureInfo()
{
	int m, n;
	int start, end;
	int stemIdx;
	int	stemnum = this->getNumstems();

	int preLeftStemPos = NEGATIVE_ONE, preRightStemPos = NEGATIVE_ONE;
	int curLeftStemPos = NEGATIVE_ONE, curRightStemPos = NEGATIVE_ONE;
	int pre_gid, cur_gid;
	pre_gid = cur_gid = 0;

	//folded structure
//	memset(strucutreAlignment, 0, MAX_STRUCTURE_ALIGN_LENGTH);
	char	tmpFoldedstrucutre[MAX_STRUCTURE_ALIGN_LENGTH];
	memset(tmpFoldedstrucutre, 0, MAX_STRUCTURE_ALIGN_LENGTH);
    memset(foldedstrucutre, 0, MAX_STRUCTURE_ALIGN_LENGTH);

	//support stemid index in two pastalines 20081015
	char *	tmp_sub_charid_idx = NULL;
	char	tmp_charid_idx[MAX_STRUCTURE_ALIGN_LENGTH];
	memset(tmp_charid_idx, 0, MAX_STRUCTURE_ALIGN_LENGTH);
	memset(foldedstrucutre_charid_idx, 0, MAX_STRUCTURE_ALIGN_LENGTH);

	const int left=0, right=1;
	char *	tmpSubFold = NULL;

	//adding the empty-stem model.  20090211
	int		gap_len = 0;
	int		empty_start = 0;

//	cout<<"----------------"<<endl;
	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;
			// 20100215
			//cout<<"["<<curLeftStemPos<<","<<curRightStemPos<<"] ";
//			cout<<cur_gid<<"["<<stemIdx<<"]\tL("<<curLeftStemPos<<","<<curRightStemPos<<")"<<endl;

			///adding the empty-stem model.  20090211
			if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
			{
				if(stemIdx == 0) {
					tmpSubFold = this->buildstemfoldedstructure(stemIdx, stemCandIdxArray[m], left);
					strcat(tmpFoldedstrucutre, tmpSubFold);

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
					{
						tmp_sub_charid_idx = strrpl(tmpSubFold, toupper(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
						strcat(tmp_charid_idx, tmp_sub_charid_idx);
						if(tmp_sub_charid_idx != NULL) {
							delete [] tmp_sub_charid_idx;
							tmp_sub_charid_idx = NULL;
						}
					}

//					//for debug
//					cout<<tmpFoldedstrucutre<<endl;

					if(tmpSubFold != NULL) {
						delete [] tmpSubFold;
						tmpSubFold = NULL;
					}
				} else {
					//adding the empty-stem model.  20090211
					if(preRightStemPos != NEGATIVE_ONE)
					{
						if(curLeftStemPos > preRightStemPos) 
						{
							if(curLeftStemPos > (preRightStemPos+1))
							{
								//to loop, just append '.'
								start	= preRightStemPos + 1;
								end		= curLeftStemPos - 1;
								for(n=start; n<=end; n++) {
									strcat(tmpFoldedstrucutre, ".");
									if(this->getNumPastaLine() == TWO_PASTALINE) {
										//support stemid index in two pastalines 20081015
										strcat(tmp_charid_idx, ".");
									}
								}
							}

							//add the next stem sequence/alignment info
							tmpSubFold = this->buildstemfoldedstructure(stemIdx, stemCandIdxArray[m], left);
							//for debug
//							cout<<tmpSubFold<<endl;

							strcat(tmpFoldedstrucutre, tmpSubFold);

							//support stemid index in two pastalines 20081015
							if(this->getNumPastaLine() == TWO_PASTALINE) 
							{
								tmp_sub_charid_idx = strrpl(tmpSubFold, toupper(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
								strcat(tmp_charid_idx, tmp_sub_charid_idx);
								if(tmp_sub_charid_idx != NULL) {
									delete [] tmp_sub_charid_idx;
									tmp_sub_charid_idx = NULL;
								}
							}
							
//							//for debug
//							cout<<tmpFoldedstrucutre<<endl;

							if(tmpSubFold != NULL) {
								delete [] tmpSubFold;
								tmpSubFold = NULL;
							}
						} else {
							//overlap
							int	iOverlap = preRightStemPos - curLeftStemPos + 1;
							int	iTmpLength = strlen(tmpFoldedstrucutre);
							for(n=0; n<iOverlap; n++) {
								tmpFoldedstrucutre[iTmpLength-1-n] = '*';
								if(this->getNumPastaLine() == TWO_PASTALINE) {
									//support stemid index in two pastalines 20081015
									tmp_charid_idx[iTmpLength-1-n] = '*';
								}
							}

							iTmpLength = strlen(candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
							char * cAliOther = new char[iTmpLength - iOverlap + 1];
							for(n=0; n<(iTmpLength - iOverlap); n++) {
								cAliOther[n] = candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment[n+iOverlap];
							}
							cAliOther[n] = '\0';
							strcat(tmpFoldedstrucutre, cAliOther);

//							//for debug
//							cout<<tmpFoldedstrucutre<<endl;

							if(this->getNumPastaLine() == TWO_PASTALINE) {
								//support stemid index in two pastalines 20081015
								tmp_sub_charid_idx = strrpl(cAliOther, toupper(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
								strcat(tmp_charid_idx, tmp_sub_charid_idx);
								if(tmp_sub_charid_idx != NULL) {
									delete [] tmp_sub_charid_idx;
									tmp_sub_charid_idx = NULL;
								}
							}
							if(cAliOther != NULL) {
								delete [] cAliOther;
								cAliOther = NULL;
							}
						}
					} 
					else	//current is not empty stem, but at least previous is empty stem
					{
						//adding the empty-stem model.  20090211
						gap_len = curLeftStemPos - empty_start;//preLeftStemPos;

						//at the same time, do it for this current iteration, that is, adding the current stem part
						//add the next stem sequence/alignment info
						tmpSubFold = this->buildstemfoldedstructure(stemIdx, stemCandIdxArray[m], left);

						if(gap_len >= 0) {
							for(n=0; n<gap_len; n++)
								strcat(tmpFoldedstrucutre, EMPTY_STEM_LABEL);
							strcat(tmpFoldedstrucutre, tmpSubFold);

							//making k-empty stem model support two pasta line format 20081221
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								for(n=0; n<gap_len; n++)
									strcat(tmp_charid_idx, EMPTY_STEM_LABEL);

								tmp_sub_charid_idx = strrpl(tmpSubFold, toupper(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
								strcat(tmp_charid_idx, tmp_sub_charid_idx);
								if(tmp_sub_charid_idx != NULL) {
									delete [] tmp_sub_charid_idx;
									tmp_sub_charid_idx = NULL;
								}
							}
						}
						else
						{	//adding allowing overlap when checking position confliction 20081221
							//overlap case
							gap_len = abs(gap_len);
							for(n=0; n<gap_len; n++)
								tmpFoldedstrucutre[strlen(tmpFoldedstrucutre)-1-n] = '*';
							string strtmp(tmpSubFold);
							strtmp = strtmp.substr(gap_len, strtmp.length() - gap_len);
							strcat(tmpFoldedstrucutre, strtmp.c_str());

							//making k-empty stem model support two pasta line format 20081221
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								for(n=0; n<gap_len; n++)
									tmp_charid_idx[strlen(tmp_charid_idx)-1-n] = '*';

								tmp_sub_charid_idx = strrpl(tmpSubFold, toupper(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
								string strcharidx(tmp_sub_charid_idx);
								strcharidx = strcharidx.substr(gap_len, strcharidx.length() - gap_len);
								strcat(tmp_charid_idx, strcharidx.c_str());
								
								if(tmp_sub_charid_idx != NULL) {
									delete [] tmp_sub_charid_idx;
									tmp_sub_charid_idx = NULL;
								}
							}
						}
/*						
						//support stemid index in two pastalines 20081015
						if(this->getNumPastaLine() == TWO_PASTALINE) 
						{
							tmp_sub_charid_idx = strrpl(tmpSubFold, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
							strcat(tmp_charid_idx, tmp_sub_charid_idx);
							if(tmp_sub_charid_idx != NULL) {
								delete [] tmp_sub_charid_idx;
								tmp_sub_charid_idx = NULL;
							}
						}
*/						

//						//for debug
//						cout<<tmpFoldedstrucutre<<endl;

						if(tmpSubFold != NULL) {
							delete [] tmpSubFold;
							tmpSubFold = NULL;
						}
					}
				}
			}//end of if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
		} 
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;
			// 20100215
//			cout<<cur_gid<<"["<<stemIdx<<"]\tR("<<curLeftStemPos<<","<<curRightStemPos<<")"<<endl;

			//adding the empty-stem model.  20090211
			if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
			{
				if(preRightStemPos != NEGATIVE_ONE)
				{
					if(curLeftStemPos > preRightStemPos) 
					{
						if(curLeftStemPos > (preRightStemPos+1))
						{
							//add loop structure alignment info here.
							start	= preRightStemPos + 1;
							end		= curLeftStemPos - 1;
							for(n=start; n<=end; n++) {
								strcat(tmpFoldedstrucutre, ".");
								if(this->getNumPastaLine() == TWO_PASTALINE) {
									//support stemid index in two pastalines 20081015
									strcat(tmp_charid_idx, ".");
								}
							}
						}

						//add the next stem sequence/alignment info
						tmpSubFold = this->buildstemfoldedstructure(stemIdx, stemCandIdxArray[m], right);
						strcat(tmpFoldedstrucutre, tmpSubFold);

						//support stemid index in two pastalines 20081015
						if(this->getNumPastaLine() == TWO_PASTALINE) 
						{
							tmp_sub_charid_idx = strrpl(tmpSubFold, tolower(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
							strcat(tmp_charid_idx, tmp_sub_charid_idx);
							if(tmp_sub_charid_idx != NULL) {
								delete [] tmp_sub_charid_idx;
								tmp_sub_charid_idx = NULL;
							}
						}
						
//						//for debug
//						cout<<tmpFoldedstrucutre<<endl;

						if(tmpSubFold != NULL) {
							delete [] tmpSubFold;
							tmpSubFold = NULL;
						}
					} else {
						//overlap
						int	iOverlap = preRightStemPos - curLeftStemPos + 1;
						int	iTmpLength = strlen(tmpFoldedstrucutre);
						for(n=0; n<iOverlap; n++) {
							tmpFoldedstrucutre[iTmpLength-1-n] = '*';
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								//support stemid index in two pastalines 20081015
								tmp_charid_idx[iTmpLength-1-n] = '*';
							}
						}

						iTmpLength = strlen(candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
						char * cAliOther = new char[iTmpLength - iOverlap + 1];
						for(n=0; n<(iTmpLength - iOverlap); n++) {
							cAliOther[n] = candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment[n+iOverlap];
						}
						cAliOther[n] = '\0';
						strcat(tmpFoldedstrucutre, cAliOther);
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							//support stemid index in two pastalines 20081015
							tmp_sub_charid_idx = strrpl(cAliOther, tolower(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
							strcat(tmp_charid_idx, tmp_sub_charid_idx);
							if(tmp_sub_charid_idx != NULL) {
								delete [] tmp_sub_charid_idx;
								tmp_sub_charid_idx = NULL;
							}
						}
						
//						//for debug
//						cout<<tmpFoldedstrucutre<<endl;

						if(cAliOther != NULL) {
							delete [] cAliOther;
							cAliOther = NULL;
						}
					}
				}
				else	//current is not empty stem, but at least previous is empty stem
				{
					//adding the empty-stem model.  20090211
					gap_len = curLeftStemPos - empty_start;//preLeftStemPos;

					//at the same time, do it for this current iteration, that is, adding the current stem part
					//add the next stem sequence/alignment info
					tmpSubFold = this->buildstemfoldedstructure(stemIdx, stemCandIdxArray[m], right);

					if(gap_len >= 0) {
						for(n=0; n<gap_len; n++)
							strcat(tmpFoldedstrucutre, EMPTY_STEM_LABEL);
						strcat(tmpFoldedstrucutre, tmpSubFold);

						//making k-empty stem model support two pasta line format 20081221
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							for(n=0; n<gap_len; n++)
								strcat(tmp_charid_idx, EMPTY_STEM_LABEL);

							tmp_sub_charid_idx = strrpl(tmpSubFold, tolower(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
							strcat(tmp_charid_idx, tmp_sub_charid_idx);
							if(tmp_sub_charid_idx != NULL) {
								delete [] tmp_sub_charid_idx;
								tmp_sub_charid_idx = NULL;
							}
						}
					}
					else
					{	//adding allowing overlap when checking position confliction 20081221
						//overlap case
						gap_len = abs(gap_len);
						for(n=0; n<gap_len; n++)
							tmpFoldedstrucutre[strlen(tmpFoldedstrucutre)-1-n] = '*';
						string strtmp(tmpSubFold);
						strtmp = strtmp.substr(gap_len, strtmp.length() - gap_len);
						strcat(tmpFoldedstrucutre, strtmp.c_str());

						//making k-empty stem model support two pasta line format 20081221
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							for(n=0; n<gap_len; n++)
								tmp_charid_idx[strlen(tmp_charid_idx)-1-n] = '*';

							tmp_sub_charid_idx = strrpl(tmpSubFold, tolower(this->pStem[stemIdx].charid), this->pStem[stemIdx].charid_idx);
							string strcharidx(tmp_sub_charid_idx);
							strcharidx = strcharidx.substr(gap_len, strcharidx.length() - gap_len);
							strcat(tmp_charid_idx, strcharidx.c_str());
												
							if(tmp_sub_charid_idx != NULL) {
								delete [] tmp_sub_charid_idx;
								tmp_sub_charid_idx = NULL;
							}
						}
					}
/*
					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE) 
					{
						tmp_sub_charid_idx = strrpl(tmpSubFold, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
						strcat(tmp_charid_idx, tmp_sub_charid_idx);
						if(tmp_sub_charid_idx != NULL) {
							delete [] tmp_sub_charid_idx;
							tmp_sub_charid_idx = NULL;
						}
					}
*/						
//					//for debug
//					cout<<tmpFoldedstrucutre<<endl;

					if(tmpSubFold != NULL) {
						delete [] tmpSubFold;
						tmpSubFold = NULL;
					}
				}
			}
		}

		if(curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE)
		{
			preLeftStemPos	= curLeftStemPos;
			preRightStemPos	= curRightStemPos;
		} else {
			if(preLeftStemPos != NEGATIVE_ONE && preRightStemPos != NEGATIVE_ONE)
			{
				empty_start = preRightStemPos + 1;
				preLeftStemPos	= curLeftStemPos;
				preRightStemPos	= curRightStemPos;
			}
			else
			{
				preLeftStemPos	= empty_start;
				preRightStemPos	= curRightStemPos;
			}
		}
		
		pre_gid	= cur_gid;
	}//for
//	cout<<tmpFoldedstrucutre<<endl
//		<<cResultOne<<endl
//		<<cResultTwo<<endl;

	//support stemid index in two pastalines 20081015
	if(this->getNumPastaLine() == TWO_PASTALINE)
	{
		int i, length = strlen(tmpFoldedstrucutre);
		if(length != strlen(tmp_charid_idx)) {
			cout<<"Something wrong here DynamicP::buildFoldedStructureInfo()"<<endl;
			exit(0);
		}
		for(i=0; i<length; i++) {
			if((tmpFoldedstrucutre[i] == 'l' && tmp_charid_idx[i]=='l')	//ascii
			||(tmpFoldedstrucutre[i] == 'r' && tmp_charid_idx[i]=='r'))
			{
				tmpFoldedstrucutre[i] = '.';
				tmp_charid_idx[i] = '.';
			}
		}
		strcpy(foldedstrucutre, tmpFoldedstrucutre);
		strcpy(foldedstrucutre_charid_idx, tmp_charid_idx);
	} else {
		char * cResultOne = strrpl(tmpFoldedstrucutre, 'l', '.');
		char * cResultTwo = strrpl(cResultOne, 'r', '.');
		strcpy(foldedstrucutre, cResultTwo);
		if(cResultOne != NULL) {
			delete [] cResultOne;
			cResultOne = NULL;
		}
		if(cResultTwo != NULL) {
			delete [] cResultTwo;
			cResultTwo = NULL;
		}

		char * cCharidOne = strrpl(tmp_charid_idx, 'l', '.');
		char * cCharidTwo = strrpl(cCharidOne, 'r', '.');
		strcpy(foldedstrucutre_charid_idx, cCharidTwo);
		if(cCharidOne != NULL) {
			delete [] cCharidOne;
			cCharidOne = NULL;
		}
		if(cCharidTwo != NULL) {
			delete [] cCharidTwo;
			cCharidTwo = NULL;
		}
	}
}

//void DynamicP::printFoldedStructureInfo(bool bDebug)
//{
//	if(bDebug)
//		cout<<"#";
//	cout<<foldedstrucutre<<endl;
//}
void DynamicP::printFoldedStructureInfo(int start, int end, int real_start, int real_end, ofstream &outfile)
{
	int i, length = strlen(foldedstrucutre);
	if(length != strlen(foldedstrucutre_charid_idx)) {
		cout<<"Something wrong here DynamicP::printFoldedStructureInfo()"<<endl;
		exit(0);
	}
	for(i=0; i<length; i++) {
		if((foldedstrucutre[i] == 'm' && foldedstrucutre_charid_idx[i]==109)	//ascii
		||(foldedstrucutre[i] == 'd' && foldedstrucutre_charid_idx[i]==100)
		||(foldedstrucutre[i] == 'i' && foldedstrucutre_charid_idx[i]==105))
		{
			foldedstrucutre[i] = '.';
			foldedstrucutre_charid_idx[i] = '.';
		}
	}
	//detect memory leak problem in rnatops  20090815
	char *c_foldseq = strextract(cWindowSequence, start, end);

	string foldedseq = c_foldseq;
	//support stemid index in two pastalines 20081015
	if(this->getNumPastaLine() == TWO_PASTALINE)
	{
		formatStringOutput(	string("Folded structure"), 
							string(foldedstrucutre),
							string(foldedstrucutre_charid_idx),		
							foldedseq, 
							real_start,
							real_end,
							1,
							outfile);
		oneHit.setFoldedStruPastaCharId(foldedstrucutre_charid_idx);
	} else {
		formatStringOutput(	string("Folded structure"), 
							string(foldedstrucutre),
							string(""),
							foldedseq, 
							real_start,
							real_end,
							0,
							outfile);
		oneHit.setFoldedStruPastaCharId("");
	}
	//support hit sorting 20081121
	oneHit.setFoldedStruSeq(foldedseq);
	oneHit.setFoldedStruAlign(foldedstrucutre);

	//detect memory leak problem in rnatops  20090815
	if(c_foldseq != NULL) {
		delete [] c_foldseq;
		c_foldseq = NULL;
	}
}
/*
 * print out some format like the following pasta structure info.
 *
 * 1........10........20........30....
 * ....|....|....|....|....|....|.....
 */
void DynamicP::printHeader(int startPos, int endPos, bool bDebug)
{
	int m;
	int iLength = endPos - startPos + 1;

	if(bDebug)
		cout<<"#";
	for(m=1; m<iLength+1; m++) 
	{
		if((m / 100) > 0)
		{
			if(m%10 == 0)
				cout<<m;
//			else if(m%10 != 1 && m%10 != 2)
//				cout<<".";
			else
			{
				if(m < 1000)
				{
					if(m%10 != 1 && m%10 != 2)
						cout<<".";
				} else if(m < 10000) {
					if(m%10 != 1 && m%10 != 2 && m%10 != 3)
						cout<<".";
				}
			}
		} else {
			if(m==1)
				cout<<m;
			else if(m%10 == 0)
				cout<<m;
			else if(m%10 != 1)
				cout<<".";
		}
	}
	cout<<endl;

	if(bDebug)
		cout<<"#";
	for(m=1; m<iLength+1; m++) {
		if(m % 5 ==0)
			cout<<"|";
		else
			cout<<".";
	}
	cout<<endl;

	if(bDebug)
		cout<<"#";
	cout<<"Folded structure"<<endl;	//

	if(bDebug)
		cout<<"#";
	for(m=startPos; m<endPos+1; m++)
		cout<<cWindowSequence[m];
	cout<<endl;
}
//build the folded structure for every stem based on the sequence/alignment of left/right arm of the stem
char * DynamicP::buildstemfoldedstructure(int stemIdx, int rank, int arm)
{
	int m, n, count, m_num;
	int left_length	= strlen(candidateStems[stemIdx].leftArm[rank].alignment);
	int right_length= strlen(candidateStems[stemIdx].rightArm[rank].alignment);
	char * left_seq					= new char[left_length + 1];;
	char * left_align				= new char[left_length + 1];
	char * right_seq				= new char[right_length + 1];
	char * right_align				= new char[right_length + 1];

	strcpy(left_seq, candidateStems[stemIdx].leftArm[rank].sequence);
	strcpy(left_align, candidateStems[stemIdx].leftArm[rank].alignment);

	strcpy(right_seq, candidateStems[stemIdx].rightArm[rank].sequence);
	strcpy(right_align, candidateStems[stemIdx].rightArm[rank].alignment);

	char * tmp_folded_structure = NULL;
	char * folded_structure		= NULL;
	m_num = 0;
	if(arm == 0)	//left arm
	{
		tmp_folded_structure= new char[left_length + 1];
		folded_structure	= new char[left_length + 1];

		strcpy(tmp_folded_structure, left_align);//candidateStems[stemIdx].leftArm[rank].alignment);
		
		//we need to build its left folded structure 
		//based on its left alignment and right sequence
		m_num = 0;
		for(n=0; n<right_length; n++)
		{
			//from right arm
			//find all the '-'. count the 'M' for every '-'
			if((right_seq[right_length-1-n] != '-') && (right_align[right_length-1-n] == tolower(this->pStem[stemIdx].charid)))//'M'))	
			{
				//count the num of matching
				m_num++;
			} 
			else
			{
//				if(right_seq[right_length-1-n] == '-')
				if(right_seq[right_length-1-n] == '-' && right_align[right_length-1-n] != 'd')	// 20100130
				{
					m_num++;	//including '-'
					//go to the left align and put '.' there
					m = 0;
					count = 0;
					while(count<m_num) {
						if(left_align[m] == toupper(this->pStem[stemIdx].charid))//'M')
							count++;
						m++;
					}
					tmp_folded_structure[m-1] = '.';
					m_num=0;	// 20100130
				}
			}
		}//for
		m = 0;
		for(n=0; n<left_length; n++)
		{
			if(left_seq[n] != '-') {
				folded_structure[m++] = tmp_folded_structure[n];
			}
		}
		folded_structure[m] = '\0';

		if(tmp_folded_structure != NULL)
		{
			delete [] tmp_folded_structure;
			tmp_folded_structure = NULL;
		}
	}
	else	//right arm
	{
		tmp_folded_structure= new char[right_length + 1];
		folded_structure	= new char[right_length + 1];
		strcpy(tmp_folded_structure, right_align);//candidateStems[stemIdx].rightArm[rank].alignment);

		//we need to build its right folded structure 
		//based on its right alignment and left sequence
		m_num = 0;
		for(n=0; n<left_length; n++)
		{
			//from right arm
			//find all the '-'. count the 'M' for every '-'
			if((left_seq[n] != '-') && (left_align[n] == toupper(this->pStem[stemIdx].charid)))//'M'))	
			{
				//count the num of matching
				m_num++;
			} 
			else
			{
//				if(left_seq[n] == '-')
				if(left_seq[n] == '-' && left_align[n] != 'd')	// 20100130
				{
					m_num++;	//including '-'
					//go to the left align and put '.' there
					m = right_length - 1;
					count = 0;
					while(count<m_num) {
						if(right_align[m] == tolower(this->pStem[stemIdx].charid))//'M')
							count++;
						m--;
					}
					tmp_folded_structure[m+1] = '.';
					m_num=0;	// 20100130
				}
			}
		}//for
		m = 0;
		for(n=0; n<right_length; n++)
		{
			if(right_seq[n] != '-') {
				folded_structure[m++] = tmp_folded_structure[n];
			}
		}
		folded_structure[m] = '\0';

		if(tmp_folded_structure != NULL)
		{
			delete [] tmp_folded_structure;
			tmp_folded_structure = NULL;
		}
	}
	//free the memory
	if(left_seq != NULL) {
		delete [] left_seq;
		left_seq = NULL;
	}
	if(left_align != NULL) {
		delete [] left_align;
		left_align = NULL;
	}
	if(right_seq != NULL) {
		delete [] right_seq;
		right_seq = NULL;
	}
	if(right_align != NULL) {
		delete [] right_align;
		right_align = NULL;
	}

	return folded_structure;
}

void DynamicP::printScoreForEveryComponent(Loop *pLoopModel, bool bDebug, ofstream &outfile)
{
	int m;
	int stemIdx, preStemIdx;
	int	stemnum = this->getNumstems();

	int preLeftStemPos = 0, preRightStemPos = 0;
	int curLeftStemPos = 0, curRightStemPos = 0;
	int pre_gid, cur_gid;
	pre_gid = cur_gid = 0;

	VTBSearch vtber(1);
	vtber.setAllowedInsNum(this->getAllowedNullLoopInsNum());
	int loopid, start, end;

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;

			if( preLeftStemPos != NEGATIVE_ONE && preRightStemPos != NEGATIVE_ONE
			&&	curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE) 
			{
				if(stemIdx == 0) {
					outfile<<"Score for Stem "<<this->pStem[stemIdx].charid<<" : "<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<endl;
				} else {
					loopid	= this->identifyLoopId(pre_gid, cur_gid);
					if(loopid != NEGATIVE_ONE)	//20100120
					{
						preStemIdx = (pre_gid - 1)/2;
						if((preRightStemPos+1) <= (curLeftStemPos-1)+1+this->getNumOfNtsOverlapBetweenStem())
						{
							if(this->getScoreOption() != SCORE_STEM)
							{
								//add loop structure alignment info here.
								start	= preRightStemPos + 1;
								end		= curLeftStemPos - 1;

								double dLoopScore = vtber.VTBSearchMax(&pLoopModel[loopid], cWindowSequence, start, end, NULL, NULL);

								//only consider the non-negative loop score
								if(this->getScoreOption() == SCORE_STEM_POSLOOP)
								{
									if(dLoopScore >= 0.0)
									{
										//support weight factor for stem and loop	20081109
										outfile<<"Score for Loop ("<<this->pStem[preStemIdx].charid<<"->"<<this->pStem[stemIdx].charid<<"): "<<this->getLoopWeight()*dLoopScore<<endl;
									}
								} else {
									outfile<<"Score for Loop ("<<this->pStem[preStemIdx].charid<<"->"<<this->pStem[stemIdx].charid<<"): "<<this->getLoopWeight()*dLoopScore<<endl;
								}
							}
						} else {
							outfile<<"Something wrong in computing the score"<<endl;
						}
					}
					//add the next stem sequence/alignment info
					outfile<<"Score for Stem "<<this->pStem[stemIdx].charid<<" : "<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<endl;
				}
			}
		} 
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;
			
			if( preLeftStemPos != NEGATIVE_ONE && preRightStemPos != NEGATIVE_ONE
			&&	curLeftStemPos != NEGATIVE_ONE && curRightStemPos != NEGATIVE_ONE) 
			{
				loopid	= this->identifyLoopId(pre_gid, cur_gid);
				if(loopid != NEGATIVE_ONE)	//20100120
				{
					preStemIdx = (pre_gid - 1)/2;
					if((preRightStemPos+1) <= (curLeftStemPos-1)+1+this->getNumOfNtsOverlapBetweenStem())
					{
						if(this->getScoreOption() != SCORE_STEM)
						{
							//add loop structure alignment info here.
							start	= preRightStemPos + 1;
							end		= curLeftStemPos - 1;

							double dLoopScore = vtber.VTBSearchMax(&pLoopModel[loopid], cWindowSequence, start, end, NULL, NULL);

							//only consider the non-negative loop score
							if(this->getScoreOption() == SCORE_STEM_POSLOOP)
							{
								if(dLoopScore >= 0.0)
								{
									//support weight factor for stem and loop	20081109
									outfile<<"Score for Loop ("<<this->pStem[preStemIdx].charid<<"->"<<this->pStem[stemIdx].charid<<"): "<<this->getLoopWeight()*dLoopScore<<endl;
								}
							} else {
								outfile<<"Score for Loop ("<<this->pStem[preStemIdx].charid<<"->"<<this->pStem[stemIdx].charid<<"): "<<this->getLoopWeight()*dLoopScore<<endl;
							}
						}
					} else {
						outfile<<"Something wrong in computing the score"<<endl;
					}
				}
			}
		}
		preLeftStemPos	= curLeftStemPos;
		preRightStemPos	= curRightStemPos;
		
		pre_gid	= cur_gid;
	}//for
}

int DynamicP::getJumpStrategy()
{
	return this->iJumpStrategy;
}
void DynamicP::setJumpStrategy(int strategy)
{
	this->iJumpStrategy = strategy;
}

int DynamicP::getStepSize()
{
	return this->iStepSize;
}
void DynamicP::setStepSize(int stepsize)
{
	this->iStepSize = stepsize;
}
double DynamicP::getCpuTimeAll()
{
	return this->cpu_time_hours_all;
}
double DynamicP::getCpuTimePreprocessing()
{
	return this->cpu_time_hours_preprocessing;
}
double DynamicP::getCpuTimeDP()
{
	return this->cpu_time_hours_dp;
}
double DynamicP::getScoreThresInJump()
{
	return this->dScoreThresInJump;
}
void DynamicP::setScoreThresInJump(double dScoreThres)
{
	this->dScoreThresInJump = dScoreThres;
}
//flag of searching reverse genome sequence
int	DynamicP::getFlagSearchReverseGenome()
{
	return this->iFlagSearchReverseGenome;
}

void DynamicP::setFlagSearchReverseGenome(int flag)
{
	this->iFlagSearchReverseGenome = flag;
//	if(flag >= 1)
//		this->iFlagSearchReverseGenome = 1;
//	else
//		this->iFlagSearchReverseGenome = 0;
}

void DynamicP::setSearchParams(int		top_k,
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
							   int		max_empty_stem_num	//20090628
							   )
{
	//detect memory leak problem in rnatops  20090815
	//limiting the number of empty-stem.  20090628
	if(max_empty_stem_num > 0)
	{
		this->setFlagSupportEmptyStem(TRUE);
		this->setMaxEmptyStemNum(max_empty_stem_num);
	} else {
		this->setFlagSupportEmptyStem(FALSE);
		this->setMaxEmptyStemNum(ZERO);	//make sure ept=0 also work.  20100312
	}

	//adding the empty-stem model.  20090211
	if(this->getFlagSupportEmptyStem() == TRUE)
		this->setLeadingNum(top_k + 1);	//adding one empty stem to every stem candidates
	else
		this->setLeadingNum(top_k);	//[k value]

	this->setThreshold(threshold);	//[threshold]
	this->setNumOfNtsOverlapBetweenStem(num_nts_overlap_between_stem);	//[number of nts in stem when overlap is allowed]
	this->setFlagMergeCandInPreprocess(iFlagMergeCandInPreprocess);	//whether taking the merge-candidate strategy in preprocessing
	this->setFlagCandwithShortestLength(iFlagCandwithShortestLength);
	this->setShiftNumMergeCand(iShiftNumMergeCand);		//[num of shift when merge candidate in preprocessing]
	this->setAllowedNullLoopInsNum(iAllowedNullLoopInsNum);//[number of insertion allowed in null loop]
	this->setPcoeff(pcoeff);		//pcoeff
	this->setJumpStrategy(iJumpStrategy);	//[skip-and-jump strategy]
	this->setStepSize(iStepSize);			//[stepsize in skip-and-jump strategy]
	this->setScoreThresInJump(dScoreThresholdInJump);			//[score_filtering threshold in jump strategy]
	this->setFlagPrintStrAlignInfo(iFlagPrintStrAlignInfo);	//[structure alignment info -s|-n]
	this->setFlagPrintScoreInfo(iFlagPrintScoreInfo);	//
	this->setFlagPrintDebugInfo(iFlagPrintDebugInfo);			//[debug info -n|-d]
	this->setFlagPrintDebugInputInfo(iFlagPrintDebugInputInfo);
	this->setFlagPrintDebugTreeInfo(iFlagPrintDebugTreeInfo);
	this->setFlagPrintDebugDPWinInfo(iFlagPrintDebugDPWinInfo);
	this->setFlagPrintDebugTDInfo(iFlagPrintDebugTDInfo);
	this->setFlagPrintDebugDPSearchInfo(iFlagPrintDebugDPSearchInfo);

	this->setNumPastaLine(PastaLineNum);
	this->setScoreOption(score_option);
	this->setStemWeight(weight_stem);
	this->setLoopWeight(weight_loop);


	//limiting the number of empty-stem.  20090628
//	if(max_empty_stem_num > 0)
//	{
//		this->setFlagSupportEmptyStem(TRUE);
//		this->setMaxEmptyStemNum(max_empty_stem_num);
//	} else {
//		this->setFlagSupportEmptyStem(FALSE);
//	}
}

int	DynamicP::getNumHit()
{
	return this->numhit;
}

int	DynamicP::getPreHitNum()
{
	return this->pre_hit_num;
}

void DynamicP::setPreHitNum(int num)
{
	this->pre_hit_num = num;
}

int	DynamicP::getNumPastaLine()
{
	return this->num_pastaline;
}
void DynamicP::setNumPastaLine(int num)
{
	this->num_pastaline = num;
}

int DynamicP::getFlagFilterSearchBased()
{
	return this->flag_filtersearch_based;
}

void DynamicP::setFlagFilterSearchBased(int bValue)
{
	if(bValue == TRUE)
		this->flag_filtersearch_based = TRUE;
	else
		this->flag_filtersearch_based = FALSE;
}

int	DynamicP::getFilterHitBeginPos()
{
	return this->filterHitBeginPos;
}

void DynamicP::setFilterHitBeginPos(int begin_pos)
{
	this->filterHitBeginPos = begin_pos;
}

int	DynamicP::getScoreOption()
{
	return this->score_option;
}

void DynamicP::setScoreOption(int option)
{
	this->score_option = option;
}

double DynamicP::getStemWeight()
{
	return this->stem_weight;
}

void DynamicP::setStemWeight(double weight)
{
	this->stem_weight = weight;
}

double DynamicP::getLoopWeight()
{
	return this->loop_weight;
}

void DynamicP::setLoopWeight(double weight)
{
	this->loop_weight = weight;
}

//support stub test 20081116
int	DynamicP::getFlagStubTest()
{
	return this->flagStubTest;
}
void DynamicP::setFlagStubTest(int flag)
{
	this->flagStubTest = flag;
}

int	DynamicP::getFlagInvalidCombination()
{
	return this->flagInvalidCombination;
}
void DynamicP::setFlagInvalidCombination(int flag)
{
	this->flagInvalidCombination = flag;
}

void DynamicP::initHitObject()
{
	oneHit.setPosStart(NEGATIVE_ONE);
	oneHit.setPosEnd(NEGATIVE_ONE);
}

SearchHit DynamicP::getHitObject()
{
	return oneHit;
}

//adding the empty-stem model.  20090211
void DynamicP::setMaxEmptyStemNum(int num)
{
	this->max_empty_stem_num = num;
}
int	DynamicP::getMaxEmptyStemNum()
{
	return this->max_empty_stem_num;
}

int DynamicP::getFlagSupportEmptyStem()
{
	return this->flagSupportEmptyStem;
}

void DynamicP::setFlagSupportEmptyStem(int flag)
{
	this->flagSupportEmptyStem = flag;
}

//do local structure alignment  20090723
int DynamicP::getFlagLocalStructureAlign()
{
	return this->flagLocalStructureAlign;
}

void DynamicP::setFlagLocalStructureAlign(int flag)
{
	this->flagLocalStructureAlign = flag;
}

//check and remove the position confliction 20081209
void DynamicP::printStemCandIdxPos()
{
	int i;
	int cur_gid;
	int cur_stem_id;
	int cur_candidate_id;
	int cur_pos_left, cur_pos_right;
	int stemnum = this->getNumstems();
	cur_gid = NEGATIVE_ONE;
//	double pre_score = 0.0;
	double cur_score = 0.0;
	for(i=0; i<2*stemnum; i++)
	{
		cur_gid = stemIdxArray[i];
		cur_stem_id = (cur_gid - 1)/2;
		cur_candidate_id = stemCandIdxArray[i];
        if(cur_gid % 2 != 0)
		{
			//left arm
			cur_pos_left = candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].a;
			cur_pos_right= candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].b;
			cur_score	 = candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].p_score;
		} else {
			//right arm
			cur_pos_left = candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].c;
			cur_pos_right= candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].d;
			cur_score	 = candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].p_score;
		}
		cout<<cur_candidate_id<<"\t["<<cur_pos_left<<","<<cur_pos_right<<"]\t"<<cur_score<<endl;
	}
}

//check and remove the position confliction 20081209
//method 1: my strategy is to always change the current candidate
//method 2: check every conflictions
//Note: needs to minus the score of removed stem
void DynamicP::checkStemCandIdxConfliction(double & score)
{
	int i=0, pre_i=0;
//	int left, right;
	int pre_gid=NEGATIVE_ONE, cur_gid=NEGATIVE_ONE;
	int pre_stem_id=NEGATIVE_ONE, cur_stem_id=NEGATIVE_ONE;
	int pre_candidate_id=0, cur_candidate_id=0;
	int pre_pos_left=NEGATIVE_ONE, pre_pos_right=NEGATIVE_ONE, cur_pos_left=NEGATIVE_ONE, cur_pos_right=NEGATIVE_ONE;
	int stemnum = this->getNumstems();
	int pair_gid = 0;
	int pair_id_candidxarray = 0;
//	pre_gid = cur_gid = NEGATIVE_ONE;
	bool changed = false;
	double pre_score = 0.0, cur_score = 0.0;
	while(true)
	{
		changed = false;
		//20081231 debug
		pre_pos_left = pre_pos_right = cur_pos_left = cur_pos_right = NEGATIVE_ONE;
		for(i=0; i<2*stemnum && !changed; i++)
		{
			pair_gid = 0;
//			left = right = 0;
			cur_gid = stemIdxArray[i];
			cur_stem_id = (cur_gid - 1)/2;
			cur_candidate_id = stemCandIdxArray[i];
			if(cur_gid % 2 != 0)
			{
				//left arm
//				left = 1;
				cur_pos_left = candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].a;
				cur_pos_right= candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].b;
			} else {
				//right arm
//				right = 1;
				cur_pos_left = candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].c;
				cur_pos_right= candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].d;
			}
			if(i != 0)
			{
				//
				if(pre_pos_left != NEGATIVE_ONE && pre_pos_right != NEGATIVE_ONE
				&& cur_pos_left != NEGATIVE_ONE && cur_pos_right != NEGATIVE_ONE)
				{
					//check the confliction with the previous pos
//					if(cur_pos_left <= pre_pos_right)
					//adding allowing overlap when checking position confliction 20081221
//					if(cur_pos_left <= (pre_pos_right - this->getNumOfNtsOverlapBetweenStem() + 1))
					if(cur_pos_left < (pre_pos_right - this->getNumOfNtsOverlapBetweenStem() + 1))
					{
/*
						//method 1
						//position confliction and let the stem related with the current gid be empty stem
						if(left)
						{
							//left and right
							stemCandIdxArray[i] = this->getLeadingNum() - 1;
							pair_gid = cur_gid + 1;
							pair_id_candidxarray = getIdxInCandIdxArray(pair_gid);
							stemCandIdxArray[pair_id_candidxarray] = this->getLeadingNum() - 1;
						} else {
							//right and left
							stemCandIdxArray[i] = this->getLeadingNum() - 1;
							pair_gid = cur_gid - 1;
							pair_id_candidxarray = getIdxInCandIdxArray(pair_gid);
							stemCandIdxArray[pair_id_candidxarray] = this->getLeadingNum() - 1;
						}
						changed = true;
*/
						//method 2
						//position confliction and let the stem with smaller score be empty stem
						//compare scores
						pre_stem_id = (pre_gid - 1)/2;
						pre_score = candidateStems[pre_stem_id].s_locinfo[pre_candidate_id].p_score;
						cur_score = candidateStems[cur_stem_id].s_locinfo[cur_candidate_id].p_score;
						if(pre_score > cur_score)
						{
							//remove current stem
							stemCandIdxArray[i] = this->getLeadingNum() - 1;
							if(cur_gid % 2 != 0)
							{
								//left and right
								pair_gid = cur_gid + 1;
							} else {
								//right and left
								pair_gid = cur_gid - 1;
							}
							score -= cur_score;
						} else {
							//remove previous stem
							stemCandIdxArray[pre_i] = this->getLeadingNum() - 1;
							if(pre_gid % 2 != 0)
							{
								//left and right
								pair_gid = pre_gid + 1;
							} else {
								//right and left
								pair_gid = pre_gid - 1;
							}
							score -= pre_score;
						}
						pair_id_candidxarray = getIdxInCandIdxArray(pair_gid);
						stemCandIdxArray[pair_id_candidxarray] = this->getLeadingNum() - 1;
						changed = true;
					}
				}
			}
			if(cur_pos_left != NEGATIVE_ONE && cur_pos_right != NEGATIVE_ONE)// && !changed)
			{
				pre_i = i;
				pre_gid = cur_gid;
				pre_candidate_id = cur_candidate_id;
				pre_pos_left = cur_pos_left;
				pre_pos_right= cur_pos_right;
			} else {
				//
			}
		}
		if(!changed) {
			break;
		}
	}
}

int DynamicP::getIdxInCandIdxArray(int gid)
{
	int i = -1;
	bool bFound = false;
	int stemnum = this->getNumstems();
	for(i=0; i<2*stemnum && !bFound; i++)
	{
		if(gid == stemIdxArray[i])
			bFound = true;
	}
	return (i-1);
}

//Sort by the nodeid in treebag 20090219
//by the sequence in the stemIdxArray
void DynamicP::arrangeNodelist(Tree_bag * root)
{
	int i;

	if(root == NULL)
		return;

	int stemIdxNum = 2 * this->getNumstems() + 2;

	Node_list * root_nhead;
	bool bFound;
	int currentNodeIdx;
	Node_list * new_child_nhead		= NULL;
	Node_list * new_child_nhead_tail= NULL;
	Node_list *new_child_nhead_head = NULL;
	
	for(i=0; i<stemIdxNum; i++)
	{
		bFound = false;
		currentNodeIdx = stemIdxFullArray[i];
		root_nhead = root->nhead;
		while(root_nhead && !bFound)
		{
			if(currentNodeIdx == root_nhead->g_node) {
				bFound = true;
			} 
			root_nhead = root_nhead->next;
		}
		if(bFound)
		{
			new_child_nhead = new Node_list[1];
			new_child_nhead->g_node	= currentNodeIdx;
			new_child_nhead->next	= NULL;
			if(new_child_nhead_tail != NULL) {
				new_child_nhead_tail->next	= new_child_nhead;
				new_child_nhead_tail		= new_child_nhead;
			} else {
				new_child_nhead_tail		= new_child_nhead;
			}
			if(new_child_nhead_head == NULL)
				new_child_nhead_head		= new_child_nhead;
		}
	}
	//free memory for nodelist before plugin the new node list.  20100125
	int nodenum = root->nodenum;
	this->free_treebag_nhead(root->nhead, nodenum);
	//
	root->nhead = new_child_nhead_head;

	int childnum = root->childnum;
	for(i=0; i<childnum; i++) {
		arrangeNodelist(root->pChild[i]);
	}
}

void DynamicP::postprocessTDDP(Tree_bag * root)
{
	this->freeTreeNodePath();
	this->freeStemIdxArray();
	this->freeStemIdxFullArray();	//Sort by the nodeid in treebag 20090219
	this->freeStemCandIdxArray();
	this->free_node_mapping();	//node_mapping[] is needed in the scan_genome function
								//opposite to the previous dp.allocate_node_mapping();
	this->postprocess_tree(root);
	this->removeTreeBagParent(root);
	this->free_tree(root);

}
