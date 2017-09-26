#pragma once

#include <iostream>
#include <vector>
#include <map>

#include <stdlib.h>	//gcc3.4

template< typename T > class TRnaGraphNode
{
public:
	TRnaGraphNode( int nID, int nNext, int nPair, T& rPayload ) : m_nID( nID ), m_nNext( nNext ),
		m_nPair( nPair ), m_nNestLevel( 0 ), m_nCrossing( 0 ), m_bRemoved( false ), m_rPayload( rPayload )
	{
	}
	
	int ID(){return m_nID;}
	int Next(){return m_nNext;}
	int Pair(){return m_nPair;}
	int NestLevel(){return m_nNestLevel;}
	int Crossing(){return m_nCrossing;}
	bool Removed(){return m_bRemoved;}

	void ID( int id ){ m_nID = id; }
	void Next( int next ){ m_nNext = next; }
	void Pair( int pair ){ m_nPair = pair; }
	void NestLevel( int nestlevel ){ m_nNestLevel = nestlevel; }
	void Crossing( int crossing ){ m_nCrossing = crossing; }
	void Removed( bool removed ){ m_bRemoved = removed; }
	T& Payload()
	{ 
		return m_rPayload; 
	}

	TRnaGraphNode<T>& operator =(const TRnaGraphNode<T>&){return *this;}

private:
	int m_nID;		// This node's identity
	int m_nNext;	// Next node in sequence order (implicit really)
	int m_nPair;	// This node's pair
	int m_nNestLevel;
	int m_nCrossing;
	bool m_bRemoved;
	T& m_rPayload;
};

template< typename T, typename B = int > struct TRnaGraphTreeBag
{
	TRnaGraphTreeBag(){}
	TRnaGraphTreeBag<T,B>& operator =(const TRnaGraphTreeBag<T,B>&){return *this;}
	int m_nID;
	int m_nParentID;
	std::vector< TRnaGraphNode<T>* > m_vecpNodes;
	std::vector< TRnaGraphTreeBag<T, B>* > m_vecChildren;
	B* m_pPayload;
};

typedef struct TKNOTRECORD
{
	int m_nStart;
	int m_nEnd;
	int m_nStartBag;
	int m_nEndBag;
}
TKnotRecord;

// TRnaGraphTd - tree decomposition class, specialized for RNA graph
template< typename T, typename B = int > class TRnaGraphTd
{
public:
	TRnaGraphTd( std::vector< TRnaGraphNode<T> >& rvecNodes ) : m_prgtbRoot( 0 ),
		m_nCurrentBag(0), m_rvecNodes(rvecNodes), m_nodeStart(-1,0,-1,m_nStartNode),
		m_nodeEnd(rvecNodes.size(),-1,rvecNodes.size(),m_nEndNode),
		m_nStartNode(0), m_nEndNode(rvecNodes.size())
	{
		earlyknot.m_nStart = -1;
		earlyknot.m_nEnd = -1;
		earlyknot.m_nStartBag = -1;
		earlyknot.m_nEndBag = -1;
	}

	~TRnaGraphTd(void){}

	// Returns the tree-decomposition, i.e. a pointer to the root treebag
	TRnaGraphTreeBag<T,B>* Decomposition(){ return m_prgtbRoot; }

	// Prints all the treebags out in order
	void Print( std::ostream & ros = std::cout )
	{
		// Iterete over bags, printing info. as we go.
		int i = 0;
		ros << std::endl;
		while( i < m_vecpBag.size() )
		{
			TRnaGraphTreeBag<T,B>* pBag = m_vecpBag[i];
			ros << std::endl;
			ros << pBag->m_nID << "\t" << pBag->m_nParentID << "\t\t" << pBag->m_vecpNodes.size() << "\t";
			int j = 0;
			while( j < pBag->m_vecpNodes.size() )
			{
				int idx = pBag->m_vecpNodes[j]->ID();
				if( idx >= 0 && idx < m_rvecNodes.size() )
//					ros << pBag->m_vecNodes[j] << ':' << m_rvecNodes[idx].Payload() << ' ';
					ros << pBag->m_vecpNodes[j]->ID() << ' ';
				else if( -1 == idx )
					ros << m_nodeStart.ID() << ' ';
				else
					ros << m_nodeEnd.ID() << ' ';
				j++;
			}
			ros << "-- ";
			j = 0;
			while( j < pBag->m_vecChildren.size() )
			{
				ros << pBag->m_vecChildren[j]->m_nID << ' ';
				j++;
			}
			i++;
		}
		ros << std::endl;
	}

	
////////////////////////////////////////////////////////////////////////////////
// Starts the construction of a tree-decomposition for the graph
	int Construct()
	{
		int success = 1;

		// First we build a list of edges with non-zero number of crossings
		typename std::vector< TRnaGraphNode<T> >::iterator nit = m_rvecNodes.begin();
//		int nCurrentEnd = m_rvecNodes.back().m_nPair;
		std::vector< int > vecCrossedEdges;
		int nTotalCrossings = 0;

		while( m_rvecNodes.end() != nit )
		{
			typename std::vector< TRnaGraphNode<T> >::iterator nit2 = nit;
			nit2++;
			while( m_rvecNodes.end() != nit2 && nit->Pair() > nit2->ID() )
			{
				if( nit2->Pair() > nit->Pair() || nit2->Pair() < nit->ID() )
				{
					nit->Crossing(nit->Crossing() + 1 );
					nTotalCrossings++;
				}
				nit2++;
			}

			if( nit->ID() < nit->Pair() )
			{
				vecCrossedEdges.push_back( nit->ID() );
//				std::cout << nit->ID() << ", " << nit->Pair() << ", " << nit->Crossing() << std::endl;
			}

			++nit;
		}

		// Iteratively find a maximal crossing edge and eliminate it by zeroing its 
		// crossing count and then reducing the crossing count for any edge with an end
		// inside its range, until all counts are zero.
		
//std::cout << nTotalCrossings << std::endl;

		//bool flag = false;
		do
		{
			int MaxID = -1;
			int MaxCross = 0;

			// Find current maximal crossing edge
			nit = m_rvecNodes.begin();
			while( m_rvecNodes.end() != nit )
			{
				if( nit->Crossing() > MaxCross )
				{
					MaxCross = nit->Crossing();
					MaxID = nit->ID();
//std::cout << MaxID << ": " << m_rvecNodes[MaxID].Crossing() << std::endl;
				}
				else
//std::cout << nit->ID() << ": " << m_rvecNodes[nit->ID()].Crossing() << std::endl;
				++nit;
			}
			
			if( nTotalCrossings )
			{
                          if( !( -1 < MaxID ) )
                          {
//                            std::cout<< "Bad crossing count" << std::endl << std::endl;
                            exit(1);
                          }
//std::cout << nTotalCrossings << " - a"  << std::endl;
				nTotalCrossings -= m_rvecNodes[MaxID].Crossing();
				m_rvecNodes[MaxID].Crossing( 0 );
				m_rvecNodes[MaxID].Removed( true );
//std::cout << nTotalCrossings << " - A"  << std::endl;
	
				// Run over nodes within current edge span, deducting one from every edge
				// that has only one end within the span.
				int id = MaxID + 1;
				while( nTotalCrossings && id < m_rvecNodes[MaxID].Pair() )
				{
					{
					if( m_rvecNodes[id].Crossing() ||  m_rvecNodes[m_rvecNodes[id].Pair()].Crossing() )
						if ( m_rvecNodes[id].Pair() > m_rvecNodes[MaxID].Pair() )
						{
							m_rvecNodes[id].Crossing( m_rvecNodes[id].Crossing() - 1 );
							nTotalCrossings--;
//std::cout << id << ": " <<nTotalCrossings << " - b" << std::endl;
						}
						else if( m_rvecNodes[id].Pair() < m_rvecNodes[MaxID].ID() )
						{
							m_rvecNodes[m_rvecNodes[id].Pair()].Crossing( m_rvecNodes[m_rvecNodes[id].Pair()].Crossing() - 1 );
							nTotalCrossings--;
//std::cout << id << ": " <<nTotalCrossings << " - c" << std::endl;
						}
					}
					id++;
				}
			}
//			if(flag || 2==nTotalCrossings)
//			{
//                          if(flag)
//                          {
//				exit(0);
//                          }
//                          flag = true;
//			}
		}
		while( nTotalCrossings );
		
		// Build the graph
		
		// Root
		delete m_prgtbRoot;
		m_prgtbRoot = MakeBag();//new TRnaGraphTreeBag;
		m_prgtbRoot->m_nID = 0;
		m_nCurrentBag = 1;
		m_prgtbRoot->m_nParentID = 0;
		int s = -1;
		int t = static_cast<int>( m_rvecNodes.size() );
		m_nodeEnd.ID(t);
		m_nEndNode = t+1;
		m_rvecNodes.push_back( m_nodeEnd );
		// Root bag gets "start" and end "nodes"
		m_prgtbRoot->m_vecpNodes.push_back(&m_nodeStart);
		m_prgtbRoot->m_vecpNodes.push_back(&m_nodeEnd);
//std::cout << m_prgtbRoot->m_nID << "\t" << m_prgtbRoot->m_nParentID << "\t\t" << s << "," << m_nodeEnd.ID() << std::endl;		
		
		BuildTreeA( s, t, m_prgtbRoot );
		
		AddCrossingPaths();
		
		return success;
	}
	
	void Release()
	{
		while( m_vecpBag.size() != 0 )
		{
			delete m_vecpBag[m_vecpBag.size()-1];
			m_vecpBag.pop_back();
		}
	}


private:

	TRnaGraphTreeBag<T,B>* MakeBag()
	{
		TRnaGraphTreeBag<T,B>* pBag = new TRnaGraphTreeBag<T,B>;
		m_vecpBag.push_back( pBag );
		return pBag;
	}

TKNOTRECORD earlyknot;

/////////////////////////////////////////////////////////////////////////////////////
// 
	void BuildTreeA( int s, int t, TRnaGraphTreeBag<T,B>* pRoot )
	{
		if( s+1 == t )
		{
//std::cout << "++ HP: "  << s+1 << ", " << m_rvecNodes[s+1].Pair() << std::endl;
			return;
		}
		int left = s, mid = s+1, right = t;
		TRnaGraphTreeBag<T,B> * pChildBag, * pNextBag;
		pChildBag = MakeBag();
		if( -1 == left )
		{
			pChildBag->m_vecpNodes.push_back(&m_nodeStart);
		}
		else
		{
		if( s == t || s+1 > m_rvecNodes.size() || t+1 > m_rvecNodes.size() ) 
			return;
			pChildBag->m_vecpNodes.push_back(&m_rvecNodes[left]);
		}
		pChildBag->m_vecpNodes.push_back(&m_rvecNodes[mid]);
		pChildBag->m_vecpNodes.push_back(&m_rvecNodes[right]);
		pRoot->m_vecChildren.push_back(pChildBag);
		pChildBag->m_nParentID = pRoot->m_nID;
		pChildBag->m_nID = m_nCurrentBag++;
//std::cout << pChildBag->m_nID << "\t" << pChildBag->m_nParentID << "\ta\t" << left << "," << mid << "," << right << std::endl;		
		
		if( m_rvecNodes[mid].Removed() || m_rvecNodes[m_rvecNodes[mid].Pair()].Removed() )
		{
//std::cout << "++ PK start = "  << mid << ", PK end = " << m_rvecNodes[mid].Pair() << std::endl;
			TKNOTRECORD kr;
			if( s+1 < m_rvecNodes[mid].Pair() )
			{
				if( earlyknot.m_nStart == mid )
				{
					earlyknot.m_nEndBag = pChildBag->m_nID;
					m_vecKnot.push_back(earlyknot);
//std::cout << "Early knot found..." << std::endl;
				}
				else
				{
				kr.m_nStart = mid;
				kr.m_nEnd = m_rvecNodes[mid].Pair();
				kr.m_nStartBag = pChildBag->m_nID;
				m_vecKnot.push_back(kr);
				}
			}
			else
			{
				int i = 0;
				while( m_vecKnot[i].m_nStart != m_rvecNodes[mid].Pair() && i < m_rvecNodes.size() )
				{
					i++;
				}
				m_vecKnot[i].m_nEndBag = pChildBag->m_nID;
			}
			if( mid+1 == t ) 
			{
				return;
			}
			else
			{
//				TRnaGraphTreeBag<T,B>* pChildBagk;
				left = mid;
				mid++;
//				pChildBagk = MakeBag();
//				pChildBagk->m_vecpNodes.push_back(&m_rvecNodes[left]);
//				pChildBagk->m_vecpNodes.push_back(&m_rvecNodes[mid]);
//				pChildBagk->m_vecpNodes.push_back(&m_rvecNodes[right]);
//				pChildBag->m_vecChildren.push_back(pChildBagk);
//				pChildBagk->m_nParentID = pChildBag->m_nID;
//				pChildBagk->m_nID = m_nCurrentBag++;
//std::cout << pChildBagk->m_nID << "\t" << pChildBagk->m_nParentID << "\tk\t" << left << "," << mid << "," << right << std::endl;		
				
				BuildTreeA( left, right, pChildBag );
				pNextBag = pChildBag;
				//return;
			}
		}
		else if( m_rvecNodes[s+1].Pair() == s )
		{
//std::cout << "++ HP: "  << s+1 << ", " << m_rvecNodes[s+1].Pair() << std::endl;
			return;
		}
		else
		{
			if( m_rvecNodes[mid].Pair() > mid )
			{
				left = mid;
				right = m_rvecNodes[mid].Pair();
			TRnaGraphTreeBag<T,B>* pChildBags;
			pChildBags = MakeBag();
			pChildBags->m_vecpNodes.push_back(&m_rvecNodes[left]);
			pChildBags->m_vecpNodes.push_back(&m_rvecNodes[right]);
			pChildBags->m_vecpNodes.push_back(&m_rvecNodes[t]);
			pChildBag->m_vecChildren.push_back(pChildBags);
			pChildBags->m_nParentID = pChildBag->m_nID;
			pChildBags->m_nID = m_nCurrentBag++;
//std::cout << pChildBags->m_nID << "\t" << pChildBags->m_nParentID << "\ts\t" << left << "," << right << "," << t << std::endl;

			BuildTreeA( left, right, pChildBags );
			pNextBag = pChildBags;
			}
			else
			{
				return;
			}
		}
		if( right < t-1 /*&& t < m_rvecNodes.size()-1 *//* && m_rvecNodes[mid].Pair() > mid && m_rvecNodes[mid].Pair()+2 < t*/ ) 
		{
			left = right;
			mid = right+1;
			right = t;
//			TRnaGraphTreeBag<T,B>* pChildBagn;
//			pChildBagn = MakeBag();
//			pChildBagn->m_vecpNodes.push_back(&m_rvecNodes[left]);
//			pChildBagn->m_vecpNodes.push_back(&m_rvecNodes[mid]);
//			pChildBagn->m_vecpNodes.push_back(&m_rvecNodes[right]);
//			pChildBag->m_vecChildren.push_back(pChildBagn);
//			pChildBagn->m_nParentID = pChildBag->m_nID;
//			pChildBagn->m_nID = m_nCurrentBag++;
//std::cout << pChildBagn->m_nID << "\t" << pChildBagn->m_nParentID << "\tn\t" << left << "," << mid << "," << right << std::endl;

			//if( m_rvecNodes[mid].Pair()+2 < m_rvecNodes.size() )
			{
				BuildTreeA( left, right, pNextBag );
			}
		}
	}


	void BuildTree( int s, int t, TRnaGraphTreeBag<T,B>* pRoot )
	{
		TRnaGraphTreeBag<T,B>* pChildBag;
		
		if( s == -1 )
		{
			TRnaGraphTreeBag<T,B>* pChildBags;
			pChildBag = MakeBag();
			pChildBag->m_vecpNodes.push_back(&m_nodeStart);
			pChildBag->m_vecpNodes.push_back(&m_rvecNodes[s+1]);
			pChildBag->m_vecpNodes.push_back(&m_nodeEnd);
			pRoot->m_vecChildren.push_back(pChildBag);
			pChildBag->m_nParentID = pRoot->m_nID;
			pChildBag->m_nID = m_nCurrentBag++;
//std::cout << pChildBag->m_nID << "\t" << pChildBag->m_nParentID << "\ts\t" << s << "," << s+1 << "," << t << std::endl;
//			BuildTree( s+1, t, pChildBag );

			pChildBags = MakeBag();
			pChildBags->m_vecpNodes.push_back(&m_rvecNodes[s+1]);
//			pChildBags->m_vecpNodes.push_back(&m_rvecNodes[t-1]);
			pChildBags->m_vecpNodes.push_back(&m_rvecNodes[m_rvecNodes[s+1].Pair()]);
			pChildBags->m_vecpNodes.push_back(&m_nodeEnd);
			pChildBag->m_vecChildren.push_back(pChildBags);
			pChildBags->m_nParentID = pChildBag->m_nID;
			pChildBags->m_nID = m_nCurrentBag++;
//std::cout << pChildBags->m_nID << "\t" << pChildBags->m_nParentID << "\ts\t" << s+1 << "," << m_rvecNodes[s+1].Pair() << "," << t << std::endl;

//			if( m_rvecNodes[t-1].Removed() || m_rvecNodes[m_rvecNodes[t-1].Pair()].Removed() )
//			{
//std::cout << "++ PK start = "  << t-1 << ", PK end = " << m_rvecNodes[t-1].Pair() << std::endl;
//			TKNOTRECORD kr;
//			if( t-1 < m_rvecNodes[t-1].Pair() )
//			{
//				kr.m_nStart = t-1;
//				kr.m_nEnd = m_rvecNodes[t-1].Pair();
//				kr.m_nStartBag = pChildBags->m_nID;
//				m_vecKnot.push_back(kr);
//			}
//			else
//			{
//				int i = 0;
//				while( i < m_vecKnot.size() && m_vecKnot[i].m_nStart != m_rvecNodes[t-1].Pair() )
//				{
//					i++;
//				}
//				if( i < m_vecKnot.size() )
//				{
//				 m_vecKnot[i].m_nEndBag = pChildBags->m_nID;
//				}
//                                else
//                                {
//				earlyknot.m_nEnd = t-1;
//				earlyknot.m_nStart = m_rvecNodes[t-1].Pair();
//				earlyknot.m_nStartBag = pChildBags->m_nID;
//				m_vecKnot.push_back(kr);
//                                }
//			}
//			}
			
			BuildTree( s+1, m_rvecNodes[s+1].Pair(), pChildBags );
//			BuildTree( s+1, t-1, pChildBags );
			if( m_rvecNodes[s+1].Pair()+1 < t) 
			{
				BuildTree( m_rvecNodes[s+1].Pair(), t, pChildBags );
			}

			return;
		}

		if( s > m_rvecNodes.size() || t > m_rvecNodes.size() ) return;
		
		if( m_rvecNodes[s+1].Removed() || m_rvecNodes[m_rvecNodes[s+1].Pair()].Removed() )
		{
			// Start pk path
//			m_nPkStartNode = s+1;
//			m_nPkEndNode = m_rvecNodes[s+1].Pair();
			pChildBag = MakeBag();
			pChildBag->m_vecpNodes.push_back(&m_rvecNodes[s]);
			pChildBag->m_vecpNodes.push_back(&m_rvecNodes[s+1]);
			pChildBag->m_vecpNodes.push_back(&m_rvecNodes[t]);
			pRoot->m_vecChildren.push_back(pChildBag);
			pChildBag->m_nParentID = pRoot->m_nID;
			pChildBag->m_nID = m_nCurrentBag++;
//std::cout << pChildBag->m_nID << "\t" << pChildBag->m_nParentID << "\te\t" << s << "," << s+1 << "," << t << std::endl;		
//std::cout << "++ PK start = "  << s+1 << ", PK end = " << m_rvecNodes[s+1].Pair() << std::endl;
			TKNOTRECORD kr;
			if( s+1 < m_rvecNodes[s+1].Pair() )
			{
				if( earlyknot.m_nStart == s+1 )
				{
					earlyknot.m_nEndBag = pChildBag->m_nID;
					m_vecKnot.push_back(earlyknot);
//std::cout << "Early knot found..." << std::endl;
				}
				else
				{
				kr.m_nStart = s+1;
				kr.m_nEnd = m_rvecNodes[s+1].Pair();
				kr.m_nStartBag = pChildBag->m_nID;
				m_vecKnot.push_back(kr);
				}
			}
			else
			{
				int i = 0;
				while( m_vecKnot[i].m_nStart != m_rvecNodes[s+1].Pair() && i < m_rvecNodes.size() )
				{
					i++;
				}
				m_vecKnot[i].m_nEndBag = pChildBag->m_nID;
			}

			if( s+2 < t) 
			{
//				BuildTree( s+1, t, pChildBag );
			}
//			BuildTree( s+1, m_rvecNodes[s+1].Pair(), pChildBag1a, m_rvecNodes );
			//return;
		}
		else if( m_rvecNodes[s].Pair() == t )
		{
			TRnaGraphTreeBag<T,B>* pChildBag1;
			if( t > s+1 )
			{
				pChildBag = MakeBag();
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[s]);
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[s+1]);
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[t]);
				pRoot->m_vecChildren.push_back(pChildBag);
				pChildBag->m_nParentID = pRoot->m_nID;
				pChildBag->m_nID = m_nCurrentBag++;
//std::cout << pChildBag->m_nID << "\t" << pChildBag->m_nParentID << "\ta1\t" << s << "," << s+1 << "," << t << std::endl;

				pChildBag1 = MakeBag();
				pChildBag1->m_vecpNodes.push_back(&m_rvecNodes[s+1]);
				pChildBag1->m_vecpNodes.push_back(&m_rvecNodes[m_rvecNodes[s+1].Pair()]);
				pChildBag1->m_vecpNodes.push_back(&m_rvecNodes[t]);
				pChildBag->m_vecChildren.push_back(pChildBag1);
				pChildBag1->m_nParentID = pChildBag->m_nID;
				pChildBag1->m_nID = m_nCurrentBag++;
//std::cout << pChildBag1->m_nID << "\t" << pChildBag1->m_nParentID << "\ta2\t" << s+1 << "," << m_rvecNodes[s+1].Pair() << "," << t << std::endl;
			
			
				BuildTree( s+1, m_rvecNodes[s+1].Pair(), pChildBag1 );
			}
			else
			{
//std::cout << "s=" << s << ", t=" << t << std::endl;
			}
			if( t > m_rvecNodes[s+1].Pair() + 1 )
			{
				BuildTree( m_rvecNodes[s+1].Pair(), t, pChildBag1 );
			}
		}
		else
		{
			if( m_rvecNodes[s].Removed() || m_rvecNodes[m_rvecNodes[s].Pair()].Removed() || m_rvecNodes[s].Pair() < s )
			{		
				TRnaGraphTreeBag<T,B>* pChildBag1;
				pChildBag = MakeBag();
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[s]);
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[s+1]);
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[t]);
				pChildBag->m_nParentID = pRoot->m_nID;
				pChildBag->m_nID = m_nCurrentBag++;
				pRoot->m_vecChildren.push_back(pChildBag);

//std::cout << pChildBag->m_nID << "\t" << pChildBag->m_nParentID << "\tc1\t" << s << "," << s+1 << "," << t << std::endl;		

				if( m_rvecNodes[s+1].Pair() < t )
				{
					pChildBag1 = MakeBag();
					pChildBag1->m_vecpNodes.push_back(&m_rvecNodes[s+1]);
					pChildBag1->m_vecpNodes.push_back(&m_rvecNodes[m_rvecNodes[s+1].Pair()]);
					pChildBag1->m_vecpNodes.push_back(&m_rvecNodes[t]);
					pChildBag->m_vecChildren.push_back(pChildBag1);
					pChildBag1->m_nParentID = pChildBag->m_nID;
					pChildBag1->m_nID = m_nCurrentBag++;
//std::cout << pChildBag1->m_nID << "\t" << pChildBag1->m_nParentID << "\tc2\t" << s+1 << "," << m_rvecNodes[s+1].Pair() << "," << t << std::endl;
				}
				else
				{
					pChildBag1 = pChildBag;
				}
				
				BuildTree( s+1, m_rvecNodes[s+1].Pair(), pChildBag1 );
				if( t > m_rvecNodes[s+1].Pair()+2 )
				{
					BuildTree( m_rvecNodes[s+1].Pair(), t, pChildBag1 );
				}
			}
			else
			{
				pChildBag = MakeBag();
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[s]);
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[m_rvecNodes[s].Pair()]);
				pChildBag->m_vecpNodes.push_back(&m_rvecNodes[t]);
				pChildBag->m_nParentID = pRoot->m_nID;
				pChildBag->m_nID = m_nCurrentBag++;

//std::cout << pChildBag->m_nID << "\t" << pChildBag->m_nParentID << "\tb\t" << s << "," << m_rvecNodes[s].Pair() << "," << t << std::endl;		

				pRoot->m_vecChildren.push_back(pChildBag);
				BuildTree( s, m_rvecNodes[s].Pair(), pChildBag );
				if( t > m_rvecNodes[s].Pair() + 1 )
				{
					BuildTree( m_rvecNodes[s].Pair(), t, pChildBag );
				}
			}
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////
// Adds the node id for the start of each crossing stem to the bags on the
// path connecting the bags in the graph containing the start and end.
	void AddCrossingPaths()
	{
//		std::cout << std::endl;
		for( int i = 0; i < m_vecKnot.size(); i++ )
		{
//			std::cout << "Knot " << i << " Nodes: " << m_vecKnot[i].m_nStart << " -> " << m_vecKnot[i].m_nEnd << ", Bags: " << m_vecKnot[i].m_nStartBag << " -> " << m_vecKnot[i].m_nEndBag << std::endl;
			
			TRnaGraphTreeBag<T,B>* pBag = m_vecpBag[m_vecKnot[i].m_nStartBag];
			TRnaGraphTreeBag<T,B>* pParent = m_vecpBag[pBag->m_nParentID];

//			for( typename std::vector< TRnaGraphTreeBag<T,B>* >::iterator it = pParent->m_vecChildren.begin(); it != pParent->m_vecChildren.end(); it++ )
//			{
//				std::cout << pBag->m_vecChildren.size() << std::endl;
//				if( pBag == *it )
//				{
//					pParent->m_vecChildren.erase(it);
//					for( int i = 0; i < pBag->m_vecChildren.size(); i++ )
//					{
//						pParent->m_vecChildren.push_back(pBag->m_vecChildren[i]);
//					}
//					break;
//				}
//			}

			bool finished = false;
			while( !finished )
			{
				if( pParent->m_vecpNodes[0]->ID() <= m_vecKnot[i].m_nStart && ( 0 == pParent->m_nID || pParent->m_vecpNodes[2]->ID() >= m_vecKnot[i].m_nEnd ) )
				{
					finished = true;
				}
				else
				{
					if( pParent->m_vecpNodes[0]->ID() != m_vecKnot[i].m_nStart )
						pParent->m_vecpNodes.push_back(&m_rvecNodes[m_vecKnot[i].m_nStart]);
					pParent = m_vecpBag[pParent->m_nParentID];
				}
			}
			
			pBag = m_vecpBag[m_vecKnot[i].m_nEndBag];
			pParent = m_vecpBag[pBag->m_nParentID];
			finished = false;
			
			pBag->m_vecpNodes.push_back(&m_rvecNodes[m_vecKnot[i].m_nStart]);
			while( !finished )
			{
				if( pParent->m_vecpNodes[0]->ID() <= m_vecKnot[i].m_nStart && pParent->m_vecpNodes[2]->ID() >= m_vecKnot[i].m_nEnd )
				{
					finished = true;
				}
				else
				{
					if( pParent->m_vecpNodes[0]->ID() != m_vecKnot[i].m_nStart )
						pParent->m_vecpNodes.push_back(&m_rvecNodes[m_vecKnot[i].m_nStart]);
					pParent = m_vecpBag[pParent->m_nParentID];
				}
			}
			if( pParent->m_vecpNodes[0]->ID() != m_vecKnot[i].m_nStart )
				pParent->m_vecpNodes.push_back(&m_rvecNodes[m_vecKnot[i].m_nStart]);
			
		}
	}

	void RemoveRedundantBags()
	{
		
	}

	TRnaGraphTd( const TRnaGraphTd& );
	TRnaGraphTd& operator=( const TRnaGraphTd& );

	TRnaGraphTreeBag<T,B>* m_prgtbRoot;
	int m_nCurrentBag;
	std::vector< TRnaGraphNode<T> >& m_rvecNodes; 
	std::vector< TRnaGraphTreeBag<T,B>* > m_vecpBag;
	std::vector< TKNOTRECORD > m_vecKnot;

	TRnaGraphTreeBag<T,B>* m_pCurrent;
	int m_nIterateChild;
	TRnaGraphNode<T> m_nodeStart;
	TRnaGraphNode<T> m_nodeEnd;
	int m_nStartNode;
	int m_nEndNode;
};
