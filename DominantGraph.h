#ifndef _DOMINANT_GRAPH_H
#define _DOMINANT_GRAPH_H

#include "Structure.h"

#define OUTSIDE_THRESHOLD -1

struct node
{
	int ID;
	//FloatVec data;
	FloatVec data;
	float score;
	int ranklowerb;
	int rankupperb;

	IntVec parrents;
	IntVec children;

	node()
	{
		ranklowerb = -1;
		rankupperb = -1;
	}
};


typedef list<int> Layer;
typedef list< Layer > DomGaph;


void ReconstructDG(DomGaph & dominantgraph, map<int, node>& realdata, bool HasLowerbound);
list<cell> TopK(map<int, node> realdata, DomGaph dominantgraph, int k, Weight w);
int ProgressiveTopK(map<int, node>& realdata, DomGaph& dominantgraph, float missing_score, Weight w, ResultPool &resultpool  , int threshold_k, int ID, bool usedrule1 = true/*,list<cell> &result*/);
vector<cell> TopKDG(map<int, node>& realdata, DomGaph& dominantgraph, int k, Weight w);
bool compfn2(node a,node b);


#endif
