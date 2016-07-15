#ifndef _RANK_LOWER_BOUND_H
#define _RANK_LOWER_BOUND_H

#include "Structure.h"

#define positive 1
#define zero 0
#define negative -1


struct node{
	int v;
	struct node * next;

	node(){
		v = -1;
		next = NULL;
	}
}Node;

typedef map<int, list<int> > Graph;

Points FindInComp(Tuples allTuples, Tuple missing_tuple, IntVec position, int & num_dominate);
bool equal(Tuple tuple, Tuple missing_tuple, IntVec position);
bool Dominate(Tuple tuple, Tuple missing_tuple, IntVec position);
bool beDominate(Tuple tuple, Tuple missing_tuple, IntVec position);

Bounds Rankbounds(Points incomp_tuples, int num_dominate, int dim);

bool CanBothWorse(DoubleVec a, DoubleVec b);
bool CanBothBetter(DoubleVec a,DoubleVec b);
int VertexCover(Graph g );


#endif
