#ifndef _PROGRESSIVETOPK_H
#define _PROGRESSIVETOPK_H

#include "Structure.h"

typedef DoubleVec Vector;

list<Cell> TATopk(vector<vector<Cell> >& sortedTuples, map<int, vector<float> >& index, int k, Weight w);

int ProgressiveTopk(vector<vector<Cell> > sortedTuples, map<int, vector<float> > index, float missing_tuple_score, Weight w, list<Vector > &resultpool);

void PreComputeTA(char * filename, vector<vector<Cell> >& storebyrow, map<int, vector<float> >& index);


#endif
