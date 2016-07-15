#ifndef _BUILD_DOMINANT_GRAPH_H
#define _BUILD_DOMINANT_GRAPH_H

#include "Structure.h"


using namespace std;

struct tuple {
	int tid;
	float key;
	DoubleVec sVal;

	unsigned long address;
};

struct node
{
	DoubleVec data;
	float score;
	int ranklowerb;
	int rankupperb;

	list<int> parrents;
	list<int> children;

	node()
	{
		ranklowerb = -1;
		rankupperb = -1;
		score = 0;
	}
};

typedef list<int> Layer;
typedef list< Layer > DomGaph;


vector<tuple> vtup;
map<int, node > realdata;

int qry_pNum;
int qry_gANum,qry_sANum;
int qry_gDOM,qry_sDOM;

// stat. variables
long long num_dom_comp=0;
long num_nonleaf_domcomp=0;
long num_leaf_domcomp=0;

int num_IO=0;
int total_rslts=0;

int BlockSize=20;		// in terms of tuples
int MemorySize=1000;	// in terms of tuples


int DiskPageSize=4096;	// in terms of bytes	[NEW]




/// *** OSP STRUCTURES
IntVec OSP_SkyAttrs,OSP_rslt;
unsigned long MAXDMNT = 1;
bool IS_StrictDominance=false;	// [22 August 2009]


/// *** EFFICIENT SORTING

struct SPairStruct {
	int id;
	float key;
};

struct SPairTupComp {
	bool operator () (SPairStruct left,SPairStruct right) const
	{
		return (left.key < right.key);
	}
};

inline unsigned long PAddress(int base,int obj) {
	//assert(base && Obj);
	unsigned long PADD = 0;
	bool PREC = true;

	tuple& tupbase=vtup[base];
	tuple& tupobj=vtup[obj];

	bool equal = true;

	for (int index = 0; index < OSP_SkyAttrs.size(); index++) {
		int j=OSP_SkyAttrs[index];
		if(tupbase.sVal[j] != tupobj.sVal[j])
			equal = false;
		if (tupbase.sVal[j] < tupobj.sVal[j]) {
			PREC = false;
			PADD |= ((unsigned long)1<<index);
		} else if ((IS_StrictDominance==false)&&(tupbase.sVal[j] == tupobj.sVal[j]))
			PADD |= ((unsigned long)1<<index);
	}

	if(equal)
		PADD = 0;

	num_dom_comp++;	// update the counter
	if(PREC && PADD!=MAXDMNT)
		PADD = 0;

	return PADD;
}


struct SNODE;	// forward declaration
int ALLOC_csize=0;
vector<SNODE*> ALLOC_nodes;

struct SNODE {
	unsigned long address;
	long size;
	bool inSkyline;
	int pObj;
	SNODE* sibling;
	SNODE* child;


	SNODE() {
		sibling=NULL;
		child=NULL;
	}

	static void Init() {
		ALLOC_csize=0;
	}

	static SNODE* CreateNode(int pobj, unsigned long padd, SNODE * sb=NULL, SNODE * chld=NULL,bool inSky=true) {
		if (ALLOC_csize>=ALLOC_nodes.size()) {	// if there is not enough space
			ALLOC_nodes.push_back(new SNODE());
		}
		// reuse storage if there is enough space!

		//assert(ALLOC_size<ALLOC_max);
		SNODE* curnode=ALLOC_nodes[ALLOC_csize];
		ALLOC_csize++;	// move the position pointer

		curnode->address = padd;
		curnode->pObj = pobj;
		curnode->size = 1;
		curnode->sibling = sb;
		curnode->child = chld;
		curnode->inSkyline = inSky;
		if (inSky)
			OSP_rslt.push_back(pobj);

		return curnode;
	}

	void Traverse(int& max_level,int& sum_level,int& count,int cur_level=0) { // new
		count++;
		sum_level+=cur_level;
		max_level=max(max_level,cur_level);

		// find out MAX level tomorrow!
		if (child!=NULL)
			child->Traverse(max_level,sum_level,count,cur_level+1);
		if (sibling!=NULL)
			sibling->Traverse(max_level,sum_level,count,cur_level);
	}


	bool Dominate(int Obj, unsigned long padd, bool isLocated=false) {
		//drill down dominating partition nodes or skip over to next sibling
		if ((address|padd) == padd){
			unsigned long PADD = PAddress(pObj,Obj);
			if( PADD == MAXDMNT)
				return true;

			bool isL = ((address==padd) && isLocated);
			if(child && child->address <= PADD){
				if(child->Dominate(Obj,PADD,isL))
					return true;
				if(isL)
					return false;
			} else {
				if(isL){
					child = CreateNode(Obj,PADD,child);
					return false;
				}
			}
		}

		//access next neighbor sibling
		if(sibling && sibling->address <= padd){
			if (sibling->Dominate(Obj,padd,isLocated))
				return true;
		} else {
			if (isLocated && address<padd) {
				sibling = CreateNode(Obj,padd,sibling);
			}
		}

		return false;
	};

	~SNODE(){
		delete sibling;
		delete child;
	};

};

	


void MemorySort(IntVec& vec);

void OSPSFS(IntVec& rslt,IntVec& vec,IntVec& SkyAttrs);

void ReadData(char* filename);
void WriteDGraph(char* filename, Layer layer);

void BuildGraph(DomGaph & dominantgraph, IntVec SkyAttrs);

#endif

