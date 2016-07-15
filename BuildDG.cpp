#include "BuildDG.h"

bool HasLowerBound = false;

inline float GetSkyKey(tuple& tup,IntVec& SkyAttrs) {
	float skey=0;
	float maxval=0,minval=qry_sDOM;

	for (int iter=0;iter<SkyAttrs.size();iter++) {
		int j=SkyAttrs[iter];
		skey+=tup.sVal[j];

		maxval=max(maxval,tup.sVal[j]);
	}

	skey=100*maxval*SkyAttrs.size()+1*skey;

	return skey;
}

// *** OPS BEGIN ***
int STAT_OSP_ssum=0,STAT_OSP_smax=0,STAT_OSP_count=0,STAT_OSP_ROUNDS=0;



// memory sorting of tuples, in ascending order of key
// no I/O cost
// *** Efficient implemention of sorting, as we don't retrieve items from "vtup" frequently  [August 13, 2009]
void MemorySort(IntVec& vec) {
	static SPairStruct* xvec=NULL;
	if (xvec==NULL) {	// array initialized once!
		assert(qry_pNum>0);
		xvec=new SPairStruct[qry_pNum];
	} else
		assert(qry_pNum>=vec.size());

	int len=vec.size();
	for (int iter=0;iter<len;iter++) {
		int i=vec[iter];
		xvec[iter].id=i;
		xvec[iter].key=vtup[i].key;
	}
	sort(xvec,xvec+len,SPairTupComp());
	//qsort (xvec, len, sizeof(SPairStruct), compare_spair);

	for (int iter=0;iter<len;iter++)
		vec[iter]=xvec[iter].id;

	// [OLD]
	//sort(vec.begin(),vec.end(),AbstractTupComp());
}


void OSPSFS(IntVec& rslt,IntVec& vec,IntVec& SkyAttrs) {
	// note, when using SFS, the original sorting key gets destroyed
	rslt.clear();
	OSP_rslt.clear();	// NEW

	// sorting step
	
	for (int i=0;i<vec.size();i++) {	// for each object
		int cur_id=vec[i];
		//vtup[cur_id].tid;// 	compute index
		tuple& curtup=vtup[cur_id];
		curtup.address=0;		// 	compute address
		curtup.key=GetSkyKey(curtup,SkyAttrs);//  compute projection
				
		//cout<<"vec value:"<<vec[i]<<",";
		//cout<<"Tid:"<<vtup[cur_id].tid<<",";
		//cout<<"Key value:"<<vtup[cur_id].key<<endl;
	}
	
	MemorySort(vec);


	// init the global parameters for OSP_SkyAttrs
	OSP_SkyAttrs.clear();
	OSP_SkyAttrs.assign(SkyAttrs.begin(),SkyAttrs.end());
	MAXDMNT = (unsigned long)pow(2,OSP_SkyAttrs.size())-1;

	SNODE::Init();	// [NEW]

	SNODE *root=NULL;
	for (int iter=0;iter<vec.size();iter++) {	// for each object
		int cur_id=vec[iter];

		if (iter==0) {
			root = SNODE::CreateNode(cur_id,vtup[cur_id].address, NULL);
		} else {
			root->Dominate(cur_id,vtup[cur_id].address,true);
		}
	}

	int cur_OSP_sum=0,cur_OSP_count=0,cur_OSP_max=0;
	root->Traverse(cur_OSP_max,cur_OSP_sum,cur_OSP_count,0);
	STAT_OSP_smax+=cur_OSP_max;
	STAT_OSP_ssum+=cur_OSP_sum;
	STAT_OSP_count+=cur_OSP_count;
	STAT_OSP_ROUNDS++;


	// NOTE: no need to deallocate space for the new technique
	//if (root!=NULL)
	//	delete root;	// dellocate the space


	rslt.assign(OSP_rslt.begin(),OSP_rslt.end());
}

void ReadData(char* filename, IntVec SkyAttrs)
{
	ifstream file(filename, ifstream::in);

	if(!file.is_open()){
		cout<<"Can not open file!!!"<<endl;
		return;
	}

	else{
		string line;
		while(getline(file, line)){
			istringstream s(line);
			bool isID = true;

			tuple record;
			string field;
			node tempnode;
			DoubleVec data;

			while(getline(s, field, ' '))
			{
			 	if(isID)
			 	{
					record.tid = atoi(field.c_str());
					isID = false;
			 	}
				else
				{
					data.push_back(0-atof(field.c_str()));
					//record.sVal.push_back(0-atoi(field.c_str()));
					//tempnode.data.push_back(atoi(field.c_str()));
				}
			}

			for(IntVec::size_type ix = 0;ix != SkyAttrs.size();ix++)
			{
			//cout<<"Attribute "<<ix<<":"<<data[SkyAttrs[ix]]<<endl;
				record.sVal.push_back(data[SkyAttrs[ix]]);
				tempnode.data.push_back(0-data[SkyAttrs[ix]]);
			}
			//cout<<"data size"<<data.size();
			//exit(-1);
			vtup.push_back(record);

			if(HasLowerBound)
			{
				IntVec::size_type upperb_position =  data.size()-1;
				tempnode.rankupperb = 0-data[upperb_position];
				IntVec::size_type lowerb_position =  data.size()-2;
				tempnode.ranklowerb = 0-data[lowerb_position];
			}
			
			realdata.insert(pair<int, node>(record.tid, tempnode));
			
		}
		file.close();
		return;
	}

}

void WriteDGraph(char* filename, Layer layer)
{
	ofstream file(filename, ofstream::out);

	if(!file.is_open())
	{
		cout<<"can not open file!!!"<<endl;
		return;
	}

	for(Layer::iterator nodeID = layer.begin();nodeID!=layer.end();nodeID++)
	{
		node tempnode = realdata.find(*nodeID)->second;

		file<<"ID "<<*nodeID<<endl;
		
		file<<"Value";
		for(DoubleVec::size_type ix = 0; ix!=tempnode.data.size();ix++)
		{
			file<<' '<<tempnode.data[ix];
		}
		file<<endl;

		//file<<"Rank-Lower-Bound ";
		//if(HasLowerBound)
		//{
		//	file<<tempnode.ranklowerb;
		//}
		//file<<endl;

		//file<<"Rank-Upper-Bound ";
		//if(HasLowerBound)
		//{
		//	file<<tempnode.rankupperb;
		//}
		//file<<endl;
		
		file<<"Parrent";
		for(list<int>::iterator iter = tempnode.parrents.begin();iter!=tempnode.parrents.end();iter++)
		{
			file<<' '<<*iter;
		}
		file<<endl;

		file<<"Children";
		for(list<int>::iterator iter = tempnode.children.begin();iter!=tempnode.children.end();iter++)
		{
			file<<' '<<*iter;
		}

		file<<endl;
	}

	file.close();
}


void BuildGraph(DomGaph & dominantgraph, IntVec SkyAttrs)
{
	IntVec vec, rslt;

	for(int i=0; i!=vtup.size();i++){
		vec.push_back(i);
		tuple temp = vtup[i];
	}

	while(!vec.empty())
	{
		cout<<"Remain size:"<<vec.size()<<endl;
		OSPSFS(rslt, vec, SkyAttrs);

        set<int> tempresult;
		Layer tempLayer;
		for(IntVec::iterator iter = rslt.begin();iter!=rslt.end();iter++)
		{
			tempLayer.push_back(*iter+1);
            tempresult.insert(*iter);
		}

		dominantgraph.push_back(tempLayer);
        //if(dominantgraph.size()>100)
           // break;

		for(IntVec::iterator iter = vec.begin();iter!=vec.end();)
		{
			bool found = false;
            /*
			for(IntVec::iterator iter1 = rslt.begin();iter1!=rslt.end();iter1++){
				if((*iter)==(*iter1))
				{
					found = true;
					rslt.erase(iter1);
					break;
				}
			}
            */
            if(tempresult.find(*iter)!=tempresult.end())
                found = true;
			if(found)
				iter = vec.erase(iter);
			else
				iter++;
		}

	}
	cout<<"OSP is end!"<<endl;
	DomGaph::iterator parrent_layer = dominantgraph.begin();
	DomGaph::iterator child_layer = parrent_layer;
	child_layer++;
	int i = 0;
	for(;child_layer!=dominantgraph.end();)
	{
		cout<<"Remain layer:"<<dominantgraph.size()-i<<endl;
		i++;
		for(Layer::iterator parrentID = parrent_layer->begin();parrentID!=parrent_layer->end();parrentID++)
		{
			node parrentnode = realdata.find(*parrentID)->second;
			for(Layer::iterator childID = child_layer->begin();childID!=child_layer->end();childID++)
			{
				node childnode = realdata.find(*childID)->second;
				bool bedominate = true;
				
				for(DoubleVec::size_type ix = 0;ix!=parrentnode.data.size();ix++)
				{
					if(parrentnode.data[ix]<childnode.data[ix])
					{
						bedominate = false;
						break;
					}
				}
				
				if(bedominate)
				{
					realdata.find(*parrentID)->second.children.push_back(*childID);
					realdata.find(*childID)->second.parrents.push_back(*parrentID);
				}
			}
			
		}
		parrent_layer++;
		child_layer++;
	}
	
	int layer_num = 0;

	ofstream report("./DGraph/report.txt",ofstream::out);
	report<<"Number-of-Layer "<<dominantgraph.size()<<endl;
	report.close();
	
	for(DomGaph::iterator iter = dominantgraph.begin();iter!=dominantgraph.end();iter++)
	{
		Layer templayer = *iter;
		char filename[64];
		sprintf(filename, "./DGraph/Layer%d", layer_num);
		WriteDGraph(filename, templayer);
		layer_num++;
	}
}

int main(int argc, char** argv)
{
	char * filename = argv[1];

	HasLowerBound = false;

    int dimension;

    cout<<"Input dimension:"<<endl;
    scanf("%d", &dimension);

	DomGaph dominantgraph;
	IntVec SkyAttrs(dimension);
	for(int i = 0;i<dimension;i++)
        SkyAttrs[i] = i;
	//SkyAttrs[3] = 3;
	//SkyAttrs[4] = 4;
	//SkyAttrs[5] = 5;
	//SkyAttrs[6] = 6;
	//SkyAttrs[7] = 7;

	ReadData(filename, SkyAttrs);
	qry_sANum = dimension;
	qry_pNum = vtup.size();

	//int num_layer = 0;
	cout<<"Constructing dominant graph"<<endl;
	BuildGraph(dominantgraph, SkyAttrs);
	cout<<"Completed!!!"<<endl;
	
	return 0;
}

