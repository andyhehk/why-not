#include "ProgressiveTopk.h"

void PreComputeTA(char * filename, vector<vector<Cell> >& storebyrow, map<int, vector<float> >& index)
{
	Tuples allTuples;
	
	ReadData(filename, allTuples);

    int dimension;
    cout<<"Input Dimension:"<<endl;
    scanf("%d", &dimension);

	vector<list<Cell> >  sortedTuples(dimension);

    int position[5] = {0,1,2, 3, 4};

	for(Tuples::iterator iter = allTuples.begin(); iter!=allTuples.end(); iter++)
	{
		Tuple tuple = *iter;
		vector<float> needed_attr;
		
		for(int i=0; i<dimension; i++)
		{
			Cell temp;
			temp.ID = tuple.ID;
			try{
				temp.data = tuple.data.at(position[i]);
				needed_attr.push_back(temp.data);
				sortedTuples.at(i).push_back(temp);
			}
			catch(out_of_range outOfRange)
			{
				cout << "\n\nException1: "<<outOfRange.what()<<endl;
				cout<< tuple.data.size()<<","<<tuple.ID;
				return;
			}
		}

		index.insert(pair<int, vector<float> >(tuple.ID, needed_attr));
	}


	for(vector<list<Cell> >::iterator iter = sortedTuples.begin(); iter!=sortedTuples.end(); iter++)
	{
		(*iter).sort(compfn);

		list<Cell>  temp = *iter;

		if(storebyrow.empty())
		{
			for(list<Cell>::iterator iter1 = temp.begin(); iter1!=temp.end();iter1++)
			{
				vector<Cell> tempvector;
				Cell tempCell;

				tempCell.ID = iter1->ID;
				tempCell.data = iter1->data;

				tempvector.push_back(tempCell);

				storebyrow.push_back(tempvector);	
			}
		}

		else
		{
			unsigned int i = 0;
			for(list<Cell>::iterator iter1 = temp.begin(); iter1!=temp.end();iter1++, i++)
			{
				Cell tempCell;

				tempCell.ID = iter1->ID;
				tempCell.data = iter1->data;

				storebyrow[i].push_back(tempCell);
			}

			assert(i == temp.size());
		}
	}	


}

list<Cell> TATopk(vector<vector<Cell> >& sortedTuples, map<int, vector<float> >& index, int k, Weight w)
{
	//map<int, float> resultset;
	//map<int, float> buffer;
    list<Cell> resultset;
    reservable_priority_queue<Cell> buffer;
    buffer.reserve(100000);
    set<int> computed;
	float threshold = INFINITE;
	list<Cell> result;
    int num = 0;	
	for(vector<vector<Cell> >::iterator iter = sortedTuples.begin(); iter != sortedTuples.end(); iter++)
	{
		vector<Cell> temp = *iter;
		float tempThreshold = 0;
		
		for(vector<Cell>::size_type ix = 0; ix!=temp.size(); ix++)
		{
			tempThreshold += temp[ix].data * w.weighting[ix];

			int tempID = temp[ix].ID;

            if(computed.find(tempID)==computed.end())
			//if(buffer.find(tempID)==buffer.end()&&resultset.find(tempID)==resultset.end())
			{
				float score = 0;

				vector<float> attrs = index.find(tempID)->second;

				for(vector<float>::size_type iy = 0; iy!=attrs.size(); iy++)
				{
					score += attrs[iy]*w.weighting[iy];
				}
                computed.insert(tempID);
				//buffer.insert(pair<int, float>(tempID, score));
                Cell tempCell;
                tempCell.ID = tempID;
                tempCell.data = score;
                //buffer.push_back(tempCell);
                buffer.push(tempCell);
                num++;
			}
		}

		if(tempThreshold < threshold)
			threshold = tempThreshold;
	/*	
        for(list<Cell>::iterator iter2 = buffer.begin();iter2!=buffer.end();)
		//for(map<int, float>::iterator iter2 = buffer.begin(); iter2!=buffer.end();)
		{
            Cell tempCell = *iter2;

			if(tempCell.data >= threshold)
			{
				resultset.push_back(tempCell);				
				buffer.erase(iter2++);
			}
			else
				iter2++;
		}
*/
        while(!buffer.empty())
        {
            if(buffer.top().data>=threshold)
            {
                resultset.push_back(buffer.top());
                buffer.pop();
            }
            
            else
                break;
        }

		if(resultset.size()>=k)
		{
            /*
			Cell temp;
			for(map<int, float>::iterator iter3 = resultset.begin(); iter3!=resultset.end(); iter3++)
			{
				temp.ID = iter3->first;
				temp.data = iter3->second;
				result.push_back(temp);
			}
			result.sort(compfn);
			return result;
            */
            cout<<"Computed Num:"<<num<<endl;
            resultset.sort(compfn);
            return resultset;
		}
		
	}
    cout<<"Computed Num:"<<num<<endl;
	resultset.sort(compfn);
	return resultset;
}

bool beInBuffer(list<Cell> buffer, int ID)
{
	for(list<Cell>::iterator iter = buffer.begin(); iter!=buffer.end(); iter++)
	{
		if(iter->ID == ID)
			return true;
	}

	return false;
}

int ProgressiveTopk(vector<vector<Cell> > sortedTuples, map<int, vector<float> > index, float missing_tuple_score, Weight w, list<Vector >& resultpool)
{
	map<int, float> resultset;
	map<int, float> buffer;
	float threshold = INFINITE;
	list<Cell> result;

	resultpool.clear();
	
	for(vector<vector<Cell> >::iterator iter = sortedTuples.begin(); iter != sortedTuples.end(); iter++)
	{
		vector<Cell> temp = *iter;
		float tempThreshold = 0;
		
		for(vector<Cell>::size_type ix = 0; ix!=temp.size(); ix++)
		{
			tempThreshold += temp[ix].data * w.weighting[ix];

			int tempID = temp[ix].ID;

			if(buffer.find(tempID)==buffer.end()&&resultset.find(tempID)==resultset.end())
			{
				float score = 0;

				vector<float> attrs = index.find(tempID)->second;

				for(vector<float>::size_type iy = 0; iy!=attrs.size(); iy++)
				{
					score += attrs[iy]*w.weighting[iy];
				}

				buffer.insert(pair<int, float>(tempID, score));
			}
		}

		if(tempThreshold < threshold)
			threshold = tempThreshold;
		
		for(map<int, float>::iterator iter2 = buffer.begin(); iter2!=buffer.end();)
		{
			int ID = iter2->first;
			float score = iter2->second;

			if(score > missing_tuple_score)
			{
				resultset.insert(pair<int, float>(ID, score));	
				resultpool.push_back(index.find(ID)->second);
				buffer.erase(iter2++);
			}
			else
				iter2++;
		}

		if(missing_tuple_score>=threshold)
		{
			Cell temp;
			for(map<int, float>::iterator iter3 = resultset.begin(); iter3!=resultset.end(); iter3++)
			{
				temp.ID = iter3->first;
				temp.data = iter3->second;
				result.push_back(temp);
			}
			return result.size()+1;
		}
		
	}

	return result.size()+1;
}

int main(int argc, char** argv)
{
    char* filename = argv[1];

    clock_t start, stop;
    vector< vector<Cell> > storebyrow;
    map<int, vector<float> > index;
    
    int dimension = 3;
    int k;
    Weight weight;

    for(int i = 0; i<dimension;i++)
        weight.weighting.push_back(1.0/dimension);

    PreComputeTA(filename, storebyrow, index);

    cout<<"Input k:"<<endl;
    scanf("%d", &k);

    start = clock();
    list<Cell> result = TATopk(storebyrow, index, k, weight);
    stop = clock();
    int i = 0;
    for(list<Cell>::iterator iter = result.begin();iter != result.end(); iter++)
    {
        if(i>k)
            break;
        cout<<"ID:"<<iter->ID<<'\t'<<iter->data<<endl;
        i++;

    }   
    cout<<"Total running time:"<<float(stop-start)/CLOCKS_PER_SEC<<endl;


    return 0;
}

