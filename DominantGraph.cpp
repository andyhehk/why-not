#include "DominantGraph.h"

void ReconstructDG(DomGaph & dominantgraph , map<int, node> &realdata, bool HasLowerbound)
{
	ifstream report(".//DGraph//report.txt", ifstream::in);

	if(!report.is_open())
	{
		cout<<"Dominant graph report file is missing, can not reconstruct!!!"<<endl;
		return;
	}

	int num_layer;
	string line;
	string field;
	
	getline(report, line);
	istringstream s(line);	

	getline(s, field, ' ');
	//cout<<field<<endl;
	getline(s, field, ' ');
	num_layer = atoi(field.c_str());
	//cout<<num_layer<<endl;
	report.close();

	for(int i=0;i<num_layer;i++)
	{
		Layer templayer;

		char filename[64];

		sprintf(filename, ".//DGraph//Layer%d",i);

		ifstream IndexFile(filename, ifstream::in);

		if(!IndexFile.is_open())
		{
			cout<<"Can not open file:"<<filename<<endl;
			cout<<"May be some index file is missing"<<endl;
			return;
		}

		while(!IndexFile.eof())
		{
			int ID;
			node tempnode;

			getline(IndexFile, line);
            if(line.empty())
                continue;
			//cout<<"Line 1 "<<line<<endl;

			s.clear();
			s.str(line.c_str());
			//cout<<"String 1 "<<s<<endl;
			getline(s, field, ' ');
			getline(s, field, ' ');
			ID = atoi(field.c_str());
			//cout<<"ID "<<ID<<endl;
			
			getline(IndexFile, line);
			//cout<<"Line 2 "<<line<<endl;
			s.clear();
			s.str(line.c_str());
			//cout<<"String 2 "<<s<<endl;
			getline(s, field, ' ');
			//cout<<"Value ";
			while(getline(s, field, ' '))
			{
				//cout<<atof(field.c_str())<<' ';
				tempnode.data.push_back(atof(field.c_str()));
			}
			//cout<<endl;
			
			//getline(IndexFile, line);
			//s.clear();
			//if(HasLowerbound)
			//{
			//	s.str(line.c_str());
			//	getline(s, field, ' ');
				//cout<<field<<' ';
			//	getline(s, field, ' ');
				//cout<<atoi(field.c_str())<<endl;
			//	tempnode.ranklowerb = atoi(field.c_str());
			//}

            //getline(IndexFile, line);
			//s.clear();
			//if(HasLowerbound)
			//{
			//	s.str(line.c_str());
			//	getline(s, field, ' ');
				//cout<<field<<' ';
			//	getline(s, field, ' ');
				//cout<<atoi(field.c_str())<<endl;
			//	tempnode.rankupperb = atoi(field.c_str());
			//}

			getline(IndexFile, line);
			s.clear();
			s.str(line.c_str());
			getline(s, field, ' ');
			//cout<<"Parrents ";
			while(getline(s, field, ' '))
			{
				tempnode.parrents.push_back(atoi(field.c_str()));
				//cout<<atof(field.c_str())<<' ';
			}
			//cout<<endl;

			getline(IndexFile, line);
			s.clear();
			s.str(line.c_str());
			getline(s, field, ' ');
			//cout<<"Children ";
			while(getline(s, field, ' '))
			{
				tempnode.children.push_back(atoi(field.c_str()));
				//cout<<atof(field.c_str())<<' ';
			}
            //cout<<endl;

			tempnode.ID = ID;
			realdata.insert(pair<int, node>(ID, tempnode));
			templayer.push_back(ID);
		}
        //exit(-1);
		IndexFile.close();
		dominantgraph.push_back(templayer);
	}
	
}

//use rank lower bound
list<cell> TopK(map<int, node> realdata, DomGaph dominantgraph, int k, Weight w)
{
	//list<cell> result;

	/*
	ofstream file("result1.log",ofstream::out);
	if(!file.is_open())
	{
		cout<<"file is not open"<<endl;
		return result;
	}*/

	
	list<cell> resultset, candidatelist;
	map<int, list<node> > waitinglist;
	set<int> computed;

	DomGaph::iterator layerone = dominantgraph.begin();
	Layer one = *layerone;
	int num = 0;
	
	for(Layer::iterator iter = one.begin();iter!=one.end();iter++)
	{
		node tempnode = realdata.find(*iter)->second;
		
		if(tempnode.ranklowerb>1)
		{
			if(waitinglist.find(tempnode.ranklowerb)==waitinglist.end())
			{
				list<node> templist;
				templist.push_back(tempnode);
				waitinglist.insert(pair<int, list<node> >(tempnode.ranklowerb, templist));
			}
		
			else
			{
				waitinglist.find(tempnode.ranklowerb)->second.push_back(tempnode);
			}
			continue;
		}
		

		float score = 0;
		for(FloatVec::size_type ix = 0; ix!=tempnode.data.size();ix++)
		{
			score += tempnode.data[ix]*w.weighting[ix];
		}
		//tempnode.score = score;
		computed.insert(*iter);

		cell tempcell;
		tempcell.ID = *iter;
		tempcell.data = score;
		candidatelist.push_back(tempcell);

		//file<<"child ID:"<<iter->ID<<" is acctually computed in first layer"<<endl;

		num++;
	}

	candidatelist.sort(compfn);
	//node R = *(candidatelist.begin());
	//resultset.push_back(R);
	//candidatelist.pop_front();

	node R = realdata.find(candidatelist.begin()->ID)->second;
	resultset.push_back(*(candidatelist.begin()));

	cell tempcell;
	//tempcell.ID = R.ID;
	//tempcell.data = candidatelist.begin()->data;
	//result.push_back(tempcell);
	
	candidatelist.erase(candidatelist.begin());
	int n=1;
	



	while(resultset.size()<k)
	{	
		//file<<"R ID:"<<R.ID<<endl;
		bool updated = false;
		for(IntVec::iterator childID = R.children.begin();childID!=R.children.end();childID++)
		{
			if(computed.count(*childID))
				continue;

			node tempnode = realdata.find(*childID)->second;

			bool allParrent_in_RS = true;

			for(IntVec::iterator parrentID = tempnode.parrents.begin(); parrentID!= tempnode.parrents.end();parrentID++)
			{
				bool inRS = false;
				for(list<cell>::iterator iter = resultset.begin();iter!=resultset.end();iter++)
				{
					if(iter->ID==*parrentID)
						inRS = true;
				}
				if(!inRS)
				{
					allParrent_in_RS = false;
					break;
				}
			}

			if(!allParrent_in_RS)
				continue;

			if(tempnode.ranklowerb>n+1)
			{
				if(waitinglist.find(tempnode.ranklowerb)==waitinglist.end())
				{
					list<node> templist;
					templist.push_back(tempnode);
					waitinglist.insert(pair<int, list<node> >(tempnode.ranklowerb, templist));
				}

				else
				{
					waitinglist.find(tempnode.ranklowerb)->second.push_back(tempnode);
				}
				continue;
			}

			float score =0;
			for(FloatVec::size_type ix = 0; ix!=tempnode.data.size();ix++)
			{
				score += tempnode.data[ix]*w.weighting[ix];
			}
			//tempnode.score = score;
			computed.insert(*childID);
			//candidatelist.push_back(tempnode);

			tempcell.ID = tempnode.ID;
			tempcell.data = score;
			candidatelist.push_back(tempcell);
			updated = true;
			num++;

		}

		if(waitinglist.find(n+1)!=waitinglist.end())
		{
			list<node> templist = waitinglist.find(n+1)->second;
			for(list<node>::iterator iter = templist.begin();iter!=templist.end();iter++)
			{
				float score =0;
				for(FloatVec::size_type ix = 0; ix!=iter->data.size();ix++)
				{
					score += iter->data[ix]*w.weighting[ix];
				}
				//iter->score = score;
				computed.insert(iter->ID);
				//candidatelist.push_back(*(iter));
				tempcell.ID = iter->ID;
				tempcell.data = score;
				candidatelist.push_back(tempcell);
				updated = true;
				num++;
				//file<<"child ID:"<<iter->ID<<" is acctually computed by ranking lower estimated"<<endl;
			}
		}

		if(updated)
			candidatelist.sort(compfn);
		//cout<<"Result set size:"<<resultset.size()<<endl;
		//cout<<"Candidate list size:"<<candidatelist.size()<<endl;
		while(candidatelist.size()>k-n)
		{
			//cout<<"Pop up something"<<endl;
			candidatelist.pop_back();
		}
		
		//R = *(candidatelist.begin());
		//resultset.push_back(R);
		//candidatelist.pop_front();
		R = realdata.find(candidatelist.begin()->ID)->second;
		resultset.push_back(*(candidatelist.begin()));

		//tempcell.ID = R.ID;
		//tempcell.data = candidatelist.begin()->data;
		//result.push_back(tempcell);
		
		candidatelist.erase(candidatelist.begin());
		n++;


	}
	
	//resultset.sort(compfn2);
	/*
	for(list<cell>::iterator iter = resultset.begin();iter!=resultset.end();iter++)
	{
		file<<"ID:"<<iter->ID<<"\t"<<"score:"<<iter->data<<endl;
	}
	
	file<<"number of computed record"<<num<<endl;
	
	file.close();
	*/

	//cout<<"number of computed record"<<num<<endl;
	return resultset;
}

int ProgressiveTopK(map<int, node>& realdata ,DomGaph& dominantgraph, float missing_score, Weight w, ResultPool & resultpool , int threshold_k, int ID, bool usedrule1/*, list<cell> &result*/)
{
    if(threshold_k<0)
        return OUTSIDE_THRESHOLD;
    //missing_score = ((int)(missing_score*1000))/1000.0;
    resultpool.ID = ID;
    resultpool.result.clear();
    resultpool.result.reserve(100000);
	//result.clear();
	/*
	ofstream file("result3.log",ofstream::out);

	if(!file.is_open())
	{
		cout<<"file is not open"<<endl;
		return -1;
	}
	*/
    clock_t start, stop;
    start =  clock();
    reservable_priority_queue<cell> candidatelist;
    candidatelist.reserve(200000);
    //priority_queue<cell, deque<cell> >candidatelist;
	//map<int, list<node> > waitinglist;
	set<int> computed;
    set<int> result;
	DomGaph::iterator layerone = dominantgraph.begin();
	Layer one = *layerone;
	int num = 0;
	
    int dimension = w.weighting.size();


	for(Layer::iterator iter = one.begin();iter!=one.end();iter++)
	{
		node* tempnode = &(realdata.find(*iter)->second);
		

		float score = 0;
		for(int ix = 0; ix<dimension;ix++)
		{
			score += tempnode->data[ix]*w.weighting[ix];
		}
        //score = ((int)(score*1000))/1000.0;
		//tempnode.score = score;
		if(score>missing_score)
		{
			computed.insert(*iter);
			//candidatelist.push_back(tempnode);
			cell tempcell;
			tempcell.ID = tempnode->ID;
			tempcell.data = score;
			//candidatelist.push_back(tempcell);
            candidatelist.push(tempcell);
			num++;
	        resultpool.result.push_back(tempnode->data);
            if(num>threshold_k&&usedrule1)
                return OUTSIDE_THRESHOLD;
		}

	    //resultpool.result.push_back(tempnode->data);

		//file<<"child ID:"<<iter->ID<<" is acctually computed in first layer"<<endl;


	}

	if(num == 0)
	{
		return num+1;
	}
	
	//candidatelist.sort(compfn);
	//node R = *(candidatelist.begin());
	//resultset.push_back(R);
	//candidatelist.pop_front();
	//node R = realdata.find(candidatelist.begin()->ID)->second;
    node* R = &(realdata.find(candidatelist.top().ID)->second);
    result.insert(R->ID);
	//resultset.push_back(*(candidatelist.begin()));
    //resultset.push_back(candidatelist.top());
	//candidatelist.erase(candidatelist.begin());
    candidatelist.pop();
	
	int n=1;

	//resultpool.result.push_back(R->data);

	while(result.size()<=threshold_k)
	//while(true)
	{	
		//file<<"R ID:"<<R.ID<<endl;
		bool update = false;
		for(IntVec::iterator childID = R->children.begin();childID!=R->children.end();childID++)
		{
			if(computed.count(*childID))
				continue;

			node *tempnode = &(realdata.find(*childID)->second);

			bool allParrent_in_RS = true;

			for(IntVec::iterator parrentID = tempnode->parrents.begin(); parrentID!= tempnode->parrents.end();parrentID++)
			{
                /*
				bool inRS = false;
				for(list<cell>::iterator iter = resultset.begin();iter!=resultset.end();iter++)
				{
					if(iter->ID==*parrentID)
						inRS = true;
				}
				if(!inRS)
				{
					allParrent_in_RS = false;
					break;
				}
                */
                if(result.find(*parrentID)==result.end())
                {
                    allParrent_in_RS = false;
                    break;
                }
			}

			if(!allParrent_in_RS)
				continue;

			float score =0;
			for(int ix = 0; ix<dimension;ix++)
			{
				score += tempnode->data[ix]*w.weighting[ix];
			}
            //score = ((int)(score*1000))/1000.0;
			//tempnode.score = score;
			computed.insert(*childID);
			if(score>missing_score)
			{
				cell tempcell;
				tempcell.ID = tempnode->ID;
				tempcell.data = score;
				//candidatelist.push_back(tempcell);
                candidatelist.push(tempcell);
				update = true;
				num++;

			    resultpool.result.push_back(tempnode->data);
                if(num>threshold_k&&usedrule1)
                    return OUTSIDE_THRESHOLD;
			}
	        //resultpool.result.push_back(tempnode->data);

		}

		if(!candidatelist.empty())
		{
			//if(update)
				//candidatelist.sort(compfn);
			//cout<<"Result set size:"<<resultset.size()<<endl;
			//cout<<"Candidate list size:"<<candidatelist.size()<<endl;
		
			//R = *(candidatelist.begin());
			//resultset.push_back(R);
			//candidatelist.pop_front();

			//R = realdata.find(candidatelist.begin()->ID)->second;
			//resultset.push_back(*(candidatelist.begin()));
			//candidatelist.erase(candidatelist.begin());
            
            R = &(realdata.find(candidatelist.top().ID)->second);
            result.insert(R->ID);
            //resultset.push_back(candidatelist.top());
            candidatelist.pop();

			//resultpool.result.push_back(R->data);
		}
		else
			break;

	}
	
	//resultset.sort(compfn2);
	
//	for(list<cell>::iterator iter = resultset.begin();iter!=resultset.end();iter++)
//	{
//		cout<<"ID:"<<iter->ID<<"\t"<<"score:"<<iter->data<<endl;
//	}
//	
	//cout<<"number of computed record "<<num<<endl;

	//file.close();
	
	//result = resultset;
    stop = clock();
    //cout<<"Time:"<<float(stop-start)/CLOCKS_PER_SEC<<endl;
    //cout<<"threshold_k:"<<threshold_k<<endl;
    //cout<<"What is the progressive num:"<<resultpool.result.size()<<endl;
    //cout<<"Progressive ranking:"<<result.size()+1<<endl;
    //if(resultpool.size()>=(w.weighting.size()-1)*2)

	if(num>=threshold_k)
		return OUTSIDE_THRESHOLD;

	else
    {  
        //cout<<"threshold_k:"<<threshold_k<<endl;
	    //cout<<"number of computed record "<<num<<endl;
        //cout<<"What is the progressive num:"<<resultpool.result.size()<<endl;
        //cout<<"Progressive ranking:"<<result.size()+1<<endl;
		
        return result.size()+1;
    }
}



//without use rank lower bound
vector<cell> TopKDG(map<int, node>& realdata, DomGaph& dominantgraph, int k, Weight w)
{
	/*
	ofstream file("result2.log",ofstream::out);

	if(!file.is_open())
	{
		cout<<"file is not open"<<endl;
		return 0;
	}
	*/
	
	//cout<<"IN"<<endl;
	vector<cell> resultset;
    resultset.reserve(k);
    //priority_queue<cell, deque<cell> > candidatelist;  
    reservable_priority_queue<cell>  candidatelist;
    candidatelist.reserve(200000);
	//map<int, list<node> > waitinglist;
	set<int> computed;
    set<int> result;

	DomGaph::iterator layerone = dominantgraph.begin();
	Layer one = *layerone;
	int num = 0;

	clock_t start_find, stop_find;
	float time_for_find = 0;

    int dimension = w.weighting.size();
	
	for(Layer::iterator iter = one.begin();iter!=one.end();iter++)
	{

		node tempnode = realdata.find(*iter)->second;		

		float score = 0;
		for(int ix = 0; ix<dimension;ix++)
		{
			score += tempnode.data[ix]*w.weighting[ix];
		}
        //score = ((int)(score*1000))/1000.0;

		//tempnode.score = score;
		computed.insert(*iter);

		cell tempcell;
		tempcell.ID = tempnode.ID;
		tempcell.data = score;
		candidatelist.push(tempcell);
        //candidatelist.push_back(tempcell);
		//file<<"child ID:"<<iter->ID<<" is acctually computed in first layer"<<endl;

		num++;
	}

	//candidatelist.sort(compfn);
	//node R = *(candidatelist.begin());
	//node R = realdata.find(candidatelist.begin()->ID)->second;
	node R = realdata.find(candidatelist.top().ID)->second;
	//resultset.push_back(R);
    result.insert(R.ID);
	resultset.push_back(candidatelist.top());
	//candidatelist.erase(candidatelist.begin());
	candidatelist.pop();
	int n=1;

	
	while(resultset.size()<k)
	{	
		bool update = false;
		//file<<"R ID:"<<R.ID<<endl;
		for(IntVec::iterator childID = R.children.begin();childID!=R.children.end();childID++)
		{
			if(computed.count(*childID))
				continue;

			node tempnode = realdata.find(*childID)->second;

			//start = clock();
			bool allParrent_in_RS = true;

			for(IntVec::iterator parrentID = tempnode.parrents.begin(); parrentID!= tempnode.parrents.end();parrentID++)
			{
				//bool inRS = false;
				//for(vector<node>::iterator iter = resultset.begin();iter!=resultset.end();iter++)
				//node temp = realdata.find(*childID)->second;
				//for(list<cell>::iterator iter = resultset.begin();iter!=resultset.end();iter++)
				//{
				//	if(iter->ID==*parrentID)
				//		inRS = true;
				//}
                if(result.find(*parrentID)==result.end())
                {
                    allParrent_in_RS = false;
                    break;
                }
				//if(!inRS)
				//{
					//allParrent_in_RS = false;
					//break;
				//}
			}
			//stop = clock();
			//cout<<"Total time:"<<float(stop- start)*20000000000.0/ CLOCKS_PER_SEC <<endl;
			//exit(-1);
			
			if(!allParrent_in_RS)
				continue;

			float score =0;
			for(int ix = 0; ix<dimension;ix++)
			{
				score += tempnode.data[ix]*w.weighting[ix];
			}
            //score = ((int)(score*1000))/1000.0;
			//tempnode.score = score;
			computed.insert(*childID);

			cell tempcell;
			tempcell.ID = tempnode.ID;
			tempcell.data = score;
			//candidatelist.push_back(tempnode);
            candidatelist.push(tempcell);
			//candidatelist.push_back(tempcell);
			num++;
			update = true;

		}
		//start_find = clock();
		//if(update)
			//candidatelist.sort(compfn);
					
				
			//stop_find = clock();

			//time_for_find += float(stop_find- start_find)/ CLOCKS_PER_SEC;
			//cout<<"Time for find:"<<time_for_find<<endl;
			//exit(-1);
		//cout<<"Result set size:"<<resultset.size()<<endl;
		//cout<<"Candidate list size:"<<candidatelist.size()<<endl;
		//while(candidatelist.size()>k-n)
		//{
			//cout<<"Pop up something"<<endl;
			//candidatelist.pop_back();
		//}
		
		//R = *(candidatelist.begin());
		R = realdata.find(candidatelist.top().ID)->second;
		//R = realdata.find(candidatelist.begin()->ID)->second;
		//resultset.push_back(R);
        result.insert(R.ID);
		resultset.push_back(candidatelist.top());
		//resultset.push_back(*(candidatelist.begin()));
		candidatelist.pop();
		//candidatelist.erase(candidatelist.begin());
		n++;
	}
	
	//resultset.sort(compfn2);
	/*
	for(list<cell>::iterator iter = resultset.begin();iter!=resultset.end();iter++)
	{
		file<<"ID:"<<iter->ID<<"\t"<<"score:"<<iter->data<<endl;
	}
	
	
	
	//cout<<"Time for find:"<<time_for_find<<endl;
	
	file.close();
	*/
    //cout<<"The final candidate size:"<<candidatelist.size()<<endl;
    //cout<<"number of computed record "<<num<<endl;
    /*  
    int count = 0;
	for(list<cell>::iterator iter = resultset.begin();iter!=resultset.end();iter++)
	{
        count++;
		cout<<"ID:"<<iter->ID<<"\t"<<"score:"<<iter->data<<endl;
   //     if(count==100)
    //        cout<<"--------------------------------"<<endl;
    }
	*/
	return resultset;
}

bool compfn2(node a,node b)
{
	if(a.score>b.score)
		return true;
	else
		return false;

}

/*
int main(int argc, char** argv)
{
    bool HasLowerbound = false;	
	DomGaph dominantgraph;
	map<int , node> realdata;

	cout<<"Constructing!!!"<<endl;
	ReconstructDG(dominantgraph, realdata, HasLowerbound);	
	cout<<"Construction Complete!!!"<<endl;

//	cout<<"Size of real data"<<realdata.size()<<endl;

	//node missing_node = realdata.find(540848)->second;
	//float missing_score = 0;

    int dimension;
    cout<<"Input dimension:"<<endl;
    scanf("%d", &dimension);
    cout<<"Dimension is:"<<dimension<<endl;
	Weight w;
    //IntVec position;
	for(int i=0; i<dimension; i++)
	{
   		w.weighting.push_back(1.0/dimension);
        //position.push_back(i);
	}
		
    char c;
	//list<FloatVec> resultpool;
    do{
	
	int topk;
	cout<<"Input k value:"<<endl;
	scanf("%d", &topk);
	cout<<"Computing top-"<<topk<<endl;

    clock_t start, stop;
    start = clock();

	TopKDG(realdata,dominantgraph, topk, w);
			
	stop = clock();
	printf("Total time: %lf\n", float(stop- start) / CLOCKS_PER_SEC);
    cout<<"Continue? Y / N"<<endl;
    scanf("%c", &c);
    }while(c!='N');
	return 0;
		
}
*/
