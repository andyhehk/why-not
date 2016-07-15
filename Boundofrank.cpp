#include "Boundofrank.h"


int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		printf("Input only allow one file name!");
		return 0;
	}

	Tuples allTuples;
	Tuple attributes;
	char * filename;
	
	filename = argv[1];
	ReadData(filename, allTuples);
 	int dim = 5;
	IntVec position(5);
	position[0] = 0;
	position[1] = 1;
	position[2] = 2;
	position[3] = 3;
	position[4] = 4;
	//position[5] = 5;
	//position[6] = 6;
	//position[7] = 7;
	//position[8] = 8;

	for(Tuples::iterator iter = allTuples.begin(); iter!=allTuples.end(); iter++)
	{
		int num_dominate = 0;
		Tuple missing_tuple = *iter;

		if(missing_tuple.ID!=2350)
			continue;

		Points incomp_tuples = FindInComp(allTuples, missing_tuple, position, num_dominate);

		cout<<"Computing ID:"<<missing_tuple.ID<<endl;
		cout<<"Num of dominate:"<<num_dominate<<endl;
		cout<<"Num of incomp:"<<incomp_tuples.size()<<endl;
		exit(-1);
		Point missing_point;
		missing_point.ID = missing_tuple.ID;
		for(IntVec::size_type ix = 0; ix != position.size(); ix++)
		{
			missing_point.coordinate.push_back(missing_tuple.data[position[ix]]);
			//cout<<"attribute"<<ix<<":"<<missing_point.coordinate[ix]<<endl;
		}

		
		Bounds bounds = Rankbounds(incomp_tuples, num_dominate, dim);

		//exit(-1);

		iter->data.push_back(bounds.lowerbound);
		iter->data.push_back(bounds.upperbound);
	
	}

	//WriteCSV( filename, allTuples, attributes);
	WriteData( filename, allTuples);
	return 0;

}

Points FindInComp(Tuples allTuples, Tuple missing_tuple, IntVec position, int & num_dominate)
{
	Points incomp_tuples;
	num_dominate = 0;
	
	for(Tuples::iterator iter = allTuples.begin(); iter!=allTuples.end(); iter++)
	{
		Tuple temp = *iter;
		DoubleVec a;
		DoubleVec data;

		bool equal = true;
		bool dominate = true; //dominate the missing tuple
		bool bedominate = true;//be dominated by the missing tuple
		
		for(DoubleVec::size_type ix = 0;ix!=position.size();ix++){
			float result = temp.data[position[ix]] - missing_tuple.data[position[ix]];
			data.push_back(0-result);
			//data.push_back(result);
			if(result!=0)
				equal = false;
			if(result>0)
				bedominate = false;
			else if(result<0)
				dominate = false;
		}

		if(!equal&&dominate)
			num_dominate++;

		else if(!equal&&!dominate&&!bedominate){
			Point temppoint;

			temppoint.ID = temp.ID;
			temppoint.coordinate.assign(data.begin(),data.end());

			incomp_tuples.push_back(temppoint);
		}
		/*
		if(!equal(temp, missing_tuple, position))
		{
			if(Dominate(temp, missing_tuple, position))
			{
				num_dominate++;
			}
			else if(beDominate(temp, missing_tuple, position))
				;
			else
			{
				Point incomp_tuple;
				incomp_tuple.ID = temp.ID;

				for(IntVec::size_type ix = 0; ix!=position.size(); ix++)
					incomp_tuple.coordinate.push_back(temp.data[position[ix]]);
				
				incomp_tuples.push_back(incomp_tuple);
			}
		}
		*/
	}
	return incomp_tuples;
}

//no use now
bool equal(Tuple tuple, Tuple missing_tuple, IntVec position)
{
	for(IntVec::size_type ix = 0; ix!=position.size(); ix++)
	{
		if(tuple.data[position[ix]]!=missing_tuple.data[position[ix]])
			return false;
	}

	return true;
}

//no use now
bool Dominate(Tuple tuple, Tuple missing_tuple, IntVec position)
{
	for(IntVec::size_type ix = 0; ix!=position.size(); ix++)
	{
		if(tuple.data[position[ix]]<missing_tuple.data[position[ix]])
			return false;
	}

	return true;
}

//no use now
bool beDominate(Tuple tuple, Tuple missing_tuple, IntVec position)
{
	for(IntVec::size_type ix = 0; ix!=position.size(); ix++)
	{
		if(tuple.data[position[ix]]>missing_tuple.data[position[ix]])
			return false;
	}
	return true;
}

Bounds Rankbounds(Points incomp_tuples, int num_dominate, int dim)
{
	//Edges edgeGraphL;
	//Edges edgeGraphU;
	Graph g;
	int num_edge = 0;

	int num_computed = 0;
	int this_num = 0;
	Bounds rankBounds;
	int k = 100;

	//srand((int)time(NULL));
	if(num_dominate<k){
		
		map<unsigned int, list<int> > vec;
		int position = 0;
		for(Points::iterator iter = incomp_tuples.begin();iter!=incomp_tuples.end();iter++)
		{
			unsigned int index = 0;
			DoubleVec data = iter->coordinate;
			for(DoubleVec::iterator iter1 = data.begin();iter1!=data.end();iter1++)
			{
				float data = *iter1;
				index = index<<1;
				if(data<=0)
					index++;
			}


			if(vec.find(index)!=vec.end()){
				vec.find(index)->second.push_back(position);
			}
			else{
				list<int> tempvec;
				tempvec.push_back(position);
				vec.insert(pair<unsigned int, list<int> >(index, tempvec));
			}
			position++;
		}

		position = 0;
		for(Points::iterator iter = incomp_tuples.begin();iter!=incomp_tuples.end();iter++)
		{
			Point a = *iter;
			unsigned int original_index = 0;
			unsigned int index = 0;
			int count = 0;
			DoubleVec data = iter->coordinate;
			for(DoubleVec::iterator iter1 = data.begin();iter1!=data.end();iter1++){
				float data = *iter1;
				original_index = original_index<<1;
				index = index<<1;
				if(data>0)
					original_index++;
				if(data>=0)
					index++;
			}
			original_index = pow(2,dim) - 1 - original_index; 
			//cout<<"Before:"<<vec.find(inverse)->second.size()<<endl;
			//cout<<"What is it?"<<*vec.find(inverse)->second.begin()<<endl;
			if(vec.find(original_index)!= vec.end())
				vec.find(original_index)->second.pop_front();
			//cout<<"After:"<<vec.find(inverse)->second.size()<<endl;
			//exit(-1);

			for(map<unsigned int, list<int> >::iterator iter2= vec.begin();iter2!=vec.end();iter2++)
			{
				unsigned int indicate =iter2->first;

				if( (indicate&index)==index)
				{
					//list<int> tempvec = iter2->second;
					
					for(list<int>::iterator iter3=iter2->second.begin();iter3!=iter2->second.end();)
					{						
							num_computed++;
							Point b = incomp_tuples[*iter3];

							float lowerbound = 0;
							float upperbound = INFINITE;
							bool canBothWorse = false;
							for(DoubleVec::size_type ix=0; ix!=b.coordinate.size(); ix++)
							{
								
								float attr1 = a.coordinate[ix];
								float attr2 = b.coordinate[ix];
					
								float templowerb;
								float tempupperb;		

								//if(attr1 == 0&& attr2>0)
								//{
									//this_num++;
									//canBothWorse = true;
									//break;
								//}
								if(attr1 == 0&& attr2==0)
									continue;
								else if(attr1 > 0)
								{
									tempupperb = -attr2/attr1;
									if(tempupperb < upperbound)
										upperbound = tempupperb;
								}
								else
								{
									templowerb = -attr2/attr1;
									if(templowerb > lowerbound)
										lowerbound = templowerb;
								}
								if(lowerbound>upperbound)
								{
									canBothWorse = true;
									break;
								}
									
							}
							if(!canBothWorse)
							{
								//cout<<"Edge (u,v):"<<a.ID<<","<<b.ID<<endl;
								num_edge++;
								if(g.find(a.ID)!=g.end())
								{
									g.find(a.ID)->second.push_back(b.ID);
								}
								else
								{
									list<int> temp;
									temp.push_back(b.ID);
									g.insert(pair<int, list<int> >(a.ID, temp));
								}

								if(g.find(b.ID)!=g.end())
								{
									g.find(b.ID)->second.push_back(a.ID);
								}
								else
								{
									list<int> temp;
									temp.push_back(a.ID);
									g.insert(pair<int, list<int> >(b.ID, temp));
								}
							}

							iter3++;
						
					}
				}
				}
			position++;
		}		                		
		cout<<"Number of computing:"<<num_computed<<endl;
		
	}

	cout<<"number of edge L:"<<num_edge<<endl;
	//cout<<"This num:"<<this_num<<endl;
	//cout<<"number of edge U:"<<edgeGraphU.size()<<endl;

	rankBounds.lowerbound = num_dominate + VertexCover(g)/2+1;
	rankBounds.upperbound = num_dominate + incomp_tuples.size()+1;
	//cout<<"lower bound:"<<rankBounds.lowerbound<<endl;
	//cout<<"upper bound:"<<rankBounds.upperbound<<endl;
	return rankBounds;
}


//no use now
bool CanBothBetter(DoubleVec a, DoubleVec b ) //check if two tuples can both rank better than the missing tuple at the same time
{
	float lowerbound = 0;
	float upperbound = INFINITE;
	
	if(a.size() == b.size())
	{
		float templowerb;
		float tempupperb;
		
		for(DoubleVec::size_type ix = 0; ix!=a.size(); ix++)
		{
			if(a[ix] == 0)
				continue;
			else if(a[ix] > 0)
			{
				tempupperb = -b[ix]/a[ix];

				if(tempupperb < upperbound)
					upperbound = tempupperb;
			}
			else
			{
				templowerb = -b[ix]/a[ix];

				if(templowerb > lowerbound)
					lowerbound = templowerb;
			}

			if(lowerbound>upperbound)
				return true;
		}
	}
	return false;
}

//no use now
bool CanBothWorse(DoubleVec a,DoubleVec b) //check if two tuples can both rank worse than the missing tuples at the same time
{
	if(a.size() == b.size())
	{
		for(DoubleVec::size_type ix =0; ix!=a.size(); ix++)
		{
			a[ix] = 0-a[ix];
			b[ix] = 0-b[ix];
		}
		return CanBothBetter(a,b);
	}
	
	return false;
}

int VertexCover(Graph g)
{
	int rank = 0;
	/*
	while(!edgeGraph.empty())
	{

		Edges::iterator iter = edgeGraph.begin();

		Edge tempEdge = *iter;

		rank +=2;
		int node1 = tempEdge.node1;
		int node2 = tempEdge.node2;
		
		while(iter!=edgeGraph.end())
		{			
			if(iter->node1 == node1 || iter->node1 == node2 || iter->node2 == node1 || iter->node2 == node2)
			{
				iter = edgeGraph.erase(iter);
			}

			else
			{
				iter++;
			}
		}
	}
*/
	while(!g.empty())
	{
		rank += 2;
		
		Graph::iterator iteru = g.begin();
		int u = iteru->first;
		list<int> edgeu = iteru->second;
		int v = *(edgeu.begin());	
		g.erase(iteru);				

		list<int>::iterator iter1 = edgeu.begin();
		
		for(iter1++;iter1!=edgeu.end();iter1++)
		{
			Graph::iterator iter = g.find(*iter1);
			for(list<int>::iterator iter2 = iter->second.begin();iter2!=iter->second.end();iter2++)
			{
				if(*iter2==u)
				{
					iter->second.erase(iter2);
					break;
				}
			}
			if(iter->second.empty())
				g.erase(iter);
		}

		Graph::iterator iterv =  g.find(v);
		list<int> edgev = iterv->second;
		g.erase(iterv);

		for(list<int>::iterator iter1 = edgev.begin();iter1!=edgev.end();iter1++)
		{
			Graph::iterator iter = g.find(*iter1);
			if(iter!=g.end())
			{
				for(list<int>::iterator iter2 = iter->second.begin();iter2!=iter->second.end();iter2++)
				{
					if(*iter2==v)
					{
						iter->second.erase(iter2);
						break;
					}
				}
				if(iter->second.empty())
					g.erase(iter);
			}
		}

	}

	return rank;
}

