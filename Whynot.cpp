#include "Whynot.h"
//#include <octave/oct.h>
//#include <octave/octave.h>
//#include <octave/parse.h>
//#include <octave/toplev.h>
//extern "C"{
//    #include "libqp.h"
//}
bool testall = true;   

typedef struct GoodCase{
	int ID;
	float winpercentage;
	
	friend bool operator < (GoodCase a, GoodCase b)
	{
		return a.winpercentage > b.winpercentage;
	}
}GoodCase;


//double *mat_H = NULL;
//uint32_t nVars;

//const double *get_col_of_mat_H( uint32_t i )
//{
//    return( &mat_H[ nVars*i ] );
//}

//static boost::mt19937 generator;
//static boost::uniform_real<> uni_dist(0,1);
//static boost::variate_generator<boost::mt19937, boost::uniform_real<> > uni(generator, uni_dist);

bool CompareSkyAns(Answer a, Answer b)
{
	return a.w.deltaw>b.w.deltaw;
}

bool Gauss(vector<DoubleVec>& matrix, DoubleVec &ans, int row, int column, double & sum2);
bool GaussianJordan(vector< DoubleVec >& matrix, DoubleVec &ans, int row, int column, double & sum2)
{
	//cout<<"inside1"<<endl;
	int i =0, j=0;

	while(i<row&&j<column)
	{
		int maxi = i;
		for(int k = i+1; k<row; k++)
		{
			if(fabs(matrix[k][j]) > fabs(matrix[maxi][j]))
				maxi = k;
		}

		if(matrix[maxi][j]!=0)
		{
			double divider = matrix[maxi][j];
			
			if(i!=maxi)
			{
				swap(matrix[i],matrix[maxi]);
			}

			for(int v = j; v<column;v++)
			{
				matrix[i][v] /= divider;
			}
		
			for(int u = i+1;u<row;u++)
			{
				for(int v = column-1; v>=0;v--)
				{
					matrix[u][v] = matrix[u][v] - matrix[u][j]*matrix[i][v];
				}
			}
			//for(int ii=0;ii<row;ii++)
				//for(int jj=0;jj<column;jj++)
					//cout<<"Matrix["<<ii<<"]"<<"["<<jj<<"]:"<<matrix[ii][jj]<<endl;
				
			i++;
		}

		else
			return false;
		j++;
	}
    
    for(int ii=0;ii<row;ii++)
    {
        for(int kk = ii+1;kk<row;kk++)
        {
            double temp = matrix[ii][kk];
            if(temp==0)
                continue;
            for(int jj=kk;jj<column;jj++)
            {
                matrix[ii][jj] -= matrix[kk][jj]*temp;
            }
        }
    }
/*
	for(int ii=0;ii<row;ii++)
    {
		for(int jj=0;jj<column;jj++)
			cout<<matrix[ii][jj]<<"\t";
        cout<<endl;
    }
*/
    if(row<column-1)
        return true; 
    
    else if(row == column -2)
    {
        if(matrix[row][column]!=0)
            return false;
    }
/*
	for(i = row-1;i>=0;i--)
	{
		
		float sum = matrix[i][column-1];

		int k = 0;
		
		for(j = column-2;j>i;j--)
		{
			//cout<<"value of k:"<<k<<endl;
			//cout<<"Size of answer:"<<ans.size()<<endl;
			sum -= matrix[i][j]*ans[k];
			k++;
		}

		float value = sum/matrix[i][i];
		if(value<=0)
			return false;
        sum2 += value;
		ans.push_back(value);
	}
*/
	
	//cout<<"inside2"<<endl;
	return true;
}


//void SampleWeightsFromP(list<Weight>& sampleweight, double quality_of_answer, double probability_guarantee, Weight w_origin, const vector< FloatVec >& incompar_points, long sizeneed = 0)
//{
//    long samplesize;
//    list<Weight> templist; 
//    int count = 0;
//    int count2 = 0;
//
//	int dimension = w_origin.weighting.size();
//	int size = incompar_points.size();
//
//	samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer))+1;
//    //samplesize = sizeneed;
//
//    if(sizeneed>0)
//        samplesize = sizeneed;
//
//    vector<Var> sortedplane;
//    sortedplane.reserve(size);
//
//    for(int i = 0;i<size;i++)
//    {
//        break;
//        FloatVec normvector = incompar_points[i];
//
//        float distance = 0;
//        float scale = 0;
//        for(int j = 0;j<dimension;j++)
//        {
//             distance += w_origin.weighting[j] * normvector[j];
//             scale += normvector[j]*normvector[j];
//        }
//        distance = fabs(distance)/sqrt(scale);
//        Var temp;
//        temp.ID = i;
//        temp.data = distance;
//        sortedplane.push_back(temp);
//    }
//
//    //cout<<"Size:"<<size<<endl;
//    //cout<<"Segment:"<<segment<<endl;
//    //cout<<"residual:"<<residual<<endl;
//
//    sort(sortedplane.begin(), sortedplane.end()); 
//
//    int segment = size/samplesize;
//    //segment = 1;
//    float min = sortedplane[0].data;
//    float max = sortedplane[size-1].data;
//    float inteval = (max-min)/segment;
//
//    vector<interv> terminatepoint;
//    terminatepoint.reserve(segment+1);
//    int jj = 0;
//    for(int i = 0;i<size;)
//    {
//        break;
//        int s = 0;
//        //cout<<"i"<<i<<endl;
//        //cout<<"Max:"<<max<<endl;
//        float terminate = min+(jj+1)*inteval;
//        interv tempinterv;
//        tempinterv.begin = i;
//        //cout<<"Terminate:"<<terminate<<endl;
//        while(sortedplane[i].data<=terminate)
//        {
//            s++;
//            i++;
//            if(!(i<size))
//                break;
//        }
//        //cout<<"Total size:"<<size<<endl;
//        //cout<<"Size:"<<s<<endl;
//        //cout<<"I:"<<i<<endl;
//        //cout<<"J:"<<jj<<endl;
//        tempinterv.end = i-1;
//        //terminatepoint[jj].end = i-1;
//        if((tempinterv.end-tempinterv.begin)!=-1)
//            terminatepoint.push_back(tempinterv);
//        jj++;
//    }
//    set<int> used;
//
//    //for(int i = 0;i<terminatepoint.size();i++)
//    //{
//        //cout<<"Begin:"<<terminatepoint[i].begin<<"\t"<<"End:"<<terminatepoint[i].end<<endl;
//    //}
//    vector< vector<int> > indexvec;
//    indexvec.resize(terminatepoint.size());
//    //indexvec.reserve(segment+1);
//    
//
//    vector< map<int, vector<int> > > mapdirecvec;
//    mapdirecvec.resize(terminatepoint.size());
//    //mapdirecvec.reserve(segment+1);
//
//    for(int i = 0;i<terminatepoint.size();i++)
//    {
//        break;
//        for(int j=terminatepoint[i].begin;j<=terminatepoint[i].end;j++)
//        {
//            int index = 0;
//
//            FloatVec temp = incompar_points[sortedplane[j].ID];
//            for(int k = 0;k<dimension;k++)
//            {
//                if(temp[k]>0)
//                {
//                    index = (index+1)<<1;
//                }
//                else
//                {
//                    index = index<<1;
//                }
//            }
//
//            if(mapdirecvec[i].find(index)==mapdirecvec[i].end())
//            {
//                indexvec[i].push_back(index);
//                vector<int> tempvec;
//                tempvec.push_back(j);
//                mapdirecvec[i].insert(pair<int, vector<int> >(index, tempvec));
//            }
//            else
//            {
//                mapdirecvec[i].find(index)->second.push_back(j);
//            }
//        }
//    }
//
//    //cout<<"Come here?"<<endl;
//    IntVec candidates;
//    candidates.reserve(dimension);
//    vector<DoubleVec> matrix;
//    matrix.reserve(dimension);
//
//    boost::uniform_int<> random_hplane(1, dimension-1);
//    boost::uniform_int<> random_candidate(0, size-1);
//    boost::variate_generator<boost::mt19937, boost::uniform_int<> > die1(generator, random_hplane);
//    boost::variate_generator<boost::mt19937, boost::uniform_int<> > die2(generator, random_candidate); 
//
//	while(samplesize>0)
//	{
//        candidates.clear();
//		matrix.clear();
//
//		//int num_hplane = rand()%(dimension-1)+1;
//        //int num_hplane = die1();
//        int num_hplane = 1;
//        int candidate_num = num_hplane;
//        //if(size>3000)
//        if(false)
//        {
//            while(candidate_num>0)
//		    {
//                int value = rand()%(terminatepoint.size());
//                {
//                    if((terminatepoint[value].end - terminatepoint[value].begin) ==-1)
//                        continue;
//                    int value1 = rand()%(indexvec[value].size());
//                    int index = indexvec[value][value1];
//                    vector<int> candidatewithdirec = mapdirecvec[value].find(index)->second;
//
//                    int value2 = rand()%(candidatewithdirec.size());
//                    int candidate = sortedplane[candidatewithdirec[value2]].ID;
//                    //cout<<"This candidate:"<<candidate<<endl;
//                    //int value2 = rand()%(terminatepoint[value].end-terminatepoint[value].begin+1);
//                    //int candidate = sortedplane[terminatepoint[value].begin+value2].ID;
//
//                    //if(used.find(candidate)==used.end())
//                    {
//                        //used.insert(candidate);
//                        candidates.push_back(candidate);
//                        candidate_num--;
//                    }
//                }
//            //cout<<"candidate: "<<i+1<<candidates[i]<<"\t";
//		    }
//        }
//
//        else
//        {
//            
//	        for(int i = 0;i<num_hplane;i++)
//		    {
//			    candidates.push_back(rand()%size);
//			    //candidates.push_back(die2());
//                //cout<<"candidate "<<i+1<<":"<<candidates[i]<<"\t";
//		    }
//        }
//        //cout<<endl;
//
//        vector<FloatVec> candidatevec;
//        candidatevec.reserve(num_hplane);
//
//        for(IntVec::iterator iter = candidates.begin();iter!=candidates.end();iter++)
//        {
//            candidatevec.push_back(incompar_points[*iter]);
//        }
//
//        for(int i = 0;i<num_hplane;i++)
//        {
//            DoubleVec rowvec;
//            rowvec.reserve(dimension+1);
//            for(int j = 0; j<dimension;j++)
//                rowvec.push_back(candidatevec[i][j]);
//            rowvec.push_back(0);
//            matrix.push_back(rowvec);
//        }
//
//
//        DoubleVec rowvec;
//        rowvec.assign(dimension+1, 1);
//        matrix.push_back(rowvec);
//
//        DoubleVec weight_i;
//        weight_i.reserve(dimension);
//
//        DoubleVec ans;
//        int row = num_hplane + 1;
//        int column = dimension + 1;
//        double sum2;
//        GaussianJordan(matrix, ans, row, column, sum2);
//
//        vector<DoubleVec> coeficient(column - row -1);
//
//        for(int i = 0;i<coeficient.size();i++)
//        {
//            for(int j =  0; j< row; j++)
//            {
//                coeficient[i].push_back(matrix[j][i+row]); 
//            }
//
//        }
//        //cout<<"---------------------------"<<endl;
//        /*
//        for(int i = 0;i<coeficient.size();i++)
//        {
//            DoubleVec temp = coeficient[i];
//            for(int j = 0;j<temp.size();j++)
//            {
//                cout<<temp[j]<<" ";
//            }
//            cout<<endl;
//        }
//        */
//        for(int i = 0;i < coeficient.size();i++)
//        {
//            rowvec.clear();
//            for(int k = 0; k<row;k++)
//                rowvec.push_back(coeficient[i][k]);
//            
//            for(int j = 0;j<dimension - row;j++)
//                rowvec.push_back(0);
//
//            if(i!=coeficient.size()-1)
//                rowvec[row+i] = 1;
//
//            double sum = 0;
//
//            for(int k = 0;k<rowvec.size();k++)
//            {
//                sum += w_origin.weighting[k]*rowvec[k];
//            }
//
//            rowvec.push_back(sum);
//            matrix.push_back(rowvec);
//
//        }
//        
//        if(Gauss(matrix, ans, dimension, column, sum2))
//        {
//            //break;
//            double sum = 0;
//
//            Weight weights;
//            weights.weighting.resize(dimension);
//
//            for(int i = 0; i<dimension;i++)
//            {
//                weights.weighting[i] = ans[i];
//            }
//
//		    double deltaW = 0;
//		
//		    for(FloatVec::size_type ix =0; ix!=weights.weighting.size(); ix++)
//		    {
//			    deltaW +=  (weights.weighting[ix]-w_origin.weighting[ix])*(weights.weighting[ix]-w_origin.weighting[ix]);
//		    }
//
//            //cout<<deltaW<<endl;
//		    weights.deltaw= sqrt(deltaW);
//		    sampleweight.push_back(weights);
//            //templist.push_back(weights);
//
//            count++;
//
//            samplesize--;
//            continue;
//        }
//        /////////////////////////////////////////////////////////////////////////////
//        ///////////////////////////
//        /////////////////
//        	
//        {   
//            IntVec defaultVec;
//            defaultVec.reserve(dimension);
//            for(int i=0;i<dimension;i++)
//                defaultVec.push_back(i);
//    
//            IntVec varIndicate, knowIndicate;
//            varIndicate.reserve(dimension);
//            knowIndicate.reserve(dimension);
//
//            varIndicate.assign(defaultVec.begin(),defaultVec.end());
//				
//            Weight weights;
//            weights.weighting.resize(dimension);
//
//            for(int i=0;i<dimension-num_hplane;i++)
//            {
//                int value = rand()%(dimension-i);
//                //int value = die1();
//                knowIndicate.push_back(varIndicate[value]);
//                varIndicate.erase(varIndicate.begin()+value);
//            }
//
//            sort(knowIndicate.begin(),knowIndicate.end());
//            sort(varIndicate.begin(),varIndicate.end());
//		    double sum_weight = 0;
//
//            for(IntVec::iterator iter = knowIndicate.begin();iter!=knowIndicate.end();iter++)
//		    {
//			    float value = (float)(rand())/RAND_MAX;
//			    //float value = uni();
//			    sum_weight += value;
//			    //ans1.push_back(value);
//                weights.weighting[*iter]=value;
//		    }
//	
//            vector<DoubleVec> matrix2;
//            matrix2.reserve(dimension);
//            DoubleVec rowvec;
//            rowvec.reserve(dimension+1);
//            for(int i = 0;i<num_hplane;i++)
//		    {
//    
//                rowvec.clear();
//			    float sum = 0;
//
//			    int indicate = 0;
//
//                for(IntVec::iterator iter = knowIndicate.begin();iter!=knowIndicate.end();iter++)
//			    {
//				    //sum += candidate_vec[indicate]*ans1[indicate];
//                    sum += incompar_points[candidates[i]][*iter] * weights.weighting[*iter];
//				    //indicate++;
//			    }
//                for(IntVec::iterator iter = varIndicate.begin();iter!=varIndicate.end();iter++)
//			    {
//				    //row.push_back(candidate_vec[indicate]);
//                    rowvec.push_back(incompar_points[candidates[i]][*iter]);
//				    //indicate++;
//			    }
//			    rowvec.push_back(-sum);
//			//row.push_back(0);
//
//			    matrix2.push_back(rowvec);
//			
//		    }
//
//
//            double sum_ans2 = 0;
//            int row = num_hplane;
//            int column = row+1;
//            DoubleVec ans2;
//
//     
//            int attempts = 0;
//            bool getonbound = false;
//            //while(attempts<10)
//            {
//                attempts++;
//                if(Gauss(matrix2, ans2, row, column, sum_ans2))
//                {
//                    sum_weight += sum_ans2;
//                    int ans2index = 0;
//			        for(DoubleVec::iterator iter = ans2.begin();iter!=ans2.end();iter++)
//			        {
//				        weights.weighting[varIndicate[ans2index]]=*iter/sum_weight;
//                        ans2index++;
//			        }
//
//
//                    for(IntVec::iterator iter = knowIndicate.begin();iter!=knowIndicate.end();iter++)
//                    {
//                        weights.weighting[*iter] = weights.weighting[*iter]/sum_weight;
//                    }
//            
//                    /* 
//                    for(int i = 0;i<candidates.size();i++)
//                    {
//			            float test_sum = 0;
//			            for(FloatVec::size_type ix = 0;ix!=incompar_points[candidates[i]].size();ix++)
//			            {
//				            test_sum +=incompar_points[candidates[i]][ix]*weights.weighting[ix];
//			            }
//
//			            cout<<"Look at this sum, is it zero? "<<test_sum<<endl;
//                    }
//                    */
//                    double deltaW = 0;
//		
//		            for(FloatVec::size_type ix =0; ix!=weights.weighting.size(); ix++)
//		            {
//			            deltaW +=  (weights.weighting[ix]-w_origin.weighting[ix])*(weights.weighting[ix]-w_origin.weighting[ix]);
//		            }
//
//		            weights.deltaw= sqrt(deltaW);
//		            sampleweight.push_back(weights);
//                    //templist.push_back(weights);
//
//                    count2++;
//                    samplesize--;
//                    continue;
//                    getonbound = true;
//                    //break;
//                }
//            }
//        }
//        /////////////////////////////////////////////////////////////////////////////////
//        continue;
//
//        {
//            Weight w;
//
//            w.weighting.reserve(dimension);
//		    float sum = 0;
//		    for(int j = 0; j < dimension; j++)
//		    {
//			    float temp = (float)(rand())/RAND_MAX;
//			    w.weighting.push_back(temp);
//			    sum += temp;
//		    }
//		
//		    float deltaW = 0;
//		
//		    for(int ix =0; ix<dimension; ix++)
//		    {
//			    w.weighting[ix]= w.weighting[ix]/sum;
//
//			    deltaW +=  (w.weighting[ix]-w_origin.weighting[ix])*(w.weighting[ix]-w_origin.weighting[ix]);
//		    }
//
//		    w.deltaw= sqrt(deltaW);
//		    //sampleweight.push_back(w);
//            templist.push_back(w);
//            samplesize--;
//        }
//    }
//   
//    //cout<<"The total attempt:"<<count<<endl;
//    //cout<<"The total attempt2:"<<count2<<endl;
//	sampleweight.sort(compWeight);
//    return; 
//    
//    list<Weight>::iterator iter1 = templist.begin();
//    sampleweight.push_back(*iter1);
//    Weight tempweight = *iter1;
//    
//    templist.erase(iter1);
//    while(!templist.empty())
//    {
//        //cout<<"come here"<<endl;
//        float min = INFINITE;
//        Weight candidate;
//        list<Weight>::iterator itertodel;
//        for(list<Weight>::iterator iter = templist.begin();iter!=templist.end();iter++)
//        {
//           float deltaw = 0;
//           
//           for(int i = 0;i<dimension;i++)
//           {
//                deltaw = (tempweight.weighting[i]-iter->weighting[i])*(tempweight.weighting[i]-iter->weighting[i]);
//           }
//
//           if(min > deltaw)
//           {
//                candidate = *iter;
//                min = deltaw;
//                itertodel = iter;
//           }
//        }
//        tempweight = candidate;
//        sampleweight.push_back(candidate);
//        templist.erase(itertodel);
//
//    }
//   
//}
//

list<Weight> SampleWeights(float quality_of_answer, float probability_guarantee, Weight w_origin)
{
	long samplesize;
	int dimension = w_origin.weighting.size();
	//srand ( time(NULL) );

	samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer/100))+1;
	list<Weight> sampleweight;
    list<Weight> templist;
	
	for(int i = 0; i < samplesize; i++)
	{
		Weight w;
        //w.weighting.reserve(dimension);
		float sum = 0;
		for(int j = 0; j < dimension; j++)
		{
			//int temp = rand()%5;
			float temp = (float)(rand())/RAND_MAX;
			//float temp = uni();
            //cout<<"temp "<<j+1<<":"<<temp<<"\t";
			w.weighting.push_back(temp);
			sum += temp;
		}
		//cout<<endl;

		float deltaW = 0;
		
		for(int ix =0; ix<dimension; ix++)
		{
			w.weighting[ix]= w.weighting[ix]/sum;

			deltaW +=  (w.weighting[ix]-w_origin.weighting[ix])*(w.weighting[ix]-w_origin.weighting[ix]);
		}

		w.deltaw= sqrt(deltaW);
		sampleweight.push_back(w);
        templist.push_back(w);
	}

	sampleweight.sort(compWeight);
    return sampleweight;
	list<Weight>::iterator iter1 = templist.begin();
    sampleweight.push_back(*iter1);
    Weight tempweight = *iter1;
    
    templist.erase(iter1);
    while(!templist.empty())
    {
        //cout<<"come here"<<endl;
        float min = INFINITE;
        Weight candidate;
        list<Weight>::iterator itertodel;
        for(list<Weight>::iterator iter = templist.begin();iter!=templist.end();iter++)
        {
           float deltaw = 0;
           
           for(int i = 0;i<dimension;i++)
           {
                deltaw = (tempweight.weighting[i]-iter->weighting[i])*(tempweight.weighting[i]-iter->weighting[i]);
           }

           if(min > deltaw)
           {
                candidate = *iter;
                min = deltaw;
                itertodel = iter;
           }
        }
        tempweight = candidate;
        sampleweight.push_back(candidate);
        templist.erase(itertodel);

    }

	//list<Weight>::iterator iter = sampleweight.begin();
	//cout<<"The first weight:"<<iter->deltaw<<endl;
	
	return sampleweight;
}

//bool Gauss(vector<DoubleVec>& matrix, DoubleVec &ans, int row, int column, double & sum2)
//{
//    double tol = 0.000001;
//
//    DoubleVec s(row);
//    ans.resize(row);
//
//    for(int i = 0;i<row;i++)
//    {
//        s[i] = fabs(matrix[i][0]);
//        for(int j = 1;j<row;j++)
//        {
//            if(fabs(matrix[i][j])>s[i])
//                s[i] = fabs(matrix[i][j]);
//        }
//    }
//
//    for(int k = 0; k<row-1;k++)
//    {
//        int p = k;
//        double big = fabs(matrix[k][k]/s[k]);
//        for(int ii = k+1; ii< row;ii++)
//        {
//            double dummy = fabs(matrix[ii][k]/s[ii]);
//            if(dummy > big)
//            {
//                big = dummy;
//                p = ii;
//            }
//        }
//        if(p!=k)
//        {
//            for(int jj = 0;jj<row;jj++)
//            {
//                double dummy = matrix[p][jj];
//                matrix[p][jj] = matrix[k][jj];
//                matrix[k][jj] = dummy;
//            }
//            double dummy = matrix[p][column-1];
//            matrix[p][column-1] = matrix[k][column-1];
//            matrix[k][column-1] = dummy;
//
//            dummy = s[p];
//            s[p] = s[k];
//            s[k] = dummy;
//        }
//
//        if(fabs(matrix[k][k]/s[k])<tol)
//        {
//            return false;
//        }
//
//        for(int i = k+1; i<row;i++)
//        {
//            double factor = matrix[i][k]/matrix[k][k];
//             
//            for(int j = k+1;j<row;j++)
//            {
//                matrix[i][j] -= factor*matrix[k][j];
//            }
//
//            matrix[i][column-1] -= factor*matrix[k][column-1];
//        }
//
//        if(fabs(matrix[row-1][row-1]/s[row-1])<tol)
//            return false;
//
//    }
//    
//    
//    if(matrix[row-1][row-1]!=0)    
//        sum2 = ans[row-1] = matrix[row-1][column-1]/matrix[row-1][row-1];
//    else
//        return false;
//    if(ans[row-1]<0)
//        return false;
//    for(int i = row-2;i>=0;i--)
//    {
//        double sum = 0;
//        for(int j = i+1; j<row;j++)
//        {
//            sum += matrix[i][j]*ans[j];
//        }
//        ans[i] = (matrix[i][column-1]-sum)/matrix[i][i];
//        sum2 += ans[i];
//        if(ans[i]<0)
//            return false;
//    }
//
//    /*
//    for(int i = 0;i<row;i++)
//    {
//        for(int j = 0;j<column;j++)
//        {
//            cout<<matrix[i][j]<<"\t";
//            
//        }
//        cout<<endl;
//    }
//    */
//
//    return true;
//}



//void SampleWeightsFromB(list<Weight>& sampleweight, double quality_of_answer, double probability_guarantee, Weight w_origin, const vector< FloatVec >& incompar_points, long sizeneed)
//{
//
//	long samplesize;
//    	int count =0;
//	int dimension = w_origin.weighting.size();
//	int size = incompar_points.size();
//
//
//    list<Weight> templist;
//	samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer))+1;
//    //samplesize = sizeneed;
//
//    if(sizeneed>0)
//        samplesize = sizeneed;
//
//    vector<Var> sortedplane;
//    sortedplane.reserve(size);
//
//    for(int i = 0;i<size;i++)
//    {
//        break;
//        FloatVec normvector = incompar_points[i];
//
//        float distance = 0;
//        float scale = 0;
//        for(int j = 0;j<dimension;j++)
//        {
//             distance += w_origin.weighting[j] * normvector[j];
//             scale += normvector[j]*normvector[j];
//        }
//        distance = fabs(distance)/sqrt(scale);
//        Var temp;
//        temp.ID = i;
//        temp.data = distance;
//        sortedplane.push_back(temp);
//    }
//
//    //cout<<"Size:"<<size<<endl;
//    int segment = size/samplesize;
//    //cout<<"Segment:"<<segment<<endl;
//    int residual = size - segment*samplesize;
//    //cout<<"residual:"<<residual<<endl;
//    int samplesizebackup = samplesize;
//
//    sort(sortedplane.begin(), sortedplane.end()); 
//
//    float min = sortedplane[0].data;
//    float max = sortedplane[size-1].data;
//    float inteval = (max-min)/samplesize;
//
//    vector<interv> terminatepoint;
//    terminatepoint.resize(samplesize+1);
//    int jj = 0;
//    for(int i = 0;i<size;)
//    {
//        break;
//        int s = 0;
//        //cout<<"Max:"<<max<<endl;
//        float terminate = min+(jj+1)*inteval;
//        //cout<<"Terminate:"<<terminate<<endl;
//        terminatepoint[jj].begin = i;
//        while(sortedplane[i].data<=terminate)
//        {
//            s++;
//            i++;
//            if(!(i<size))
//                break;
//        }
//        //cout<<"Total size:"<<size<<endl;
//        //cout<<"Size:"<<s<<endl;
//        //cout<<"I:"<<i<<endl;
//        //cout<<"J:"<<jj<<endl;
//        terminatepoint[jj].end = i-1;
//        jj++;
//    }
//
//    DoubleVec row;
//
//    row.reserve(dimension+1);
//
//    IntVec candidates;
//    candidates.reserve(dimension);
//	DoubleVec ans1;
//	DoubleVec ans2;
//    ans1.reserve(dimension);
//    ans2.reserve(dimension);
//
//    boost::uniform_int<> random_hplane(1, dimension-1);
//    boost::uniform_int<> random_candidate(0, size-1);
//    boost::variate_generator<boost::mt19937, boost::uniform_int<> > die1(generator, random_hplane);
//    boost::variate_generator<boost::mt19937, boost::uniform_int<> > die2(generator, random_candidate); 
//
//	vector< DoubleVec > matrix;
//    matrix.reserve(dimension);
//	
//    IntVec defaultVec;
//    defaultVec.reserve(dimension);
//    for(int i=0;i<dimension;i++)
//        defaultVec.push_back(i);
//    
//    IntVec varIndicate, knowIndicate;
//    varIndicate.reserve(dimension);
//    knowIndicate.reserve(dimension);
//
//    set<int> used;
//	while(samplesize>0)
//	{
//        candidates.clear();
//        ans1.clear();
//        ans2.clear();
//        matrix.clear();
//        knowIndicate.clear();
//        varIndicate.assign(defaultVec.begin(),defaultVec.end());
//				
//        Weight weights;
//        weights.weighting.resize(dimension);
//		//int num_equa = rand()%(dimension-1)+1;
//        //int num_equa = die1();
//        int num_equa = 1;
//		//int num_equa = rand()%(dimension-3)+3;
//        //int num_equa = dimension-1;
//
//        //cout<<"number of equation:"<<num_equa<<endl;
//        int candidate_num = num_equa;
//
//        if(false)
//        while(candidate_num>0)
//		{
//            int value = rand()%(terminatepoint.size());
//            {
//                if((terminatepoint[value].end - terminatepoint[value].begin) ==-1)
//                    continue;
//                int value2 = rand()%(terminatepoint[value].end-terminatepoint[value].begin+1);
//                int candidate = sortedplane[terminatepoint[value].begin+value2].ID;
//
//                if(used.find(candidate)==used.end())
//                {
//                    used.insert(candidate);
//                    candidates.push_back(candidate);
//                    candidate_num--;
//                }
//            }
//            //cout<<"candidate: "<<i+1<<candidates[i]<<"\t";
//		}
//
//     	else
//        {
//            
//	        for(int i = 0;i<num_equa;i++)
//		    {
//			    candidates.push_back(rand()%size);
//			    //candidates.push_back(die2());
//                //cout<<"candidate "<<i+1<<":"<<candidates[i]<<"\t";
//		    }
//        }
//
//        for(int i=0;i<dimension-num_equa;i++)
//        {
//            int value = rand()%(dimension-i);
//            knowIndicate.push_back(varIndicate[value]);
//            varIndicate.erase(varIndicate.begin()+value);
//        }
//
//        sort(knowIndicate.begin(),knowIndicate.end());
//        sort(varIndicate.begin(),varIndicate.end());
//        //for(IntVec::iterator iter=knowIndicate.begin();iter!=knowIndicate.end();iter++)
//          //cout<<"what is know:"<<*iter<<endl;
//        
//        //for(IntVec::iterator iter=varIndicate.begin();iter!=varIndicate.end();iter++)
//          //cout<<"what is var:"<<*iter<<endl;
//        //exit(-1);
//		double sum_weight = 0;
//
//        for(IntVec::iterator iter = knowIndicate.begin();iter!=knowIndicate.end();iter++)
//		{
//			//float value = (float)(rand())/RAND_MAX;
//            float value = uni();
//            //cout<<"Know value:"<<value<<"\t";
//			sum_weight += value;
//			//ans1.push_back(value);
//            weights.weighting[*iter]=value;
//		}
//        //cout<<endl;
//		//cout<<"Num of varable now:"<<num_var<<endl;
//		for(int i = 0;i<num_equa;i++)
//		{
//			//FloatVec row;
//            //row.reserve(num_equa+1);
//			//FloatVec candidate_vec = incompar_points[candidates[i]];
//            row.clear();
//			float sum = 0;
//
//			int indicate = 0;
//
//            for(IntVec::iterator iter = knowIndicate.begin();iter!=knowIndicate.end();iter++)
//			{
//				//sum += candidate_vec[indicate]*ans1[indicate];
//                sum += incompar_points[candidates[i]][*iter] * weights.weighting[*iter];
//				//indicate++;
//			}
//            //cout<<"what about here????"<<endl;
//			//for(int j=0;j<num_equa;j++)
//			//for(int j=0;j<dimension;j++)
//            for(IntVec::iterator iter = varIndicate.begin();iter!=varIndicate.end();iter++)
//			{
//				//row.push_back(candidate_vec[indicate]);
//                row.push_back(incompar_points[candidates[i]][*iter]);
//				//indicate++;
//            }
//
//
//			row.push_back(-sum);
//			//row.push_back(0);
//
//			//cout<<endl;
//            //for(int i = 0;i<row.size();i++)
//                    //cout<<row[i]<<"\t";
//
//			//exit(-1);
//
//            //cout<<endl;
//			matrix.push_back(row);
//			
//		}
//
//        //row.assign(dimension+1, 1);
//        //matrix.push_back(row);
//			//FloatVec row;
//			//for(int j=0;j<dimension;j++)
//			//{
//			//	row.push_back(1);
//			//}
//			//row.push_back(1);
//			//matrix.push_back(row);
//        double sum_ans2 = 0;
//        int row = num_equa;
//        int column = row+1;
//
//        if(Gauss(matrix, ans2, row, column, sum_ans2))
//        {
//            sum_weight += sum_ans2;
//            int ans2index = 0;
//			for(DoubleVec::iterator iter = ans2.begin();iter!=ans2.end();iter++)
//			{
//				weights.weighting[varIndicate[ans2index]]=*iter/sum_weight;
//                ans2index++;
//			}
//
//
//            for(IntVec::iterator iter = knowIndicate.begin();iter!=knowIndicate.end();iter++)
//            {
//                weights.weighting[*iter] = weights.weighting[*iter]/sum_weight;
//            }
//            
//            /* 
//            for(int i = 0;i<candidates.size();i++)
//            {
//			    float test_sum = 0;
//			    for(FloatVec::size_type ix = 0;ix!=incompar_points[candidates[i]].size();ix++)
//			    {
//				    test_sum +=incompar_points[candidates[i]][ix]*weights.weighting[ix];
//			    }
//
//			    cout<<"Look at this sum, is it zero? "<<test_sum<<endl;
//            }
//            */
//			float deltaW = 0;
//			for(int ix =0; ix<dimension; ix++)
//			{
//				deltaW +=  (weights.weighting[ix]-w_origin.weighting[ix])*(weights.weighting[ix]-w_origin.weighting[ix]);
//			}
//			weights.deltaw= sqrt(deltaW);
//			//templist.push_back(weights);
//            sampleweight.push_back(weights);
//            samplesize--;
//            count++;
//        }
//        
//        else
//        {
//            continue;
//            Weight w;
//
//            w.weighting.reserve(dimension);
//		    float sum = 0;
//		    for(int j = 0; j < dimension; j++)
//		    {
//			//int temp = rand()%5;
//			float temp = (float)(rand())/RAND_MAX;
//			w.weighting.push_back(temp);
//			sum += temp;
//		    }
//		
//		    float deltaW = 0;
//		
//		    for(int ix =0; ix<dimension; ix++)
//		    {
//			    w.weighting[ix]= w.weighting[ix]/sum;
//
//			    deltaW +=  (w.weighting[ix]-w_origin.weighting[ix])*(w.weighting[ix]-w_origin.weighting[ix]);
//		    }
//
//		    w.deltaw= sqrt(deltaW);
//		    sampleweight.push_back(w);
//            //templist.push_back(w);
//            samplesize--;
//        }
//    
//   
// 
//          
//	}//loop for samplesize
//
//    cout<<"Total attempt:"<<count<<endl;
//	//cout<<"Total time:"<<sumtime / CLOCKS_PER_SEC <<endl;
//	sampleweight.sort(compWeight);
//
//    return;
//    
//	list<Weight>::iterator iter1 = templist.begin();
//    sampleweight.push_back(*iter1);
//    Weight tempweight = *iter1;
//    
//    templist.erase(iter1);
//    while(!templist.empty())
//    {
//        //cout<<"come here"<<endl;
//        float min = INFINITE;
//        Weight candidate;
//        list<Weight>::iterator itertodel;
//        for(list<Weight>::iterator iter = templist.begin();iter!=templist.end();iter++)
//        {
//           float deltaw = 0;
//           
//           for(int i = 0;i<dimension;i++)
//           {
//                deltaw = (tempweight.weighting[i]-iter->weighting[i])*(tempweight.weighting[i]-iter->weighting[i]);
//           }
//
//           if(min > deltaw)
//           {
//                candidate = *iter;
//                min = deltaw;
//                itertodel = iter;
//           }
//        }
//        tempweight = candidate;
//        sampleweight.push_back(candidate);
//        templist.erase(itertodel);
//
//    }
//    
//    //for(FloatVec::size_type ix=0;ix!=iter->weighting.size();ix++)
//	//	cout<<"Weight "<<ix+1<<" is"<<iter->weighting[ix]<<endl;
//
//}
/*
void SampleWeightsFromQP(list<Weight>& sampleweight, double quality_of_answer, double probability_guarantee, Weight w_origin, const vector< FloatVec >& incompar_points, long sizeneed = 0)
{
    long samplesize;
    int count =0;
	int dimension = w_origin.weighting.size();
	int size = incompar_points.size();
	srand( time(NULL) );

	samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer))+1;

    if(sizeneed>0)
        samplesize = sizeneed;
 
    //--------begin qp-------------------------//
    int nMatrices;
    int i, j, idx;
    //mtf_matrix_t *matrix_array = NULL;
    uint32_t nConstr;
    double *diag_H = NULL;
    double *vec_f = NULL;
    double *vec_b = NULL;
    uint32_t *vec_I = NULL;
    uint8_t *vec_S = NULL;
    double *vec_x = NULL;
    libqp_state_T exitflag;
    double TolAbs = 0.0;
    double TolRel = 1e-9;
    uint32_t MaxIter = 0xFFFFFFFF;
    double QP_TH = -LIBQP_PLUS_INF;

    exitflag = libqp_splx_solver(get_col_of_mat_H, diag_H, vec_f, vec_b,vec_I, vec_S, vec_x, nVars, MaxIter, TolAbs, TolRel, QP_TH, NULL);
}
*/
/*
void SampleWeightsFromQP(Engine * ep,list<Weight>& sampleweight, double quality_of_answer, double probability_guarantee, Weight w_origin, const vector< FloatVec >& incompar_points, long sizeneed = 0)
{
    long samplesize;
    int count =0;
	int dimension = w_origin.weighting.size();
	int size = incompar_points.size();

	samplesize = (long)(log(1-probability_guarantee)/log(1-quality_of_answer))+1;

    if(sizeneed>0)
        samplesize = sizeneed;
   
    while(samplesize>0)
    {
        
        int  num_hplane = rand()%(dimension-1)+1;
        //int  num_hplane = 1;
        IntVec candidates;
        candidates.reserve(dimension);

        for(int i = 0; i< num_hplane;i++)
        {
            int value = rand()%size;
            candidates.push_back(value);
            //cout<<"candidate"<<i+1<<":"<<value<<endl;
        }


    //---------------begin matlab quadratic programming------------------//
    
        mxArray *f = NULL, *Aeq = NULL, *d = NULL, *num = NULL, *result = NULL;

        double *w = NULL, *A = NULL;
        w = new double[dimension];

        for(int i = 0;i<dimension;i++)
        {
            w[i] = -w_origin.weighting[i];
        }
        
        A = new double[(num_hplane+1)*dimension];
    
        for(int i = 0;i<num_hplane;i++)
        {
            for(int j = 0;j<dimension;j++)
            {
                A[i*dimension+j] = incompar_points[candidates[i]][j];
                //cout<<A[i*dimension+j]<<"\t";
            }
            //cout<<endl;
        }

        for(int i = 0;i<dimension;i++)
        {
            A[num_hplane*dimension+i] = 1;
            //cout<<A[num_hplane*dimension+i]<<"\t";
        }
        
        //cout<<endl;
        int *dim = NULL;
        dim =  new int[1];
        dim[0] = dimension;
        int *num_equa = NULL;
        num_equa = new int[1];
        num_equa[0] = num_hplane;

        d = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
        num = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
        f = mxCreateDoubleMatrix(dimension, 1, mxREAL);
        Aeq = mxCreateDoubleMatrix(dimension, num_hplane+1, mxREAL);

        memcpy((void *)mxGetPr(d), (void *)dim,sizeof(dim));
        memcpy((void *)mxGetPr(num), (void *)num_equa,sizeof(num_equa));
        memcpy((void *)mxGetPr(Aeq), (void*)A,sizeof(A)*dimension*(num_hplane+1));
        memcpy((void *)mxGetPr(f), (void *)w,sizeof(w)*dimension);

        if(ep==NULL)
            cout<<"What the fuck"<<endl;
        engPutVariable(ep, "d", d);
        engPutVariable(ep, "num", num);
        engEvalString(ep, "H = eye(d)");
        engPutVariable(ep, "f", f);
        engPutVariable(ep, "Aeq", Aeq);
        //engEvalString(ep, "f = [0.2;0.2;0.2;0.2;0.2]");
        //engEvalString(ep, "Aeq = [1 2 3 -4 -5;1 1 1 1 1]");
        engEvalString(ep, "Aeq = Aeq'");
        engEvalString(ep, "beq = [zeros(num,1);1]");
        engEvalString(ep, "lb = zeros(d,1)");
        engEvalString(ep, "ub = ones(d,1)");
        engEvalString(ep, "options = optimset('LargeScale','off','TolX', 1e-2, 'TolFun', 1e-2)");

        engEvalString(ep, "x = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options)");

        //mxArray *result = NULL;

        //result = mxCreateDoubleMatrix(dimension, 1, mxREAL);

        Weight weights;
        float deltaW = 0;
        weights.weighting.resize(dimension);
        bool allpositive = true;    
        if((result = engGetVariable(ep, "x"))==NULL)
        {
            cout<<"Oops! you didn't create a variable X"<<endl;
        }
        else{
            //cout<<"M:"<<mxGetM(result)<<endl;
            //cout<<"N:"<<mxGetN(result)<<endl;
            double * weight = (double *)mxGetPr(result);
            for(int i = 0;i<dimension;i++)
            {
                if(fabs(weight[i])<0.0001)
                {
                    weight[i] = 0;
                }
                else if(weight[i]<0)
                {
                    allpositive = false;
                    break;
                }

				deltaW +=  (weights.weighting[i]-w_origin.weighting[i])*(weights.weighting[i]-w_origin.weighting[i]);
                weights.weighting[i] = weight[i];
            }
        }

        if(allpositive)
        {
			float sum = 0;
            for(int i = 0;i<dimension;i++)
            {
                sum += weights.weighting[i];
            }

            if(sum==0)
                break;
            weights.deltaw = sqrt(deltaW);
            sampleweight.push_back(weights);
            samplesize--;


   
             
            cout<<"sum:"<<sum<<endl;
            for(int i = 0;i<candidates.size();i++)
            {
			    float test_sum = 0;
			    for(FloatVec::size_type ix = 0;ix!=incompar_points[candidates[i]].size();ix++)
			    {
				    test_sum +=incompar_points[candidates[i]][ix]*weights.weighting[ix];
			    }

			    cout<<"Look at this sum, is it zero? "<<test_sum<<endl;
            }
            
        }
        
        //double *weight =(double *) mxGetPr(result);

        //for(int i = 0;i<dimension;i++)
        //{
            //cout<<"Weighting "<<i+1<<":"<<weight[i]<<endl;
        //}
        //cout<<"----------------"<<endl;
        delete [] w;
        delete [] dim;
        delete [] num_equa;
        delete [] A;
        //delete [] weight;
        mxDestroyArray(f);
        mxDestroyArray(Aeq);
        mxDestroyArray(d);
        mxDestroyArray(num);
        mxDestroyArray(result);

        count++;
        if(count%1000==0)
            cout<<"Total Attemp"<<endl;
    }

	sampleweight.sort(compWeight);
    return;
}
*/
//bool GaussianElim(vector< DoubleVec >& matrix, DoubleVec &ans, int row, int column, double & sum2)
//{
//	//cout<<"inside1"<<endl;
//	int i =0, j=0;
//
//	while(i<row&&j<column)
//	{
//		int maxi = i;
//		for(int k = i+1; k<row; k++)
//		{
//			if(fabs(matrix[k][j]) > fabs(matrix[maxi][j]))
//				maxi = k;
//		}
//
//		if(matrix[maxi][j]!=0)
//		{
//			double divider = matrix[maxi][j];
//			
//			if(i!=maxi)
//			{
//				swap(matrix[i],matrix[maxi]);
//			}
//
//			for(int v = j; v<column;v++)
//			{
//				matrix[i][v] /= divider;
//			}
//		
//			for(int u = i+1;u<row;u++)
//			{
//				for(int v = column-1; v>=0;v--)
//				{
//					matrix[u][v] = matrix[u][v] - matrix[u][j]*matrix[i][v];
//				}
//			}
//			//for(int ii=0;ii<row;ii++)
//				//for(int jj=0;jj<column;jj++)
//					//cout<<"Matrix["<<ii<<"]"<<"["<<jj<<"]:"<<matrix[ii][jj]<<endl;
//				
//			i++;
//		}
//
//		else
//			return false;
//		j++;
//	}
//   /* 
//    */
//	for(i = row-1;i>=0;i--)
//	{
//		
//		double sum = matrix[i][column-1];
//
//		int k = 0;
//		
//		for(j = column-2;j>i;j--)
//		{
//			//cout<<"value of k:"<<k<<endl;
//			//cout<<"Size of answer:"<<ans.size()<<endl;
//			sum -= matrix[i][j]*ans[k];
//			k++;
//		}
//
//		double value = sum/matrix[i][i];
//	    if(fabs(value)<=0.000001)
//            value = 0;
//        //else if(value <0)
//            //return false;
////			return false;
//        sum2 += value;
//		ans.push_back(value);
//	}
//
//    for(int i = 0;i<ans.size();i++)
//    {
//        cout<<"The "<<i+1<<" value:"<<ans[i]<<endl;
//    }
//	
//	//cout<<"inside2"<<endl;
//	return true;
//}

bool compWeight(Weight a, Weight b)
{
	if(a.deltaw<=b.deltaw)
		return true;
	else
		return false;
}
/*
void updateSkyline(deque<Answer> & skylineAnswer, int K, Weight W, ResultPool resultpool, float score)
{
	if(!skylineAnswer.empty())
	{
		Answer ans= skylineAnswer.back();

		if(ans.k<K||(ans.k == K&&ans.w.deltaw<W.deltaw))
		return;
	}
	
	Answer temp;
	temp.k = K;
	temp.w = W;
	temp.result = resultpool;
	temp.score = score;

	skylineAnswer.push_back(temp);
	return;	
}
*/
bool beDominate(deque<Answer>& skylineAnswer, Weight W, float missingscore, int kmin, int ID)
{
    //return false;
    //list<cell> sortAnswer;
    list<cell> sortAnswer;
    int dimension = W.weighting.size();
    int position = 0;
	for(deque<Answer>::iterator iter = skylineAnswer.begin(); iter!=skylineAnswer.end();iter++)
    {
        float deltaw = 0;
        for(int ix = 0;ix<dimension;ix++)
        {
            deltaw += (iter->w.weighting[ix]-W.weighting[ix])*(iter->w.weighting[ix]-W.weighting[ix]); 
        }
        cell tempcell;
        tempcell.ID = position;
        tempcell.data = deltaw;
        position++;
		//iter->w.deltaw = deltaw;
        sortAnswer.push_back(tempcell);
    }

	//sort(skylineAnswer.begin(),skylineAnswer.end(),CompareSkyAns);
    sortAnswer.sort();
	//for(deque<Answer>::reverse_iterator iter = skylineAnswer.rbegin(); iter!=skylineAnswer.rend();iter++)
	
	/*
	for(deque<Answer>::iterator iter = skylineAnswer.begin(); iter!=skylineAnswer.end();iter++)
	{
		if(Pruningbyrule2(iter->result.result, W, kmin, missingscore))
			return true;
	}
	*/

	
    for(list<cell>::iterator iter = sortAnswer.begin();iter!=sortAnswer.end();iter++)
	{
        deque<Answer>::iterator iterAns = skylineAnswer.begin()+iter->ID;
        if(iterAns->result.ID == ID)
            continue;
        //cout<<"skyline answer id:"<<iter->result.ID<<endl;
		if(Pruningbyrule2(iterAns->result.result, W, kmin, missingscore))
			return true;
        //if(count>0)
            //break;
	}
	
	return false;
}

bool Pruningbyrule2(const vector<FloatVec>& result, Weight w , int kmin, float missingscore)
{
    //return false;
	FloatVec score;
    int k = 0;
    //missingscore = ((int)(missingscore*1000))/1000.0;
    //cout<<"Missing score:"<<missingscore<<endl;
	for(vector<FloatVec >::const_iterator iter = result.begin(); iter!=result.end(); iter++)
	{
		float tempscore =0;
		//FloatVec temp = *iter;
		for(int ix=0; ix!=w.weighting.size(); ix++)
			tempscore+=(*(iter->begin()+ix))*(*(w.weighting.begin()+ix));


        //tempscore = ((int)(tempscore*1000))/1000.0;
		score.push_back(tempscore);
        if(tempscore>missingscore)
            k++;
        if(k>=kmin-1)
            return true;
	}
	//sort(score.begin(),score.end(), Cmp);
    //if(score[kmin-2]>missingscore)
		//return true;
	//else
	return false;
	//return score[kmin-2];
}

int Cmp(float a, float b)
{
	if(a>b)
		return true;
	else
		return false;
}

float Score(vector<Tuple>& missing_node, Weight& w)
{
    float score = INFINITE;
    for(vector<Tuple>::iterator iter = missing_node.begin();iter!=missing_node.end();iter++)
    {
        float minscore = 0;
	
	    for(FloatVec::size_type ix=0; ix!=w.weighting.size(); ix++)
            minscore += iter->data[ix]*w.weighting[ix];

        //minscore = ((int)(minscore*1000))/1000.0;
        if(minscore<score)
            score = minscore;
    }
    
    //score = ((int)(score*1000))/1000.0;
	return score;
}

int Estimatebounds(const Tuples& allTuples, vector<Tuple>& missing_node, IntVec position, vector< FloatVec >& incompar_points , const Index & index)
{  
    int size = missing_node.size();
    //set<int> num_dominate;
    int num_dominate = 0;
	
	float result;
	bool equal, dominate, bedominate;

	FloatVec data(position.size());
	data.reserve(position.size());
	//incompar_points.reserve(allTuples.size()*3);
    for(vector<Tuple>::iterator iter2 = missing_node.begin();iter2!=missing_node.end();iter2++)
    {
        float min = INFINITE;
	    for(FloatVec::size_type ix = 0;ix!=position.size();ix++){
                if(iter2->data[position[ix]]<min)
                    min = iter2->data[position[ix]];
        }

	//for(Tuples::iterator iter = allTuples.begin(); iter!=allTuples.end(); iter++)
	   for(Index::const_iterator iterindex = index.begin();iterindex!=index.end();iterindex++)
	   {
        if(min>=iterindex->key)
            break;
		//Tuple temp = *iter;		
        //Tuple temp = allTuples[iterindex->ID-1];
		equal = true;
		dominate = true; //dominate the missing tuple
		bedominate = true;//be dominated by the missing tuple
		//data.clear();
//cout<<"Tuple's ID:"<<temp.ID<<endl;
		    for(FloatVec::size_type ix = 0;ix!=position.size();++ix){
                //if(iter2->data[position[ix]]<min)
                    //min = iter2->data[position[ix]];
			    //float result = temp.data[position[ix]] - iter2->data[position[ix]];
			    //data[ix] = allTuples[iterindex->ID-1].data[position[ix]] - iter2->data[position[ix]];
			    data[ix] = allTuples[iterindex->ID-1].data[ix] - iter2->data[ix];
			//cout<<"result "<<ix+1<<":"<<result<<"\t";
			//}

			//for(FloatVec::size_type ix = 0;ix!=data.size();++ix){
			    if(data[ix]!=0)
				    equal = false;
				if(data[ix]>0)
				    bedominate = false;
			    else if(data[ix]<0)
				    dominate = false;
				
			}

			if(!equal&&dominate)
			{
                	//num_dominate.insert(allTuples[iterindex->ID-1].ID);
				++num_dominate;
					//break;
			}
			else if(!equal&&!dominate&&!bedominate){
			    incompar_points.push_back(data);
			}

			//next:;
		   //cout<<"THis ID:"<<temp.ID<<endl;
        }
	}
//exit(-1);
	//bounds bound;
	//bound.lowerbound = num_dominate + 1;
	//bound.upperbound = num_dominate + num_incomp_tuples + 1;

    return num_dominate + missing_node.size();
}

bool Equal(Tuple a, Tuple b)
{
    FloatVec dataA = a.data;
    FloatVec dataB = b.data;

    for(FloatVec::size_type ix = 0;ix!=dataA.size();ix++)
    {
        if(dataA[ix]!=dataB[ix])
            return false;
    }
    return true;
}

bool Dominate(Tuple a, Tuple b)
{
    FloatVec dataA = a.data;
    FloatVec dataB = b.data;

    for(FloatVec::size_type ix = 0;ix!=dataA.size();ix++)
    {
        if(dataA[ix]<dataB[ix])
            return false;
    }
    return true;
}

bool BeDominate(Tuple a, Tuple b)
{
    FloatVec dataA = a.data;
    FloatVec dataB = b.data;

    for(FloatVec::size_type ix = 0;ix!=dataA.size();ix++)
    {
        if(dataA[ix]>dataB[ix])
            return false;
    }
    return true;
}

void ReadIndex(char* filename, Index & index)
{
    ifstream file(filename, ifstream::in);

    if(!file.is_open())
    {
        cout<<"Can't open file!!"<<endl;
        return;
    }
    
    while(!file.eof())
    {
        Indexcell tempcell;
        file>>tempcell.ID;
        file>>tempcell.key;
        index.push_back(tempcell);
    }
    file.close();
    return;
}

int main(int argc, char** argv)
{
	char* filename = argv[1];
    char* indexfile = argv[2];
    int dimension = atoi(argv[3]);

	//cout<<"Dimension you issue:"<<dimension<<endl;
    //cout<<"Dimension:"<<dimension<<endl;

	MyConfigType config;
	map<string, string> configmap;

	if(!readConfig("whynot.conf", config, configmap)) {
		return -1;
	}

	float vary_preference[3] = {config.pma, config.nm, config.pmk};
	listConfig(configmap);
	
	string tvstime = "Tvstime.dat";
    string tvspenalty = "Tvspenalty.dat";
    string pvstime = "Pvstime.dat";
    string pvspenalty = "Pvspenalty.dat";
    string missingnum = "MissingNum.dat";
    string korigin = "Korigin.dat";
    string thresholdk = "Thresholdk.dat";
	string missingpos = "MissingPosvstime.dat";

	string path;
	if(config.datatype == "i") {
		path = "../paper/Graph/I/";
	}

	else if(config.datatype == "c") {
		path = "../paper/Graph/C/";
	}

	else if(config.datatype == "ac") {
		path = "../paper/Graph/AC/";
	}

	tvstime = path + tvstime;
	tvspenalty = path + tvspenalty;
	pvstime = path + pvstime;
	pvspenalty = path + pvspenalty;
	missingnum = path + missingnum;
	korigin = path + korigin;
	thresholdk = path + thresholdk;
	missingpos = path + missingpos;
    
	Index index;
    ReadIndex(indexfile, index);
    cout<<"Index size:"<<index.size()<<endl;
    
	srand( time(NULL) );
	//vector<vector<cell> > storebyrow;
    //map<int, FloatVec > index;

	//PreComputeTA(filename, storebyrow, index);	
        //--------------------------------//
    //
    //
    //const char* argvv [] = {"", "--silent"};
    //octave_main(2, (char **) argvv, true);

    //----------------------------------//
    //    Engine *ep;
    //    if (!(ep = engOpen("\0"))) {
    //        fprintf(stderr, "\nCan't start MATLAB engine\n");

    //    }



    //--------------------------------//
    DomGaph dominantgraph;
	map<int, node> realdata;
	bool HasLowerbound = false;

	cout<<"Constructing!!!"<<endl;
	ReconstructDG(dominantgraph, realdata, HasLowerbound);	
	cout<<"Construction completes!!!"<<endl;

	//---------------------------------//
	Tuples allTuples;
    ReadData(filename, allTuples);
    //cout<<"Size of realdata:"<<realdata.size()<<endl;
	//cout<<"Size of tuples:"<<allTuples.size()<<endl;
	//--------test time------------------//
	clock_t start_sample, stop_sample;
	clock_t start, stop;
	clock_t start_query, stop_query;
    clock_t start_pruning, stop_pruning;

    //cout<<"Input dimension:"<<endl;
    //scanf("%d", &dimension);
	srand(time(NULL));
//-----------------initialize original weighting vector--------------------------//
    
	float deltak;
	float deltaw;

    int default_k_origin = config.korg;
    float default_quality = config.t;
    float default_prob = config.prob;
    int repeat = config.repeat;

    IntVec position;
    position.reserve(dimension);

	Weight weight_origin;
	weight_origin.deltaw = 0;
	
    float deltaWU = 0;
	for(int i=0; i<dimension; i++)
	{
		//weight_origin.weighting.push_back(value[i]/sum);
        float value = 1.0/dimension;
		weight_origin.weighting.push_back(value);
        position.push_back(i);
        deltaWU +=value*value;
		//cout<<"Weight "<<i<<":"<<weight_origin.weighting[i]<<endl;
    }

    deltaWU = sqrt(1+deltaWU);
    
    //deltaWU = (dimension-1)*1.0/(dimension*dimension);
    //deltaWU += (dimension-1)*1.0*(dimension-1)/(dimension*dimension);
    //deltaWU = sqrt(deltaWU);
    cout<<"Delta weighting upper:"<<deltaWU<<endl;
    

if(testall)
{
    FloatVec qualityvec = config.tvary;
	
	int quality_size = qualityvec.size();

    FloatVec timeresult(quality_size*3, 0);
    FloatVec penaltyresult(quality_size*3, 0);
    IntVec numtopk(quality_size*3, 0), numbuffer(quality_size*3, 0);

    IntVec samplesizevec;
    for(int i = 1;i<=5;i++)
    {
        samplesizevec.push_back(1000*i);
    }
    //timeresult.resize(24);
    //penaltyresult.reserve(24);
    ///FloatVec pvspresult(12);
//////////////////////////////////////////////////////////////////////////////////////
//
// The code below is the main code for answering why-not question. I repeat the main code 
// many times to test the efficiency of the algorithm while varying different parameters.
// you will see that I'm vary lazy, because I just copy the same code many times.
//
//////////////////////////////////////////////////////////////////////////////////////
for(int num_of_test = 0;num_of_test<repeat;num_of_test++)
{
	//break;
    for(DoubleVec::size_type qualityindex = 0;qualityindex != qualityvec.size();qualityindex++)
    {
	 
			int k_origin = default_k_origin;                 //original value of k
            int missing_num = 1;
		    //cout<<"Missing Num:"<<missing_num<<endl;
			vector<Tuple> missing_node;
			IntVec missing_ID;
		    
		    vector<cell> gen_missing_point = TopKDG(realdata, dominantgraph, 10*k_origin+1, weight_origin);
		    missing_ID.push_back(gen_missing_point.rbegin()->ID);
		    for(IntVec::iterator iter=missing_ID.begin();iter!=missing_ID.end();iter++)
		    {
			missing_node.push_back(allTuples[*iter-1]); 
		    }
		    //exit(-1); 
		    for(IntVec::iterator iter = missing_ID.begin();iter!=missing_ID.end();iter++)
			    cout<<"Missing ID:"<<*iter<<endl;

		    double quality_of_answer = qualityvec[qualityindex];             //quallity of the answer
        
			float probability_guarantee = default_prob;          //probability to get such an answer  
			vector< FloatVec > incompar_points;
		    incompar_points.reserve(allTuples.size()*missing_num);	
	

		    float sampletime = 0;
		    float estimateboundtime = 0;
		    long sample_size = 0;
	        sample_size = (long)(log(1-probability_guarantee)/log(1-quality_of_answer/100))+1;
            //cout<<"Total Sample size:"<<sample_size<<endl;
            FloatVec BestScoreVec(3, INFINITE);
            FloatVec TotalTimeVec(3, 0);
    
		    start = clock();
			//-----------------------rank lower and upper bounds---------------------------//
			int lowerbound = Estimatebounds(allTuples,missing_node, position, incompar_points, index);
		    stop = clock();
		    estimateboundtime = float(stop-start)/CLOCKS_PER_SEC;
      while(sample_size>0)
      {
			list<Weight> weights;
            long sampleforeach = 100000;
            if(sample_size<sampleforeach)
            {
                sampleforeach = sample_size;
                sample_size = 0;
            }
            else
                sample_size -= sampleforeach;
			//cout<<"Total time:"<<float(stop- start) / CLOCKS_PER_SEC <<endl;
		    //cout<<"What is the lower bound:"<<lowerbound<<endl;
		    //cout<<"What is the num of incomparable:"<<incompar_points.size()<<endl;
		    //exit(-1);
	
		    start_sample=clock();
		    weights = SampleWeights(quality_of_answer,probability_guarantee, weight_origin);
            //weights.clear();
		    //SampleWeightsFromB(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points, sampleforeach);
		    //SampleWeightsFromP(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);

			stop_sample=clock();
		    sampletime = float(stop_sample - start_sample) / CLOCKS_PER_SEC;
			cout<<"Sample size:"<<weights.size()<<endl;
			cout<<"Sample time:"<<float(stop_sample- start_sample) / CLOCKS_PER_SEC <<endl;
		    //exit(-1);


		    FloatVec penalty(3);
		    penalty[0] = 0.9;
		    penalty[1] = 0.5;
		    penalty[2] = 0.1;

		    for(FloatVec::size_type ix = 0;ix!=penalty.size();ix++)
		    {    
			    float penalty_k = penalty[ix];          //ipenalty for changing k
			    float penalty_w = 1 - penalty_k;             //penalty for changing w

			cout<<"penalty_k:"<<penalty_k<<endl;
		 
			int threshold_k = allTuples.size();
			float totaltime = 0;
			float pruningtime = 0;
            float topk_time = 0;
			float bestquality = INFINITE;
			float worstquality = 0;
			float sumquality = 0;
            clock_t prune_start, prune_stop;
            clock_t topk_start, topk_stop;
		  
			int num_bufferused = 0;
		 
			    int ranklowerb;
			    int rankupperb;	
			    float tmpthreshold = INFINITE;
			float threshold_w = INFINITE;
	
			if(lowerbound>threshold_k)
			{
			    cout<<"The rank lower bound:"<<lowerbound<<endl;
			    cout<<"Can't rank into the threshold!"<<endl;
			    exit(-1);
			}

			int deltaKlowerbound = 0;
			    if(lowerbound - k_origin > deltaKlowerbound)
				    deltaKlowerbound = lowerbound - k_origin;

			    Answer BestAns;
			    list<cell> resultset;
			    deque<Answer> skylineAnswer;
	
			start = clock();
			    float missing_score = Score(missing_node, weight_origin);

			    int k_old, kmin, kprune;
                int k_delta;
			int topk_id = 1;
			kmin = allTuples.size();
			ResultPool resultpool;

                topk_start = clock(); 
		    	k_old = ProgressiveTopK(realdata, dominantgraph, missing_score, weight_origin, resultpool, kmin, topk_id);
                topk_stop = clock();
                topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
			    kmin = k_old;
				BestAns.k = k_old;
                BestAns.result.ID = resultpool.ID;
				BestAns.result.result.resize(1000000);
			    BestAns.result.result.assign(resultpool.result.begin(),resultpool.result.end());
				BestAns.w = weight_origin;
				//BestAns.score = penalty_k*(BestAns.k-k_origin) + penalty_w* BestAns.w.deltaw;
				BestAns.score = penalty_k*(BestAns.k-k_origin)/(k_old-k_origin) + penalty_w* BestAns.w.deltaw/deltaWU;
				skylineAnswer.push_back(BestAns);
		
				tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
				//tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;

			    if(threshold_w > tmpthreshold)
				    threshold_w = tmpthreshold;

             

			    unsigned int i = 0, j = 0, num = 1;
			    int test_count = 0;
				
				/////////////////////////////////////////////////////////////////////////
				//
				// The main loop for testing the candidate weigthing vectors
				//
				/////////////////////////////////////////////////////////////////////////
			    for(list<Weight>::iterator iter = weights.begin(); iter!=weights.end(); iter++)
			    {
				    Weight w = *iter;
				    if(w.deltaw>threshold_w)
                    {
					    break;
                        //continue;
                    }
                    k_delta = (BestAns.score-penalty_w*w.deltaw/deltaWU)*(k_old-k_origin)/penalty_k;

                    if(k_delta<deltaKlowerbound)
                    {
                        //num++;
                        //continue;
                        break;
                    }
                    if((k_delta+k_origin)<kmin)
                    {
                        kprune = k_delta+k_origin+1;
                        //cout<<"kmin"<<kmin<<endl;
                        //cout<<"Kprune"<<kprune<<endl;
                    }
                    else
                        kprune = kmin;

                    //kprune = kmin;
			        bool prunebyrule2 = false;
			    
				    missing_score = Score(missing_node, w);

                    //cout<<"Missing score:"<<missing_score<<endl;
                    
				    if(Pruningbyrule2(resultpool.result, w, kprune, missing_score))
				    {
					    ++i;
			        	prunebyrule2 = true;
			        }
			        
		
			        if(prunebyrule2)
                    {
			            continue;
                    }
                    bool prunebyrule3 = false;
				    if(beDominate(skylineAnswer, w, missing_score, kprune, topk_id))
				    {
					    ++j;
					    prunebyrule3 = true;
				    }
		
		
			        if(prunebyrule3)
                    {
			            continue;
                    }

			        topk_id++;
                    topk_start = clock();
		 		    int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, kprune-1, topk_id);
		 		    //int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, allTuples.size(), topk_id, false);
                    topk_stop = clock();
                    topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
				    ++num;

                    //cout<<"k:"<<k<<endl;
				    if(k == OUTSIDE_THRESHOLD)		
				    {
					    //continue;
                        Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            skylineAnswer.push_back(tempAnswer);

				    }
	
				    //if(k<kmin)
                     else
    			     {
                        if(k < k_origin)
                            k = k_origin;
                        if(kmin > k)
                            kmin = k;
				        float score = penalty_k*(k-k_origin)/(k_old-k_origin) + penalty_w*w.deltaw/deltaWU;
				        //float score = penalty_k*(kmin-k_origin) + penalty_w*w.deltaw;
                        //cout<<"score:"<<score<<endl;
                        //cout<<"Bestscore:"<<BestAns.score<<endl;

			            Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            tempAnswer.score = score;
			            skylineAnswer.push_back(tempAnswer);
		
				        if(score < BestAns.score)
				        {
					        BestAns = tempAnswer;

					        //float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;
					        float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
			
		        			if(threshold_w > tmpthreshold)
				        		threshold_w = tmpthreshold;
				        }

			        }
			    }

			//stop_query=clock();
			stop=clock();
	
			//int size = weights.size();
			totaltime += float(stop- start) / CLOCKS_PER_SEC;
            TotalTimeVec[ix] += sampletime+totaltime;
            if(BestAns.score<BestScoreVec[ix])
                BestScoreVec[ix] = BestAns.score;
		    //cout<<"Avarage bufferresult used:"<<num_bufferused<<endl;
            //cout<<"Number of Topk:"<<num+1<<endl;
			//cout<<"------------------------------------"<<endl;
		    //cout<<"Sample size:"<<sample_size<<endl;
		    //cout<<"Sample Time:"<<sampletime<<endl;
		    //cout<<"Pruning Time:"<<pruningtime<<endl;
		    //cout<<"Topk time:"<<topk_time/num_of_test<<endl;
		    //cout<<"Assign time"<<assign_time/num_of_test<<endl;
            //

			cout<<"Average Time:"<<sampletime+totaltime+estimateboundtime<<endl;
		    //cout<<"Average quality"<<sumquality<<endl;	
			
            cout<<"------------------------------------"<<endl;
            //numtopk[qualityindex*3+ix] += num_bufferused;
            //numbuffer[qualityindex*3+ix] +=num+1;

		}//loop penalty
            
      }//while sample size decrease
      for(int i = 0;i<3;i++) 
      {
		    timeresult[qualityindex*3+i] += TotalTimeVec[i]+estimateboundtime;
		    penaltyresult[qualityindex*3+i] += BestScoreVec[i];
      }
    }//loop for quality increase 
  }//num of test    


/*
    ofstream usedfile("Statis.dat", ofstream::out);

    for(int i = 0;i<qualityvec.size();i++)
    {
        usedfile<<i;
        for(int j = 0;j<3;j++)
        {
            usedfile<<' '<<numtopk[i*3+j]/repeat<<","<<numbuffer[i*3+j]/repeat;
        
        }
        usedfile<<endl;
    }
usedfile.close();
*/

    ofstream tvstfile(tvstime.c_str(), ofstream::out);
    ofstream tvspfile(tvspenalty.c_str(), ofstream::out);
    
    for(int i = 0;i<qualityvec.size();i++)
    {
        tvstfile<<qualityvec[i]<<"\%";
        for(int j=0;j<3;j++)
        {
            tvstfile<<' '<<timeresult[i*3+j]/repeat;
        }
        tvstfile<<endl;
    }
    tvstfile.close();

    FloatVec maxpenalty(3, 0);

    for(int i = 0;i<qualityvec.size();i++)
    {
        for(int j = 0;j<3;j++)
        {
            if(maxpenalty[j]<penaltyresult[i*3+j])
                maxpenalty[j] = penaltyresult[i*3+j];
        }
    }

    for(int i = 0;i<qualityvec.size();i++)
    {
        tvspfile<<qualityvec[i]<<"\%";
        for(int j=0;j<3;j++)
        {
            //tvspfile<<' '<<penaltyresult[i*3+j]/(maxpenalty[j]);
            tvspfile<<' '<<penaltyresult[i*3+j]/repeat;
        }
        tvspfile<<endl;
    }
    tvspfile.close();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The code below is the same as the code above section, just vary the Pr this time
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
   FloatVec probvec = config.probvary;
   int num_prob = 5;
    timeresult.assign(num_prob*3, 0);
    penaltyresult.assign(num_prob*3, 0);
for(int num_of_test = 0; num_of_test <repeat;num_of_test++)
{
    //break; 
    for(FloatVec::size_type probix = 0;probix!=probvec.size();probix++)
    {
	 
		    int k_origin = default_k_origin;                 //original value of k
		    int missing_num = 1;
		  //  cout<<"Missing Num:"<<missing_num<<endl;
			vector<Tuple> missing_node;
            IntVec missing_ID;
		    vector<cell> gen_missing_point = TopKDG(realdata, dominantgraph, 10*k_origin+1, weight_origin);
		    missing_ID.push_back(gen_missing_point.rbegin()->ID);
		    for(IntVec::iterator iter=missing_ID.begin();iter!=missing_ID.end();iter++)
		    {
			missing_node.push_back(allTuples[*iter-1]); 
		    }
		    //exit(-1); 
		    for(IntVec::iterator iter = missing_ID.begin();iter!=missing_ID.end();iter++)
			    cout<<"Missing ID:"<<*iter<<endl;

		    float quality_of_answer = default_quality;             //quallity of the answer
			float probability_guarantee = probvec[probix];          //probability to get such an answer  
			vector< FloatVec > incompar_points;
		    incompar_points.reserve(allTuples.size()*missing_num);	
	

		    float sampletime = 0;
		    float estimateboundtime = 0;
		    int sample_size = 0;
			FloatVec BestScoreVec(3, INFINITE);
            FloatVec TotalTimeVec(3, 0);
            list<Weight> weights;

		    start = clock();
			//-----------------------rank lower and upper bounds---------------------------//
			int lowerbound = Estimatebounds(allTuples,missing_node, position, incompar_points, index);
		    stop = clock();
			//cout<<"Total time:"<<float(stop- start) / CLOCKS_PER_SEC <<endl;
		    //cout<<"What is the lower bound:"<<lowerbound<<endl;
		    //cout<<"What is the num of incomparable:"<<incompar_points.size()<<endl;
		    //exit(-1);
		    estimateboundtime = float(stop-start)/CLOCKS_PER_SEC;
	
		    start_sample=clock();
            weights.clear();
		    weights = SampleWeights(quality_of_answer,probability_guarantee, weight_origin);
		    //SampleWeightsFromB(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);
		    //SampleWeightsFromP(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);

			stop_sample=clock();
			sample_size = weights.size();
			cout<<"Sample size:"<<weights.size()<<endl;
			//cout<<"Sample time:"<<float(stop_sample- start_sample) / CLOCKS_PER_SEC <<endl;
		    sampletime = float(stop_sample - start_sample) / CLOCKS_PER_SEC;
		    //exit(-1);


		    FloatVec penalty(3);
		    penalty[0] = 0.9;
		    penalty[1] = 0.5;
		    penalty[2] = 0.1;
		 
		    for(FloatVec::size_type ix = 0;ix!=penalty.size();ix++)
		    {    
			    float penalty_k = penalty[ix];          //ipenalty for changing k
			    float penalty_w = 1 - penalty_k;             //penalty for changing w

			    //cout<<"penalty_k:"<<penalty_k<<endl;
		 
			int threshold_k = allTuples.size();
			float totaltime = 0;
			float pruningtime = 0;
            float topk_time = 0;
			float bestquality = INFINITE;
			float worstquality = 0;
			float sumquality = 0;
            clock_t prune_start, prune_stop;
            clock_t topk_start, topk_stop;
		  
			int num_bufferused = 0;
		 
			    int ranklowerb;
			    int rankupperb;	
			    float tmpthreshold = INFINITE;
			float threshold_w = INFINITE;
	
			if(lowerbound>threshold_k)
			{
			    cout<<"The rank lower bound:"<<lowerbound<<endl;
			    cout<<"Can't rank into the threshold!"<<endl;
			    exit(-1);
			}

			int deltaKlowerbound = 0;
			    if(lowerbound - k_origin > deltaKlowerbound)
				    deltaKlowerbound = lowerbound - k_origin;

			    Answer BestAns;
			    list<cell> resultset;
			    deque<Answer> skylineAnswer;
	
			start = clock();
			    float missing_score = Score(missing_node, weight_origin);

			    int k_old, kmin, kprune;
                int k_delta;
			int topk_id = 10;
			kmin = allTuples.size();
			ResultPool resultpool;

                topk_start = clock(); 
		    	k_old = ProgressiveTopK(realdata, dominantgraph, missing_score, weight_origin, resultpool, kmin, topk_id);
                topk_stop = clock();
                topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
			    kmin = k_old;
				BestAns.k = k_old;
                BestAns.result.ID = resultpool.ID;
				BestAns.result.result.resize(1000000);
			    BestAns.result.result.assign(resultpool.result.begin(),resultpool.result.end());
				BestAns.w = weight_origin;
				//BestAns.score = penalty_k*(BestAns.k-k_origin) + penalty_w* BestAns.w.deltaw;
				BestAns.score = penalty_k*(BestAns.k-k_origin)/(k_old-k_origin) + penalty_w* BestAns.w.deltaw/deltaWU;
				skylineAnswer.push_back(BestAns);
		
				tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
				//tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;

			    if(threshold_w > tmpthreshold)
				    threshold_w = tmpthreshold;

             

			    unsigned int i = 0, j = 0, num = 1;
			    int test_count = 0;
			    for(list<Weight>::iterator iter = weights.begin(); iter!=weights.end(); iter++)
			    {
				    Weight w = *iter;
				    if(w.deltaw>threshold_w)
                    {
					    break;
                        //continue;
                    }
                    k_delta = (BestAns.score-penalty_w*w.deltaw/deltaWU)*(k_old-k_origin)/penalty_k;

                    if(k_delta<deltaKlowerbound)
                    {
                        //num++;
                        //continue;
                        break;
                    }
                    if((k_delta+k_origin)<kmin)
                    {
                        kprune = k_delta+k_origin+1;
                        //cout<<"kmin"<<kmin<<endl;
                        //cout<<"Kprune"<<kprune<<endl;
                    }
                    else
                        kprune = kmin;

                    //kprune = kmin;
			        bool prunebyrule2 = false;
			    
				    missing_score = Score(missing_node, w);

                    //cout<<"Missing score:"<<missing_score<<endl;
                    
				    if(Pruningbyrule2(resultpool.result, w, kprune, missing_score))
				    {
					    ++i;
			        	prunebyrule2 = true;
			        }
			        
		
			        if(prunebyrule2)
                    {
			            continue;
                    }
                    bool prunebyrule3 = false;
				    if(beDominate(skylineAnswer, w, missing_score, kprune, topk_id))
				    {
					    ++j;
					    prunebyrule3 = true;
				    }
		
		
			        if(prunebyrule3)
                    {
			            continue;
                    }

			        topk_id++;
                    topk_start = clock();
		 		    int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, kprune-1, topk_id);
		 		    //int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, allTuples.size(), topk_id, false);
                    topk_stop = clock();
                    topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
				    ++num;

                    //cout<<"k:"<<k<<endl;
				    if(k == OUTSIDE_THRESHOLD)		
				    {
					    //continue;
                        Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            skylineAnswer.push_back(tempAnswer);

				    }
	
				    //if(k<kmin)
                     else
    			     {
                        if(k < k_origin)
                            k = k_origin;
                        if(kmin > k)
                            kmin = k;
				        float score = penalty_k*(k-k_origin)/(k_old-k_origin) + penalty_w*w.deltaw/deltaWU;
				        //float score = penalty_k*(kmin-k_origin) + penalty_w*w.deltaw;
                        //cout<<"score:"<<score<<endl;
                        //cout<<"Bestscore:"<<BestAns.score<<endl;

			            Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            tempAnswer.score = score;
			            skylineAnswer.push_back(tempAnswer);
		
				        if(score < BestAns.score)
				        {
					        BestAns = tempAnswer;

					        //float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;
					        float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
			
		        			if(threshold_w > tmpthreshold)
				        		threshold_w = tmpthreshold;
				        }

			        }
			    }

			//stop_query=clock();
			stop=clock();
	
			//int size = weights.size();
			totaltime += float(stop- start) / CLOCKS_PER_SEC;
            TotalTimeVec[ix] += sampletime+totaltime;
            if(BestAns.score<BestScoreVec[ix])
                BestScoreVec[ix] = BestAns.score;
		    //cout<<"Avarage bufferresult used:"<<num_bufferused<<endl;
            //cout<<"Number of Topk:"<<num+1<<endl;
			//cout<<"------------------------------------"<<endl;
		    //cout<<"Sample size:"<<sample_size<<endl;
		    //cout<<"Sample Time:"<<sampletime<<endl;
		    //cout<<"Pruning Time:"<<pruningtime<<endl;
		    //cout<<"Topk time:"<<topk_time/num_of_test<<endl;
		    //cout<<"Assign time"<<assign_time/num_of_test<<endl;
            //

			cout<<"Average Time:"<<sampletime+totaltime+estimateboundtime<<endl;
		    //cout<<"Average quality"<<sumquality<<endl;	
			
            cout<<"------------------------------------"<<endl;
            //numtopk[qualityindex*3+ix] += num_bufferused;
            //numbuffer[qualityindex*3+ix] +=num+1;

		}//loop penalty
            for(int i = 0;i<3;i++) 
            {
		        timeresult[probix*3+i] += TotalTimeVec[i]+estimateboundtime;
		        penaltyresult[probix*3+i] += BestScoreVec[i];
            }
        }//loop P
 } //num of test  
    ofstream pvstfile(pvstime.c_str(), ofstream::out);
    ofstream pvspfile(pvspenalty.c_str(), ofstream::out);
    
    for(int i = 0;i<probvec.size();i++)
    {
        pvstfile<<probvec[i];
        for(int j=0;j<3;j++)
        {
            pvstfile<<' '<<timeresult[i*3+j]/repeat;
        }
        pvstfile<<endl;
    }
    pvstfile.close();
    
    FloatVec max(3, 0) ;
    for(int j = 0;j<3;j++)
        for(int i = 0;i<probvec.size();i++)
        {
            if(penaltyresult[i*3+j]>max[j])
                max[j] = penaltyresult[i*3+j];
        }

    for(int i = 0;i<probvec.size();i++)
    {
        pvspfile<<probvec[i];
        for(int j=0;j<3;j++)
        {
            //pvspfile<<' '<<penaltyresult[i*3+j]/max[j];
            pvspfile<<' '<<penaltyresult[i*3+j]/repeat;
        }
        pvspfile<<endl;
    }
    pvspfile.close();


    //return 0;

    IntVec koriginvec = config.kvary;
	int kvarysize = koriginvec.size();
    FloatVec koresult;
    //koresult.reserve(12);
    koresult.assign(kvarysize*3, 0);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// The code is the same as above, just varrying the k value of the top-k query
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout<<"Varying K origin:"<<endl;
    for(int num_of_test = 0; num_of_test < repeat; num_of_test++)
    {
       for(IntVec::size_type koindex = 0;koindex!=koriginvec.size();koindex++)
	   {
			//int k_origin = default_k_origin;                 //original value of k
            int k_origin = koriginvec[koindex];
		    //cout<<"Missing Num:"<<missing_num<<endl;
			vector<Tuple> missing_node;
			IntVec missing_ID;
			vector<cell> gen_missing_point;
			if(koindex==2)
				gen_missing_point = TopKDG(realdata, dominantgraph, 10*k_origin, weight_origin);
			else
		    	gen_missing_point = TopKDG(realdata, dominantgraph, 10*k_origin+1, weight_origin);

		    missing_ID.push_back(gen_missing_point.rbegin()->ID);
            /*
		    vector<cell> gen_missing_point = TopKDG(realdata, dominantgraph, 6000, weight_origin);
            for(int i = 0;i<missing_num;i++)
            {
		        missing_ID.push_back(gen_missing_point[1000*(i+1)].ID);
            }
            */
		    for(IntVec::iterator iter=missing_ID.begin();iter!=missing_ID.end();iter++)
		    {
				missing_node.push_back(allTuples[*iter-1]); 
		    }
		    //exit(-1); 
		    for(IntVec::iterator iter = missing_ID.begin();iter!=missing_ID.end();iter++)
			    cout<<"Missing ID:"<<*iter<<endl;

		    float quality_of_answer = default_quality;             //quallity of the answer
			float probability_guarantee = default_prob;          //probability to get such an answer  
			vector< FloatVec > incompar_points;
		    incompar_points.reserve(allTuples.size());	
	

		    float sampletime = 0;
		    float estimateboundtime = 0;
		    int sample_size = 0;
			list<Weight> weights;

		    start = clock();
			//-----------------------rank lower and upper bounds---------------------------//
			int lowerbound = Estimatebounds(allTuples,missing_node, position, incompar_points, index);
		    stop = clock();
			//cout<<"Total time:"<<float(stop- start) / CLOCKS_PER_SEC <<endl;
		    //cout<<"What is the lower bound:"<<lowerbound<<endl;
		    //cout<<"What is the num of incomparable:"<<incompar_points.size()<<endl;
		    //exit(-1);
		    estimateboundtime = float(stop-start)/CLOCKS_PER_SEC;
	
		    start_sample=clock();
		    weights = SampleWeights(quality_of_answer,probability_guarantee, weight_origin);
            //weights.clear(); 
		    //SampleWeightsFromB(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);
		    //SampleWeightsFromP(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);

			stop_sample=clock();
			sample_size = weights.size();
			//cout<<"Sample size:"<<weights.size()<<endl;
			//cout<<"Sample time:"<<float(stop_sample- start_sample) / CLOCKS_PER_SEC <<endl;
		    sampletime = float(stop_sample - start_sample) / CLOCKS_PER_SEC;
		    //exit(-1);


		    FloatVec penalty(3);
		    penalty[0] = 0.9;
		    penalty[1] = 0.5;
		    penalty[2] = 0.1;
		 
		    for(FloatVec::size_type ix = 0;ix!=penalty.size();ix++)
		    {    
			    float penalty_k = penalty[ix];          //ipenalty for changing k
			    float penalty_w = 1 - penalty_k;             //penalty for changing w

			    //cout<<"penalty_k:"<<penalty_k<<endl;
		 
			int threshold_k = allTuples.size();
			float totaltime = 0;
			float pruningtime = 0;
            float topk_time = 0;
			float bestquality = INFINITE;
			float worstquality = 0;
			float sumquality = 0;
            clock_t prune_start, prune_stop;
            clock_t topk_start, topk_stop;
		  
			int num_bufferused = 0;
		 
			    int ranklowerb;
			    int rankupperb;	
			    float tmpthreshold = INFINITE;
			float threshold_w = INFINITE;
	
			if(lowerbound>threshold_k)
			{
			    cout<<"The rank lower bound:"<<lowerbound<<endl;
			    cout<<"Can't rank into the threshold!"<<endl;
			    exit(-1);
			}

			int deltaKlowerbound = 0;
			    if(lowerbound - k_origin > deltaKlowerbound)
				    deltaKlowerbound = lowerbound - k_origin;

			    Answer BestAns;
			    list<cell> resultset;
			    deque<Answer> skylineAnswer;
	
			start = clock();
			    float missing_score = Score(missing_node, weight_origin);

			    int k_old, kmin, kprune;
                int k_delta;
			int topk_id = 10;
			kmin = allTuples.size();
			ResultPool resultpool;

                topk_start = clock(); 
		    	k_old = ProgressiveTopK(realdata, dominantgraph, missing_score, weight_origin, resultpool, kmin, topk_id);
                topk_stop = clock();
                topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
			    kmin = k_old;
				BestAns.k = k_old;
                BestAns.result.ID = resultpool.ID;
				BestAns.result.result.resize(1000000);
			    BestAns.result.result.assign(resultpool.result.begin(),resultpool.result.end());
				BestAns.w = weight_origin;
				//BestAns.score = penalty_k*(BestAns.k-k_origin) + penalty_w* BestAns.w.deltaw;
				BestAns.score = penalty_k*(BestAns.k-k_origin)/(k_old-k_origin) + penalty_w* BestAns.w.deltaw/deltaWU;
				skylineAnswer.push_back(BestAns);
		
				tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
				//tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;

			    if(threshold_w > tmpthreshold)
				    threshold_w = tmpthreshold;

             

			    unsigned int i = 0, j = 0, num = 1;
			    int test_count = 0;
			    for(list<Weight>::iterator iter = weights.begin(); iter!=weights.end(); iter++)
			    {
				    Weight w = *iter;
				    if(w.deltaw>threshold_w)
                    {
					    break;
                        //continue;
                    }
                    k_delta = (BestAns.score-penalty_w*w.deltaw/deltaWU)*(k_old-k_origin)/penalty_k;

                    if(k_delta<deltaKlowerbound)
                    {
                        //num++;
                        //continue;
                        break;
                    }
                    if((k_delta+k_origin)<kmin)
                    {
                        kprune = k_delta+k_origin+1;
                        //cout<<"kmin"<<kmin<<endl;
                        //cout<<"Kprune"<<kprune<<endl;
                    }
                    else
                        kprune = kmin;

                    //kprune = kmin;
			        bool prunebyrule2 = false;
			    
				    missing_score = Score(missing_node, w);

                    //cout<<"Missing score:"<<missing_score<<endl;
                    
				    if(Pruningbyrule2(resultpool.result, w, kprune, missing_score))
				    {
					    ++i;
			        	prunebyrule2 = true;
			        }
			        
		
			        if(prunebyrule2)
                    {
			            continue;
                    }
                    bool prunebyrule3 = false;
				    if(beDominate(skylineAnswer, w, missing_score, kprune, topk_id))
				    {
					    ++j;
					    prunebyrule3 = true;
				    }
		
		
			        if(prunebyrule3)
                    {
			            continue;
                    }

			        topk_id++;
                    topk_start = clock();
		 		    int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, kprune-1, topk_id);
		 		    //int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, allTuples.size(), topk_id, false);
                    topk_stop = clock();
                    topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
				    ++num;

                    //cout<<"k:"<<k<<endl;
				    if(k == OUTSIDE_THRESHOLD)		
				    {
					    //continue;
                        Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            skylineAnswer.push_back(tempAnswer);

				    }
	
				    //if(k<kmin)
                     else
    			     {
                        if(k < k_origin)
                            k = k_origin;
                        if(kmin > k)
                            kmin = k;
				        float score = penalty_k*(k-k_origin)/(k_old-k_origin) + penalty_w*w.deltaw/deltaWU;
				        //float score = penalty_k*(kmin-k_origin) + penalty_w*w.deltaw;
                        //cout<<"score:"<<score<<endl;
                        //cout<<"Bestscore:"<<BestAns.score<<endl;

			            Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            tempAnswer.score = score;
			            skylineAnswer.push_back(tempAnswer);
		
				        if(score < BestAns.score)
				        {
					        BestAns = tempAnswer;

					        //float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;
					        float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
			
		        			if(threshold_w > tmpthreshold)
				        		threshold_w = tmpthreshold;
				        }

			        }
			    }

			//stop_query=clock();
			stop=clock();
	
			//int size = weights.size();
			//cout<<"Weight size:"<<size<<endl;
		    
            /*
			cout<<"Pruned by rule 1:"<<(weights.size()+1-i-j-num)<<endl;
			cout<<"Pruned by rule 2:"<<i<<endl;
			cout<<"Pruned by rule 3:"<<j<<endl;
		    
                                      
			cout<<"Actual top-k:"<<num<<endl;
		    */
            /*            		     	
			vector<cell> result = TopKDG(realdata, dominantgraph, k_origin, weight_origin);
	
			for(vector<cell>::iterator iter = result.begin(); iter!=result.end(); iter++)
			{
				cout<<"ID:"<<iter->ID<<"----"<<"Score:"<<iter->data<<endl;
			}
	
			if(BestAns.k!=OUTSIDE_THRESHOLD)
				result = TopKDG(realdata, dominantgraph, (BestAns.k<k_origin?k_origin:BestAns.k+10), BestAns.w);
			cout<<"After Why-not:"<<endl;
			cout<<"Original ranking:"<<k_old<<endl;
			cout<<"New Ranking:"<<BestAns.k<<endl;
	
			for(vector<cell>::iterator iter = result.begin(); iter!=result.end(); iter++)
			{
				cout<<"ID:"<<iter->ID<<"----"<<"Score:"<<iter->data<<endl;
			}
		    cout<<"New weighting:"<<endl; 
            for(FloatVec::size_type ii = 0; ii != BestAns.w.weighting.size();ii++)
            {
                cout<<"Weighting "<<ii+1<<":"<<BestAns.w.weighting[ii]<<endl;
            }
            */
    
			totaltime += float(stop- start) / CLOCKS_PER_SEC;
		    sumquality +=BestAns.score;
		    //cout<<"Avarage bufferresult used:"<<num_bufferused<<endl;
			//cout<<"------------------------------------"<<endl;
		    //cout<<"Sample size:"<<sample_size<<endl;
		    //cout<<"Sample Time:"<<sampletime<<endl;
		    //cout<<"Pruning Time:"<<pruning_time<<endl;
		    //cout<<"Topk time:"<<topk_time/num_of_test<<endl;
		    //cout<<"Assign time"<<assign_time/num_of_test<<endl;
			cout<<"Average Time:"<<(totaltime+sampletime+estimateboundtime)<<endl;
		    //cout<<"Average quality"<<sumquality<<endl;	
			//cout<<"------------------------------------"<<endl;

            koresult[koindex*3+ix] +=(totaltime+sampletime+estimateboundtime);
		}//loop penalty
      }//loop for k_origin
    }//loop for num of test 
    ofstream kofile(korigin.c_str(), ofstream::out);

    for(int i = 0;i<koriginvec.size();i++)
    {
        kofile<<koriginvec[i];
        for(int j = 0;j<3;j++)
        {
            kofile<<' '<<koresult[i*3+j]/repeat;
        }
        kofile<<endl;
    }
    kofile.close();
//////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The code is the same as above, just varying Mising number this time.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
cout<<"Varying Missing Num:"<<endl;

       FloatVec missingnumresult;
       missingnumresult.assign(15, 0);

       int missing_num = 1; 
    for(int num_of_test = 0; num_of_test < repeat; num_of_test++)
	{
       for(missing_num = 1; missing_num<=5; missing_num++)
	   {
           //break;
			int k_origin = default_k_origin;                 //original value of k
            //int k_origin = koriginvec[koindex];
		    //cout<<"Missing Num:"<<missing_num<<endl;
			vector<Tuple> missing_node;
			IntVec missing_ID;
            
		    vector<cell> gen_missing_point = TopKDG(realdata, dominantgraph, 600, weight_origin);
            for(int i = 0;i<missing_num;i++)
            {
		        missing_ID.push_back(gen_missing_point[100*(i+1)].ID);
            }
            
		    for(IntVec::iterator iter=missing_ID.begin();iter!=missing_ID.end();iter++)
		    {
			    missing_node.push_back(allTuples[*iter-1]); 
		    }

			/*
		    for(IntVec::iterator iter = missing_ID.begin();iter!=missing_ID.end();iter++)
            {
               cout<<"Missing ID:"<<*iter<<endl;
               for(int i = 0; i <dimension;i++)
               {
                   cout<<"Value "<<i+1<<":"<<allTuples[*iter-1].data[i]<<endl;
               }
            }
			*/

		    //exit(-1); 

		    //for(IntVec::iterator iter = missing_ID.begin();iter!=missing_ID.end();iter++)
			    //cout<<"Missing ID:"<<*iter<<endl;


		    float quality_of_answer = default_quality;             //quallity of the answer
			float probability_guarantee = default_prob;          //probability to get such an answer  
			vector< FloatVec > incompar_points;
		    incompar_points.reserve(allTuples.size()*missing_num);	
	

		    float sampletime = 0;
		    float estimateboundtime = 0;
		    int sample_size = 0;
			list<Weight> weights;

		    start = clock();
			//-----------------------rank lower and upper bounds---------------------------//
			int lowerbound = Estimatebounds(allTuples,missing_node, position, incompar_points, index);
		    stop = clock();
			//cout<<"Total time:"<<float(stop- start) / CLOCKS_PER_SEC <<endl;
		    //cout<<"What is the lower bound:"<<lowerbound<<endl;
		    //cout<<"What is the num of incomparable:"<<incompar_points.size()<<endl;
		    //exit(-1);
		    estimateboundtime = float(stop-start)/CLOCKS_PER_SEC;
	
		    start_sample=clock();
		    weights = SampleWeights(quality_of_answer,probability_guarantee, weight_origin);
            //weights.clear();
		    //SampleWeightsFromB(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);
		    //SampleWeightsFromP(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);

			stop_sample=clock();
			sample_size = weights.size();
			//cout<<"Sample size:"<<weights.size()<<endl;
			//cout<<"Sample time:"<<float(stop_sample- start_sample) / CLOCKS_PER_SEC <<endl;
		    sampletime = float(stop_sample - start_sample) / CLOCKS_PER_SEC;
		    //exit(-1);


		    FloatVec penalty(3);
		    penalty[0] = 0.9;
		    penalty[1] = 0.5;
		    penalty[2] = 0.1;
		

		    for(FloatVec::size_type ix = 0;ix!=penalty.size();ix++)
		    {    
			    float penalty_k = penalty[ix];          //ipenalty for changing k
			    float penalty_w = 1 - penalty_k;             //penalty for changing w

			    //cout<<"penalty_k:"<<penalty_k<<endl;
		 
			int threshold_k = allTuples.size();
			float totaltime = 0;
			float pruningtime = 0;
            float topk_time = 0;
			float bestquality = INFINITE;
			float worstquality = 0;
			float sumquality = 0;
            clock_t prune_start, prune_stop;
            clock_t topk_start, topk_stop;
		  
			int num_bufferused = 0;
		 
			    int ranklowerb;
			    int rankupperb;	
			    float tmpthreshold = INFINITE;
			float threshold_w = INFINITE;
	
			if(lowerbound>threshold_k)
			{
			    cout<<"The rank lower bound:"<<lowerbound<<endl;
			    cout<<"Can't rank into the threshold!"<<endl;
			    exit(-1);
			}

			int deltaKlowerbound = 0;
			    if(lowerbound - k_origin > deltaKlowerbound)
				    deltaKlowerbound = lowerbound - k_origin;

			    Answer BestAns;
			    list<cell> resultset;
			    deque<Answer> skylineAnswer;
	
			start = clock();
			    float missing_score = Score(missing_node, weight_origin);

			    int k_old, kmin, kprune;
                int k_delta;
			int topk_id = 10;
			kmin = allTuples.size();
			ResultPool resultpool;

                topk_start = clock(); 
		    	k_old = ProgressiveTopK(realdata, dominantgraph, missing_score, weight_origin, resultpool, kmin, topk_id);
                topk_stop = clock();
                topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
			    kmin = k_old;
				BestAns.k = k_old;
                BestAns.result.ID = resultpool.ID;
				BestAns.result.result.resize(1000000);
			    BestAns.result.result.assign(resultpool.result.begin(),resultpool.result.end());
				BestAns.w = weight_origin;
				//BestAns.score = penalty_k*(BestAns.k-k_origin) + penalty_w* BestAns.w.deltaw;
				BestAns.score = penalty_k*(BestAns.k-k_origin)/(k_old-k_origin) + penalty_w* BestAns.w.deltaw/deltaWU;
				skylineAnswer.push_back(BestAns);
		
				tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
				//tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;

			    if(threshold_w > tmpthreshold)
				    threshold_w = tmpthreshold;

             

			    unsigned int i = 0, j = 0, num = 1;
			    int test_count = 0;
			    for(list<Weight>::iterator iter = weights.begin(); iter!=weights.end(); iter++)
			    {
				    Weight w = *iter;
				    if(w.deltaw>threshold_w)
                    {
					    break;
                        //continue;
                    }
                    k_delta = (BestAns.score-penalty_w*w.deltaw/deltaWU)*(k_old-k_origin)/penalty_k;

                    if(k_delta<deltaKlowerbound)
                    {
                        //num++;
                        //continue;
                        break;
                    }
                    if((k_delta+k_origin)<kmin)
                    {
                        kprune = k_delta+k_origin+1;
                        //cout<<"kmin"<<kmin<<endl;
                        //cout<<"Kprune"<<kprune<<endl;
                    }
                    else
                        kprune = kmin;

                    //kprune = kmin;
			        bool prunebyrule2 = false;
			    
				    missing_score = Score(missing_node, w);

                    //cout<<"Missing score:"<<missing_score<<endl;
                    
				    if(Pruningbyrule2(resultpool.result, w, kprune, missing_score))
				    {
					    ++i;
			        	prunebyrule2 = true;
			        }
			        
		
			        if(prunebyrule2)
                    {
			            continue;
                    }
                    bool prunebyrule3 = false;
				    if(beDominate(skylineAnswer, w, missing_score, kprune, topk_id))
				    {
					    ++j;
					    prunebyrule3 = true;
				    }
		
		
			        if(prunebyrule3)
                    {
			            continue;
                    }

			        topk_id++;
                    topk_start = clock();
		 		    int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, kprune-1, topk_id);
		 		    //int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, allTuples.size(), topk_id, false);
                    topk_stop = clock();
                    topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
				    ++num;

                    //cout<<"k:"<<k<<endl;
				    if(k == OUTSIDE_THRESHOLD)		
				    {
					    //continue;
                        Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            skylineAnswer.push_back(tempAnswer);

				    }
	
				    //if(k<kmin)
                     else
    			     {
                        if(k < k_origin)
                            k = k_origin;
                        if(kmin > k)
                            kmin = k;
				        float score = penalty_k*(k-k_origin)/(k_old-k_origin) + penalty_w*w.deltaw/deltaWU;
				        //float score = penalty_k*(kmin-k_origin) + penalty_w*w.deltaw;
                        //cout<<"score:"<<score<<endl;
                        //cout<<"Bestscore:"<<BestAns.score<<endl;

			            Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            tempAnswer.score = score;
			            skylineAnswer.push_back(tempAnswer);
		
				        if(score < BestAns.score)
				        {
					        BestAns = tempAnswer;

					        //float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;
					        float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
			
		        			if(threshold_w > tmpthreshold)
				        		threshold_w = tmpthreshold;
				        }

			        }
			    }

			//stop_query=clock();
			stop=clock();
	
			//int size = weights.size();
			//cout<<"Weight size:"<<size<<endl;
		    
                        		     	
			totaltime += float(stop- start) / CLOCKS_PER_SEC;
		    sumquality +=BestAns.score;
            /*
		    cout<<"Avarage bufferresult used:"<<num_bufferused<<endl;
			cout<<"------------------------------------"<<endl;
		    cout<<"Sample size:"<<sample_size<<endl;
		    cout<<"Sample Time:"<<sampletime<<endl;
		    cout<<"Pruning Time:"<<pruning_time<<endl;
		    //cout<<"Topk time:"<<topk_time/num_of_test<<endl;
		    //cout<<"Assign time"<<assign_time/num_of_test<<endl;
			cout<<"Average Time:"<<(totaltime+sampletime+estimateboundtime)<<endl;
		    cout<<"Average quality"<<sumquality<<endl;	
			cout<<"------------------------------------"<<endl;
            */
            missingnumresult[(missing_num-1)*3+ix] += (totaltime+sampletime+estimateboundtime);
		}//loop penalty
      }//loop for missing num
	}//loop for repeat
    ofstream missingnumfile(missingnum.c_str(), ofstream::out);

    for(int i =0;i<5;i++)
    {
        missingnumfile<<i+1;
        for(int j = 0;j<3;j++)
        {
            missingnumfile<<' '<<missingnumresult[i*3+j]/repeat;
        
        }
        missingnumfile<<endl;
    }
    
    missingnumfile.close();
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The code below is the same as above, it is just convinient for me to vary the missing position of the 
// missing tuples than output the data to its file.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
cout<<"Varying Missing Position:"<<endl;


       IntVec missingposvec = config.rankvary;
	   int rankvarysize = missingposvec.size();
       
       FloatVec missingposresult;
       //missingposresult.reserve(15);
       missingposresult.assign(rankvarysize*3, 0);

       for(int i = 0; i<repeat;i++)
       {
       for(IntVec::size_type missingposindex = 0; missingposindex!=missingposvec.size();missingposindex++)
	   {
			int k_origin = default_k_origin;                 //original value of k
            //int k_origin = koriginvec[koindex];
		    cout<<"Missing position"<<missingposvec[missingposindex]<<endl;

			vector<Tuple> missing_node;
			IntVec missing_ID;
		    //cout<<"Korigin:"<<k_origin<<endl; 
		    vector<cell> gen_missing_point = TopKDG(realdata, dominantgraph, missingposvec[missingposindex], weight_origin);
		    missing_ID.push_back(gen_missing_point.rbegin()->ID);
            
            
		    for(IntVec::iterator iter=missing_ID.begin();iter!=missing_ID.end();iter++)
		    {
			missing_node.push_back(allTuples[*iter-1]); 
		    }
		    //exit(-1); 
		    //for(IntVec::iterator iter = missing_ID.begin();iter!=missing_ID.end();iter++)
			    //cout<<"Missing ID:"<<*iter<<endl;

		    float quality_of_answer = default_quality;             //quallity of the answer
			float probability_guarantee = default_prob;          //probability to get such an answer  
			vector< FloatVec > incompar_points;
		    incompar_points.reserve(allTuples.size());	
	

		    float sampletime = 0;
		    float estimateboundtime = 0;
		    int sample_size = 0;
			list<Weight> weights;

		    start = clock();
			//-----------------------rank lower and upper bounds---------------------------//
			int lowerbound = Estimatebounds(allTuples,missing_node, position, incompar_points, index);
		    stop = clock();
			//cout<<"Total time:"<<float(stop- start) / CLOCKS_PER_SEC <<endl;
		    //cout<<"What is the lower bound:"<<lowerbound<<endl;
		    //cout<<"What is the num of incomparable:"<<incompar_points.size()<<endl;
		    //exit(-1);
		    estimateboundtime = float(stop-start)/CLOCKS_PER_SEC;
	
		    start_sample=clock();
		    weights = SampleWeights(quality_of_answer,probability_guarantee, weight_origin);
            //weights.clear();
		    //SampleWeightsFromB(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);
		    //SampleWeightsFromP(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);
			stop_sample=clock();
			sample_size = weights.size();
			//cout<<"Sample size:"<<weights.size()<<endl;
			//cout<<"Sample time:"<<float(stop_sample- start_sample) / CLOCKS_PER_SEC <<endl;
		    sampletime = float(stop_sample - start_sample) / CLOCKS_PER_SEC;
		    //exit(-1);


		    FloatVec penalty(3);
		    penalty[0] = 0.9;
		    penalty[1] = 0.5;
		    penalty[2] = 0.1;

		    for(FloatVec::size_type ix = 0;ix!=penalty.size();ix++)
		    {    
			    float penalty_k = penalty[ix];          //ipenalty for changing k
			    float penalty_w = 1 - penalty_k;             //penalty for changing w

			    cout<<"penalty_k:"<<penalty_k<<endl;
		 
			int threshold_k = allTuples.size();
			float totaltime = 0;
			float pruningtime = 0;
            float topk_time = 0;
			float bestquality = INFINITE;
			float worstquality = 0;
			float sumquality = 0;
            clock_t prune_start, prune_stop;
            clock_t topk_start, topk_stop;
		  
			int num_bufferused = 0;
		 
			    int ranklowerb;
			    int rankupperb;	
			    float tmpthreshold = INFINITE;
			float threshold_w = INFINITE;
	
			if(lowerbound>threshold_k)
			{
			    cout<<"The rank lower bound:"<<lowerbound<<endl;
			    cout<<"Can't rank into the threshold!"<<endl;
			    exit(-1);
			}

			int deltaKlowerbound = 0;
			    if(lowerbound - k_origin > deltaKlowerbound)
				    deltaKlowerbound = lowerbound - k_origin;

			    Answer BestAns;
			    list<cell> resultset;
			    deque<Answer> skylineAnswer;
	
			start = clock();
			    float missing_score = Score(missing_node, weight_origin);

			    int k_old, kmin, kprune;
                int k_delta;
			int topk_id = 10;
			kmin = allTuples.size();
			ResultPool resultpool;

                topk_start = clock(); 
		    	k_old = ProgressiveTopK(realdata, dominantgraph, missing_score, weight_origin, resultpool, kmin, topk_id);
                topk_stop = clock();
                topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
			    kmin = k_old;
				BestAns.k = k_old;
                BestAns.result.ID = resultpool.ID;
				BestAns.result.result.resize(1000000);
			    BestAns.result.result.assign(resultpool.result.begin(),resultpool.result.end());
				BestAns.w = weight_origin;
				//BestAns.score = penalty_k*(BestAns.k-k_origin) + penalty_w* BestAns.w.deltaw;
				BestAns.score = penalty_k*(BestAns.k-k_origin)/(k_old-k_origin) + penalty_w* BestAns.w.deltaw/deltaWU;
				skylineAnswer.push_back(BestAns);
		
				tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
				//tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;

			    if(threshold_w > tmpthreshold)
				    threshold_w = tmpthreshold;

             

			    unsigned int i = 0, j = 0, num = 1;
			    int test_count = 0;
			    for(list<Weight>::iterator iter = weights.begin(); iter!=weights.end(); iter++)
			    {
				    Weight w = *iter;
				    if(w.deltaw>threshold_w)
                    {
					    break;
                        //continue;
                    }
                    k_delta = (BestAns.score-penalty_w*w.deltaw/deltaWU)*(k_old-k_origin)/penalty_k;

                   
                    if(k_delta<deltaKlowerbound)
                    {
                        //num++;
                        //continue;
                        break;
                    }
                    if((k_delta+k_origin)<kmin)
                    {
                        kprune = k_delta+k_origin+1;
                        //cout<<"kmin"<<kmin<<endl;
                        //cout<<"Kprune"<<kprune<<endl;
                    }
                    else
                        kprune = kmin;

                    //kprune = kmin;
			        bool prunebyrule2 = false;
			    
				    missing_score = Score(missing_node, w);

                    //cout<<"Missing score:"<<missing_score<<endl;
				    if(Pruningbyrule2(resultpool.result, w, kprune, missing_score))
				    {
					    ++i;
			        	prunebyrule2 = true;
			        }
			     
		
			        if(prunebyrule2)
                    {
			            continue;
                    }
                    bool prunebyrule3 = false;
				    if(beDominate(skylineAnswer, w, missing_score, kprune, topk_id))
				    {
					    ++j;
					    prunebyrule3 = true;
				    }
		
		
			        if(prunebyrule3)
                    {
			            continue;
                    }

			        topk_id++;
                    topk_start = clock();
		 		    int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, kprune-1, topk_id);
		 		    //int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, allTuples.size(), topk_id, false);
                    topk_stop = clock();
                    topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
				    ++num;

                    //cout<<"k:"<<k<<endl;
				    if(k == OUTSIDE_THRESHOLD)		
				    {
					    //continue;
                        Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            skylineAnswer.push_back(tempAnswer);

				    }
	
				    //if(k<kmin)
                     else
    			     {
                        if(k < k_origin)
                            k = k_origin;
                        if(kmin > k)
                            kmin = k;
				        float score = penalty_k*(k-k_origin)/(k_old-k_origin) + penalty_w*w.deltaw/deltaWU;
				        //float score = penalty_k*(kmin-k_origin) + penalty_w*w.deltaw;
                        //cout<<"score:"<<score<<endl;
                        //cout<<"Bestscore:"<<BestAns.score<<endl;

			            Answer tempAnswer;
			            tempAnswer.k = k;
			            tempAnswer.w = w;
                        tempAnswer.result.ID = resultpool.ID;
			            tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
			            tempAnswer.score = score;
			            skylineAnswer.push_back(tempAnswer);
		
				        if(score < BestAns.score)
				        {
					        BestAns = tempAnswer;

					        //float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;
					        float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
			
		        			if(threshold_w > tmpthreshold)
				        		threshold_w = tmpthreshold;
				        }

			        }
			    }

			//stop_query=clock();
			stop=clock();
	
			//int size = weights.size();
			//cout<<"Weight size:"<<size<<endl;
		    
            
			cout<<"Pruned by rule 1:"<<(weights.size()+1-i-j-num)<<endl;
			cout<<"Pruned by rule 2:"<<i<<endl;
			cout<<"Pruned by rule 3:"<<j<<endl;
		    
                                      
			cout<<"Actual top-k:"<<num<<endl;
		    
            /*            		     	
			vector<cell> result = TopKDG(realdata, dominantgraph, k_origin, weight_origin);
	
			for(vector<cell>::iterator iter = result.begin(); iter!=result.end(); iter++)
			{
				cout<<"ID:"<<iter->ID<<"----"<<"Score:"<<iter->data<<endl;
			}
	
			if(BestAns.k!=OUTSIDE_THRESHOLD)
				result = TopKDG(realdata, dominantgraph, (BestAns.k<k_origin?k_origin:BestAns.k+10), BestAns.w);
			cout<<"After Why-not:"<<endl;
			cout<<"Original ranking:"<<k_old<<endl;
			cout<<"New Ranking:"<<BestAns.k<<endl;
	
			for(vector<cell>::iterator iter = result.begin(); iter!=result.end(); iter++)
			{
				cout<<"ID:"<<iter->ID<<"----"<<"Score:"<<iter->data<<endl;
			}
		    cout<<"New weighting:"<<endl; 
            for(FloatVec::size_type ii = 0; ii != BestAns.w.weighting.size();ii++)
            {
                cout<<"Weighting "<<ii+1<<":"<<BestAns.w.weighting[ii]<<endl;
            }
            */
    
			totaltime += float(stop- start) / CLOCKS_PER_SEC;
		    sumquality +=BestAns.score;
            /*
		    cout<<"Avarage bufferresult used:"<<num_bufferused<<endl;
			cout<<"------------------------------------"<<endl;
		    cout<<"Sample size:"<<sample_size<<endl;
		    cout<<"Sample Time:"<<sampletime<<endl;
		    cout<<"Pruning Time:"<<pruning_time<<endl;
		    //cout<<"Topk time:"<<topk_time/num_of_test<<endl;
		    //cout<<"Assign time"<<assign_time/num_of_test<<endl;
			cout<<"Average Time:"<<(totaltime+sampletime+estimateboundtime)<<endl;
		    cout<<"Average quality"<<sumquality<<endl;	
			cout<<"------------------------------------"<<endl;
            */
            missingposresult[missingposindex*3+ix] += (totaltime+sampletime+estimateboundtime);
		}//loop penalty
      }//loop for missing position 
    }//loop for repeat times
      ofstream missingposfile(missingpos.c_str(), ofstream::out);

      for(int i =0; i<missingposvec.size();i++)
      {
        missingposfile<<missingposvec[i];
        for(int j = 0;j<3;j++)
        {
            missingposfile<<' '<<missingposresult[i*3+j]/repeat;
        }
        missingposfile<<endl;
      }
      missingposfile.close();

}//condition for testall 
/////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//
//
// The code below is used for test the case study, the code is the same as above.
//
//
//
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////

//	priority_queue<GoodCase> collect;
//	int datasize =  allTuples.size();
//	int testid = 0;
//	while(testid<datasize){
//		
//		vector<float> uniform_sample_score;
//		vector<float> our_sample_score;
//
//	 	testid++;
//		for(int useuniform = 0;useuniform<2;useuniform++){
//
//			if(useuniform)
//				cout<<"Using Naive Sampling Method"<<endl;
//
//			else
//				cout<<"Using Our Sampling Method"<<endl;
//
//		if(!testall)
//		{
//			int k_origin = 3;                 //original value of k
//            int missing_num = 1;
//		    //cout<<"Missing Num:"<<missing_num<<endl;
//			vector<Tuple> missing_node;
//			IntVec missing_ID;
//            //missing_ID.push_back(2);
//            //missing_ID.push_back(2564); //O'neal in all player
//            //missing_ID.push_back(2876);
//            //missing_ID.push_back(263);
//		    //missing_ID.push_back(1766);    //jordan in all player
//		    //missing_ID.push_back(428);      //Kobe in scorer
//		    //missing_ID.push_back(1706);      //Magic Johnson in scorer
//            //missing_ID.push_back(3532);     
//            //missing_ID.push_back(3634);     
//		    //missing_ID.push_back(1605);     //Iverson in scorer
//		    //missing_ID.push_back(389);// yao ming in Center
//            //missing_ID.push_back(570);// wilt chamberlain in scorer
//            //missing_ID.push_back(263);//larry bird in scorer
//        
//            //missing_ID.push_back(666);//Iverson in Guard
//            //missing_ID.push_back(176);//Kobe in Guard
//
//		    //vector<cell> gen_missing_point = TopKDG(realdata, dominantgraph, 10*k_origin+1, weight_origin);
//            //missing_ID.push_back(gen_missing_point.rbegin()->ID);
//		    //vector<cell> gen_missing_point = TopKDG(realdata, dominantgraph, 6000, weight_origin);
//		    //missing_ID.push_back(gen_missing_point[50].ID);
//		    //missing_ID.push_back(gen_missing_point[100].ID);
//		    //missing_ID.push_back(gen_missing_point[200].ID);
//		    //missing_ID.push_back(gen_missing_point[300].ID);
//		    //missing_ID.push_back(gen_missing_point[400].ID);
//		    //missing_ID.push_back(gen_missing_point[500].ID);
//			
//			missing_ID.push_back(testid);
//
//		    for(IntVec::iterator iter=missing_ID.begin();iter!=missing_ID.end();iter++)
//		    {
//			    missing_node.push_back(allTuples[*iter-1]); 
//		    }
//		    //exit(-1); 
//		    for(IntVec::iterator iter = missing_ID.begin();iter!=missing_ID.end();iter++)
//            {
//               cout<<"Missing ID:"<<*iter<<endl;
//               for(int i = 0; i <dimension;i++)
//               {
//                   cout<<"Value "<<i+1<<":"<<allTuples[*iter-1].data[i]<<endl;
//               }
//            }
//
//            //exit(-1);
//		    float quality_of_answer = default_quality;             //quallity of the answer
//		    //float quality_of_answer = 0.001;             //quallity of the ainswer
//			float probability_guarantee = default_prob;          //probability to get such an answer  
//			//float probability_guarantee = 0.9;          //probability to get such an answer  
//
//			vector< FloatVec > incompar_points;
//		    incompar_points.reserve(allTuples.size()*missing_num);	
//	
//
//		    float sampletime = 0;
//		    float estimateboundtime = 0;
//		    int sample_size = 0;
//
//		    start = clock();
//			//-----------------------rank lower and upper bounds---------------------------//
//			int lowerbound = Estimatebounds(allTuples,missing_node, position, incompar_points, index);
//		    stop = clock();
//			cout<<"Total time:"<<float(stop- start) / CLOCKS_PER_SEC <<endl;
//		    cout<<"What is the lower bound:"<<lowerbound<<endl;
//		    //cout<<"What is the num of incomparable:"<<incompar_points.size()<<endl;
//		    //exit(-1);
//		    estimateboundtime = float(stop-start)/CLOCKS_PER_SEC;
//	
//		    FloatVec penalty(3);
//		    penalty[0] = 0.9;
//		    penalty[1] = 0.5;
//		    penalty[2] = 0.1;
//            
//            FloatVec deltakvec, deltawvec, scorevec;
//            deltakvec.assign(3, 0);
//            deltawvec.assign(3, 0);
//            scorevec.assign(3 ,0);
//            int repeat = 1;
//            FloatVec timevec;
//            timevec.assign(3, 0);
//
//            float sampletime_count = 0;
//            float estimate_time_cout = 0;
//        	for(int i = 0; i<repeat;i++)
//        	{
//
//				list<Weight> weights;
//				//list<Weight> uniformweights;
//		    	start_sample=clock();
//
//				if(incompar_points.size()==0)
//					break;
//
//				//if(useuniform)
//		    		weights = SampleWeights(quality_of_answer,probability_guarantee, weight_origin);
//				//else
//		    		//SampleWeightsFromB(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);
//		    	//SampleWeightsFromQP( ep, weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);
//		    	//SampleWeightsFromP(weights, quality_of_answer,probability_guarantee, weight_origin, incompar_points);
//				stop_sample=clock();
//				sample_size = weights.size();
//				cout<<"Sample size:"<<weights.size()<<endl;
//				cout<<"Sample time:"<<float(stop_sample- start_sample) / CLOCKS_PER_SEC <<endl;
//		    	sampletime = float(stop_sample - start_sample) / CLOCKS_PER_SEC;
//            	sampletime_count += sampletime;
//		    	//exit(-1);
//
//		    	for(FloatVec::size_type ix = 0;ix!=penalty.size();ix++)
//		    	{    
//			    	float penalty_k = penalty[ix];          //ipenalty for changing k
//			    	float penalty_w = 1 - penalty_k;             //penalty for changing w
//
//					//penalty_k = 0.1;
//			    	cout<<"penalty_k:"<<penalty_k<<endl;
//		 
//					int threshold_k = allTuples.size();
//					float totaltime = 0;
//					float pruningtime = 0;
//            		float topk_time = 0;
//					float bestquality = INFINITE;
//					float worstquality = 0;
//					float sumquality = 0;
//            		clock_t prune_start, prune_stop;
//            		clock_t topk_start, topk_stop;
//		  
//					int num_bufferused = 0;
//		 
//			    	int ranklowerb;
//			    	int rankupperb;	
//			    	float tmpthreshold = INFINITE;
//					float threshold_w = INFINITE;
//	
//					if(lowerbound>threshold_k)
//					{
//			    		cout<<"The rank lower bound:"<<lowerbound<<endl;
//			    		cout<<"Can't rank into the threshold!"<<endl;
//			    		exit(-1);
//					}
//
//					int deltaKlowerbound = 0;
//			    	if(lowerbound - k_origin > deltaKlowerbound)
//				    	deltaKlowerbound = lowerbound - k_origin;
//
//			    	Answer BestAns;
//			    	list<cell> resultset;
//			    	deque<Answer> skylineAnswer;
//	
//					start = clock();
//			    	float missing_score = Score(missing_node, weight_origin);
//
//			    	int k_old, kmin, kprune;
//                	int k_delta;
//					int topk_id = 1;
//					kmin = allTuples.size();
//					ResultPool resultpool;
//
//                	topk_start = clock(); 
//		    		k_old = ProgressiveTopK(realdata, dominantgraph, missing_score, weight_origin, resultpool, kmin, topk_id);
//                	topk_stop = clock();
//                	topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
//			    	kmin = k_old;
//					BestAns.k = k_old;
//                	BestAns.result.ID = resultpool.ID;
//					BestAns.result.result.resize(1000000);
//			    	BestAns.result.result.assign(resultpool.result.begin(),resultpool.result.end());
//					//for(vector<FloatVec>::size_type ii = 0;ii!=resultpool.result.size();ii++)
//					//{
//					//	cout<<"ii:"<<ii<<endl;
//					//	BestAns.result.result[ii] = resultpool.result[ii];
//					//}
//					//cout<<"Get here?"<<endl;
//					BestAns.w = weight_origin;
//					//BestAns.score = penalty_k*(BestAns.k-k_origin) + penalty_w* BestAns.w.deltaw;
//					BestAns.score = penalty_k*(BestAns.k-k_origin)/(k_old-k_origin) + penalty_w* BestAns.w.deltaw/deltaWU;
//					skylineAnswer.push_back(BestAns);
//		
//					tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
//					//tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;
//
//			    	if(threshold_w > tmpthreshold)
//				    	threshold_w = tmpthreshold;
//
//             
//
//			    	unsigned int i = 0, j = 0, num = 1;
//			    	int test_count = 0;
//			    	for(list<Weight>::iterator iter = weights.begin(); iter!=weights.end(); iter++)
//			    	{
//				    	Weight w = *iter;
//				    	if(w.deltaw>threshold_w)
//                    	{
//					    	//break;
//                        	//continue;
//                    	}
//                    	k_delta = (BestAns.score-penalty_w*w.deltaw/deltaWU)*(k_old-k_origin)/penalty_k;
//
//                    	if(k_delta<deltaKlowerbound)
//                    	{
//                        	//num++;
//                        	//continue;
//                        	break;
//                    	}
//                    	if((k_delta+k_origin)<kmin)
//                    	{
//                        	kprune = k_delta+k_origin+1;
//                        	//cout<<"kmin"<<kmin<<endl;
//                        	//cout<<"Kprune"<<kprune<<endl;
//                    	}
//                    	else
//                        	kprune = kmin;
//
//                    	//kprune = kmin;
//			        	bool prunebyrule2 = false;
//			    	
//				    	missing_score = Score(missing_node, w);
//
//                    	//cout<<"Missing score:"<<missing_score<<endl;
//					
//						
//				    	if(Pruningbyrule2(resultpool.result, w, kprune, missing_score))
//				    	{
//					    	++i;
//			        		prunebyrule2 = true;
//			        	}
//			     
//		
//			        	if(prunebyrule2)
//                    	{
//			            	continue;
//                    	}
//                    	bool prunebyrule3 = false;
//				    	if(beDominate(skylineAnswer, w, missing_score, kprune, topk_id))
//				    	{
//					    	++j;
//					    	prunebyrule3 = true;
//				    	}
//						
//		
//			        	if(prunebyrule3)
//                    	{
//			            	continue;
//                    	}
//						
//			        	topk_id++;
//                    	topk_start = clock();
//		 		    	int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, kprune-1, topk_id);
//		 		    	//int k = ProgressiveTopK(realdata, dominantgraph, missing_score, w, resultpool, allTuples.size(), topk_id, false);
//                    	topk_stop = clock();
//                    	topk_time += (float)(topk_stop-topk_start)/CLOCKS_PER_SEC; 
//				    	++num;
//
//                    	//cout<<"k:"<<k<<endl;
//				    	if(k == OUTSIDE_THRESHOLD)		
//				    	{
//                        	//continue;
//			            	Answer tempAnswer;
//			            	tempAnswer.k = k;
//			            	tempAnswer.w = w;
//                        	tempAnswer.result.ID = resultpool.ID;
//			            	tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
//			            	skylineAnswer.push_back(tempAnswer);
//				    	}
//	
//				    	//if(k<kmin)
//                     	else
//    			     	{
//                        	if(k < k_origin)
//                            	k = k_origin;
//                        	if(kmin > k)
//                            	kmin = k;
//				        	float score = penalty_k*(k-k_origin)/(k_old-k_origin) + penalty_w*w.deltaw/deltaWU;
//				        	//float score = penalty_k*(kmin-k_origin) + penalty_w*w.deltaw;
//                        	//cout<<"score:"<<score<<endl;
//                        	//cout<<"Bestscore:"<<BestAns.score<<endl;
//
//			            	Answer tempAnswer;
//			            	tempAnswer.k = k;
//			            	tempAnswer.w = w;
//                        	tempAnswer.result.ID = resultpool.ID;
//			            	tempAnswer.result.result.assign(resultpool.result.begin(),resultpool.result.end());
//			            	tempAnswer.score = score;
//			            	skylineAnswer.push_back(tempAnswer);
//
//							//cout<<"what happen 1?"<<endl;
//		
//				        	if(score < BestAns.score)
//				        	{
//								//cout<<"what happen 2?"<<endl;
//					        	BestAns = tempAnswer;
//
//					        	//float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound)/penalty_w;
//					        	float tmpthreshold = (BestAns.score - penalty_k*deltaKlowerbound/(k_old-k_origin))*deltaWU/penalty_w;
//			
//		        				if(threshold_w > tmpthreshold)
//				        			threshold_w = tmpthreshold;
//				        	}
//
//			        	}
//			    	}
//
//				//stop_query=clock();
//				stop=clock();
//				//int size = weights.size();
//				//cout<<"Weight size:"<<size<<endl;
//		    
//            	/*            
//				cout<<"Pruned by rule 1:"<<(weights.size()+1-i-j-num)<<endl;
//				cout<<"Pruned by rule 2:"<<i<<endl;
//				cout<<"Pruned by rule 3:"<<j<<endl;
//		    
//                                      
//				cout<<"Actual top-k:"<<num<<endl;
//		    
//                          		     	
//				vector<cell> result = TopKDG(realdata, dominantgraph, k_origin, weight_origin);
//	
//            
//				for(vector<cell>::iterator iter = result.begin(); iter!=result.end(); iter++)
//				{
//					cout<<"ID:"<<iter->ID<<"----"<<"Score:"<<iter->data<<endl;
//				}
//			
//				if(BestAns.k!=OUTSIDE_THRESHOLD)
//					result = TopKDG(realdata, dominantgraph, (BestAns.k<k_origin?k_origin:BestAns.k+10), BestAns.w);
//				cout<<"After Why-not:"<<endl;
//				cout<<"Original ranking:"<<k_old<<endl;
//				cout<<"New Ranking:"<<BestAns.k<<endl;
//	
//				for(vector<cell>::iterator iter = result.begin(); iter!=result.end(); iter++)
//				{
//					cout<<"ID:"<<iter->ID<<"----"<<"Score:"<<iter->data<<endl;
//				}
//		    	cout<<"New weighting:"<<endl; 
//            	for(FloatVec::size_type ii = 0; ii != BestAns.w.weighting.size();ii++)
//            	{
//                	cout<<"Weighting "<<ii+1<<":"<<BestAns.w.weighting[ii]<<endl;
//            	}
//            	*/
//            
//				totaltime += float(stop- start) / CLOCKS_PER_SEC;
//		    	sumquality +=BestAns.score;
//            
//				//cout<<"Avarage bufferresult used:"<<num_bufferused/num_of_test<<endl;
//				//cout<<"------------------------------------"<<endl;
//		    	//cout<<"Sample size:"<<sample_size<<endl;
//		    	//cout<<"Sample Time:"<<sampletime<<endl;
//		    	//cout<<"Pruning Time:"<<pruningtime<<endl;
//		    	//cout<<"Topk time:"<<topk_time<<endl;
//		    	//cout<<"Assign time"<<assign_time/num_of_test<<endl;
//				//cout<<"Average Time:"<<(totaltime+sampletime+estimateboundtime)<<endl;
//		    	//cout<<"Average quality"<<sumquality<<endl;	
//			
//            	//cout<<"------------------------------------"<<endl;
//            	timevec[ix] += totaltime+sampletime+estimateboundtime;
//            
//            	deltakvec[ix] += BestAns.k - k_origin;
//            	deltawvec[ix] += BestAns.w.deltaw;
//            	scorevec[ix] += BestAns.score;
//
//			}//loop penalty
//    	}//loop for repeat;
//         
//        cout<<"----------------------------------"<<endl;
//        for(int i = 0; i<3; i++)
//        {
//            //cout<<"Estimate Bound time:"<<estimateboundtime<<endl;
//            //cout<<"Sample Time:"<<sampletime_count/repeat<<endl;
//            cout<<"Average Time:"<<timevec[i]/repeat<<endl;
//            //cout<<"Query Time:"<<timevec[i]/repeat-sampletime_count/repeat-estimateboundtime<<endl;
//            cout<<"Penalty value k:"<<penalty[i]<<endl;
//            cout<<"Delta K:"<<deltakvec[i]/repeat<<endl;
//            cout<<"Delta W:"<<deltawvec[i]/repeat<<endl;
//            cout<<"Penalty:"<<scorevec[i]/repeat<<endl;
//
//			if(useuniform)
//				uniform_sample_score.push_back(scorevec[i]/repeat);				
//			else
//				our_sample_score.push_back(scorevec[i]/repeat);
//
//            cout<<"--------------------------------"<<endl;
//        }
//        
//		}// if test all
//
//		
//
//		}//useuniform
//		
//			if(our_sample_score[0]<uniform_sample_score[0]&&our_sample_score[0]!=0)
//			{
//				float winpercentage = uniform_sample_score[0]/our_sample_score[0];
//
//				if(collect.size()<10)
//				{
//					GoodCase temp;
//					temp.ID = testid;
//					temp.winpercentage = winpercentage;
//					collect.push(temp);
//				}
//				else
//				{
//					if(collect.top().winpercentage<winpercentage)
//					{
//						GoodCase temp;
//						temp.ID = testid;
//						temp.winpercentage = winpercentage;
//						collect.pop();
//						collect.push(temp);
//					}
//				}
//			}
//
//	}//while collection not full
//
//
//	cout<<"The set of id I want to get:"<<endl;
//
//	while(!collect.empty())
//	{
//		cout<<"ID:"<<collect.top().ID<<"\t"<<"Win Percentage:"<<collect.top().winpercentage<<endl;
//		collect.pop();
//	}
    //do_octave_atexit();
    
    //engClose(ep);
    return 0;
}
