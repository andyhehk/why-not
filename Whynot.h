#ifndef _WHY_NOT_TOP_K_H
#define _WHY_NOT_TOP_K_H

#include "DominantGraph.h"
#include "config.h"
//#include "engine.h"
#define  BUFSIZE 256
//#include <boost/random/linear_congruential.hpp>
//#include <boost/random/uniform_int.hpp>
//#include <boost/random/uniform_real.hpp>
//#include <boost/random/variate_generator.hpp>
//#include <boost/random.hpp>


// Sun CC doesn't handle boost::iterator_adaptor yet
//#if !defined(__SUNPRO_CC) || (__SUNPRO_CC > 0x530)
//#include <boost/generator_iterator.hpp>
//#endif

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std {
    using ::time;
}
#endif



struct Answer
{
	int k;
	Weight w;
    ResultPool result;
	float score;

	Answer(){
		score = 1;
		k = OUTSIDE_THRESHOLD;
	}
    friend bool operator <(Answer a, Answer b)
    {
        return a.w.deltaw<b.w.deltaw;
    }
};

struct SubEquation
{
    float coeficient;
    int varid;
};

struct Var
{
    int ID;
    
    float data;
    
    friend bool operator <(Var a, Var b)
    {
        return a.data<b.data;
    }
    
};

struct interv
{
    int begin;
    int end;
};


list<Weight> SampleWeights(float quality_of_answer, float probability_guarantee, Weight w_origin);
void SampleWeightsFromB(list<Weight>& sampleweight, double quality_of_answer, double probability_guarantee, Weight w_origin, const vector< FloatVec >& incompar_points, long sizeneed = 0);
bool GaussianElim(vector< DoubleVec >& matrix, DoubleVec &ans, int row, int column, double& sum2);
bool compWeight(Weight a, Weight b);

bool Pruningbyrule2(const vector<FloatVec>& result, Weight w, int kmin, float missingscore);
int Cmp(float a, float b);

float Score(const vector<Tuple>& missing_node, Weight & w);

//void updateSkyline(list<Answer> & skylineAnswer, int K, Weight W, ResultPool resultpool, float score);

//bool beDominate(list<Answer>& skylineAnswer, Weight W, float missingscore, int kmin, int ID);
int Estimatebounds(const Tuples& allTuples, const vector<Tuple>& missing_node, IntVec position, vector< FloatVec >& incompar_points, const Index & index);

#endif
