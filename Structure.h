#ifndef _STRUCTURE_H
#define _STRUCTURE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <queue>
#include <deque>
#include <limits>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <cstdlib>
#include <time.h>
#include <algorithm>

#include <cassert>

using namespace std;

#define INFINITE numeric_limits<float>::max()
#define random(x) (rand()%100)

//typedef vector<Fraction> FloatVec;

#define IntVec	vector<int>
#define FloatVec    vector<float>
#define DoubleVec   vector<double>




struct Point
{
	int ID;
	FloatVec coordinate;
};


struct Edge
{
	int node1;
	int node2;
};


struct bounds
{
	int lowerbound;
	int upperbound;
};

struct Tuple
{
	int ID;
	FloatVec data;
//	bool beskyline;

//	Tuple(){beskyline = false;}
};

struct Weight
{
    DoubleVec weighting;
	float deltaw;
};

typedef vector < Tuple > Tuples;

struct cell
{
	int ID;
	float data;

	friend bool operator <(cell c1,cell c2)
	{
		return c1.data<c2.data;
	}
};

struct Indexcell
{
    int ID;
    float key;
    
    friend bool operator <(Indexcell c1, Indexcell c2)
    {
        return c1.key>c2.key;
    }
    
};

struct ResultPool
{
    int ID;
    vector<FloatVec> result;
};

typedef list<Indexcell> Index;

typedef vector< Point > Points;
typedef vector<string> Attribute;
typedef list <Edge> Edges;

bool compfn(cell a, cell b);

void QuickSort(vector<cell> &  arr, int left, int right) ;

void ReadCSVbeta(char * filename, Tuples& allTuples, Tuple& attributes);

void ReadData(char * filename, Tuples& allTuples);
void WriteData(char * filename, Tuples& allTuples);

void ReadCSV(char * filename, Tuples& allTuples);

void WriteCSV(char *filename, Tuples allTuples);


class Fraction
{
private:
   int numerator;
   int denominator; 

friend ostream &operator<<(ostream &output,const Fraction &fraction)
{
   if(fraction.numerator%fraction.denominator!=0)
    output<<fraction.numerator<<"/"<<fraction.denominator;
   else
    output<<fraction.numerator/fraction.denominator;
   return output;
}

    friend istream &operator>>(istream &input,Fraction &fraction) 
    { 
        input>>fraction.numerator>>fraction.denominator; 
        return input; 
    }

public:
	   Fraction(int numerator=0,int denominator=1);
	   Fraction(char *str);
	   const Fraction operator+();
	   const Fraction operator+(const Fraction &);
	   const Fraction operator-();
	   const Fraction operator-(const Fraction &);
	   const Fraction operator*(const Fraction &);
	   const Fraction operator/(const Fraction &);
	   bool operator<(const Fraction &); 
	      bool operator>(const Fraction &); 
	      bool operator==(const Fraction &); 
	      bool operator!=(const Fraction &);
	      Fraction Easiest(Fraction &); 
	      int gcd(int a,int b); 
};

template <class T>
class reservable_priority_queue: public std::priority_queue<T>
{
public:
    typedef typename std::priority_queue<T>::size_type size_type;
    reservable_priority_queue(size_type capacity = 0) { reserve(capacity); };
    void reserve(size_type capacity) { this->c.reserve(capacity); } 
    size_type capacity() const { return this->c.capacity(); } 
};

#endif
