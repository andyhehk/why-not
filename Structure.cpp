#include "Structure.h"

bool compfn(cell a, cell b)
{
	if(a.data>b.data||a.data==b.data)
		return true;
	else
		return false;
	
}

void ReadData(char * filename, Tuples& allTuples)
{
	ifstream file(filename, ifstream::in);
	string line;
	size_t found;
	allTuples.clear();
	
	if(!file.is_open())
	{
		printf("Can not open file!\n");
		return;
	}

	while(getline(file, line))
	{
		istringstream s(line);
		bool isID = true;

		Tuple record;
		string field;

		while(getline(s, field, ' '))
		{
		 	if(isID)
		 	{
				record.ID = atoi(field.c_str());
				isID = false;
		 	}
			else
			{
				if(!strcmp(field.c_str(), "\r"))
					continue;
				record.data.push_back(atof(field.c_str()));
				//record.data.push_back(ato(field.c_str()));
			}
		}

		allTuples.push_back(record);
		
	}
	file.close();
	return;

}

void WriteData(char * filename, Tuples& allTuples)
{
	ofstream file(filename, ofstream::out);
	string line;

	if(!file.is_open())
	{
		printf("Can not open file!\n");
		return;
	}

	for(Tuples::iterator iter = allTuples.begin();iter!=allTuples.end();iter++)
	{
		Tuple temp = *iter;
		FloatVec data = temp.data;
		file<<temp.ID;
		for(FloatVec::iterator iter1 = data.begin();iter1!=data.end();iter1++)
		{
			file<<' '<<*iter1;
		}
		file<<endl;
	}
	file.close();
	return;
}


void ReadCSV(char * filename, Tuples& allTuples)
{	
	ifstream file(filename, ifstream::in);
	string line;

	if(!file.is_open())
	{
		printf("Can not open file!\n");
		return;
	}

	allTuples.clear();
	int id = 0;
	while(getline( file, line))
	{
				
		istringstream s(line);
		size_t found;
		
		Tuple record;
		string field;
		bool isID = true;
		record.ID = ++id;
		while(getline(s, field, ','))
		{

			//if(isID){
				//record.ID = atoi(field.c_str());
				//isID = false;
			//}
			//else
				record.data.push_back(atof(field.c_str()));
		}
		allTuples.push_back(record);		
		
	}

	file.close();	
	return;
}


void WriteCSV(char *filename, Tuples allTuples)
{
	ofstream file(filename, ofstream::out);

	if(!file.is_open())
	{
		printf("Can not open file!\n");
		return;
	}

	for(Tuples::iterator iter = allTuples.begin(); iter!=allTuples.end(); iter++)
	{
		file<<iter->ID;
		FloatVec data = iter->data;

		for(FloatVec::iterator iter1 = data.begin(); iter1!=data.end(); iter1++)
			file<<","<<*iter1;

		file<<endl;
		
	}

	file.close();
	return;
	
}

Fraction::Fraction(int num,int den)
{
	numerator=num;
	denominator=den; 
}

Fraction::Fraction(char *str)
{
	sscanf(str,"%d/%d",&numerator,&denominator);
}

const Fraction Fraction::operator+()
{
	Fraction Value;
	Value.numerator=numerator;
	Value.denominator=denominator;
	return Value;
}


const Fraction Fraction::operator-()
{
	Fraction Value;
	Value.numerator=-numerator;
	Value.denominator=denominator;
	return Value;
}


const Fraction Fraction::operator+(const Fraction &right)
{
	Fraction resultValue;
	resultValue.numerator=numerator*right.denominator+denominator*right.numerator;
	resultValue.denominator=denominator*right.denominator;
	return Easiest(resultValue); 
}

const Fraction Fraction::operator-(const Fraction &right)
{
	Fraction resultValue;
	resultValue.numerator=numerator*right.denominator-denominator*right.numerator;
	resultValue.denominator=denominator*right.denominator;
	return Easiest(resultValue);   
}

const Fraction Fraction::operator*(const Fraction &right)
{
	Fraction resultValue;
	resultValue.numerator=numerator*right.numerator;
	resultValue.denominator=denominator*right.denominator;
	return Easiest(resultValue);
}

const Fraction Fraction::operator/(const Fraction &right)
{
	Fraction resultValue;
	resultValue.numerator=numerator*right.denominator;
	resultValue.denominator=denominator*right.numerator;
	return Easiest(resultValue);
}

bool Fraction::operator<(const Fraction &right)
{
	int leftValue,rightValue;
	leftValue=numerator*right.denominator;
	rightValue=right.numerator*denominator;
	return leftValue<rightValue;
}

bool Fraction::operator>(const Fraction &right)
{
	int leftValue,rightValue;
	leftValue=numerator*right.denominator;
	rightValue=right.numerator*denominator;
	return leftValue>rightValue;
}

bool Fraction::operator==(const Fraction &right)
{
	int leftValue,rightValue;
	leftValue=numerator*right.denominator;
	rightValue=right.numerator*denominator;
	return leftValue==rightValue;
}

bool Fraction::operator!=(const Fraction &right)
{
	int leftValue,rightValue;
	leftValue=numerator*right.denominator;
	rightValue=right.numerator*denominator;
	return leftValue!=rightValue;
}

Fraction Fraction::Easiest(Fraction &fraction) 
{ 
   int getGcd=gcd(fraction.numerator,fraction.denominator);
   Fraction value(fraction.numerator/getGcd,fraction.denominator/getGcd);
   return value;
}

int Fraction::gcd(int a,int b)
{
	int temp;
	while(b)
	{
	   temp=b;
	   b=a%b;
	   a=temp;
	}
	return a;
}

