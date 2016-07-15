#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <map>
#include <cassert>

using namespace std;

//The structure to store the configuration
struct MyConfigType {
	float pmk;
	float pma;
	float nm;
	int korg;
	int rankorg;
	float t;
	float prob;
	int repeat;
	string datatype; //i | c | ac

	bool EARLYSTOP;
	bool USEBUFFERING;
	bool USENEWRANKFUNCTION;

	vector<int> kvary;
	vector<int> rankvary;
	vector<float> tvary;
	vector<float> probvary;
};

bool isPropertySet(map<string,string> configmap, string properties);
bool readConfig(string filename, MyConfigType & config, map<string,string> & configmap); 
void listConfig(map<string,string> & configmap);
#endif
