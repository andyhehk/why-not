#include "config.h"

//test if the property is set or not
bool isPropertySet(map<string,string> configmap, string properties) {
	vector<string> values;

	istringstream s(properties);
	while(!s.eof()) {
		string field;
		getline(s,field,'|');
		values.push_back(field);
	}

	for(vector<string>::iterator iter = values.begin(); iter != values.end(); iter++) {
		if(configmap.find(*iter)==configmap.end()) {
			return false;
		}
	}
	return true;
}

//Read the parameter from the configuration file
bool readConfig(string filename, MyConfigType & config, map<string,string> & configmap) {
	ifstream input(filename.c_str(), ifstream::in);

	if(!input.is_open()) {
		cout<<"The configuration file is missing!!!"<<endl;
		return false;
	}

	while(!input.eof()) {
		string line;
		string propertyname;
		string propertyvalue;

		getline(input, line);
		if(line.empty()||line.at(0)=='#') 
			continue;

		istringstream s(line);
		getline(s, propertyname, '=');
		getline(s, propertyvalue, '=');

		remove(propertyname.begin(),propertyname.end(),' ');
		remove(propertyvalue.begin(),propertyvalue.end(),' ');

		configmap.insert(pair<string,string>(propertyname,propertyvalue));
	}

	//Check if the default parameters are all set

	string dflt_properties = "nm|korg|rankorg|t|prob|repeat";
	if(!isPropertySet(configmap,dflt_properties)) {
		cout<<"At least one of the default parameters: korg, rankorg, t, prob, nm, repeat are not set!!!"<<endl;
		return false;
	}

	for(map<string,string>::iterator iter = configmap.begin();iter!=configmap.end();iter++) {
		string name = (*iter).first;
		string value = (*iter).second;

		istringstream s(value);
		
		if(name == "pmk") {
			config.pmk = atof(value.c_str());
		}

		else if(name == "pma") {
			config.pma = atof(value.c_str());
		}

		else if(name == "nm") {
			config.nm = atof(value.c_str());
		}

		else if(name=="korg") {
			config.korg = atoi(value.c_str());
		}
				
		else if(name=="rankorg") {
			config.rankorg = atoi(value.c_str());
		}

	 	else if(name=="t") {
			config.t = atof(value.c_str());
		}

		else if(name=="prob") {
			config.prob = atof(value.c_str());
		}
		
		else if(name=="repeat") {
			config.repeat = atoi(value.c_str());
		}

		else if(name=="datatype") {
			config.datatype = value;
		}


		else if(name=="kvary") {
			while(!s.eof()) {
				string field;
				getline(s,field,',');
				config.kvary.push_back(atoi(field.c_str()));
			}
		}

		else if(name=="rankvary") {
			while(!s.eof()) {
				string field;
				getline(s,field,',');
				config.rankvary.push_back(atoi(field.c_str()));
			}
		}
			
		else if(name=="tvary") {
			while(!s.eof()) {
				string field;
				getline(s,field,',');
				config.tvary.push_back(atof(field.c_str()));
			}
		}

		else if(name=="probvary") {
			while(!s.eof()) {
				string field;
				getline(s,field,',');
				config.probvary.push_back(atof(field.c_str()));
			}
		}

		else if(name=="EARLYSTOP") {
			int indicator = atoi(value.c_str());

			if(indicator==1)
				config.EARLYSTOP = true;
			else if(indicator==0)
				config.EARLYSTOP = false; 
			else{
				cout<<"We expect 1 or 0 for "<<name<<endl;
				exit(-1);
			}
		}
		
		else if(name=="USEBUFFERING") {
			int indicator = atoi(value.c_str());
		
			if(indicator==1)
				config.USEBUFFERING = true;
			else if(indicator==0)
				config.USEBUFFERING = false; 
			else{
				cout<<"We expect 1 or 0 for "<<name<<endl;
				exit(-1);
			}
		}

		else if(name=="USENEWRANKFUNCTION") {
			int indicator = atoi(value.c_str());
			if(indicator==1) {
				config.USENEWRANKFUNCTION = true;
			}
			else if(indicator==0)
				config.USENEWRANKFUNCTION = false; 
			else{
				cout<<"We expect 1 or 0 for "<<name<<endl;
				exit(-1);
			}
		}
	}

	return true;
}

//Display the configure parameter 
void listConfig(map<string,string> & configmap) {
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout<<"The configuration fot the why-not program::"<<endl;
	for(map<string,string>::iterator iter = configmap.begin(); iter != configmap.end(); iter++) {
		cout<<(*iter).first<<"="<<(*iter).second<<endl;
	}
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
}
