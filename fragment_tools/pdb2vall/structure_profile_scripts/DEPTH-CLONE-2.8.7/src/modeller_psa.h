# include <iostream>
# include <fstream>
# include <math.h>
# include <string>
# include <map>
# include <cmath>
# include "easystring.h"
using namespace std;


class residue_asa{ // class object to holder modeller accessibility information
	private:
	float _all_sum, _all_per, _nonP_sum, _nonP_per, _polar_sum, _polar_per, _side_sum, _side_per, _main_sum, _main_per;
	string _resname, _resSeq;

	public:
	residue_asa(){}
	~residue_asa(){}

	// set and get functions
	void set_resSeq(string s){
		_resSeq = s;
	} // end set_resSeq
	void set_resname(string s){
		_resname = s;
	} // end set_resname
	void set_all_sum(float s){
		_all_sum = s;
	} // end set_		
	void set_all_per(float s){
		_all_per = s;
	} // end set_		
	void set_nonP_sum(float s){
		_nonP_sum = s;
	} // end set_		
	void set_nonP_per(float s){
		_nonP_per = s;
	} // end set_		
	void set_polar_sum(float s){
		_polar_sum = s;
	} // end set_		
	void set_polar_per(float s){
		_polar_per = s;
	} // end set_		
	void set_side_sum(float s){
		_side_sum = s;
	} // end set_		
	void set_side_per(float s){
		_side_per = s;
	} // end set_		
	void set_main_sum(float s){
		_main_sum = s;
	} // end set_		
	void set_main_per(float s){
		_main_per = s;
	} // end set_

	float all_sum(){
		return _all_sum;
	} // end 
	float all_per(){
		return _all_per;
	} // end 
	float nonP_sum(){
		return _nonP_sum;
	} // end 
	float nonP_per(){
		return _nonP_per;
	} // end 
	float polar_sum(){
		return _polar_sum;
	} // end 
	float polar_per(){
		return _polar_per;
	} // end 
	float side_sum(){
		return _side_sum;
	} // end 
	float side_per(){
		return _side_per;
	} // end 
	float main_sum(){
		return _main_sum;
	} // end 
	float main_per(){
		return _main_per;
	} // end
	string resname(){
		return _resname;
	} // end resname
	string resSeq(){
		return _resSeq;
	} // end resSeq

}; // end class

class ASA{
	private:
	vector <residue_asa> data;
	string line, resSeq_tmp;
	int d;
	map <string, int> finder;

	public:
	ASA(string fname){
		d = -1;
		ifstream fin; fin.open(fname.c_str());
		while (!fin.eof()){
			getline(fin, line);
			if (line.size() == 0){
				continue;
			} // end if
	
			// read data
			if (line.substr(0,6) == "ACCESS"){
				d = d + 1;
				data.push_back(residue_asa());

				resSeq_tmp = strip(line.substr(17,3), " ") + ":" + strip(line.substr(6,5), " ");
				data[d].set_resSeq(resSeq_tmp);
				data[d].set_resname(strip(line.substr(11,6), " "));

				finder[resSeq_tmp] = d;

				data[d].set_all_sum(atof(strip(line.substr(20,6), " ").c_str()));
				data[d].set_all_per(atof(strip(line.substr(26,6), " ").c_str()));
				data[d].set_nonP_sum(atof(strip(line.substr(34,6), " ").c_str()));
				data[d].set_nonP_per(atof(strip(line.substr(40,6), " ").c_str()));
				data[d].set_polar_sum(atof(strip(line.substr(48,6), " ").c_str()));
				data[d].set_polar_per(atof(strip(line.substr(54,6), " ").c_str()));
				data[d].set_side_sum(atof(strip(line.substr(62,6), " ").c_str()));
				data[d].set_side_per(atof(strip(line.substr(68,6), " ").c_str()));
				data[d].set_main_sum(atof(strip(line.substr(76,6), " ").c_str()));
				data[d].set_main_per(atof(strip(line.substr(82,6), " ").c_str()));
			} // end if
		} // end while
		fin.close(); fin.clear();
	} // end ASA
	~ASA(){}

	residue_asa residue(string resSeq){
		return data[finder[resSeq]];
	} // end residue_asa

	residue_asa residue(int k){
		return data[k];
	} // end residue_asa

}; // end class
