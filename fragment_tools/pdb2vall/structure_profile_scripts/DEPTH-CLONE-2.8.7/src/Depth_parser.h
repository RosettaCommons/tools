# include <iostream>
# include <fstream>
# include <math.h>
# include <string>
# include <map>
# include <cmath>
# include "easystring.h"
using namespace std;

class residue_depth{ // class object to holder modeller accessibility information
	private:
	float _depth_all, _depth_MC, _depth_SC, _depth_SCP, _depth_SCNP;
	string _restype, _resSeq;

	public:
	residue_depth(){}
	~residue_depth(){}

	// set and get functions
	void set_resSeq(string s){
		_resSeq = s;
	} // end set_resSeq
	void set_restype(string s){
		_restype = s;
	} // end set_restype
	void set_depth_all(float s){
		_depth_all = s;
	} // end set_		
	void set_depth_SCNP(float s){
		_depth_SCNP = s;
	} // end set_		
	void set_depth_SCP(float s){
		_depth_SCP = s;
	} // end set_		
	void set_depth_SC(float s){
		_depth_SC = s;
	} // end set_		
	void set_depth_MC(float s){
		_depth_MC = s;
	} // end set_

	float depth_all(){
		return _depth_all;
	} // end 
	float depth_SCNP(){
		return _depth_SCNP;
	} // end 
	float depth_SCP(){
		return _depth_SCP;
	} // end 
	float depth_SC(){
		return _depth_SC;
	} // end 
	float depth_MC(){
		return _depth_MC;
	} // end 
	string restype(){
		return _restype;
	} // end restype
	string resSeq(){
		return _resSeq;
	} // end resSeq

}; // end class


class Depth{
	private:
	vector <residue_depth> data;
	vector <string> bufferline, _all_residues;
	string line, resSeq, restype;
	int d;
	map <string, int> finder;

	public:
	Depth(string fname){
		d = -1;
		ifstream fdepth; fdepth.open(fname.c_str());
		while (!fdepth.eof()){
			getline(fdepth, line);
			if (line.size() == 0){
				continue;
			} // end if
			if (line.substr(0, 1) == "#"){
				continue;
			} // end if
	
			// read data
			bufferline = split(line, "\t");
			resSeq = bufferline[0];
			restype = bufferline[1];

			data.push_back(residue_depth());
			d = d + 1;
			finder[resSeq] = d;
			_all_residues.push_back(resSeq);

			if (isNumeric(bufferline[2])){
				data[d].set_depth_all(atof(bufferline[2].c_str()));
			} else {
				data[d].set_depth_all(-99);
			} // end if

			if (isNumeric(bufferline[4])){
				data[d].set_depth_MC(atof(bufferline[4].c_str()));
			} else {
				data[d].set_depth_MC(-99);
			} // end if

			if (isNumeric(bufferline[6])){
				data[d].set_depth_SC(atof(bufferline[6].c_str()));
			} else {
				data[d].set_depth_SC(-99);
			} // end if

			if (isNumeric(bufferline[8])){
				data[d].set_depth_SCP(atof(bufferline[8].c_str()));
			} else {
				data[d].set_depth_SCP(-99);
			} // end if

			if (isNumeric(bufferline[10])){
				data[d].set_depth_SCNP(atof(bufferline[10].c_str()));
			} else {
				data[d].set_depth_SCNP(-99);
			} // end if
		} // end while
		fdepth.close(); fdepth.clear();
	} // end Depth
	~Depth(){}

	residue_depth residue(string resSeq){
		return data[finder[resSeq]];
	} // end residue_depth

	residue_depth residue(int k){
		return data[k];
	} // end residue_depth

	vector <string> all_residues(){
		return _all_residues;
	} // end all_residues

}; // end class


