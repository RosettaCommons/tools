# ifndef PDB_H
# define PDB_H

# include <string>
# include <vector>
# include "easystring.h"
# include "shorthand.h"
# include <cstdlib>
# include <set>
# include <fstream>
# include <sstream>
# include <stdio.h>
# include <iostream>
using namespace std;


class PDB{
	private:
	vector <string> _atom, _name, _altLoc, _resName, _chainID, _resSeq, _iCode, _element, _charge;
	set < pair <string, string> > _residue_list;
	vector <float> _x, _y, _z, _occupancy, _T;
	vector <int> _serial;
	map < pair <string, string>, vector <int> > _residue;

	public:
	PDB(string filename, bool read_hetero, bool read_modified){
		string residue_20 [20] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
		vector <string> normal_residues(&residue_20[0], &residue_20[19]);
		string line, __atom, __element, __iCode, __altLoc, __resName;
		bool read = 0;
		int missed_n;
		string spaces = "                              ";

		ifstream fpdb; fpdb.open(filename.c_str());
		while (!fpdb.eof()){
			getline(fpdb, line);
			if (line.size() >= 54){ // at least contains coordinates
				if (line.size() < 80){
					line = line + spaces;
				} // end if
				try{
					__atom = strip(line.substr(0,6-0));
					if (__atom == "ATOM"){
						read = 1;
					} else if (__atom == "HETATM"){
						if (read_hetero == 1){
							read = 1;
						} else {
							read = 0;
						} // end if
					} else {
							read = 0;
					} // end if
					if (read == 1) {
						__element = strip(line.substr(76, 78-76));
						if (__element != "H") {
							__iCode = strip(line.substr(26, 27-26));
							if (__iCode == "") {
								__altLoc = line.substr(16, 17-16);
								if ( ( __altLoc == " ") || (__altLoc == "A") || (__altLoc == "1") ) {
									__resName = strip(line.substr(17, 20-17));
									read = 0;
									if (is_element(__resName, normal_residues)){
										read = 1;
									} else {
										if (read_modified == 1){
											read = 1;
										} // end if
									} // end if
									if (read == 1){
										_atom.push_back(strip(__atom));
										_serial.push_back(int(atof(line.substr(6, 11-6).c_str())));
										_name.push_back(line.substr(12, 16-12));
										_altLoc.push_back(strip(__altLoc));
										_resName.push_back(strip(__resName));
										_chainID.push_back(strip(line.substr(21, 22-21)));
										_resSeq.push_back(strip(line.substr(22, 26-22)));
										_iCode.push_back(strip(line.substr(26, 27-26)));
										_x.push_back(atof(line.substr(30, 38-30).c_str()));
										_y.push_back(atof(line.substr(38, 46-38).c_str()));
										_z.push_back(atof(line.substr(46, 54-46).c_str()));
										_occupancy.push_back(atof((line.substr(54, 60-54).c_str())));
										_T.push_back(atof(line.substr(60, 66-60).c_str()));
										__element = strip(line.substr(76, 78-76));
										_charge.push_back(strip(line.substr(78, 80-78)));
									} // end if
								} // end if
							} // end if
						} // end if
					} // end if
				} catch (exception& out_of_range){
					// do nothing
				} // end try
			} // end if
		} // end while
		fpdb.close(); fpdb.clear();

		// fill in void for possible missing values
		missed_n = _x.size() - _occupancy.size();
		for (int i = 0; i < missed_n; ++i){
			_occupancy.push_back(0);
		} // end for
		missed_n = _x.size() - _T.size();
		for (int i = 0; i < missed_n; ++i){
			_T.push_back(0);
		} // end for
		missed_n = _x.size() - _element.size();
		for (int i = 0; i < missed_n; ++i){
			_element.push_back("");
		} // end for
		missed_n = _x.size() - _charge.size();
		for (int i = 0; i < missed_n; ++i){
			_charge.push_back(string(""));
		} // end for
	} // end PDB

	PDB(string filename, vector <string> keywords){ // overload
		string line, __atom, __element, __iCode, __altLoc, __resName;
		int missed_n;
		string spaces = "                              ";

		ifstream fpdb; fpdb.open(filename.c_str());
		while (!fpdb.eof()){
			getline(fpdb, line);
			if (line.size() >= 54){ // at least contains coordinates
				if (line.size() < 80){
					line = line + spaces;
				} // end if
				try{
					__atom = strip(line.substr(0,6-0));
					if (is_element(__atom, keywords)){
						_atom.push_back(strip(__atom));
						_serial.push_back(int(atof(line.substr(6, 11-6).c_str())));
						_name.push_back(strip(line.substr(13, 16-13)));
						_altLoc.push_back(line.substr(16, 17-16));
						_resName.push_back(strip(line.substr(17, 20-17)));
						_chainID.push_back(strip(line.substr(21, 22-21)));
						_resSeq.push_back(strip(line.substr(22, 26-22)));
						_iCode.push_back(strip(line.substr(26, 27-26)));
						_x.push_back(atof(line.substr(30, 38-30).c_str()));
						_y.push_back(atof(line.substr(38, 46-38).c_str()));
						_z.push_back(atof(line.substr(46, 54-46).c_str()));
						_occupancy.push_back(atof((line.substr(54, 60-54).c_str())));
						_T.push_back(atof(line.substr(60, 66-60).c_str()));
						_element.push_back(strip(line.substr(76, 78-76)));
						_charge.push_back(strip(line.substr(78, 80-78)));
						_residue_list.insert( make_pair(strip(line.substr(21, 22-21)), strip(line.substr(22, 26-22))));
					} // end if
				} catch (exception& out_of_range){
					// do nothing
				} // end try
			} // end if
		} // end while
		fpdb.close(); fpdb.clear();

		// fill in void for possible missing values
		missed_n = _x.size() - _occupancy.size();
		for (int i = 0; i < missed_n; ++i){
			_occupancy.push_back(0);
		} // end for
		missed_n = _x.size() - _T.size();
		for (int i = 0; i < missed_n; ++i){
			_T.push_back(0);
		} // end for
		missed_n = _x.size() - _element.size();
		for (int i = 0; i < missed_n; ++i){
			_element.push_back("");
		} // end for
		missed_n = _x.size() - _charge.size();
		for (int i = 0; i < missed_n; ++i){
			_charge.push_back(string(""));
		} // end for

		// group atom into residue
		for (unsigned int i = 0; i < _x.size(); ++i){
			_residue[make_pair(_chainID[i], _resSeq[i])].push_back(i);
		} // end for

	} // end PDB

	~PDB(){}

	// size
	unsigned int size(){
		return _x.size();
	} // end size

	// set and get functions
	vector <string> atom(){
		return _atom;
	} // end atom()
	string atom(int i){
		return _atom[i];
	} // end atom()
	void set_atom(int i, string s){
		_atom[i] = s;
	} // end set

	vector <int> residue(pair <string, string> resname){
		return _residue[resname];
	} // end reside

	set < pair <string, string> > residue_list(){
		return _residue_list;
	} // end residue_list

	vector <string> name(){
		return _name;
	} // end atom()
	string name(int i){
		return _name[i];
	} // end atom()
	void set_name(int i, string s){
		_name[i] = s;
	} // end set

	vector <string> altLoc(){
		return _altLoc;
	} // end atom()
	string altLoc(int i){
		return _altLoc[i];
	} // end atom()
	void set_altLoc(int i, string s){
		_altLoc[i] = s;
	} // end set

	vector <string> resName(){
		return _resName;
	} // end atom()
	string resName(int i){
		return _resName[i];
	} // end atom()
	void set_resName(int i, string s){
		_resName[i] = s;
	} // end set

	vector <string> chainID(){
		return _chainID;
	} // end atom()
	string chainID(int i){
		return _chainID[i];
	} // end atom()
	void set_chainID(int i, string s){
		_chainID[i] = s;
	} // end set

	vector <string> resSeq(){
		return _resSeq;
	} // end atom()
	string resSeq(int i){
		return _resSeq[i];
	} // end atom()
	void set_resSeq(int i, string s){
		_resSeq[i] = s;
	} // end set

	vector <string> iCode(){
		return _iCode;
	} // end atom()
	string iCode(int i){
		return _iCode[i];
	} // end atom()
	void set_iCode(int i, string s){
		_iCode[i] = s;
	} // end set

	vector <string> element(){
		return _element;
	} // end atom()
	string element(int i){
		return _element[i];
	} // end atom()
	void set_element(int i, string s){
		_element[i] = s;
	} // end set

	vector <string> charge(){
		return _charge;
	} // end atom()
	string charge(int i){
		return _charge[i];
	} // end atom()
	void set_charge(int i, string s){
		_charge[i] = s;
	} // end set

	vector <float> x(){
            return _x;
	} // end atom()
	float x(int i){
		return _x[i];
	} // end atom()
	void set_x(int i, float s){
		_x[i] = s;
	} // end set

	vector <float> y(){
            return _y;
	} // end atom()
	float y(int i){
		return _y[i];
	} // end atom()
	void set_y(int i, float s){
		_y[i] = s;
	} // end set

	vector <float> z(){
            return _z;
	} // end atom()
	float z(int i){
		return _z[i];
	} // end atom()
	void set_z(int i, float s){
		_z[i] = s;
	} // end set

	vector <float> T(){
            return _T;
	} // end atom()
	float T(int i){
		return _T[i];
	} // end atom()
	void set_T(int i, float s){
		_T[i] = s;
	} // end set

	vector <float> occupancy(){
            return _occupancy;
	} // end atom()
	float occupancy(int i){
		return _occupancy[i];
	} // end atom()
	void set_occupancy(int i, float s){
		_occupancy[i] = s;
	} // end set

	vector <int> serial(){
            return _serial;
	} // end atom()
	float serial(int i){
		return _serial[i];
	} // end atom()
	void set_serial(int i, int s){
		_occupancy[i] = s;
	} // end set

	// write to file
	void write(string outfilename){
		ofstream fout; fout.open(outfilename.c_str());
		for (unsigned int i = 0; i < _x.size(); ++i){
			fout << pdb_line(_atom[i], _serial[i], _name[i], _altLoc[i], _resName[i], _chainID[i], _resSeq[i], _iCode[i], _x[i], _y[i], _z[i], _occupancy[i], _T[i], _element[i], _charge[i]);
		} // end for
		fout.close(); fout.clear();
	} // end write
}; // end class

# endif
