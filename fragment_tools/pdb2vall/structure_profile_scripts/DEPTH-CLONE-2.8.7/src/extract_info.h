# include <iostream>
# include <string>
# include <fstream>
# include <map>
# include <stdlib.h>
# include "easystring.h"
# include "shorthand.h"
# include "PDB.h"
using namespace std;

class depth_recorder{
	private:
	vector <double> recorder;
	double mean, stdev;
	string name, restype;

	public:
	depth_recorder(){} // end depth_recorder
	~depth_recorder(){}

	void set_name(string _name){
		name = _name;
	} // end set_name

	void insert(double d){
		recorder.push_back(d);
	} // end insert

	void consolidate(){
		mean = 0; stdev = 0;
		for (unsigned int i = 0; i < recorder.size(); ++i){
			mean = mean + recorder[i]; stdev = stdev + recorder[i]*recorder[i];
		} // end for
		if (recorder.size() != 0){
			mean = mean / recorder.size();
			stdev = sqrt(stdev / recorder.size() - mean*mean);
		} else {
			mean = -99; stdev = -99;
		} // end if
	} // end consolidate
	
	void set_restype(string _restype_in){
		restype = _restype_in;
	} // end set_restype

	string res_type(){
		return restype;
	}// end res_type

	double get_mean(){
		return mean;
	} // end get_mean

	double get_stdev(){
		return stdev;
	} // end get_stdev

	string get_name(){
		return name;
	} // end get_name
}; // end class

string check_value(double s){
	if (s != -99){
		return align(strf(num2string(s),2,0), 6, "r");
	} else {
		return align("nan", 6, "r");
	} // end if
} // end check_value


class depth_assigner{ // Usage: consolidate atomic depth into residue depth
	private:
	string bufferline, resSeq, name, resname, restype, chain_i, seq_i;
	char element;
	vector <string> seq, residue_list, keywords;
	vector <int> atoms;
	map <string, depth_recorder> res_record, mc_record, sc_record, sc_polar_record, sc_nonpolar_record;

	public:
	depth_assigner(){
		keywords.push_back("ATOM");
	};
	~depth_assigner(){};

	map <string, float> assign(string input_pdb, string output_fname){
		PDB mdl = PDB(input_pdb, keywords);
		for (unsigned int i = 0; i < mdl.size(); ++i){
			resname = mdl.chainID(i)+":"+mdl.resSeq(i);
			if (is_element(resname, residue_list) == 0){
				residue_list.push_back(resname);
			} // end if
		} // end for

		for (unsigned int i = 0; i < residue_list.size(); ++i){
			resname = residue_list[i];
			chain_i = split(residue_list[i], ":")[0];
			seq_i = split(residue_list[i], ":")[1];
			atoms = mdl.residue(make_pair(chain_i, seq_i));
			restype = mdl.resName(atoms[0]);

			// begin updating //
			res_record.insert(make_pair(resname, depth_recorder()));
			res_record[resname].set_restype(restype);

			mc_record.insert(make_pair(resname, depth_recorder()));
			mc_record[resname].set_restype(restype);

			sc_record.insert(make_pair(resname, depth_recorder()));
			sc_record[resname].set_restype(restype);

			sc_polar_record.insert(make_pair(resname, depth_recorder()));
			sc_polar_record[resname].set_restype(restype);

			sc_nonpolar_record.insert(make_pair(resname, depth_recorder()));
			sc_nonpolar_record[resname].set_restype(restype);

			for (unsigned int j = 0; j < atoms.size(); ++j){
				name = mdl.name(atoms[j]); element = name[0];
				res_record[resname].insert(mdl.T(atoms[j]));
				if ((name=="N") || (name=="CA") || (name=="C") || (name=="O")){ // main chain atom
					mc_record[resname].insert(mdl.T(atoms[j]));
				} else {
					sc_record[resname].insert(mdl.T(atoms[j]));
					if (element == 'C'){
						sc_nonpolar_record[resname].insert(mdl.T(atoms[j]));
					} else {
						sc_polar_record[resname].insert(mdl.T(atoms[j]));
					} // end if
				} // end if
			} // end for
		} // end for

		ofstream fout; fout.open(output_fname.c_str());
		fout << "# chain:residue\tall-atom\tall-atom(stdev)\tMC-atom\tMC-atom(stdev)\tSC-atom\tSC-atom(stdev)\tSC-polar-atom\tSC-polar-atom(stdev)\tSC-nonpolar-atom\tSC-nonpolar-atom(stdev)" << endl;
		for (unsigned int i = 0; i < residue_list.size(); ++i){
			resname = residue_list[i];
			res_record[resname].consolidate(); 
			mc_record[resname].consolidate(); 
			sc_record[resname].consolidate(); 
			sc_polar_record[resname].consolidate();  
			sc_nonpolar_record[resname].consolidate();

			fout << resname << "\t" << res_record[resname].res_type() << "\t" << check_value(round_double(res_record[resname].get_mean(),2)) << "\t" << check_value(round_double(res_record[resname].get_stdev(),2)) << "\t" << check_value(round_double(mc_record[resname].get_mean(),2)) << "\t" << check_value(round_double(mc_record[resname].get_stdev(),2)) << "\t" << check_value(round_double(sc_record[resname].get_mean(),2)) << "\t" << check_value(round_double(sc_record[resname].get_stdev(),2)) << "\t" << check_value(round_double(sc_polar_record[resname].get_mean(),2)) << "\t" << check_value(round_double(sc_polar_record[resname].get_stdev(),2)) << "\t" << check_value(round_double(sc_nonpolar_record[resname].get_mean(),2)) << "\t" << check_value(round_double(sc_nonpolar_record[resname].get_stdev(),2)) << endl;
		} // end for
		fout.close(); fout.clear();

		string res;
		map <string, float> output;
		for (unsigned int i = 0; i < residue_list.size(); ++i){
			res = residue_list[i];
			output[res] = atof(check_value(round_double(res_record[res].get_mean(),2)).c_str());
		} // end for
		return output;

	} // end assign

};
