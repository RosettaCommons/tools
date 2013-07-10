# include <string>
# include <vector>
# include <cstdlib>
# include <fstream>
# include <sstream>
# include <stdio.h>
# include <iostream>
# include "PDB.h"
using namespace std;

class substitutor{
	private:
	vector <string> data, bufferline;
	string line, sol_name, element;
	int serial, resSeq, solvent_index;
	float x, y, z;

	public:
	substitutor(){}
	~substitutor(){}

	void substitute(string pdb_fname, string datafile, string outfilename){
		data.clear(); 
		PDB mdl = PDB(pdb_fname, 0, 1);
		ifstream fdata; fdata.open(datafile.c_str());
		while (!fdata.eof()){
			getline(fdata, line);
			if (line.size() == 0){
				continue;
			} // end if
			data.push_back(line);
		} // end while
		fdata.close(); fdata.clear();

		for (unsigned int i = 0; i < mdl.size(); ++i){
			bufferline = split(data[i], "\t");
			mdl.set_x(i, atof(bufferline[1].c_str()));
			mdl.set_y(i, atof(bufferline[2].c_str()));
			mdl.set_z(i, atof(bufferline[3].c_str()));
		} // end for
		mdl.write(outfilename);

		serial = int(mdl.serial(mdl.size()-1));
		resSeq = atoi(mdl.resSeq(mdl.size()-1).c_str());
		solvent_index = -1;


		ofstream fout; fout.open(outfilename.c_str(), ios::app);
		for (unsigned int i = mdl.size(); i < data.size(); ++i){
			serial = serial + 1;
			solvent_index = solvent_index + 1;
			if (solvent_index % 3 == 0){
				resSeq = resSeq + 1;
				sol_name = "OW";
			} else if (solvent_index % 3 == 1){
				sol_name = "HW1";
			} else if (solvent_index % 3 == 2){
				sol_name = "HW2";
			} // end if
			element = sol_name.substr(0,1);

			bufferline = split(data[i], "\t");
			x = atof(bufferline[1].c_str());
			y = atof(bufferline[2].c_str());
			z = atof(bufferline[3].c_str());

			fout << pdb_line("HETATM", serial, sol_name, " ", "SOL", " ", num2string(resSeq), " ", x, y, z, 1.00, 2.8, element, " ");
		} // end for
		fout.close(); fout.clear();
		remove(datafile.c_str());
	} // end 
	
}; // end class
