# include <iostream>
# include <fstream>
# include <stdlib.h>
# include <stdio.h>
# include <vector>
# include <set>
# include <string>
# include "easystring.h"
# include "easyoption.h"
# include "consolidator.h"
# include <algorithm>
# include "install_dir.h" // generated on compilation to specify install directory
using namespace std;

int main(int argc,char* argv[]){
	// check inputs
	if (argc == 1){
		cout << "Compare depth of chains before and after forming complexes" << endl;
		cout << "OPTION\n-i input.pdb \n-o output \n-n iteration-number (Default 25) \n-survive minimum-number-of-neighbouring-solvent (Default 2) \n-keep minimally-solvated-models-filename (Default NULL) \n-cavity cavity-size (Default 4.2)  \n-w solvent_box.pdb (Default "+LIB_DIR+"/water_models/spc216.pdb) \n-box_x length of solvation box in x-dimenson (Default 19)\n-box_y length of solvation box in y-dimension (Default 19)\n-box_z length of solvation box in z-dimension (Default 19)\n-sol_diameter solvent-diameter (Default 2.8) \n-clash_dist vdw-distance(solvent, solute) (Default 2.6)" << endl;
		exit(1);
	} // end if

	// define inputs
	option parser;
	parser.set_arg("-i","input","");
	parser.set_arg("-o","output","");
	parser.set_arg("-n","iteration","25");
	parser.set_arg("-w","water_model", LIB_DIR+"/water_models/spc216.pdb");
	parser.set_arg("-box_x","box_x","19");
	parser.set_arg("-box_y","box_y","19");
	parser.set_arg("-box_z","box_z","19");
	parser.set_arg("-sol_diameter","sol_diameter","2.8");
	parser.set_arg("-clash_dist","clash_dist","2.6");
	parser.set_arg("-cavity","cavity","4.2");
	parser.set_arg("-survive","survival_n","2");
	parser.set_arg("-keep","min_sol_name","");


	// parse the command line
	parser.parse(argc,argv);

	// get input filename
	string infilename = parser.value("input");
	string outfilename = parser.value("output");
	if (infilename == ""){
		cerr << "please provide input file (pdb file)" << endl;
		exit(1);
	} else if (outfilename == ""){
		cerr << "please provide output file name" << endl;
		exit(1);
	} // end if

	// parse input parameters
	int iteration_n = atoi(parser.value("iteration").c_str());
	int survival_n = atoi(parser.value("survival_n").c_str());
	string water_model = parser.value("water_model");
	float radius_w = atof(parser.value("sol_diameter").c_str())/2;
	float clash_dist = atof(parser.value("clash_dist").c_str());
	float cavity = atof(parser.value("cavity").c_str());
	string min_sol_name = parser.value("min_sol_name");
	float unit_box_x = atof(parser.value("box_x").c_str());
	float unit_box_y = atof(parser.value("box_y").c_str());
	float unit_box_z = atof(parser.value("box_z").c_str());

	// deconvolute complex into different chains
	string line, newfilename;
	char chainID_new;
	set <char> chains;
	ifstream fin; fin.open(infilename.c_str());
	while (!fin.eof()){
		getline(fin, line);
		if (line.size() != 0){
			if (line.substr(0,4) == "ATOM"){
				chainID_new = line.substr(21,1)[0];
				chains.insert(chainID_new);
			} // end if
		} // end if
	} // end while
	fin.close(); fin.clear();

	vector <string> chain_fnames;
	ofstream fout[chains.size()];
	map <char, int> chain_map;
	unsigned int tmp = -1;
	for (set <char>::iterator iter = chains.begin(); iter != chains.end(); ++iter){
		tmp = tmp + 1;
		chain_map[*iter] = tmp;
		fout[tmp].open((infilename+(*iter)).c_str());
		chain_fnames.push_back((infilename+(*iter)).c_str());
	} // end for


	ofstream fcomplex; fcomplex.open((infilename+"-complex").c_str());

	fin.open(infilename.c_str());
	while (!fin.eof()){
		getline(fin, line);
		if (line.size() != 0){
			if (line.substr(0,4) == "ATOM"){
				chainID_new = line.substr(21,1)[0];
				fout[chain_map[chainID_new]] << line << endl;
				fcomplex << line << endl;
			} // end if
		} // end if
	} // end while
	fin.close(); fin.clear();
	for (unsigned int i = 0; i < chains.size(); ++i){
		fout[i].close(); fout[i].clear();
	} // end for
	fcomplex.close(); fcomplex.clear();

	// compute depths for each chain/domain
	string cmd;
	for (unsigned int i = 0; i < chain_fnames.size(); ++i){
		// compute depth value
		if (min_sol_name == ""){
			cmd = EXE_DIR+"/DEPTH -i "+chain_fnames[i]+" -o "+chain_fnames[i]+" -n "+num2string(iteration_n)+" -survive "+num2string(survival_n)+" -cavity "+num2string(cavity)+" -w "+water_model+" -box_x "+num2string(unit_box_x)+" -box_y "+num2string(unit_box_y)+" -box_z "+num2string(unit_box_z)+" -sol_diameter "+num2string(radius_w)+" -clash_dist "+num2string(clash_dist);
		} else {
			cmd = EXE_DIR+"/DEPTH -i "+chain_fnames[i]+" -o "+chain_fnames[i]+" -n "+num2string(iteration_n)+" -survive "+num2string(survival_n)+" -cavity "+num2string(cavity)+" -w "+water_model+" -box_x "+num2string(unit_box_x)+" -box_y "+num2string(unit_box_y)+" -box_z "+num2string(unit_box_z)+" -sol_diameter "+num2string(radius_w)+" -clash_dist "+num2string(clash_dist)+" -keep "+min_sol_name;
		} // end if
		cout << "processing " << chain_fnames[i] << endl;
		system(cmd.c_str());
	} // end for

	// consolidate information
	ofstream fconsolidate; fconsolidate.open((infilename+"-consolidate").c_str());
	fconsolidate << "# chain residue\tall-atom\tMC-atom\tSC-atom\tSC-polar-atom\tSC-nonpolar-atom" << endl;
	ifstream fstream;
	for (unsigned int i = 0; i < chain_fnames.size(); ++i){
		newfilename = chain_fnames[i] + "-residue.depth";
		fstream.open(newfilename.c_str());
		getline(fstream, line);
		while (!fstream.eof()){
			getline(fstream, line);
			if (line.size() == 0){
				continue;
			} // end if
			if (line.substr(0,1) == "#"){
				continue;
			} // end if
			fconsolidate << line << endl;
		} // end while
		fstream.close(); fstream.clear();
		remove(newfilename.c_str());
	} // end for
	fconsolidate.close(); fconsolidate.clear();

	// compute depth for complex
	if (min_sol_name == ""){
		cmd = EXE_DIR+"/DEPTH -i "+infilename+" -o "+infilename+"-complex -n "+num2string(iteration_n)+" -survive "+num2string(survival_n)+" -cavity "+num2string(cavity)+" -w "+water_model+" -box_x "+num2string(unit_box_x)+" -box_y "+num2string(unit_box_y)+" -box_z "+num2string(unit_box_z)+" -sol_diameter "+num2string(radius_w)+" -clash_dist "+num2string(clash_dist);
	} else {
		cmd = EXE_DIR+"/DEPTH -i "+infilename+" -o "+infilename+"-complex -n "+num2string(iteration_n)+" -survive "+num2string(survival_n)+" -cavity "+num2string(cavity)+" -w "+water_model+" -box_x "+num2string(unit_box_x)+" -box_y "+num2string(unit_box_y)+" -box_z "+num2string(unit_box_z)+" -sol_diameter "+num2string(radius_w)+" -clash_dist "+num2string(clash_dist)+" -keep "+min_sol_name;
	} // end if
	cout << "processing " << infilename << "-complex " << endl;
	system(cmd.c_str());

	// compute difference between two data-set
	cout << "consolidating ..." << endl;
	consolidator C;
	C.consolidate(infilename+"-complex-residue.depth", infilename+"-consolidate", outfilename);
	remove((infilename+"-complex-residue.depth").c_str());
	remove((infilename+"-consolidate").c_str());
	remove((infilename+"-complex").c_str());
	for (unsigned int i = 0; i < chain_fnames.size(); ++i){
		remove(chain_fnames[i].c_str());
	} // end for

	return 0;
} // end main
