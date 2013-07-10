// C++ version of residue depth
// rewritten by KuanPern-Tan

// depth is the distance of an atom to its nearest surface solvent molecule atom
// main program to link and run depth program
// RESIDUE-DEPTH(CLONE)
// version 2.8
// Author: KuanPern Tan
// Date: 2 Nov 2010
// Institute: Bioinformatics Institute A*STAR, Singapore

# include <cstdlib>
# include <iostream>
# include <string>
# include <vector>
# include <map>
# include <fstream>
# include <sstream>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include "omp.h"
using namespace std;

// user defined libraries
# include "install_dir.h"
# include "PDB.h"
# include "easystring.h"
# include "easyoption.h"
# include "shorthand.h"
# include "depth.h"
# include "extract_info.h"
# include "substitutor.h"


int main(int argc,char* argv[]){
	// --- GET INPUT PARAMETERS ---
	// step 0: use current time as job identifier
	time_t seconds = time (NULL); stringstream ss; ss << seconds; string job_name_tmp = ss.str();
	string job_name = job_name_tmp.substr(job_name_tmp.size()-6,6);
	ofstream flog; flog.open((job_name+".log").c_str()); // log file
	flog << "job id: " << job_name << endl;
	// step 1: get input from user
	// step 1.1: print instruction
	if (argc == 1){
		cout << "Calculate atomic depth for an input pdb file." << endl;
		cout << "OPTION\n-i input.pdb \n-o output \n-n iteration-number (Default 25) \n-survive minimum-number-of-neighbouring-solvent (Default 2) \n-keep minimally-solvated-models-filename (Default NULL) \n-cavity cavity-size (Default 4.2)  \n-w solvent_box.pdb (Default "+LIB_DIR+"/water_models/spc216.pdb) \n-box_x length of solvation box in x-dimenson (Default 19)\n-box_y length of solvation box in y-dimension (Default 19)\n-box_z length of solvation box in z-dimension (Default 19)\n-sol_diameter solvent-diameter (Default 2.8) \n-clash_dist vdw-distance(solvent, solute) (Default 2.6)\n-thread number-of-threads" << endl;
		exit(1);
	} // end if
	// step 1.2: define parser for parameters
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
	parser.set_arg("-thread","thread","2");
	// step 1.3: parse the command line input
	parser.parse(argc,argv);
	// step 1.4: assign input parameters
	string infilename = parser.value("input");
	string outfilename = parser.value("output");
	if (infilename == ""){
		cerr << "please provide input file (pdb file)" << endl;
		exit(1);
	} else if (outfilename == ""){
		cerr << "please provide output file name" << endl;
		exit(1);
	} // end if
	static unsigned int iteration_n = atoi(parser.value("iteration").c_str());
	static unsigned int survival_n = atoi(parser.value("survival_n").c_str());
	string water_model = parser.value("water_model");
	static float radius_w = atof(parser.value("sol_diameter").c_str())/2;
	static float clash_dist = atof(parser.value("clash_dist").c_str());
	static float cavity = atof(parser.value("cavity").c_str());
	string min_sol_name = parser.value("min_sol_name");
	static float unit_box_x = atof(parser.value("box_x").c_str());
	static float unit_box_y = atof(parser.value("box_y").c_str());
	static float unit_box_z = atof(parser.value("box_z").c_str());
	static unsigned int thread = atoi(parser.value("thread").c_str());

	// step 2: read pdb for solute atom
	string line, bufferfname1, bufferfname2;
	vector <string> keywords; // decide what entry to read
	keywords.push_back("ATOM"); // protein atom included by default
	PDB mdl_in = PDB(infilename, keywords);
	unsigned int atom_number = mdl_in.size();
	float** x = new float* [iteration_n]; // note: dynamic memory allocation is used here to support large protein file
	float** y = new float* [iteration_n];
	float** z = new float* [iteration_n];
	float* x_static = new float [atom_number]; // static copy of PDB coordinates
	float* y_static = new float [atom_number];
	float* z_static = new float [atom_number];
	for (unsigned int i = 0; i < atom_number; ++i){
		x_static[i] = mdl_in.x(i); y_static[i] = mdl_in.y(i); z_static[i] = mdl_in.z(i);
	} // end for
	flog << "read solute atoms, totaling " << atom_number << " atoms" << endl;

	// step 3: read solvent model
	vector <string> solvent_keywords; solvent_keywords.push_back("ATOM");
	PDB sol_mdl = PDB(water_model, keywords);
	unsigned int solvent_number = sol_mdl.size();
	float xw[iteration_n][solvent_number], yw[iteration_n][solvent_number], zw[iteration_n][solvent_number];
	for (unsigned int i = 0; i < solvent_number; ++i){
		xw[0][i] = sol_mdl.x(i); yw[0][i] = sol_mdl.y(i); zw[0][i] = sol_mdl.z(i);
	} // end for
	for (unsigned int i = 1; i < iteration_n; ++i){
		memcpy(xw[i], xw[0], solvent_number*sizeof(float)); 
		memcpy(yw[i], yw[0], solvent_number*sizeof(float)); 
		memcpy(zw[i], zw[0], solvent_number*sizeof(float)); 
	} // end for

	// --- RUN DEPTH ALGORITHM ---
	// step 0: prepare parameters and memory block for calculation
	map <int,string> dimension;
	dimension[0] = "xy"; dimension[1] = "xz"; dimension[2] = "yz";
	static float sol_diameter = 2*radius_w;
	float random_translation[iteration_n], random_rotation[iteration_n]; int random_axis[iteration_n];
	float** depth_one = new float* [iteration_n]; // working memory for atomic depth
	float*  ave_depth = new float  [atom_number]; // average of depth
	float*  std_depth = new float  [atom_number]; // std of depth
	depth_run instance[iteration_n]; // chare array
	// step 1: assign parameter to each chare
	for (unsigned int i = 0; i < iteration_n; ++i){
		random_translation[i] = rand()/float(RAND_MAX)*sol_diameter*2 - sol_diameter;
		random_rotation[i] = rand()/float(RAND_MAX)*360;
		random_axis[i] = int(rand()/float(RAND_MAX)*3);
		instance[i].set_atom_number(atom_number);
		instance[i].set_solvent_number(solvent_number);
		instance[i].set_label(i, job_name);
		instance[i].initialize_box(6); // cut-off for cell list algorithm
		instance[i].set_parameter(survival_n, cavity, dimension[random_axis[i]], random_rotation[i], random_translation[i], radius_w, clash_dist, min_sol_name, unit_box_x, unit_box_y, unit_box_z);
	} // end for
	// step 2: run simulations

	omp_set_num_threads(thread);
	# pragma omp parallel for
	for (unsigned int i = 0; i < iteration_n; ++i){ // every char
		depth_one[i] = new float [atom_number]; // get working memory
		x[i] = new float [atom_number]; memcpy(x[i], x_static, atom_number*sizeof(float)); // coordinates memory
		y[i] = new float [atom_number]; memcpy(y[i], y_static, atom_number*sizeof(float));
		z[i] = new float [atom_number]; memcpy(z[i], z_static, atom_number*sizeof(float));
		instance[i].get_depth(depth_one[i], x[i],y[i],z[i],xw[i],yw[i],zw[i]); // and run
		for (unsigned int a = 0; a < atom_number; ++a){ // update output from char
			ave_depth[a] = ave_depth[a] + depth_one[i][a];
			std_depth[a] = std_depth[a] + depth_one[i][a]*depth_one[i][a];
		} // end for
		delete[] x[i]; delete[] y[i]; delete[] z[i]; delete[] depth_one[i]; // free memory
	} // end for
	delete[] depth_one; // free memory
	delete[] x; delete[] y; delete[] z;
	delete[] x_static; delete[] y_static; delete[] z_static;
	// step 2.1: consolidate data
	for (unsigned int a = 0; a < atom_number; ++a){
		ave_depth[a] = ave_depth[a] / iteration_n;
		std_depth[a] = sqrt(((std_depth[a]/iteration_n) - ave_depth[a]*ave_depth[a]));
	} // end for

	// --- OUTPUT ---
	// step 1: output colorized (atomic) PDB and atomic depth
	bufferfname1 = outfilename+"-atomic_depth.pdb"; bufferfname2 = outfilename+"-atom.depth"; // filename
	PDB mdl = PDB(infilename, keywords);
	ofstream fatom ; fatom.open(bufferfname2.c_str());
	fatom << "index\tresidue\tatom-type\tmean(depth)\tstdev(depth)" << endl;
	for (unsigned int i = 0; i < mdl.size(); ++i){
		mdl.set_T(i, ave_depth[i]);
		fatom << mdl.serial(i) << "\t" << mdl.chainID(i) << ":" << mdl.resSeq(i) << "\t" << mdl.name(i) << "\t" << num2string(round_double(ave_depth[i],2)) << "\t" << num2string(round_double(std_depth[i],2)) << endl;
	} // end for
	mdl.write(bufferfname1);
	fatom.close(); fatom.clear();
	// step 2: output colorized (residue) PDB and residue depth
	bufferfname1 = outfilename+"-residue_depth.pdb"; bufferfname2 = outfilename+"-residue.depth";
	depth_assigner assigner; // assign atoms to its residue
	map <string, float> residue_depth = assigner.assign(outfilename+"-atomic_depth.pdb", bufferfname2);	
	PDB mdl_depth = PDB(infilename, keywords);
	for (unsigned int i = 0; i < mdl_depth.size(); ++i){
		mdl.set_T(i, residue_depth[mdl.chainID(i)+":"+mdl.resSeq(i)]);
	} // end for
	mdl.write(bufferfname1);

	// --- LOG FILE ---
	// step 1: combine log files from iterations
	ifstream bitflog;
	for (unsigned int i = 0; i < iteration_n; ++i){
		bitflog.open((job_name+"-iter"+num2string(i)).c_str());
		flog << endl << "iteration: " << i+1 << endl;
		while (!bitflog.eof()){
			getline(bitflog, line);
			flog << line << endl;
		} // end while
		bitflog.close(); bitflog.clear();
		remove((job_name+"-iter"+num2string(i)).c_str()); // remove log files from individual simulation
	} // end for
	flog.close(); flog.clear();
	// step 2: keep minimal solvated models if requested
	if (min_sol_name != ""){
		substitutor S;
		for (unsigned int i = 0; i < iteration_n; ++i){
			S.substitute(bufferfname1, min_sol_name+"-"+num2string(i)+".coordinates", min_sol_name+"-"+num2string(i)+".pdb");
		} // end for
	} // end if
	return 0;

} // end main
