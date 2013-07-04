# include <cmath>
# include <string>
# include <map>
# include <vector>
# include <fstream>
# include <cstdlib>
# include "shorthand.h"
# include "easyoption.h"
# include "easystring.h"
# include "modeller_psa.h"
# include "Depth_parser.h"
# include "PDB.h"
# include "grid_map.h"
# include "install_dir.h" // generated on compilation to specify install directory
using namespace std;

int main(int argc,char* argv[]){
	if (argc ==1){
		cout << "Usage: Predict probability of involvement of residues in a binding cavity" << endl;
		cout << "Option:\n-d depth-profile \n-a asa-profile \n-p parameter-file \n-o output \n-c input.pdb \n-y output.pdb \n-e output.depth-asa (Default NULL)" << endl;

	} // end if

	option parser;
	parser.set_arg("-d","depth_profile","");
	parser.set_arg("-a","asa_profile","");
	parser.set_arg("-p","parameter_file",LIB_DIR+"/par/survive2.odd.map");
	parser.set_arg("-o","output_file","");
	parser.set_arg("-c","in_pdb", "");
	parser.set_arg("-e","out_dev", "");
	parser.set_arg("-y","out_struct", "");

	parser.parse(argc, argv);

	string depth_profile = parser.value("depth_profile");
	string asa_profile = parser.value("asa_profile");
	string parameter = parser.value("parameter_file");
	string outfile = parser.value("output_file");
	string in_pdb = parser.value("in_pdb");
	string out_depth_asa = parser.value("out_dev");
	string out_pdb = parser.value("out_struct");

	if ((depth_profile == "") || (asa_profile == "") || (parameter == "") || (in_pdb == "") || (outfile == "")){
		cerr << "error: please check input files" << endl;
		exit(1);
	} // end if

	vector <string> par = split(parameter,":");
	string amino_acids [20] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

	string line, resname;
	vector <string> bufferline;
	float d, a;
	// load density maps
	map <string, Grid_Map> U;
	for (unsigned int i = 0; i < 20; ++i){
		U[amino_acids[i]].read(par[0]+amino_acids[i]+"_survive"+par[1]+".odd.map");
	} // end for
	cout << "read maps" << endl;

	// read residue-depth data
	Depth D = Depth(depth_profile.c_str());
	cout << "read depth profile" << endl;
	ASA A = ASA(asa_profile.c_str());
	cout << "read asa profile" << endl;

	// get all residue names (according to depth)
	vector <string> prediction;
	vector <string> all_residues = D.all_residues();
	map <string, int> finder;

	// read residue type
	map <string, string> residue_type;
	vector <string> keyword; keyword.push_back("ATOM");
	PDB mdl = PDB(in_pdb, keyword);
	for (unsigned int i = 0; i < mdl.size(); ++i){
		residue_type[mdl.chainID(i)+":"+mdl.resSeq(i)] = mdl.resName(i);
	} // end for
	cout << "read PDB" << endl;


	float o; string restype;
	for (unsigned int i = 0; i < all_residues.size(); ++i){
		resname = all_residues[i];
		finder[resname] = i;
		restype = residue_type[resname];
		d = D.residue(resname).depth_all(); a = A.residue(resname).all_per();
		o = U[restype].value(d, a);
		prediction.push_back(strf(o));
	} // end for
	cout << "finished computing prediction" << endl;

	ofstream fout; fout.open(outfile.c_str());
	fout << "# RES Prob." << endl;
	for (unsigned int i = 0; i < all_residues.size(); ++i){
		fout << all_residues[i] << " " << prediction[i] << endl;
	} // end for
	fout.close(); fout.clear();
	cout << "finished writing output" << endl;

	// color pdb if requested
	if (out_pdb != ""){
		string resSeq;
		float value;
		PDB mdl = PDB(in_pdb, 0, 1);
		for (unsigned int i = 0; i < mdl.size(); ++i){
			resSeq = mdl.chainID(i)+":"+mdl.resSeq(i);
			value = atof(prediction[finder[resSeq]].c_str());
			mdl.set_T(i, value);
		} // end for
		mdl.write(out_pdb);
	} // end if

	if (out_depth_asa != ""){
		ofstream fdev; fdev.open(out_depth_asa.c_str());
		fdev << "# RES | ALL | MC | SC | SCP | SCNP" << endl;
		for (unsigned int i = 0; i < all_residues.size(); ++i){
			resname = all_residues[i];
			fdev << resname;
			fdev << " | " << D.residue(resname).depth_all() << " " << A.residue(resname).all_per();
			fdev << " | " << D.residue(resname).depth_MC() << " " << A.residue(resname).main_per();
			fdev << " | " << D.residue(resname).depth_SC() << " " << A.residue(resname).side_per();
			fdev << " | " << D.residue(resname).depth_SCP() << " " << A.residue(resname).polar_per();
			fdev << " | " << D.residue(resname).depth_SCNP() << " " << A.residue(resname).nonP_per();
			fdev << endl;
		} // end for
		fdev.close(); fdev.clear();
	} // end if

	return 0;
} // end main
