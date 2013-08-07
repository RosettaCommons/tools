// Shrake-Rupley algorithm to calculate accessible surface area
// version 1.0
// Author: KuanPern Tan
// Date: 28 Dec 2010
// Institute: Bioinformatics Institute A*STAR, Singapore

# include <iostream>
# include <math.h>
# include <string>
# include <fstream>
# include <set>
# include <map>
# include "easystring.h"
# include "PDB.h"
using namespace std;

# include "self_distmatrix.h"
# define PI 3.14159265359

/*
float vdw(string s){	// vdw radii of biological atom, assume all-atom model. data extracted from MODELLER::radii14.lib
	if (s == "H"){
		 return 0.80;
	} else if (s == "C"){
		 return 1.70;
	} else if (s == "N"){
		 return 1.60;
	} else if (s == "O"){
		 return 1.55;
	} else if (s == "P"){
		 return 1.80;
	} else if (s == "S"){
		 return 1.90;
	} else if (s == "ZN"){
		 return 0.74;
	} else if (s == "FE"){
		 return 0.65;
	} else if (s == "MG"){
		 return 0.72;
	} else if (s == "NE"){
		 return 1.60;
	} else {
		return 1.80; // assume Carbon atom
	} // end if
} // end vdw
*/

float vdw(string s){	// vdw radii of biological atom, assume all-atom model. data extracted from Cambridge Structural Database
	if (s == "H"){
		return 1.09;
	} else if (s == "N"){
		return 1.55;
	} else if (s == "CU"){
		return 1.40;
	} else if (s == "CL"){
		return 1.75;
	} else if (s == "C"){
		return 1.70;
	} else if (s == "O"){
		return 1.52;
	} else if (s == "I"){	
		return 1.98;
	} else if (s == "P"){
		return 1.80;
	} else if (s == "B"){
		return 2.00;
	} else if (s == "BR"){
		return 1.85;
	} else if (s == "S"){
		return 1.80;
	} else if (s == "SE"){
		return 1.90;
	} else if (s == "F"){
		return 1.47;
	} else if (s == "FE"){
		return 1.80;
	} else if (s == "K"){
		return 2.75;
	} else if (s == "MN"){
		return 2.00;
	} else if (s == "MG"){
		return 1.73;
	} else if (s == "ZN"){
		return 1.39;
	} else if (s == "HG"){
		return 1.55;
	} else if (s == "XE"){
		return 2.16;
	} else if (s == "AU"){
		return 1.66;
	} else if (s == "LI"){
		return 1.82;
	} else {
		return 1.8;
	} // end if
} // end vdw


class Point{	// a point in which collective of them form the surface
	private:
	bool removal;
	float _x, _y, _z;

	public:
	Point(){} // end Point
	~Point(){}

	void initialize(float x_in, float y_in, float z_in){
		_x = x_in; _y = y_in; _z = z_in;
		removal = 0;
	} // end initialize()

	void remove(){
		removal = 1;
	} // end remove

	bool status(){
		return !(removal);
	} // end removal

	float x(){
		return _x;
	} // end x
	float y(){
		return _y;
	} // end x
	float z(){
		return _z;
	} // end x

}; // end class

float point_distance(Point A, Point B){	// distance between two point
	float dx = (A.x() - B.x());
	float dy = (A.y() - B.y());
	float dz = (A.z() - B.z());
	return sqrt(dx*dx + dy*dy + dz*dz);
} // end distance

class ATOM{		// atom
	private:
	float x, y, z, d_area, area, clash_dist, probe_radius, d;
	bool is_mainchain, next_point;
	unsigned int _serial, survival;
	string _element;
	string _residue;
	set <int> neighbour;
	Point _position;
	vector <Point> Points; // position of atom itself

	public:
	ATOM(){
		is_mainchain = 0; survival = 0; area = 0;
	} // end ATOM
	~ATOM(){}

	// set functions
	void initialize(unsigned int serial_in, string element, string name, string residue_in, float _probe_radius, float _x, float _y, float _z){
		x = _x; y = _y; z =_z;
		_serial = serial_in;
		probe_radius = _probe_radius;
		_position.initialize(x, y ,z);
		_residue = residue_in;
		_element = element;
		if ((name=="N") || (name=="CA") || (name=="C") || (name=="O")){
			is_mainchain = 1;
		} // end if
	} // end initialize

	void add_neighbour(unsigned int s){
		neighbour.insert(s);
	} // end add_neighbour

	// get functions
	string residue(){
		return _residue;
	} // end residue

	string element(){
		return _element;
	} // end element

	Point position(){
		return _position;
	} // end position

	int serial(){
		return _serial;
	} // end serial

	// Core computation starts here
	void build_sphere(vector <float> xs, vector <float> ys, vector <float> zs){
		d = probe_radius + vdw(_element);
		d_area = 4*PI*(d*d) / xs.size(); // fractional area represented by each point
		for (unsigned int i = 0; i < xs.size(); ++i){
			Points.push_back(Point());
			Points[Points.size()-1].initialize(x + xs[i]*d, y + ys[i]*d, z + zs[i]*d);
		} // end for
	} // end build_sphere

	int n_point(){
		return Points.size();
	} // end int

	void remove_overlap(map <int, ATOM> &L){
		for (unsigned int i = 0; i < Points.size(); ++i){ // for all point
			for (set <int>::iterator iter = neighbour.begin(); iter != neighbour.end(); ++iter){ // search all neighbours
				clash_dist = vdw(L[*iter].element()) + probe_radius; // clash distance is dependent on neighbour type
				if (point_distance(Points[i], L[*iter].position()) < clash_dist){ // if clashed
					Points[i].remove(); // remove point
					break; // next point
				} // end if
			} // end for
		} // end for
	} // end remove overlap

	// update atomic area to residue area
	void update_residue(float &ASA, float &ASA_MC, float &ASA_SC, float &ASA_SCP, float &ASA_SCNP){
		survival = 0;
		for (unsigned int i = 0; i < Points.size(); ++i){
			if (Points[i].status() == 1){
				survival = survival + 1;
			} // end if
		} // end for
		area = survival * d_area;
		ASA = ASA + area;
		if (is_mainchain == 1){
			ASA_MC = ASA_MC + area;
		} else {
			ASA_SC = ASA_SC + area;
			if (_element == "C"){
				ASA_SCNP = ASA_SCNP + area;
			} else {
				ASA_SCP = ASA_SCP + area;
			} // end if
		} // end if
	} // end update_residue

	float access_area(){
		return area;
	} // end access

}; //end class

float div(float a, float b){
	if (b == 0){
		return 0;
	} else {
		return 100 * a / b;
	} // end if
} // end div


int main(int argc,char* argv[]){
	// step 0: membership, build sphere
	if (argc != 5){
		cerr << "syntax: ./ASA input.pdb resolution probe-size output" << endl;
		exit(1);
	} // end if
	vector <string> keywords; keywords.push_back("ATOM");
	PDB mdl = PDB(argv[1], keywords); // read PDB (need to add hydrogen before computation, use 'babel -h' to preprocess the structure first)
	cout << "structure: " << argv[1] << endl;
	unsigned int n = atoi(argv[2]); // determine resolution (number of point in sphere)
	float probe_radius = atof(argv[3]); // 1.4 for solvent type = water
	cout << "probe_radius of probe: " << probe_radius << endl;
	cout << "output file: " << argv[4] << endl;

	// generate a sample sphere surface defined with points from center a the fixed probe_radius (Golden Section Spiral Algorithm)
	vector <float> xs, ys, zs;
	vector <float> sphere_pts;
	vector < vector <float> > sphere;

	float y, r, phi;
	float inc = 2.4;
	float offset = 2.0 / n;
	for (unsigned int k = 0; k < n; k++){
		y = k * offset - 1 + (offset / 2.0);
		r = sqrt(1 - y*y);
		phi = k * inc;
		xs.push_back(cos(phi)*r);
		ys.push_back(y);
		zs.push_back(sin(phi)*r);
	} // end for

	// variable declaration
	map <int, string> atom2residue;
	map <string, string> residue2type;
	map <int, ATOM> L;
	map <string, float> ASA, ASA_MC, ASA_SC, ASA_SCP, ASA_SCNP; // output containers
	unsigned int this_atom; string residue;
	vector <string> sorted_residue_list;
	vector <int> sorted_atom_list;

	// read atomic data, build residue and build spheres around atoms
	for (unsigned int i = 0; i < mdl.size(); ++i){
		this_atom = mdl.serial(i);
		sorted_atom_list.push_back(this_atom);
		residue = mdl.chainID(i)+":"+mdl.resSeq(i);
		atom2residue[this_atom] = residue;
		residue2type[residue] = mdl.resName(i);

		if (is_element(residue, sorted_residue_list) == 0){ // build new residue as appropriate
			sorted_residue_list.push_back(residue);
			ASA[residue] = 0; ASA_MC[residue] = 0; ASA_SC[residue] = 0; ASA_SCP[residue] = 0; ASA_SCNP[residue] = 0;
		} // end if
		
		// build an atom
		L[this_atom] = ATOM();
		L[this_atom].initialize(mdl.serial(i), mdl.element(i), mdl.name(i), residue, probe_radius, mdl.x(i), mdl.y(i), mdl.z(i));
		// build sphere around the atom
		L[this_atom].build_sphere(xs, ys, zs);
	} // end for

	// step 1: find neighbouring atom for every atom
	float neighbour_distance = 2*(probe_radius + 2.75); // largest vdw of biomolecule = K, R(K) = 2.75 (most likely an overkill for most structure though)
	string distmap_file = string(argv[1])+".distmap";
	self_distmatrix C = self_distmatrix();
	C.compute(argv[1], neighbour_distance, distmap_file);
	ifstream fdist; fdist.open(distmap_file.c_str());
	unsigned int atom1, atom2;
	string line;
	vector <string> bufferline;
	while (!fdist.eof()){
		getline(fdist, line);
		if (line.size() == 0){
			continue;
		} // end if

		bufferline = split(line, " ");
		atom1 = atoi(bufferline[5].c_str());
		atom2 = atoi(bufferline[7].c_str());
		if ( atof(bufferline[10].c_str()) < neighbour_distance){
			if (atom1 != atom2){ // self-neighbour disallowed
				L[atom1].add_neighbour(atom2);
				L[atom2].add_neighbour(atom1);
			} // end if
		} // end if
	} // end while
	fdist.close(); fdist.clear();

	// step 2: main computation, remove overlapping region
	for(map <int, ATOM>::iterator iter = L.begin(); iter != L.end(); ++iter){
		(*iter).second.remove_overlap(L);
	} // end for

	// step 3: update atomic ASA to residue ASA
	string res;
	for(map <int, ATOM>::iterator iter = L.begin(); iter != L.end(); ++iter){
		res = (*iter).second.residue();
		(*iter).second.update_residue(ASA[res], ASA_MC[res], ASA_SC[res], ASA_SCP[res], ASA_SCNP[res]);
	} // end for

	// percentage conversion
	map <string, float> std_ALL, std_SCNP, std_SCP, std_SC, std_MC; // value extracted (statistically) from MODELLER 9v7
	std_ALL["ALA"] = 109.18 ; std_SCNP["ALA"] =  70.80 ; std_SCP["ALA"] =   0.00 ; std_SC["ALA"] =  70.80 ; std_MC["ALA"] = 38.36 ;
	std_ALL["ARG"] = 243.03 ; std_SCNP["ARG"] =  79.70 ; std_SCP["ARG"] = 125.53 ; std_SC["ARG"] = 205.22 ; std_MC["ARG"] = 37.81 ;
	std_ALL["ASN"] = 144.93 ; std_SCNP["ASN"] =  45.01 ; std_SCP["ASN"] =  61.88 ; std_SC["ASN"] = 106.90 ; std_MC["ASN"] = 38.02 ;
	std_ALL["ASP"] = 141.31 ; std_SCNP["ASP"] =  49.50 ; std_SCP["ASP"] =  53.78 ; std_SC["ASP"] = 103.28 ; std_MC["ASP"] = 38.02 ;
	std_ALL["CYS"] = 136.15 ; std_SCNP["CYS"] =  98.25 ; std_SCP["CYS"] =   0.00 ; std_SC["CYS"] =  98.25 ; std_MC["CYS"] = 37.78 ;
	std_ALL["GLN"] = 180.60 ; std_SCNP["GLN"] =  52.53 ; std_SCP["GLN"] =  90.25 ; std_SC["GLN"] = 142.78 ; std_MC["GLN"] = 37.82 ;
	std_ALL["GLU"] = 174.46 ; std_SCNP["GLU"] =  61.20 ; std_SCP["GLU"] =  75.42 ; std_SC["GLU"] = 136.63 ; std_MC["GLU"] = 37.82 ;
	std_ALL["GLY"] =  81.20 ; std_SCNP["GLY"] =  33.14 ; std_SCP["GLY"] =   0.00 ; std_SC["GLY"] =  33.14 ; std_MC["GLY"] = 48.05 ;
	std_ALL["HIS"] = 183.87 ; std_SCNP["HIS"] =  96.38 ; std_SCP["HIS"] =  51.88 ; std_SC["HIS"] = 148.24 ; std_MC["HIS"] = 35.59 ;
	std_ALL["ILE"] = 177.21 ; std_SCNP["ILE"] = 139.74 ; std_SCP["ILE"] =   0.00 ; std_SC["ILE"] = 139.74 ; std_MC["ILE"] = 36.70 ;
	std_ALL["LEU"] = 180.34 ; std_SCNP["LEU"] = 142.51 ; std_SCP["LEU"] =   0.00 ; std_SC["LEU"] = 142.51 ; std_MC["LEU"] = 37.73 ;
	std_ALL["LYS"] = 204.09 ; std_SCNP["LYS"] = 118.39 ; std_SCP["LYS"] =  47.88 ; std_SC["LYS"] = 166.28 ; std_MC["LYS"] = 37.82 ;
	std_ALL["MET"] = 197.05 ; std_SCNP["MET"] = 159.22 ; std_SCP["MET"] =   0.00 ; std_SC["MET"] = 159.22 ; std_MC["MET"] = 37.78 ;
	std_ALL["PHE"] = 199.72 ; std_SCNP["PHE"] = 164.48 ; std_SCP["PHE"] =   0.00 ; std_SC["PHE"] = 164.48 ; std_MC["PHE"] = 35.10 ;
	std_ALL["PRO"] = 140.92 ; std_SCNP["PRO"] = 120.12 ; std_SCP["PRO"] =   0.00 ; std_SC["PRO"] = 120.12 ; std_MC["PRO"] = 20.80 ;
	std_ALL["SER"] = 118.23 ; std_SCNP["SER"] =  48.19 ; std_SCP["SER"] =  31.44 ; std_SC["SER"] =  79.63 ; std_MC["SER"] = 38.60 ;
	std_ALL["THR"] = 142.04 ; std_SCNP["THR"] =  76.21 ; std_SCP["THR"] =  28.04 ; std_SC["THR"] = 104.24 ; std_MC["THR"] = 37.78 ;
	std_ALL["TRP"] = 248.07 ; std_SCNP["TRP"] = 185.87 ; std_SCP["TRP"] =  23.65 ; std_SC["TRP"] = 209.55 ; std_MC["TRP"] = 38.28 ;
	std_ALL["TYR"] = 212.28 ; std_SCNP["TYR"] = 134.88 ; std_SCP["TYR"] =  42.21 ; std_SC["TYR"] = 177.11 ; std_MC["TYR"] = 35.02 ;
	std_ALL["VAL"] = 153.32 ; std_SCNP["VAL"] = 115.85 ; std_SCP["VAL"] =   0.00 ; std_SC["VAL"] = 115.85 ; std_MC["VAL"] = 37.24 ;


	// step 4: print output
	ofstream fpre; fpre.open(argv[4]);
	fpre << "#       Res   Res      All atoms    Non P side    Polar Side    Total Side    Main Chain" << endl;
	fpre << "#       Num  type      Sum  Per.     Sum  Per.     Sum  Per.     Sum  Per.     Sum  Per." << endl;
	fpre.close(); fpre.clear();
	FILE* fout; fout = fopen(argv[4], "a");
	for (unsigned int i = 0; i < sorted_residue_list.size(); ++i){
		res = sorted_residue_list[i];
		bufferline = split(res,":");
		fprintf (fout, "%6s %4s %5s %1s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", "ACCESS", bufferline[1].c_str(), residue2type[res].c_str(), bufferline[0].c_str(), ASA[res], div(ASA[res] , std_ALL[residue2type[res]]), ASA_SCNP[res], div(ASA_SCNP[res] , std_SCNP[residue2type[res]]), ASA_SCP[res], div(ASA_SCP[res] , std_SCP[residue2type[res]]), ASA_SC[res], div(ASA_SC[res] , std_SC[residue2type[res]]), ASA_MC[res], div(ASA_MC[res] , std_MC[residue2type[res]]));
	} // end for
	fclose(fout);

	remove(distmap_file.c_str()); // remove temporary file
	return 0;
} // end main 
