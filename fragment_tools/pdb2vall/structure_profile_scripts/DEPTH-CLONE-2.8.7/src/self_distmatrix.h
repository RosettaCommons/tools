// c++ version of distance matrix script
// Author: Tan Kuan Pern

# include <iostream>
# include <string>
# include <map>
# include <math.h>
# include <fstream>
# include "PDB.h"
using namespace std;

std::vector <int> build_neighbour_valid(unsigned int i, int x_box, int y_box, int box_no){
	vector <int> output;
	int xy_box = x_box*y_box;
	int neighbour_buffer[27] = {i, i+1, i+x_box-1, i+x_box, i+x_box+1, i-1, i-x_box-1, i-x_box, i-x_box+1, i-xy_box-x_box-1, i-xy_box-x_box, i-xy_box-x_box+1, i-xy_box-1, i-xy_box, i-xy_box+1, i-xy_box+x_box-1, i-xy_box+x_box, i-xy_box+x_box+1, i+xy_box-x_box-1, i+xy_box-x_box, i+xy_box-x_box+1, i+xy_box-1, i+xy_box, i+xy_box+1, i+xy_box+x_box-1, i+xy_box+x_box, i+xy_box+x_box+1};
	for (int j = 0; j < 27; ++j){
		if (neighbour_buffer[j] >= 0){
			if (neighbour_buffer[j] < box_no+1){
				output.push_back(neighbour_buffer[j]);
			} // end if
		} // end if
	} // end for
	return output;
} // end neighbour


class self_distmatrix{
	private:

	public:
	self_distmatrix(){}
	~self_distmatrix(){}

	void compute(string input_fname1, float cut_off, string output_fname){
		// read pdb file	
		vector <string> keywords;
		keywords.push_back("ATOM"); keywords.push_back("HETATM");
		PDB pdb = PDB(input_fname1, keywords);

		// declare output file
		ofstream outfile;
		outfile.open(output_fname.c_str());

		// parameters for boxes
		int x_box, y_box, z_box, box_no;

		// max and min of each dimension
		double max_x = pdb.x(0); double max_y = pdb.y(0); double max_z = pdb.z(0);
		double min_x = pdb.x(0); double min_y = pdb.y(0); double min_z = pdb.z(0);

		for (unsigned int i = 0; i < pdb.size(); i++){
			// x dimension
			if (pdb.x(i) > max_x){
				max_x = pdb.x(i);
			} else if (pdb.x(i) < min_x){
				min_x = pdb.x(i);
			} // end if
			if (pdb.y(i) > max_y){
				max_y = pdb.y(i);
			} else if (pdb.y(i) < min_y){
				min_y = pdb.y(i);
			} // end if
			if (pdb.z(i) > max_z){
				max_z = pdb.z(i);
			} else if (pdb.z(i) < min_z){
				min_z = pdb.z(i);
			} // end if
		} // end for

		x_box = int((max_x - min_x)/cut_off) + 1;
		y_box = int((max_y - min_y)/cut_off) + 1;
		z_box = int((max_z - min_z)/cut_off) + 1;
		box_no = x_box*y_box*z_box;

		vector <int> box_content[box_no+1];
		// compute cell membership
		int x_index, y_index, z_index, cell_index;
		for (unsigned int i = 0; i < pdb.size(); i++){
			x_index = int((pdb.x(i) - min_x) / cut_off) + 1;
			y_index = int((pdb.y(i) - min_y) / cut_off) + 1;
			z_index = int((pdb.z(i) - min_z) / cut_off) + 1;
			cell_index = x_index + (y_index-1)*(x_box) + (z_index-1)*(x_box*y_box);
			box_content[cell_index].push_back(i);
		} // end for

		// declare neighbourhoods
		vector <int> neighbour[box_no+1];
		for (int i = 0; i < box_no+1; i++){
			neighbour[i] = build_neighbour_valid(i, x_box, y_box, box_no);
		} // end for

		// compute inter-distances within cell
		string compare1, compare2;
		double x1, y1, z1, x2, y2, z2, distance;
		string resName1 = ""; string resName2 = "";

		// compute inter-distances between cells
		for (int i = 1; i < box_no + 1; i++) { // for every box
			for (unsigned int n = 0; n < neighbour[i].size(); n++){ // select its neighbour
				if (neighbour[i][n] < box_no + 1 ) {
					for (unsigned int atom1 = 0; atom1 < box_content[i].size(); atom1++) {	//every atom from box
						for (unsigned int atom2 = 0; atom2 < box_content[neighbour[i][n]].size(); atom2++){ // to atom from neighbour box
							if (1){
								x1 = pdb.x(box_content[i][atom1]);
								y1 = pdb.y(box_content[i][atom1]);
								z1 = pdb.z(box_content[i][atom1]);
								x2 = pdb.x(box_content[neighbour[i][n]][atom2]);
								y2 = pdb.y(box_content[neighbour[i][n]][atom2]);
								z2 = pdb.z(box_content[neighbour[i][n]][atom2]);
								distance = dist(x1,y1,z1,x2,y2,z2);

								// write output
								if (distance <= cut_off){
									resName1 = pdb.resName(box_content[i][atom1]);
									resName2 = pdb.resName(box_content[neighbour[i][n]][atom2]);
	
									compare1 = pdb.chainID(box_content[i][atom1])+":"+pdb.resSeq(box_content[i][atom1]);
									compare2 = pdb.chainID(box_content[neighbour[i][n]][atom2])+":"+pdb.resSeq(box_content[neighbour[i][n]][atom2]);
	
									if (compare1 <= compare2) {
										outfile << compare1 << ' ' << compare2 << ' ' << resName1 << ' ' << resName2 << " ( " << pdb.serial(box_content[i][atom1]) << " , "  << pdb.serial(box_content[neighbour[i][n]][atom2]) << " ) : " << distance << endl;
									} else if (compare1 > compare2) {
										outfile << compare2 << ' ' << compare1 << ' ' << resName2 << ' ' << resName1 << " ( " << pdb.serial(box_content[neighbour[i][n]][atom2]) << " , "  << pdb.serial(box_content[i][atom1]) << " ) : " << distance << endl;
									} // end if
								} // end if
							} // end if
						} // end for
					} // end for
				} // end if
			} // end for
		} // end for		
		outfile.close();
	} // end compute
}; // end class
