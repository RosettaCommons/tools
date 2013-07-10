// core computing unit of depth program

std::pair <float,float> rotate(float x, float y, float theta){ //rotation
	float xr, yr;
	xr = x*cos(theta) - y*sin(theta);
	yr = x*sin(theta) + y*cos(theta);
	return std::make_pair(xr,yr);
} // end rotate

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

std::vector <int> build_neighbour_2(unsigned int i, int x_box, int y_box, int box_no){
	vector <int> output;
	int xy_box = x_box*y_box; 
	int neighbour_buffer[98] = {i+1+2*x_box, i+1+2*x_box+2*xy_box, i+1+2*x_box+xy_box, i+1+2*x_box-2*xy_box, i+1+2*x_box-xy_box, i+1+2*xy_box, i+1+x_box+2*xy_box, i+1+x_box-2*xy_box, i+1-2*x_box, i+1-2*x_box+2*xy_box, i+1-2*x_box+xy_box, i+1-2*x_box-2*xy_box, i+1-2*x_box-xy_box, i+1-2*xy_box, i+1-x_box+2*xy_box, i+1-x_box-2*xy_box, i+2, i+2*x_box, i+2*x_box+2*xy_box, i+2*x_box+xy_box, i+2*x_box-2*xy_box, i+2*x_box-xy_box, i+2*xy_box, i+2+2*x_box, i+2+2*x_box+2*xy_box, i+2+2*x_box+xy_box, i+2+2*x_box-2*xy_box, i+2+2*x_box-xy_box, i+2+2*xy_box, i+2+x_box, i+2+x_box+2*xy_box, i+2+x_box+xy_box, i+2+x_box-2*xy_box, i+2+x_box-xy_box, i+2+xy_box, i+2-2*x_box, i+2-2*x_box+2*xy_box, i+2-2*x_box+xy_box, i+2-2*x_box-2*xy_box, i+2-2*x_box-xy_box, i+2-2*xy_box, i+2-x_box, i+2-x_box+2*xy_box, i+2-x_box+xy_box, i+2-x_box-2*xy_box, i+2-x_box-xy_box, i+2-xy_box, i+x_box+2*xy_box, i+x_box-2*xy_box, i-1+2*x_box, i-1+2*x_box+2*xy_box, i-1+2*x_box+xy_box, i-1+2*x_box-2*xy_box, i-1+2*x_box-xy_box, i-1+2*xy_box, i-1+x_box+2*xy_box, i-1+x_box-2*xy_box, i-1-2*x_box, i-1-2*x_box+2*xy_box, i-1-2*x_box+xy_box, i-1-2*x_box-2*xy_box, i-1-2*x_box-xy_box, i-1-2*xy_box, i-1-x_box+2*xy_box, i-1-x_box-2*xy_box, i-2, i-2*x_box, i-2*x_box+2*xy_box, i-2*x_box+xy_box, i-2*x_box-2*xy_box, i-2*x_box-xy_box, i-2*xy_box, i-2+2*x_box, i-2+2*x_box+2*xy_box, i-2+2*x_box+xy_box, i-2+2*x_box-2*xy_box, i-2+2*x_box-xy_box, i-2+2*xy_box, i-2+x_box, i-2+x_box+2*xy_box, i-2+x_box+xy_box, i-2+x_box-2*xy_box, i-2+x_box-xy_box, i-2+xy_box, i-2-2*x_box, i-2-2*x_box+2*xy_box, i-2-2*x_box+xy_box, i-2-2*x_box-2*xy_box, i-2-2*x_box-xy_box, i-2-2*xy_box, i-2-x_box, i-2-x_box+2*xy_box, i-2-x_box+xy_box, i-2-x_box-2*xy_box, i-2-x_box-xy_box, i-2-xy_box, i-x_box+2*xy_box, i-x_box-2*xy_box};
	for (int j = 0; j < 98; ++j){
		if (neighbour_buffer[j] >= 0){
			if (neighbour_buffer[j] < box_no+1){
				output.push_back(neighbour_buffer[j]);
			} // end if
		} // end if
	} // end for
	return output;
} // end neighbour

class depth_run{
	private:
	float radius_w, cavity, degree, displacement, cell_spacing, clash_dist, unit_box_x, unit_box_y, unit_box_z, x_cog, y_cog, z_cog, xw_cog, yw_cog, zw_cog;
	unsigned int survival_n, index, atom_number, solvent_number, sol_number;
	int rot1, rot2, modulus, sol_index;
	float x_tmp, y_tmp, z_tmp;
	float xsol_one[3], ysol_one[3], zsol_one[3];
	int cell_index_one[3];
	std::pair<float,float> rotated; float r[3]; // transformation buffer memory
	map <string, int> dimension;
	string rot_axis, job_id, min_sol_name;

	public:
	depth_run(){ // constructor
		dimension["x"] = 0; dimension["y"] = 1; dimension["z"] = 2;
	} // end depth_run
	~depth_run(){ //destructor
	} // end ~depth_run

	void set_label(unsigned int index_in, string job_name_in){
		job_id = job_name_in;
		index = index_in;
	} // end set

	void initialize_box(float cell_spacing_in){
		cell_spacing = cell_spacing_in;
	} // end set

	void set_atom_number(unsigned int n){
		atom_number = n;
	} // end set

	void set_solvent_number(unsigned int n){
		solvent_number = n;
	} // end set

	void set_parameter(unsigned int survival_n_in, float cavity_in, string rot_axis_in, float degree_in, float displacement_in, float radius_w_in, float clash_dist_in, string min_sol_name_in, float unit_box_x_in, float unit_box_y_in, float unit_box_z_in){
		survival_n = survival_n_in;
		cavity = cavity_in;
		rot_axis = rot_axis_in; degree = degree_in; displacement = displacement_in;
		radius_w = radius_w_in; clash_dist = clash_dist_in;
		min_sol_name = min_sol_name_in;
		unit_box_x = unit_box_x_in; unit_box_y = unit_box_y_in; unit_box_z = unit_box_z_in;
	} // end set

	void get_depth(float depth[], float x[], float y[], float z[], float xw[], float yw[], float zw[]){

		// declare memory for depth value
		for (unsigned int i = 0; i < atom_number; ++i){
			*(depth+i) = 9999;
		} // end for

		ofstream flog; flog.open((job_id+"-iter"+num2string(index)).c_str()); // open log file		

		// --- ROTATION ---
		// step 1: convert coordinate relative to cog at origin (solute)
		x_cog = 0; y_cog = 0; z_cog = 0;
		for (unsigned int i = 0; i < atom_number; ++i){
			x_cog = x_cog + x[i]; y_cog = y_cog + y[i]; z_cog = z_cog + z[i];
		} // end for
		x_cog = x_cog / atom_number; y_cog = y_cog / atom_number; z_cog = z_cog / atom_number;
		for (unsigned int i = 0; i < atom_number; ++i){
			x[i] = x[i] - x_cog; y[i] = y[i] - y_cog; z[i] = z[i] - z_cog;
		} // end for
		x_cog = 0; y_cog = 0; z_cog = 0;
		// step 2: rotation
		rot1 = dimension[rot_axis.substr(0,1)]; rot2 = dimension[rot_axis.substr(1,1)];
		flog << "cavity size: " << cavity << endl;
		flog << "action: rotate along " << rot_axis << "-axis by " << degree << " degree" << endl;
		degree = degree*(3.14159/180.0); flog << "action: translate along x-axis by " << displacement << "Å" << endl;
		// report
		flog << "rotate solute atoms ";
		for (unsigned int i = 0; i < atom_number; ++i){
			r[0] = x[i]; r[1] = y[i]; r[2] = z[i];
			rotated = rotate(r[rot1],r[rot2],degree);
			r[rot1] = rotated.first; r[rot2] = rotated.second;
			x[i] = r[0]; y[i] = r[1]; z[i] = r[2];
		} // end for
		flog << "done" << endl;

		// --- BUILD SOLVENT BOX ---
		// step 1: inspect size of solute molecule
		float x_max = -9999, y_max = -9999, z_max = -9999;
		float x_min =  9999, y_min =  9999, z_min =  9999;
		float x_mid, y_mid, z_mid;
		for (unsigned int i = 0; i < atom_number; ++i){
			x_max = max(x_max, x[i]); y_max = max(y_max, y[i]); z_max = max(z_max, z[i]);
			x_min = min(x_min, x[i]); y_min = min(y_min, y[i]); z_min = min(z_min, z[i]);
		} // end for
		x_mid = (x_max + x_min)/2.0; y_mid = (y_max + y_min)/2.0; z_mid = (z_max + z_min)/2.0;
		// step 2.1 : calculate the require size of solvent box
		float bufferspace = 4*(radius_w*1.2) + clash_dist; // space for two solvents + worst-case-of-clash-buffer, 20% buffer space for displacement
		float x_length = x_max - x_min + bufferspace + 2*fabs(displacement)*1.2; // x-axis is going to move
		float y_length = y_max - y_min + bufferspace;
		float z_length = z_max - z_min + bufferspace;
		flog << "bufferspace: " << bufferspace << "Å " << endl;
		flog << "lengths: " << x_length << "Å " << y_length << "Å " << z_length << "Å " << endl;
		// step 2.2 : calculate number of block of unit solvent box needed
		int block_x = int(x_length/unit_box_x) + 1;
		int block_y = int(y_length/unit_box_y) + 1;
		int block_z = int(z_length/unit_box_z) + 1;
		flog << "solvent box dimension: " << block_x << " " << block_y << " " << block_z << endl; 
		// step 2.3 : calculate the size of the solvent box
		// step 2.3.1 : all position in unit-solvent box > 0
		float max_xw = -9999; float max_yw = -9999; float max_zw = -9999;
		float min_xw =  9999; float min_yw =  9999; float min_zw =  9999;
		for (unsigned int i = 0; i < solvent_number; ++i){
			min_xw = min(min_xw, xw[i]); min_yw = min(min_yw, yw[i]); min_zw = min(min_zw, zw[i]);
		} // end for
		for (unsigned int i = 0; i < solvent_number; ++i){
			xw[i] = xw[i] + min_xw; yw[i] = yw[i] + min_yw; zw[i] = zw[i] + min_zw;  
			max_xw = max(max_xw, xw[i]); max_yw = max(max_yw, yw[i]); max_zw = max(max_zw, zw[i]);
			min_xw = min(min_xw, xw[i]); min_yw = min(min_yw, yw[i]); min_zw = min(min_zw, zw[i]);
		} // end for
		float max_x = (block_x-1)*unit_box_x + max_xw;
		float max_y = (block_y-1)*unit_box_y + max_yw;
		float max_z = (block_z-1)*unit_box_z + max_zw;
		float min_x = min_xw;
		float min_y = min_yw;
		float min_z = min_zw;
		flog << "x max,min: " << max_x << " " << min_x << endl;
		flog << "y max,min: " << max_y << " " << min_y << endl;
		flog << "z max,min: " << max_z << " " << min_z << endl;
		// step 3 : build cell-list
		int x_index, y_index, z_index, cell_index;
		int x_box = int((max_x - min_x)/cell_spacing) + 1;
		int y_box = int((max_y - min_y)/cell_spacing) + 1;
		int z_box = int((max_z - min_z)/cell_spacing) + 1;
		int box_no = x_box*y_box*z_box;
		flog << "cell-length: " << cell_spacing << " A" << endl;
		flog << "cell dimension: " << x_box << " " << y_box << " " << z_box << endl;
		vector <int> solute_box[box_no+1], solvent_box[box_no+1];
		for (int i = 0; i < box_no+1; i++){
			solute_box[i].clear(); solvent_box[i].clear();
		} // end for
		// step 2.4 : add solvents into box
		sol_index = -1;
		sol_number = block_x*block_y*block_z*solvent_number;
		float* xsol = new float[sol_number];
		float* ysol = new float[sol_number];
		float* zsol = new float[sol_number];
		flog << "add solvent ... ";
		for (int i = 0; i < block_x; ++i){ // x-
			for (int j = 0; j < block_y; ++j){ // y-
				for (int k = 0; k < block_z; ++k){ // z-
					for (unsigned int p = 0; p < solvent_number; ++p){
						modulus = p % 3; // count atom from solvent molecule 
						x_tmp = xw[p] + i*unit_box_x;
						y_tmp = yw[p] + j*unit_box_y;
						z_tmp = zw[p] + k*unit_box_z;
						xsol_one[modulus] = x_tmp;
						ysol_one[modulus] = y_tmp; 
						zsol_one[modulus] = z_tmp;
						// step 3.1 : assign cell index to solvents
						x_index = int((x_tmp - min_x) / cell_spacing) + 1;
						y_index = int((y_tmp - min_y) / cell_spacing) + 1;
						z_index = int((z_tmp - min_z) / cell_spacing) + 1;
						cell_index = x_index + (y_index-1)*(x_box) + (z_index-1)*(x_box*y_box);
						cell_index_one[modulus] = cell_index;

						if (modulus == 2){ // all three atoms
							for (int r = 0; r < 3; ++r){
								sol_index = sol_index + 1;
								solvent_box[cell_index_one[0]].push_back(sol_index); // hydrogens follow oxygen
								*(xsol+sol_index) = xsol_one[r]; 
								*(ysol+sol_index) = ysol_one[r];
								*(zsol+sol_index) = zsol_one[r];
							} // end for
						} // end if
					} // end for
				} // end for
			} // end for
		} // end for
		flog << "done" << endl << "total number of surface solvent molecule added: " << sol_number/3 << endl;

		// ADD SOLUTE TO BOX	
		// step 1: mid-point of solvent box
		float sys_mid_x = (max_x + min_x) / 2.0;
		float sys_mid_y = (max_y + min_y) / 2.0;
		float sys_mid_z = (max_z + min_z) / 2.0;
		// step 2: put in solute
		float orient_x = - x_mid + sys_mid_x;
		float orient_y = - y_mid + sys_mid_y;
		float orient_z = - z_mid + sys_mid_z;
		for (unsigned int i = 0; i < atom_number; ++i){
			x[i] = x[i] + orient_x; y[i] = y[i] + orient_y; z[i] = z[i] + orient_z;
			x[i] = x[i] + displacement; // displacement
		} // end for
		// step 3: assign cell index to solute (need verification)
		flog << "assign cell-indices to solute solvent... " << endl;
		for (unsigned int i = 0; i < atom_number; ++i){
			x_index = int(( x[i] - min_x) / cell_spacing) + 1;
			y_index = int(( y[i] - min_y) / cell_spacing) + 1;
			z_index = int(( z[i] - min_z) / cell_spacing) + 1;
			cell_index = x_index + (y_index-1)*(x_box) + (z_index-1)*(x_box*y_box);
			solute_box[cell_index].push_back(i);
		} // end for
		flog << "assigned, build neighbour list" << endl;

		// --- DEPTH ALGORITHM ---
		// step 0: prepare memory for neighbour list 
		int this_neighbour;
		vector <int> neighbour_valid[box_no+1];
		set < std::pair <int, int> > inter_cells;
		set <int>* contact = new set <int> [sol_number];
		for (int i = 0; i < box_no + 1; ++i) {
			neighbour_valid[i].clear(); neighbour_valid[i] = build_neighbour_valid(i, x_box, y_box, box_no);
		} // end for
		for (unsigned int i = 0; i < sol_number; ++i){
			(*(contact + i)).clear();
		} // end for
		// step 1: build neighbour list
		for (int i = 1; i < box_no + 1; ++i){
			for (unsigned int n = 0; n < neighbour_valid[i].size(); ++n){
				if (i <= neighbour_valid[i][n]){
					inter_cells.insert(make_pair(i,neighbour_valid[i][n])); // choose inter-cells, prevent double-counting
				} else {
					inter_cells.insert(make_pair(neighbour_valid[i][n],i));
				} // end if
			} // end for
		} // end for
		// step 1.1: define neighbour atom for every atom (uses cell-list algorithm here)
		flog << "removing clashing atoms ... ";
		int sol_1, sol_2, solvent_box_1, solvent_box_2;
		float x1, y1, z1, x2, y2, z2, d;
		for (set < std::pair <int, int> >::iterator iter = inter_cells.begin(); iter != inter_cells.end(); ++iter){
			solvent_box_1 = (*iter).first; solvent_box_2 = (*iter).second; // betwee two cells
			if (solvent_box[solvent_box_1].size() > 0){	// cell 1 has atom(s)
				for (unsigned int n = 0; n < solvent_box[solvent_box_1].size(); ++n){
					sol_1 = solvent_box[solvent_box_1][n];
					x1 = *(xsol+sol_1); y1 = *(ysol+sol_1); z1 = *(zsol+sol_1); // every atom in cell 1
					if (solvent_box[solvent_box_2].size() > 0){ // cell 2 has atom(s)
						for (unsigned int m = 0; m < solvent_box[solvent_box_2].size(); ++m){
							sol_2 = solvent_box[solvent_box_2][m];
							x2 = *(xsol+sol_2); y2 = *(ysol+sol_2); z2 = *(zsol+sol_2); // every atom in cell 2
							if ((sol_1 % 3 == 0) && (sol_2 % 3 == 0) && (sol_1 != sol_2)){
								d = dist(x1,y1,z1,x2,y2,z2); 
								if (d <= cavity){ // compute distance, see if within neighbourhood
									(*(contact + sol_1)).insert(int(sol_2/3)*3); // record neighbour
									(*(contact + sol_2)).insert(int(sol_1/3)*3);
								} // end if
							} // end if
						} // end for
					} // end if
				} // end for
			} // end if
		} // end for
		inter_cells.clear();

		// step 2: remove clasing atom 
		int mol_index;
		vector <int> clash_atoms, clash_sols;
		vector <int> possible_cavity_sols;
		bool inclusion[sol_number];
		for (unsigned int i = 0; i < sol_number; ++i){
			inclusion[i] = 1;
		} // end for
		clash_atoms.clear();
		// step 2.1: compute solute to solvent distance (cell list algorithm)
		for (int i = 1; i < box_no + 1; ++i){ // for every box
			if (solute_box[i].size() > 0){ // that contains solute(s)
				for (unsigned int n = 0; n < neighbour_valid[i].size(); ++n){ // for every neighbouring 
					this_neighbour = neighbour_valid[i][n];
					if (solvent_box[this_neighbour].size() > 0){ // box that contain solvent(s)
						inter_cells.insert(make_pair(i,this_neighbour)); // record this pair of boxes (prevent double-counting)
					} // end if
				} // end for
			} // end if
		} // end for
		int solute_box_index, solvent_box_index;
		for (set < std::pair <int, int> >::iterator iter = inter_cells.begin(); iter != inter_cells.end(); ++iter){ // for
			solute_box_index = (*iter).first; solvent_box_index = (*iter).second; // the recorded pairs of boxes
			for (unsigned int n = 0; n < solute_box[solute_box_index].size(); ++n){ // every solute
				mol_index = solute_box[solute_box_index][n];
				x1 = x[mol_index]; y1 = y[mol_index]; z1 = z[mol_index];
				for (unsigned int solv = 0; solv < solvent_box[solvent_box_index].size(); ++solv){ // to every solvent
					sol_index = solvent_box[solvent_box_index][solv];
					x2 = *(xsol+sol_index); y2 = *(ysol+sol_index); z2 = *(zsol+sol_index);
					d = dist(x1,y1,z1,x2,y2,z2);
					if (d <= clash_dist){ // find distance and 
						clash_atoms.push_back(sol_index); // record solvent that clashes
					} // end if
				} // end for
			} // end for
		} // end for
		// step 2.2: remove clashing solvent (molecule-based (3 atom per molecule (assume H2O), remove if any one of them clashes)
		vector <int> buf;
		for (unsigned int i = 0; i < clash_atoms.size(); ++i){ // for clashing atoms, 
			buf.push_back(int(clash_atoms[i]/3)); // convert atom index to molecule index
		} // end for
		buf = non_redundant(buf); // non-redundant
		for (unsigned int j = 0; j < buf.size(); ++j){ // convert back to atom index
			clash_sols.push_back(buf[j]*3); // all atom in the same molecule recorded to be removed
			clash_sols.push_back(buf[j]*3+1);
			clash_sols.push_back(buf[j]*3+2);
		} // end for
		buf.clear();
		for (unsigned int i = 0; i < clash_sols.size(); ++i){
			inclusion[clash_sols[i]] = 0; // removal
		} // end for
		flog << "done" << endl << "removed " << clash_sols.size()/3 << " clashing solvent molecules" << endl;

		// step 3: remove cavity solvent
		flog << "computing cavity solvents ... " << endl;
		bool cavity_sols [sol_number];
		int possible_sol, cavity_count = 0, cavity_count_old = 0, layer = 0;
		for (unsigned int i = 0; i < sol_number; ++i){
			cavity_sols[i] = 0;
		} // end for
		while (1){ // iteratively search 
			layer = layer + 1;
			possible_cavity_sols.clear(); // clear memory for next round
			flog << "\tscanning layer " << layer << " ... ";
			for (unsigned int i = 0; i < clash_sols.size(); ++i){ // for contact neighbours of clashed atom (they are candidtates to be removed as cavity solvent)
				for (set <int> ::iterator iter = (*(contact + clash_sols[i])).begin(); iter != (*(contact + clash_sols[i])).end(); ++iter){
					possible_cavity_sols.push_back(int((*iter)/3)*3); // record possible cavity solvent molecule
					(*(contact + *iter)).erase(clash_sols[i]); // atom that was removed is excluded from neighbourhood
				} // end for
			} // end for
			possible_cavity_sols = non_redundant(possible_cavity_sols);
			flog << "possible cavity solvent number: " << possible_cavity_sols.size() << endl;

			clash_sols.clear(); // clear memory
			for (unsigned int i = 0; i < possible_cavity_sols.size(); ++i){ // for
				possible_sol = possible_cavity_sols[i]; // every possible cavity solvent
				if ( (*(contact + possible_sol)).size() < survival_n){ // if neighbour < survival number
					clash_sols.push_back(possible_sol); // record it as solvent to be removed
					cavity_sols[possible_sol  ] = 1; // and all its atom removed
					cavity_sols[possible_sol+1] = 1;
					cavity_sols[possible_sol+2] = 1;
				} // end if
			} // end for

			// count cavity solvents
			for (unsigned int i = 0; i < sol_number; ++i){
				if (cavity_sols[i] == 1){
					cavity_count = cavity_count + 1;
				} // end if
			} // end for

			if (cavity_count == cavity_count_old){ // if no more new cavity solvent
				break; // quit
			} else { // else
				cavity_count_old = cavity_count; // iterate
				cavity_count = 0;
			} // end if
		} // end while
		for (unsigned int i = 0; i < sol_number; ++i){ // remove all cavity atoms previously found
			if (cavity_sols[i] == 1){
				inclusion[i] = 0;
			} // end if
		} // end for
		flog << "done " << endl << "removed " << cavity_count/3 << " solvent molecules in cavities or distant solvent molecule" << endl;

		// step 4: computing atomic depth value
		flog << "computing depth value ... " << endl;
		// compute adaptive cell-list, from experience most atom will have depth value lower than 5 = first trial. 
		// atom with no solvent in its first layer of cell-list will progessively advance to second, third etc.

		// actual depth value computation, solute to solvent distances
		// first step:: compute only first shell, using formal cell-list (6A) (nearest solvent cannot be 6A more than protein surface (2.6 + 3.3), thus 6A defines the surface solvents
		// step 4.1: computing surface solvent
		flog << "computing surface solvent" << endl;
		map <int, int> contact_solvent;
		vector <int> surface_solvent;
		for (int i = 1; i < box_no + 1; ++i) { // for every cell
			if (solute_box[i].size() == 0){ // that has solute
				continue;
			} // end if 
			for (unsigned int p = 0; p < solute_box[i].size(); ++p){ // for 
				mol_index = solute_box[i][p];
				x1 = x[mol_index]; y1 = y[mol_index]; z1 = z[mol_index]; //every contained solute atom
				for (unsigned int n = 0; n < neighbour_valid[i].size(); ++n){ // for every neighbouring 
					this_neighbour = neighbour_valid[i][n]; // solvent cell
					for (unsigned int q = 0; q < solvent_box[this_neighbour].size(); ++q){ // every its
						sol_index = solvent_box[this_neighbour][q]; // solvent atom
						if (inclusion[sol_index] == 1){ // if not removed clash atom
							if (cavity_sols[sol_index] == 0){ // or cavity atom
								x2 = *(xsol+sol_index); y2 = *(ysol+sol_index); z2 = *(zsol+sol_index);
								d = dist(x1,y1,z1,x2,y2,z2); // compute their distance = a depth value
								if (d < depth[mol_index]){ // choose the smallest depth value
									depth[mol_index] = d;
									contact_solvent[mol_index] = sol_index; // record this solvent atom as surface solvent atom
								} else if (d < cavity){
									contact_solvent[mol_index] = sol_index; // record this solvent atom as surface solvent atom
								} // end if
							}// end if
						} // end if
					} // end for
				} // end for
			} // end for
		} // end for
		for (map <int, int>::iterator iter = contact_solvent.begin(); iter != contact_solvent.end(); ++iter){ // consolidate surface solvent atom recorded
			surface_solvent.push_back( int((*iter).second/3)*3   );
			surface_solvent.push_back( int((*iter).second/3)*3+1 );
			surface_solvent.push_back( int((*iter).second/3)*3+2 );
		} // end for
		surface_solvent = non_redundant(surface_solvent); // to non-redundant surface solvent molecule
		flog << "finished computing non-redundant surface solvent: " << surface_solvent.size() << endl;
		// !!! ASSUMPTION: the solvent atoms which define depth value for protein atom completely wraps the protein, i.e. they completely define the first surface layer of solvent
		// step 4.2: depth value for protein atom not yet updated
		for (unsigned int i = 0; i < atom_number; ++i){ // for protein atom
			if (depth[i] == 9999){ // not updated
				x1 = x[i]; y1 = y[i]; z1 = z[i];
				for (unsigned int j = 0; j < surface_solvent.size(); ++j){ // wrt to surface solvents
					sol_index = surface_solvent[j];
					x2 = *(xsol+sol_index); y2 = *(ysol+sol_index); z2 = *(zsol+sol_index);
					d = dist(x1,y1,z1,x2,y2,z2); // distance
					depth[i] = min(depth[i],d); // that is smallest define its depth
				} // end for
			} // end if
		} // end for

		// --- OUTPUT ---
		// step 1: write coordinates of minimally solvated structures, if requested
		if (min_sol_name != ""){
			ofstream fout;
			fout.open((min_sol_name+"-"+num2string(index)+".coordinates").c_str()); // write header
			for (unsigned int i = 0; i < atom_number; ++i){
				fout << "SOLUTE\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;
			} // end for
//			for (unsigned int i = 0; i < sol_number; ++i){
//				fout << "SOLVENT\t" << *(xsol + i) << "\t" << *(ysol + i) << "\t" << *(zsol + i) << endl;
//			} // end for
			for (unsigned int i = 0; i < surface_solvent.size(); ++i){
				fout << "SOLVENT\t" << *(xsol + surface_solvent[i]) << "\t" << *(ysol + surface_solvent[i]) << "\t" << *(zsol + surface_solvent[i]) << endl;
			} // end for
			fout.close(); fout.clear();
		} // end if
		delete[] xsol; delete[] ysol; delete[] zsol; delete[] contact;

		flog << "finished iteration " << index+1 << endl;
		flog.close(); flog.clear();
	} // end get_depths

}; // end depth_run
