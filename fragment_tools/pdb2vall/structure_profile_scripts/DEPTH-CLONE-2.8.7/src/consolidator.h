// C++ code to combine and compare depth result from complex and individual chains
# include <cstdlib>
# include <iostream>
# include <fstream>
# include <string>
# include <map>
# include "easystring.h"
using namespace std;

class recorder{
	private:
	string name, restype;
	float S_ALL, S_MC, S_SC, S_SCP, S_SCNP;
	float X_ALL, X_MC, X_SC, X_SCP, X_SCNP, D_ALL, D_MC, D_SC, D_SCP, D_SCNP;

	public:
	recorder(string name_in){
		name = name_in;
		S_ALL = -99; S_MC = -99; S_SC = -99; S_SCP = -99; S_SCNP = -99;		
	} // end recorder
	~recorder(){}

	void set_restype(string restype_in){
		restype = restype_in;
	} // end 

	void S_ALL_update(float s){
		if (S_ALL == -99){
			S_ALL = s;
		} else if (s > S_ALL) {
			S_ALL = s;
		} // end if
	} // end 
	void S_MC_update(float s){
		if (S_MC == -99){
			S_MC = s;
		} else if (s > S_MC) {
			S_MC = s;
		} // end if
	} // end 
	void S_SC_update(float s){
		if (S_SC == -99){
			S_SC = s;
		} else if (s > S_SC) {
			S_SC = s;
		} // end if
	} // end 
	void S_SCP_update(float s){
		if (S_SCP == -99){
			S_SCP = s;
		} else if (s > S_SCP) {
			S_SCP = s;
		} // end if
	} // end 
	void S_SCNP_update(float s){
		if (S_SCNP == -99){
			S_SCNP = s;
		} else if (s > S_SCNP) {
			S_SCNP = s;
		} // end if
	} // end 

	void set_X_ALL(float s){
		X_ALL = s;
	} // end set_X_ALL

	void set_X_MC(float s){
		X_MC = s;
	} // end set_X_ALL

	void set_X_SC(float s){
		X_SC = s;
	} // end set_X_ALL

	void set_X_SCP(float s){
		X_SCP = s;
	} // end set_X_ALL

	void set_X_SCNP(float s){
		X_SCNP = s;
	} // end set_X_ALL

	void set_D_ALL(float s){
		D_ALL = s;
	} // end set_X_ALL

	void set_D_MC(float s){
		D_MC = s;
	} // end set_X_ALL

	void set_D_SC(float s){
		D_SC = s;
	} // end set_X_ALL

	void set_D_SCP(float s){
		D_SCP = s;
	} // end set_X_ALL

	void set_D_SCNP(float s){
		D_SCNP = s;
	} // end set_X_ALL

	string retrieve(){
		string out_string = "";
		out_string = out_string + name + "\t" + restype + "\t";
		out_string = out_string + strf(X_ALL - D_ALL) + "\t";
		out_string = out_string + "\t"+strf(S_ALL) + "\t";
		out_string = out_string + "\t"+strf(X_ALL)+ "\t"+strf(D_ALL)+ "\t";

		out_string = out_string + strf(X_MC - D_MC) + "\t";
		out_string = out_string + "\t"+strf(S_MC) + "\t";
		out_string = out_string + "\t"+strf(X_MC)+ "\t"+strf(D_MC)+ "\t";

		out_string = out_string + strf(X_SC - D_SC) + "\t";
		out_string = out_string + "\t"+strf(S_SC) + "\t";
		out_string = out_string + "\t"+strf(X_SC)+ "\t"+strf(D_SC)+ "\t";

		out_string = out_string + strf(X_SCP - D_SCP) + "\t";
		out_string = out_string + "\t"+strf(S_SCP) + "\t";
		out_string = out_string + "\t"+strf(X_SCP)+ "\t"+strf(D_SCP)+ "\t";

		out_string = out_string + strf(X_SCNP - D_SCNP) + "\t";
		out_string = out_string + "\t"+strf(S_SCNP) + "\t";
		out_string = out_string + "\t"+strf(X_SCNP)+ "\t"+strf(D_SCNP)+ "\t";

		out_string = out_string + "\n";
		
		return out_string;
	} // end 
}; // end class

class consolidator{
	private:
	int last;
	map <string, int> finder;
	vector <recorder> data;
	vector <string> bufferline;
	string line, restype, name;


	public:
	consolidator(){}
	~consolidator(){}

	void consolidate(string complex_depth, string consolidate_depth, string out_fname){
		ifstream fcomplex; fcomplex.open(complex_depth.c_str());

		while (!fcomplex.eof()){
			getline(fcomplex, line);
			if (line.size() == 0){
				continue;
			} else if (line.substr(0,1) == "#"){
				continue;
			} // end if

			bufferline = split(line, "\t");
			name = bufferline[0];
			restype = bufferline[1];

			data.push_back(recorder(name));
			last = data.size()-1;

			data[last].set_restype(restype);
			finder[name] = last;
			data[last].set_X_ALL(atof(bufferline[2].c_str()));
			data[last].S_ALL_update(atof(bufferline[3].c_str()));

			data[last].set_X_MC(atof(bufferline[4].c_str()));
			data[last].S_MC_update(atof(bufferline[5].c_str()));

			data[last].set_X_SC(atof(bufferline[6].c_str()));
			data[last].S_SC_update(atof(bufferline[7].c_str()));

			data[last].set_X_SCP(atof(bufferline[8].c_str()));
			data[last].S_SCP_update(atof(bufferline[9].c_str()));

			data[last].set_X_SCNP(atof(bufferline[10].c_str()));
			data[last].S_SCNP_update(atof(bufferline[11].c_str()));
		} // end while
		fcomplex.close(); fcomplex.clear();

		ifstream fchains; fchains.open(consolidate_depth.c_str());
		while (!fchains.eof()){
			getline(fchains, line);
			if (line.size() == 0){
				continue;
			} else if (line.substr(0,1) == "#"){
				continue;
			} // end if

			bufferline = split(line, "\t");
			name = bufferline[0];

			data[finder[name]].set_D_ALL(atof(bufferline[2].c_str()));
			data[finder[name]].S_ALL_update(atof(bufferline[3].c_str()));

			data[finder[name]].set_D_MC(atof(bufferline[4].c_str()));
			data[finder[name]].S_MC_update(atof(bufferline[5].c_str()));

			data[finder[name]].set_D_SC(atof(bufferline[6].c_str()));
			data[finder[name]].S_SC_update(atof(bufferline[7].c_str()));

			data[finder[name]].set_D_SCP(atof(bufferline[8].c_str()));
			data[finder[name]].S_SCP_update(atof(bufferline[9].c_str()));

			data[finder[name]].set_D_SCNP(atof(bufferline[10].c_str()));
			data[finder[name]].S_SCNP_update(atof(bufferline[11].c_str()));
		} // end while
		fchains.close(); fchains.clear();

		ofstream fout; fout.open(out_fname.c_str());
		fout << "# Residue\ttype\t";
		fout << "diff(depth)\terror\tcomplex\tchain\t";
		fout << "diff(depth(main-chain))\terror\tcomplex\tchain\t";
		fout << "diff(depth(side-chain))\terror\tcomplex\tchain\t";
		fout << "diff(depth(side-chain,polar))\terror\tcomplex\tchain\t";
		fout << "diff(depth(side-chain,nonpolar))\terror\tcomplex\tchain\t";
		fout << endl;
		for (unsigned int i = 0; i < data.size(); ++i){
			fout << data[i].retrieve();
		} // end for
		fout.close(); fout.clear();
	} // end consolidate

};
