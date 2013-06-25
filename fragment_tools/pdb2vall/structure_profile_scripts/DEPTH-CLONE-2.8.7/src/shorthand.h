# ifndef SHORTHAND_H
# define SHORTHAND_H

# include <cstdlib>
# include <iostream>
# include <string>
# include <vector>
# include <map>
# include <set>
# include <algorithm>
# include <fstream>
# include <sstream>
# include <stdio.h>
# include <math.h>
# include "easystring.h"

using namespace std;


float line(float x1, float y1, float x2, float y2, float x){
	return y2 - ((y2-y1) / (x2-x1)) * (x2-x);
} // end line


double round_double(double x, int dec){
	if (dec <= 0){
		return double(int(x));
	} else {
		int w = int(x);
		int p = 10;
		for (int i = 1; i < dec; ++i){
			p = p*10;
		} // end for
		double decimal = p*(x - w);
		double decimal_new;
		if ( decimal - int(decimal) > 0.5){
			decimal_new = int(decimal + 1);
		} else {
			decimal_new = int(decimal);
		} // end if
		decimal_new = decimal_new/p;
		return w + decimal_new;
	} // end if
} // end round_double

vector <int> set2vec(set <int> S){
	vector <int> output;
	for (set <int>::iterator iter = S.begin(); iter != S.end(); ++iter){
		output.push_back(*iter);
	} // end for
	return output;
} // end set2vec

vector <string> set2vec (set <string> s){
	vector <string> output;
	for (set <string>::iterator iter = s.begin(); iter != s.end(); ++iter){
		output.push_back(*iter);
	} // end for
	return output;
} // end set2vec

vector <pair <string, string> > set2vec(set <pair <string, string> > s){
	vector <pair <string, string> > output;
	for (set <pair <string, string> >::iterator iter = s.begin(); iter != s.end(); ++iter){
		output.push_back(*iter);
	} // end for
	return output;
} // end set2vec

set <int> vec2set(vector <int> input){
	set <int> non_redun_list;
	for (unsigned int i = 0; i < input.size(); ++i){
		non_redun_list.insert(input[i]);
	} // end for
	return non_redun_list;
} // end set2vec



vector <int> non_redundant(vector <int> input){
	sort(input.begin(), input.end());
	vector<int>::iterator it;

	// using default comparison:
	it = unique(input.begin(), input.end());
	input.resize(it - input.begin() );

	return input;
} // end for

vector <string> non_redundant(vector <string> input){
	sort(input.begin(), input.end());
	vector<string>::iterator it;

	// using default comparison:
	it = unique(input.begin(), input.end());
	input.resize(it - input.begin() );

	return input;
} // end for


bool is_element(int i, set <int> box){
	if (box.find(i) != box.end()){
		return 1;
	} else {
		return 0;
	} // end if
} // end is_element


bool is_element(int p, vector <int> box){
	for (unsigned int i = 0; i < box.size(); ++i){
		if (box[i] == p){
			return 1;
		} // end if
	} // end for
	return 0;
} // end is_element

bool is_element(string p, vector <string> box){
	for (unsigned int i = 0; i < box.size(); ++i){
		if (box[i] == p){
			return 1;
		} // end if
	} // end for
	return 0;
} // end is_element

bool is_element(float i, set <float> box){
	if (box.find(i) != box.end()){
		return 1;
	} else {
		return 0;
	} // end if
} // end is_element

// round to nearest integer
int near_round(double s){
	double dec = s - int(s);
	if (dec >= 0.5){
		return int(s) + 1;
	} else {
		return int(s);
	} // end if
} // end near_round

// decimal round
double dec_round(double s, int dec){
	double up = 10;
	for (int i = 0; i < dec - 1; ++i){
		up = up*10;
	} // end for
	int m = near_round(s*up);
	return m / up;
} // end dec_round

// pseudo-pdb writer
string pdb_line(string atom, int serial, string name, string altLoc, string resName, string chainID, string resSeq, string iCode, double x, double y, double z, double occupancy, double T, string element = "", string charge = ""){
	string line = "";
	line = align(atom,6,"l") + align(num2string(serial),5,"r") + "  " + align(name,3,"l") + align(altLoc,1,"l") + align(resName,3,"l") + " " + align(chainID,1,"l") + align(resSeq,4,"r") + align(iCode,1,"l") + "   " + align(strf(dec_round(x,3),3,0),8,"r") + align(strf(dec_round(y,3),3,0),8,"r") + align(strf(dec_round(z,3),3,0),8,"r") + align(strf(dec_round(occupancy,2),2,0),6,"r") + align(strf(dec_round(T,2),2,0),6,"r") + "          " + align(element,2,"l") + align(charge,2,"r") + "\n";
	return line;
} // end pdb_line

double dist(double x1, double y1, double z1, double x2, double y2, double z2){
	double dx, dy, dz;
	dx = (x1-x2); dy = (y1-y2); dz = (z1-z2);
	return sqrt(dx*dx + dy*dy + dz*dz);
} // end dist

set <int> set_substract(set <int> A, set <int> B){
	vector <int> pre_output;
	vector <int> adapt_B = set2vec(B);
	for ( set <int> ::iterator iter = A.begin(); iter != A.end(); ++iter){
		if (!(is_element(*iter, adapt_B))){
			pre_output.push_back(*iter);
		} // end if
	} // end for
	return vec2set(pre_output);
} // end set_substract

# endif
