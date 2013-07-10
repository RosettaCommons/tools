# ifndef EASYSTRING
# define EASYSTRING

# include <string>
# include <vector>
# include <sstream>
# include <cstdlib>
using namespace std;

bool isNumeric(string s){
	float x;
	istringstream ss(s);
	if (!(ss >> x)){
		return 0;
	} else {
		return 1;
	} // end if
} // end isNumeric

string num2string(double x){
	stringstream ss;
	ss << x;	
	string output = ss.str();
	return output;
} // end double2string

string num2string(int x){
	stringstream ss;
	ss << x;
	string output = ss.str();
	return output;
} // end double2string

string num2string(unsigned int x){
	stringstream ss;
	ss << x;
	string output = ss.str();
	return output;
} // end double2string

string strip(string s, string FS = " "){
	if (s.size() > FS.size()){
		int start = s.find_first_not_of(FS);
		int end = s.find_last_not_of(FS);
		if (start == -1){
			return "";
		} // end if
		return s.substr(start,end-start+1);
	} else if (s.size() < FS.size()){
		return s;
	} else {
		if (s == FS){
			return "";
		} else {
			return s;
		} // end if
	} // end if
} // end function strip

std::vector<string> split(string s, string fs){
	string bufferline = string(s);
	string FS = string(fs);
	string FS2 = FS + FS;
	size_t finder;
	while (1){
		finder = bufferline.find(FS2);
		if (finder == string::npos){
			break;
		} // end if
		bufferline.replace(finder, FS2.size(), FS);
	} // end while

	bufferline = bufferline.substr(0,bufferline.find_last_not_of(FS)+1);
	vector <string> output_array;
	unsigned int current_pos = 0;
	int next_pos;
	while (current_pos < bufferline.length()-FS.size()){
		next_pos = bufferline.find(FS,current_pos);
		if (next_pos == -1){
			break;
		} // end if
		output_array.push_back(bufferline.substr(current_pos,next_pos-current_pos));
		current_pos = next_pos+FS.length();
	} // end while
	output_array.push_back(bufferline.substr(current_pos));
	return output_array;
} // end function split

int word_count(string s, string w){
	int occurrence = -1;
	int pos = -1;
	do {
		pos = s.find(w,pos+1);
		occurrence = occurrence + 1;
	} while(pos != -1); // end do

	return occurrence;
} // end function word_count


// string alignmnet
string align(string s, unsigned int l, string alignment = "l"){
	string blank;
	if (s.length() >= l){
		return s.substr(0,l);
	} else {
		for (unsigned int i = 0; i < l - s.length(); i++){
			blank = blank + " ";
		} // end for
	} // end if
	if (alignment == "l"){
		return s + blank;
	} else {
		return blank + s;
	} // end if
} // end lalign

string strf(string x, int n = 2, bool sign = 0){
	if (!(isNumeric(x))){
		return x;
	} else {
		string output = x;
		float s = atof(x.c_str());
		if (s >= 0){
			output = "+"+output;
		} // end if
		if (n == 0){
			return num2string(int(s));
		} else {
			size_t dec = output.find(".");
			if (dec != string::npos){
				for (int i = 0; i < n; ++i){
					output = output + "0";
				} // end for
				string output2 = output.substr(0, dec) + ".";
				output2 = output2 + output.substr(dec+1, n);
				if (sign == 1){
					return output2;
				} else {
					if (output2.substr(0,1) == "-"){
						return output2;
					} else {
						return output2.substr(1, output2.size()-1);
					} // end if
				} // end if
			} else {
				output = output + ".";
				for (int i = 0; i < n; ++i){
					output = output + "0";
				} // end for
				if (sign == 1){
					return output;
				} else {
					if (output.substr(0,1) == "-"){
						return output;
					} else {
						return output.substr(1, output.size()-1);
					} // end if
				} // end if
			} // end if
		} // end if
	} // end if
} // end strf II

string strf(float x, int n = 2, bool sign = 0){
	return strf(num2string(x), n, sign);
} // end strf I

# endif
