// header file for buiding grid map

# include <cstdlib>
# include <vector>
# include <string>
# include <stdio.h>
# include <iostream>
# include <fstream>
using namespace std;
# include "easystring.h"

class Grid_Map{
	private:
	string line;
	float x_start, y_start, x_step, y_step, _yscale;
	unsigned int c;
	int x_index, y_index, xL, yL;
	vector <string> bufferline;
	float** M;

	public:
	Grid_Map(){} // end Grid_Map
	~Grid_Map(){}

	void read(string filename){
		ifstream fin;
		fin.open(filename.c_str());
		getline(fin, line);
		bufferline = split(line, " ");
		x_start = atof(bufferline[1].c_str());
		x_step  = atof(bufferline[2].c_str());
		y_start = atof(bufferline[3].c_str());
		y_step  = atof(bufferline[4].c_str());
		xL = atoi(bufferline[5].c_str());
		yL = atoi(bufferline[6].c_str());
		_yscale = atof(bufferline[7].c_str());
		M = new float* [yL];
		for (int i = 0; i < yL; ++i){
			M[i] = new float [xL];
		} // end for
		c = -1;
		while (!fin.eof()){
			getline(fin, line);
			if (line.size() == 0){
				continue;
			} // end if
			c = c + 1;
			bufferline = split(line, "\t");
			for (unsigned int u = 0; u < bufferline.size(); ++u){
				M[c][u] = atof(bufferline[u].c_str());
			} // end for
		} // end while
		fin.close(); fin.clear();
	} // end read

	float value(float x, float y){
		y_index = int( (y*_yscale - y_start) / y_step );
		x_index = int( (x - x_start) / x_step );
		if (y_index >= yL){
			y_index = yL - 1;
			cout << "warning: map " << y << " onto " << ((yL-1)*y_step + y_start)/_yscale << endl;
		} // end if
		if (x_index >= xL){
			x_index = xL - 1;
			cout << "warning: map " << x << " onto " << (xL-1)*x_step + x_start << endl;
		} // end if
		return M[y_index][x_index];
	} // end value

	float xstart(){
		return x_start;
	} // end fx
	float ystart(){
		return y_start;
	} // end fx
	float xstep(){
		return x_step;
	} // end fx
	float ystep(){
		return y_step;
	} // end fx
	float yscale(){
		return _yscale;
	} // end fx
	int xrange(){
		return xL;
	} // end fx
	int yrange(){
		return yL;
	} // end fx

}; // end Grid_Map
