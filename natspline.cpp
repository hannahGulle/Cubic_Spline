#include<iostream>
#include<fstream>
#include<string>
#include<math.h>

using namespace std;

int main(int argc, char* argv[]){

	// Total Number of Points in NMR File
	int numpoints = 0;

	// NMR Line Variables
	string dataFile, outFile;
	double baseline, tolerance;
	int filterType, filterSize, filterPass, integration;

	// Open and Read NMR.in file
	ifstream infile;
	ofstream ofile;
	
	if(argc == 1){
		cerr << "Error: No File Given" << endl;
		exit(1);
	}
	
	string filename = argv[1];
	infile.open(filename.c_str());

	// File is Known to have 8 Defined Lines
	infile >> dataFile;
	infile >> baseline;
	infile >> tolerance;
	infile >> filterType;	
	infile >> filterSize;
	infile >> filterPass;
	infile >> integration;
	infile >> outFile;
	infile.close();


	// data.dat Variables
	string line;

	infile.open(dataFile.c_str());
	while(getline(infile,line)){
		numpoints++;
	}


	// 1D Array Containing X Values
	double X[numpoints];
	// 1D Array Containing Y Values
	double Y[numpoints];

	infile.close();
	infile.open(dataFile.c_str());
	
	for(int i = 0; i < numpoints; i++){
		double x, y;
		infile >> x >> y;
		X[i] = x;
		Y[i] = y;
	}
	infile.close();

	// TMS Variables
	int tmsIndex = -1;

	// Finds TMS Peak
	for(int i = 0; i < numpoints; i++){
		if( (Y[i] > baseline) ){
			if ( ((tmsIndex != -1) && (X[i] > X[tmsIndex])) || (tmsIndex == -1) ){
				tmsIndex = i;
			}
		}
	}	
	double tms = X[tmsIndex];
	cout << "tms= " << tms << endl;

	ofile.open("nattmsdat");
	for(int i = 0; i < numpoints; i++){
		X[i] = X[i] - tms;
		ofile << X[i] << " " << Y[i] << endl;
	}
	ofile.close();

	//cout << "filter type: " << filterType << endl;
	//TODO: Input Filtering Code from spline.cpp
	
	int m = numpoints - 1;

	double a[numpoints], b[numpoints], c[numpoints], d[numpoints];
	double h[numpoints], alpha[m], l[numpoints], u[numpoints], z[numpoints];

  	int i, J;
  	for(i = 0; i < numpoints; i++){
    		a[i] = Y[i];
		cout << "a= " << a[i] << endl;
	}
  
    	for (i = 0; i <= m; i++){
		h[i] = X[i+1] - X[i];
	}

	alpha[0] = 0.0;
	for (i = 1; i < m-1; i++) {
        	alpha[i] = (a[i+1] - a[i]) * 3.0 / h[i] - (a[i] - a[i-1]) * 3.0 / h[i-1];
	}
      	
	l[0] = 1.0;
	u[0] = 0.0;
	z[0] = 0.0;
      	
	for (i = 1; i <= m; i++) {
        	l[i] = 2.0 * (X[i+1] - X[i-1]) - h[i-1] * u[i-1];
       		u[i] = h[i] / l[i];
         	z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
	}
      	
	l[numpoints-1] = 1.0;
	z[numpoints-1] = 0.0;
	c[numpoints-1] = 0.0;

      	for (i = 0; i <= m; i++) {
        	J = m - i;
         	c[J] = z[J] - u[J] * c[J+1];
         	b[J] = (a[J+1] - a[J]) / h[J] - h[J] * (c[J+1] + 2.0 * c[J]) / 3.0;
         	d[J] = (c[J+1] - c[J]) / (3.0 * h[J]);
      	} 


	return 0;
}
