#include<fstream>
#include<iostream>
#include<cstring>
#include<string>
#include<vector>
#include<cstdlib>
#include<stdlib.h>

using namespace std;

void openFile(ifstream& infile, char* argv[], int argc);
float tmsPeak(float[][]  spectrum, float baseline);
vector< vector<float> > boxcar(vector< vector<float> > spectrum, int filterSize);
vector< vector<float> > savitzkyGolay(vector< vector<float> > spectrum, int filterSize);

int main(int argc, char* argv[]){

	// Rank 2 Vector
	vector< vector<float> > spectrum;
		
	// NMR.in Line Variables
	string dataFile, outFile;
	float baseline, tolerance;
	int filterType, filterSize, filterPass, intgrtn;

	// Open and Read in NMR.in file
	ifstream infile;
	ofstream ofile;
	openFile(infile, argv, argc);

	// File is known to have 8 defined lines
	infile >> dataFile;
	infile >> baseline;	
	infile >> tolerance;
	infile >> filterType;
	infile >> filterSize;
	infile >> filterPass;
	infile >> intgrtn;
	infile >> outFile;
	infile.close();

	// data.dat Variables
	string line;
	int currRow = 0;
	string delim =', ';

	// Open and Read in data file
	infile.open(dataFile.c_str());
	while(infile.getline(line)){
		int delimIndex = line.find(delim);
		float x = line.substr(0,delimIndex);
		float y = line.subsr(delimIndex+1);

		spectrum.at(currRow).at(0) = x;
		spectrum.at(currRow).at(1) = y;
		
		currRow++;
	}
	infile.close();
	
	// TMS Variables
	float tms;

	// Find TMS Peak in the 2D Vector
	tms = tmsPeak(spectrum, baseline);

	// Subtract TMS value from all other values
	for(int i = 0; i < spectrum.size(); i++){
		for(int j = 0; j < spectrum.size(); j++){
			spectrum.at(i).at(j) = (spectrum.at(i).at(j) - tms);
		}
	}

	// Filter Switch
	switch(filterType){
		case 0:{
				// No Filter
				if( filterSize != 0 )
					cout << "Filter: None Specified. No Filter Size Needed" << endl;
				break;
			}
		case 1:{
				// Overlay a Boxcar Filter on the 2D Vector
				if( (filterSize > 0 ) && (filterSize % 2 != 0) )
					spectrum = boxcar(spectrum, filterSize);
				else
					cerr << "Invalid Boxcar Filter Size: Must be Positive Odd" << endl;
				break;
			}
		case 2:{
				// Overlay a Savitzky-Golay Filter on the 2D Vector
				if( (filterSize == 5) || (filterSize == 11) || (filterSize == 17))
					spectrum = savitzkyGolay(spectrum, filterSize);
				else
					cerr << "Invalid Savitzky Filter Size: 5, 11, OR 17" << endl;
				break;
			}
		default:
			break;
	}		
	

	// TODO:Find All Roots
	


	// TODO:Find the Peak Manifold and Hydrogen Ratio
	


	// Integration Switch -- 4 Different Integration Methods
	
		// 1.) TODO:Composite Simpson Pg. 205

		
		// 2.) TODO:Romberg Pg. 216
		
		
		// 3.) TODO:Adaptive Quadrature Pg. 224
		

		// 4.) TODO:Guassian Quadrature
		
	return 0;
}




void openFile(ifstream& infile, char* argv[], int argc){
	if(argc == 1){
		cerr << "Error: No File Given" << endl;
		exit(1);
	}
	string filename = argv[1];
	infile.open(filename.c_str());	
}

float tmsPeak(vector< vector<float> > spectrum, float baseline){
	int i = spectrum.size();
	int j = spectrum.size();

	while(spectrum.at(i).at(j) <= baseline){
		i--;
		j--;
	}
	return spectrum.at(i).at(j);
}

// Boxcar Filter changes average of current value extending one below and 
// filterSize-2 above
vector< vector<float> > boxcar(vector< vector<float> > spectrum, int filterSize){
	float sum = 0.0;
	float avg;
	bool done = false;

	int endIndex = filterSize;
	int startIndex = 0;
	for(int j = 0; j < 2; j++){
		while(!done){
			for(int i = startIndex; i < endIndex; i++){

				if(i > spectrum.size()){
					sum += spectrum.at(abs(spectrum.size() - i).at(j);
				}
				sum += spectrum.at(i).at(j);
			}
			avg = sum / float(filterSize);
			spectrum.at(startIndex+1).at(j) = avg;
			startIndex++;
			endIndex++;
		
			if(endIndex > spectrum.size()){
				done = true;
			}
		}
	}
	return spectrum;
}

// SavitzkyGolay Filter changes average of exact middle of odd numbered collection
vector< vector<float> > savitzkyGolay(vector< vector<float> > spectrum, int filterSize){
	int mid = (filterSize/2)+1;
	int startIndex = mid;
	int endIndex = mid + (filterSize/2);

	float sum = 0.0;
	float avg;
	bool done = false;

	for(int j = 0; j<2; j++){
		while(!done){
			for(int i = startIndex; i < endIndex; i++){
				if(i > spectrum.size()){
					sum += spectrum.at(abs(spectrum.size() - i)).at(j);
				}
				sum += spectrum.at(i).at(j);
			}
			avg = sum / float(filterSize);
			spectrum.at(startIndex+1).at(j) = avg;
			startIndex++;
			endIndex++;

			if(endIndex > spectrum.size()){
				done = true;
			}
		}
	}
	return spectrum;	
}
























