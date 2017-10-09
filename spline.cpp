#include<iostream>
#include<fstream>
#include<iostream>
#include<cstring>
#include<string>

using namespace std;

void openFile(ifstream& infile, char* argv[], int argc);
int main(int argc, char* argv[]){
	
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

	// Open and Read in data file
	infile.open(dataFile);
	string line;
	while(infile.getline(line)){
		
	}

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
