#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<vector>
#include<array>

using namespace std;

double findtmsPeak(double baseline, double X[], double Y[], int numpoints);
double eqn(double x, double xi, double a, double b, double c, double d);
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

	// ********************** TMS STEP ************************************************

	double tms = findtmsPeak(baseline, X, Y, numpoints);
	ofile.open("nattmsdat");
	for(int i = 0; i < numpoints; i++){
		X[i] = X[i] - tms;
		ofile << X[i] << " " << Y[i] << endl;
	}
	ofile.close();

	// *********************** APPLY FILTERS STEP **************************************

	cout << "filter type: " << filterType << endl;

	switch(filterType){

		// NO FILTER
		case 0:{
			       cout << "Filter: None Specified. No Filter Size Needed" << endl;
			       break;
		       }

		       // BOXCAR FILTER
		case 1:{
			       // Overlay a Boxcar Filter on the 2D Vector
			       if( (filterSize > 0 ) && (filterSize % 2 != 0) ){
				       for(int z = 0; z < filterPass; z++){
					       float sum = 0.0;	// Sum of the current points in the filter
					       float avg;		// Average of the current points in the filter
					       bool done = false;	// Process End boolean

					       int endIndex = filterSize;
					       int startIndex = numpoints-1;

					       while(startIndex < numpoints){
						       for(int i = startIndex; i < endIndex; i++){

							       // Treats the vector as a circular structure
							       // Overflow returns to the beginning of the vector
							       if(i > numpoints){
								       sum += Y[i - numpoints];
							       }
							       else{
								       sum += Y[i];
							       }
							       // No underflow occurs because the starting value is 1	
						       }
						       avg = sum / float(filterSize);

						       // Overflow applies to the top limit of the filter range as well
						       if(startIndex+1 > numpoints){
							       Y[startIndex - numpoints] = avg;
						       }
						       else{
							       Y[startIndex+1] = avg;
						       }
						       if(startIndex+1 > numpoints){
							       startIndex = 0;
						       }
						       else{
							       startIndex++;
						       }
						       endIndex++;
						       sum = 0.0;	// Sum is restart for each filter pass
					       }
				       }
			       }
			       else{
				       cerr << "Invalid Boxcar Filter Size: " << filterSize << ". Must be Positive Odd" << endl;
				       exit(1);
			       }
			       break;

		       }

		       // SAVITSKY GOLAY FILTER
		case 2:{
			       // Overlay a Savitzky-Golay Filter on the 2D Vector
			       if( (filterSize == 5) || (filterSize == 11) || (filterSize == 17)){

				       int m = filterSize/2;
				       int j;
				       double norm;
				       double sg[filterSize] = {};

				       if( filterSize == 5){  
					       double tmp[filterSize] = {-3.0,12.0,17.0,12.0,-3.0};
					       for(int i = 0; i < filterSize; i++){
						       sg[i] = tmp[i];
					       }
					       norm = 35;
				       }
				       if( filterSize == 7){  
					       double tmp[filterSize] = {-2.0,3.0,6.0,7.0,6.0,3.0,-2.0};

					       for(int i = 0; i < filterSize; i++){
						       sg[i] = tmp[i];
					       }
					       norm = 21;
				       }
				       if( filterSize == 11){  
					       double tmp[filterSize] = {-36.0,9.0,44.0,69.0,84.0,89.0,84.0,69.0,44.0,9.0,-36.0};

					       for(int i = 0; i < filterSize; i++){
						       sg[i] = tmp[i];
					       }
					       norm = 429;
				       }

				       for(int z = 0; z < filterPass; z++){
					       for(int j = 0; j < numpoints; j++){
						       double ci, next, y;
						       for(int i = 0; i < filterSize; i++){
							       if( j+i > numpoints){
								       next = Y[j + i - numpoints];
							       }
							       else{
								       next = Y[j+i];
							       }

							       ci = sg[i];
							       y += (ci * next);
						       }
						       if(j+m > numpoints){
							       Y[j+m - numpoints] = y / norm;
						       }
						       else{
							       Y[j+m] = y / norm;
						       }
						       y = 0;
					       }
						cout << "filterpass= " << z << endl;
				       }
			       }
			       else{
				       cerr << "Invalid Savitzky Filter Size: " << filterSize << ". 5, 11, OR 17" << endl;
				       exit(1);
			       }
			       break;
		       }

		default:{
				break;

			}

	}

	ofile.open("sgdat");
	for(int i = 0; i < numpoints; i++){
		ofile << X[i] << " " << Y[i] << endl;
	}
	ofile.close();


	// ************************ Fit the Natural Cubic Spline Step **************************
	int m = numpoints - 1;

	double a[numpoints], b[numpoints], c[numpoints], d[numpoints];
	double h[numpoints], alpha[m], l[numpoints], u[numpoints], z[numpoints];

	int i, J;
	for(i = 0; i < numpoints; i++){
		a[i] = Y[i];
	}

	for (i = 0; i <= m; i++){
		h[i] = X[i+1] - X[i];
	}

	alpha[0] = 0.0;
	for (i = 1; i < m; i++) {
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
	// Cubic Spline Made of Coefficients in Arrays a, b, c, d


	// ***************************** INTERPOLATION STEP ************************************

	vector<double> interpX;
	vector<double> interpY;

	vector<double> roots;
	vector<int> rootIndex;

	ofile.open("testSpline");

	// Interpolate on the Natural Cubic Spline
	double k;
	double prev = Y[0];
	double curr;
	double val;
	for(i = 0; i < numpoints; i++){

		// Interpolate given value increase OR Decrease
		if( X[i+1] < X[i]){
			for(k = X[i]; k > X[i+1]; k = k - 0.001){
				interpX.push_back(k);
				val = eqn(k, X[i], a[i], b[i], c[i], d[i]);
				interpY.push_back(val);
			}
		}
		else{
			for(k = X[i]; k < X[i+1]; k = k + 0.001){
				interpX.push_back(k);
				val = eqn(k, X[i], a[i], b[i], c[i], d[i]);
				interpY.push_back(val);
			}
		}

		curr = val;

		// Root Finding
		if ( ((prev > baseline) && (curr < baseline) ) ||
				((prev < baseline) && (curr > baseline) )) {

			cout << "prev= " << prev << " curr= " << curr << endl;
			roots.push_back(prev);
			roots.push_back(curr);
			rootIndex.push_back(i);
		}
		prev = curr;
	}

	// ***************************** ROOT FINDING STEP ******************************************

	//TODO:Root Finding at Baseline using Steffenson
	// Need to hold the indexes of the roots after finding the roots based on the
	// endpoints (one above baseline and one below baseline)	 
	// INPUT: initial approximation (p0), tolerance, maximum iterations of N0

	// For each equation given each coefficient of a, b, c, d, 
	// calculate the root(s) of the equation


	// GIVEN: roots vector (double) with prev, curr, prev, curr, ...
	// and rootIndex (int)
	double rootVals[rootIndex.size()] = {};

	i = 0;
	int j = 0;
	int numroots = 0;
	cout << "rootIndex.size()= " << rootIndex.size() << endl;
	while ( (i < roots.size()) && (j < rootIndex.size()) ){

		int rootI = rootIndex[j];		

		double endpta = roots[i];
		double endptb = roots[i+1];

		int biIndex = 1;
		double fa = eqn( endpta, X[rootI], a[rootI], b[rootI], c[rootI], d[rootI]);

		while( biIndex <= 1000){


			double p = endpta + (endptb - endpta)/2;
			double fp = eqn(p, X[rootI], a[rootI], b[rootI], c[rootI], d[rootI]);

			if( (fp == baseline) || ( (endptb - endpta)/2 < tolerance) ){
				cout << "root= " << X[rootI] << endl;
				rootVals[i] = X[rootI];
				numroots++;
				break;
			}
			biIndex++;

			if( fa * fp > 0){
				endpta = p;
			}
			else{
				endptb = p;
			}
		}

		i++;
		j++;
	}

	//TODO: 4 Integration Techniques

	// Find Root Peaks

/*
	double rootpeaks[numroots/2] = {};
	double peak;

	for(int i = 0; i < numroots; i++){
		int prevRoot = rootIndex[i];
		int currRoot = rootIndex[i+1];
		peak = interpY[prevRoot];

		if( prevRoot > currRoot ){
			int index = currRoot;
			while ( index < prevRoot ){
				if( (interpY[index] > peak) && (interpX[index] < interpX[prevRoot]) && (interpX[index] > interpX[currRoot])){
					peak = interpY[index];
				}
				index++;
			}
		}
		else{
			int index = prevRoot;
			while (index < currRoot){
				if( (interpY[index] > peak) && (interpX[index] > interpX[prevRoot]) && (interpX[index] < interpX[currRoot])){
					peak = interpY[index];
				}
				index++;
			}
		}
		rootpeaks[i] = peak;
	}
	
	for(int i = 0; i < numroots/2; i++){
		cout << "peak= " << rootpeaks[i] << endl;
	}
*/

	// Composite Simpson OR Trapezoidal
	double x, xi;
	i = 0;
	while( i+1 < numroots){

		double endpointA = rootVals[i];
		double endpointB = rootVals[i+1];

		double h = (endpointB - endpointA)/50;
		double xi0 = eqn( a[rootIndex[i]], b[rootIndex[i]], c[rootIndex[i]], d[rootIndex[i]], endpointA, interpX[rootIndex[i]]) +
				eqn( a[rootIndex[i]], b[rootIndex[i]], c[rootIndex[i]], d[rootIndex[i]], endpointB, interpX[rootIndex[i+1]]);
		double xi1 = 0.0;
		double xi2 = 0.0;
		cout << "h= " << h << " xi0= " << xi0 << endl;

		for(int j = 1; j < 50; j++){

			x = endpointA + j*h;
			if( j%2 == 0){
				xi2 = xi2 + eqn( a[rootIndex[i]], b[rootIndex[i]], c[rootIndex[i]], d[rootIndex[i]], x, interpX[rootIndex[i]]);
			}
			else{
				xi1 = xi1 + eqn( a[rootIndex[i]], b[rootIndex[i]], c[rootIndex[i]], d[rootIndex[i]], x, interpX[rootIndex[i]]);
			}
		}
		i = i + 2;
		xi = h * (xi0 + (2*xi2) + (4*xi1))/3.0;
		cout << "area= " << xi << endl;
	}

	//TODO: Romberg
	
	

	//TODO: Gaussian Quadrature


	//TODO: Adaptive Quadrature


	//TODO: Compare the Areas as a Ratio of Hydrogen Atoms


	//TODO: Output the Results
	return 0;
}

double eqn(double x, double xi, double a, double b, double c, double d){

	return( a + b * (x - xi) + c * (x - xi) * (x - xi)  +  d * (x - xi) * (x - xi) * (x - xi) );
}

double findtmsPeak(double baseline, double X[], double Y[], int numpoints){
	int tmsIndex = -1;

	for(int i = 0; i < numpoints; i++){
		if( (Y[i] > baseline) ){
			if ( ((tmsIndex != -1) && (X[i] > X[tmsIndex])) || (tmsIndex == -1) ){
				tmsIndex = i;
			}
		}
	}	
	return X[tmsIndex];
}
