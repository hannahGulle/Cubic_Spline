#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<vector>
#include<array>


struct root{
	double x;
	double y;
	int index;
};

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

	vector<root> rootEndpts;

	ofile.open("testSpline");

	// Interpolate on the Natural Cubic Spline
	double k;
	for(i = 0; i < numpoints; i++){
		// Interpolate given value increase OR Decrease
		if( X[i+1] < X[i]){
			for(k = X[i]; k > X[i+1]; k = k - 0.001){
				interpX.push_back(k);
				interpY.push_back( eqn(k, X[i], a[i], b[i], c[i], d[i]) );
			}
		}
		else{
			for(k = X[i]; k < X[i+1]; k = k + 0.001){
				interpX.push_back(k);
				interpY.push_back( eqn(k, X[i], a[i], b[i], c[i], d[i]) );
			}
		}
	}

	// Find the Root Endpoints
	for(int i = 0; i < interpY.size(); i++){

		if( (interpY[i] > baseline && interpY[i+1] < baseline) ||
			(interpY[i] < baseline && interpY[i+1] > baseline)){
			root prev = {interpX[i], interpY[i], i};
			root curr = {interpX[i+1], interpY[i+1], i+1};
			rootEndpts.push_back(prev);
			rootEndpts.push_back(curr);
		}
	}

	// Print root endpoints
	for(int i = 0; i < rootEndpts.size(); i+=2){
		cout << "prev= " << rootEndpts[i].x << " curr= " << rootEndpts[i+1].x << endl;
	}

	i = 0;
	int numroots = 0;
	vector<root> roots;

	while ( i < rootEndpts.size() ){

		int INDEXA = rootEndpts[i].index;		
		double endpta = rootEndpts[i].x;
		double endptb = rootEndpts[i+1].x;

		int biIndex = 1;
		double fa = rootEndpts[i].y;
		double p;
		while( biIndex <= 1000){

			if( endpta > endptb ){
				p = endpta - fabs(endptb - endpta)/2.0;
			}
			else{
				p = endpta + fabs(endptb - endpta)/2.0;
			}
			double fp = eqn(p, interpX[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);

			if( (fp == baseline) || ( fabs(endptb - endpta)/2 < tolerance) ){
				root r;
				r.x = p;
				cout << "root= " << p << endl;
				roots.push_back(r);
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

		i+=2;
		
	}
	cout << endl;

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

	// Array X -- Initial X values with TMS Subtracted
	// Array rootVals --
	// Array rootIndex --
	// Array interpX --


	// XXX:Composite Simpson OR Trapezoidal
	double x, xi;
	i = 0;
	while( i+1 < numroots){

		double endpointA = roots[i].x;
		double endpointB = roots[i+1].x;
		int INDEXA = roots[i].index;
		int INDEXB = roots[i+1].index;

		double h = fabs(endpointB - endpointA)/50;
		double xi0 = eqn(endpointA, X[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]) +
			eqn( endpointB, X[INDEXB], a[INDEXB], b[INDEXB], c[INDEXB], d[INDEXB]);
		double xi1 = 0.0;
		double xi2 = 0.0;

		for(int j = 1; j < 50; j++){

			x = endpointA + j*h;
			if( j%2 == 0){
				xi2 = xi2 + eqn( x, X[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
			}
			else{
				xi1 = xi1 + eqn( x, X[INDEXA],  a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
			}
		}
		i = i + 2;
		xi = h * (xi0 + (2*xi2) + (4*xi1))/3.0;
		cout << "area= " << xi << endl;
	}
	cout << "end compsimp" << endl; 
	cout << endl;

	//XXX: Romberg

	i = 0;
	while( i+1 < numroots){	

		double A = roots[i].x;
		double B = roots[i+1].x;
		int INDEXA = roots[i].index;
		int INDEXB = roots[i+1].index;

		int n = 50;

		double sum, part;
		int m, tmp;

		double R[2][n];
		double h = fabs(B - A);

		R[0][0] = (h / 2.0) * (eqn(A, X[INDEXA],  a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]) +
				eqn( B, X[INDEXB], a[INDEXB], b[INDEXB], c[INDEXB], d[INDEXB]));

		for(int j = 2; j <= n; j++){
			sum = 0.0;
			m = exp((i-2) * log(2.0)) + 0.5;
			for(double k = 1; k <= m; k++){
				part = A + (h * (k - 0.5));
				sum += eqn( part, X[INDEXA],  a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
			}
			R[1][0] = 0.5 * ( R[0][0] + (h * sum) );

			for( int l = 2; l <= i; l++){
				tmp = exp(2 * (j-1) * log(2.0)) + 0.5;
				R[1][j-1] = R[1][j-2] + (R[1][j-2] - R[0][j-2]) / ( tmp - 1.0);
			}
			h = h/2.0;
			for(int m = 1; m <= i; m++){

				R[0][j-1] = R[1][j-1];
			}
		}
		cout << "area= " << R[1][n-1] << endl;
		i+=2;
	}
	cout << "end romberg" << endl;
	cout << endl;

	//TODO: Gaussian Quadrature

	// Step 1: Weights and Root
	
	// Step 2: Integrate


	// XXX: Adaptive Quadrature

	int n = 50;
	double TOL[n], A[n], H[n], FA[n], FC[n], FB[n], S[n], V[7];
	int L[n];
	double AA, BB, EPS, APP, FD, FE, S1, S2;
	int I, LEV, INDEXA, INDEXB;
	int count = 0;

	EPS = tolerance;

	while( count+1 < numroots){
		AA = roots[count].x;
		BB = roots[count+1].x;
		INDEXA = roots[count].index;
		INDEXB = roots[count+1].index;	

		APP = 0.0;
		I = 1;
		TOL[I] = 10.0 * EPS;
		A[I] = AA;
		H[I] = 0.5 * (BB - AA);
		FA[I] = eqn(AA, X[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
		FA[I] = eqn((AA + H[I]), X[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
		FB[I] = eqn(BB, X[INDEXB], a[INDEXB], b[INDEXB], c[INDEXB], d[INDEXB]);
		S[I] = H[I] * (FA[I] + 4.0 * FC[I] + FB[I]) / 3.0;
		L[I] = 1;

		while ( I > 0){

			FD = eqn( (A[I] + 0.5 + H[I]), X[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
			FE = eqn( (A[I] + 1.5 + H[I]), X[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
			S1 = H[I] * (FA[I] + 4.0 * FD + FC[I]) / 6.0;
			S2 = H[I] * (FC[I] + 4.0 * FE + FB[I]) / 6.0;

			V[0] = A[I];
			V[1] = FA[I];
			V[2] = FC[I];
			V[3] = FB[I];
			V[4] = H[I];
			V[5] = TOL[I];
			V[6] = S[I];
			LEV = L[I];
			I--;

			if( S1 + S2 - V[6] < V[5]){
				APP += S1 + S2;
			}
			else{
				I++;
				A[I] = V[0] + V[4];
				FA[I] = V[2];
				FC[I] = FE;
				FB[I] = V[3];
				H[I] = 0.5 * V[4];
				TOL[I] = 0.5 * V[5];
				S[I] = S2;
				L[I] = LEV + 1;

				I++;
				A[I] = V[0];
				FA[I] = V[1];
				FC[I] = FD;
				FB[I] = V[2];
				H[I] = H[I-1];

				S[I] = S1;
				L[I] = L[I-1];
			}
		}
		cout << "area= " << APP << endl;
		count+= 2;
	}
	cout << "end adaptive" << endl;
	cout << endl;
	

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
