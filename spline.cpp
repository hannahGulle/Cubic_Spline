// Hannah Gulle
// CSC 335 Project 1 -- Cubic Spline
// Due 11/15/17

#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<vector>
#include<algorithm>
#include<cmath>
#include<ctime>

// Holds the double X and Y values of a point
struct point{
	double x;
	double y;
};

// Holds the coefficients of a cubic spline 
// segment (a, b, c, d) and the corresponding
// double x as the lower range x value
struct spline{
	double a;
	double b;
	double c;
	double d;
	double x;
};

// Holds the 3 possible savitzky golay filter
// size arrays
struct sg{
	double five[5];
	double seven[7];
	double eleven[11];
};

using namespace std;

// AREA INTEGRATION METHODS
double romberg( double a, double b, spline spline[], int numpoints, double baseline, double tolerance);
double adaptive( double AA, double BB, spline spline[], int numpoints, double tolerance, double baseline);
double compositeSimpson( double A, double B, spline spline[], int numpoints, double baseline);
double guassian( double a, double b, spline spline[], int numpoints, double baseline);

// CUBIC SPLINE METHODS
vector<point> interpolateSpline(point points[], int numpoints, spline spline[], double baseline);
double eqn(double x, spline s, double baseline);

// ROOT FINDING, SORTING, AND TMS PEAK METHODS
double findRoot( point r1, point r2, spline s[], int numpoints, double baseline, double tolerance);
bool sortIncrX( point p1, point p2);
double findtmsPeak(double baseline, point points[], int numpoints);
double findTop( double x1, double x2, spline s[], int numpoints, double baseline);

int main(int argc, char* argv[]){

	cout.precision(12);

	// Possible savitzky golay filter sizes with corresponding weights
	sg sgfilters = { {-3.0,12.0,17.0,12.0,-3.0}, {-2.0,3.0,6.0,7.0,6.0,3.0,-2.0}, {-36.0,9.0,44.0,69.0,84.0,89.0,84.0,69.0,44.0,9.0,-36.0} };

	// Total Number of Points in NMR File
	int numpoints = 0;

	// NMR Line Variables
	string dataFile, outFile, filterName, integrationName;
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

	// Data File Variables
	string line;
	// Read Data File to Find the Number of Points
	infile.open(dataFile.c_str());
	while(getline(infile,line)){
		numpoints++;
	}
	infile.close();
	
	// Begin Finding the Time in Seconds for Analysis
	clock_t start = clock();


	// 1D Array Containing Point Structs Containing X and Y Values
	// From the Data File
	point points[numpoints];
	infile.open(dataFile.c_str());
	for(int i = 0; i < numpoints; i++){
		point p;
		infile >> p.x >> p.y;
		points[i] = p;
	}
	infile.close();

	// Sorts the points array by increasing X value
	sort( points, points + numpoints, sortIncrX );

	// ********************** TMS STEP ************************************************
	// Finds the highest Y value above the baseline with the highest x value
	// then subtracts the tms point x value from all x values.
	ofile.open("tms");
	double tms = findtmsPeak(baseline, points, numpoints);
	for(int i = 0; i < numpoints; i++){
		points[i].x = points[i].x - tms;
		ofile << points[i].x << " " << points[i].y << endl;
	}
	ofile.close();

	// *********************** APPLY FILTERS STEP **************************************
	// Three filters available: None, Boxcar, Savitzky Golay (5, 7, 11)
	switch(filterType){
		// NO FILTER
		case 0:{
			filterName = "No Filter";
			 cout << "Filter: None Specified. No Filter Size Needed" << endl;
			 break;
		       }

		       // BOXCAR FILTER
		case 1:{
			filterName = "Boxcar Filter";
			       // Overlay a Boxcar Filter on the 2D Vector
			       if( (filterSize > 0 ) && (filterSize % 2 != 0) ){
				       for(int z = 0; z < filterPass; z++){
					       double sum = 0.0;	// Sum of the current points in the filter
					       double avg;		// Average of the current points in the filter
					       bool done = false;	// Process End boolean

					       int endIndex = filterSize;
					       int startIndex = numpoints-1;

					       while(startIndex < numpoints){
						       for(int i = startIndex; i < endIndex; i++){

							       // Treats the vector as a circular structure
							       // Overflow returns to the beginning of the vector
							       if(i > numpoints){
								       sum += points[i - numpoints].y;
							       }
							       else{
								       sum += points[i].y;
							       }
							       // No underflow occurs because the starting value is 1	
						       }
						       avg = sum / double(filterSize);

						       // Overflow applies to the top limit of the filter range as well
						       if(startIndex+1 > numpoints){
							       points[startIndex - numpoints].y = avg;
						       		startIndex = 0;
							}
						       else{
							       points[startIndex+1].y = avg;
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
				filterName = "Savitzky Golay Filter";
			       // Overlay a Savitzky-Golay Filter on the 2D Vector
			       if( (filterSize == 5) || (filterSize == 11) || (filterSize == 17)){

				       int m = filterSize/2;		// Filtered Element at Middle Index
				       double norm;			// Normalization Factor
				       double sg[filterSize] = {};	// Filter Weights Array

					// Weight Arrays and Normalization Factors for the Corresponding Savitzky Filter Size
				       if( filterSize == 5){  
					       for(int i = 0; i < 5; i++) { sg[i] = sgfilters.five[i]; } 
					       norm = 35;
				       }
				       if( filterSize == 7){  
					       for(int i = 0; i < 7; i++) { sg[i] = sgfilters.seven[i]; }
					       norm = 21;
				       }
				       if( filterSize == 11){  
					       for(int i = 0; i < 11; i++) { sg[i] = sgfilters.eleven[i]; }
					       norm = 429;
				       }

					// Run the filter filterPass number of times over the entire array
					// Filtering each point as the middle index of the filter size
				       for(int z = 0; z < filterPass; z++){
					       for(int j = 0; j < numpoints; j++){
						       double ci, next, y;
						       for(int i = 0; i < filterSize; i++){
							       if( j+i > numpoints){
								       next = points[j + i - numpoints].y;
							       }
							       else{
								       next = points[j+i].y;
							       }

							       ci = sg[i];	// Filter Weight Value
							       y += (ci * next);
						       }
						       if(j+m > numpoints){
							       points[j+m - numpoints].y = y / norm;
						       }
						       else{
							       points[j+m].y = y / norm;
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

	// ************************ Fit the Natural Cubic Spline **************************
	// Coefficients of the Natural Cubic Spline stored in Spline array "Spline" with 
	// each coefficient a double property of the spline struct

	int m = numpoints - 1;

	spline spline[numpoints];
	double h[m], alpha[m], l[numpoints], u[numpoints], z[numpoints];
	double a[numpoints], b[numpoints], c[numpoints], d[numpoints];

	int i, J;
	for(i = 0; i < numpoints; i++){
		a[i] = points[i].y;
	}

	for (i = 1; i <= m; i++){
		h[i-1] = (points[i].x - points[i-1].x);
	}

	alpha[0] = 0.0;
	for (i = 1; i < m; i++) {
		alpha[i] =  (a[i+1] - a[i]) * 3.0 / h[i] - (a[i] - a[i-1]) * 3.0 / h[i-1];
	}

	// L, U, and Z matrices for Tridiagonal Matrix Solving
	l[0] = 1.0;
	u[0] = 0.0;
	z[0] = 0.0;

	for (i = 1; i < m; i++) {
		l[i] = 2.0 * (points[i+1].x - points[i-1].x) - h[i-1] * u[i-1];
		u[i] = h[i] / l[i];
		z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
	}

	l[m] = 1.0;
	u[m] = 0.0;
	spline[m].c = 0.0;

	for (J = numpoints - 2; J >= 0; J--) {
		c[J] = z[J] - u[J] * c[J+1];
		b[J] = (a[J+1] - a[J]) / h[J] - h[J] * (c[J+1] + 2.0 * c[J]) / 3.0;
		d[J] = (c[J+1] - c[J]) / (3.0 * h[J]);
	} 

	for(int i = 0; i < numpoints; i++){
		spline[i] = {a[i], b[i], c[i], d[i], points[i].x};
	}

/*	
	ofile.open("splinePts");
	vector<point> interp = interpolateSpline(points, numpoints, spline, baseline);
	for(int i = 0; i < interp.size(); i++){
		ofile << interp[i].x << " " << interp[i].y << endl;
	}
*/
	// ************************** Find the Roots ************************************
	// Uses Bisection method on the prev and curr root endpoints to find the exact value
	// of the root. If no root is found (NAN), no root is added to the roots point vector.
	
	vector<double> roots;
	for(int i = 0; i < numpoints-1; i++){
		double root1 = NAN;
		// One point must be above and the other below in either variation
		// (first below and second above OR first above and second below)
		if(points[i].y < baseline && points[i+1].y > baseline){
			// Using Bisection Method
			root1 = findRoot( points[i], points[i+1],spline, numpoints, baseline, tolerance);
		}

		if(points[i+1].y < baseline && points[i].y > baseline){
			// Using Bisection Method
			root1 = findRoot( points[i+1], points[i], spline, numpoints, baseline, tolerance);
		}

		// Push all real roots
		if( !isnan(root1)){
			roots.push_back(root1);
		}	
	}

	// *************************** INTEGRATE FOR AREA ***************************************
	// Four Integration Types are Offered
	// Composite Simpson, Romberg, Adaptive Quadrature, Gaussian Quadrature
	
	switch(integration){
		case 0:{
			integrationName = "Composite Simpson Integration";
			break;
		}

		case 1:{
			integrationName = "Romberg Integration";
			break;
		}

		case 2:{
			integrationName = "Adaptive Quadrature Integration";
			break;
		}

		case 3:{
			integrationName = "Gaussian Quadrature Integration";
			break;
		}

		default:{
			integrationName = "Invalid Input; No Integration";
			cout << "Invalid Integration Input Type. Choose from the following options:" << endl;
			cout << "Composite Simpson (0)" << endl;
			cout << "Romberg (1)" << endl;
			cout << "Adaptive Quadrature (2)" << endl;
			cout << "Gaussian Quadrature (3)" << endl;
			break;
		}
	}
	double area;
	// The roots point vector holds pairs of roots at either side of a peak
	// Each pair of roots is used to find the area of the peak
	// Composite Simspon
	for( int i = 0; i < roots.size(); i+=2){
		area = compositeSimpson(roots[i], roots[i+1], spline, numpoints, baseline);
		cout << "composite area= " << area << endl;
	} cout << endl;


/*
	// Romberg
	for( int i = 0; i < roots.size(); i+=2){
		area = romberg(roots[i], roots[i+1], spline, numpoints, baseline, tolerance);
		cout << "romberg area= " << area << endl;
	} cout << endl;


	// Adaptive Quadrature
	for( int i = 0; i < roots.size(); i+=2){
		area = adaptive( roots[i], roots[i+1], spline, numpoints, tolerance, baseline);
		cout << "adaptive area= " << area << endl;
	} cout << endl;
*/

	// Guassian Quadrature
	for( int i = 0; i < roots.size(); i+=2){
		area = guassian(roots[i], roots[i+1], spline, numpoints, baseline);
		cout << "guassian area= " << area << endl;
	} cout << endl;

/*
	
	double top;
	for(int i = 0; i < roots.size(); i+=2){
		point p1 = roots[i];
		point p2 = roots[i+1];	
		top = findTop(p1, p2, spline, numpoints, baseline);
		cout << "top= " << top << endl;
	}
*/

	clock_t end = clock();

	// ************************ OUTPUT RESULTS TO ANALYSIS.DAT ************************************
	
	// Program Options as Outlined in the nmr.in file:
	// Baseline Adjustment | Tolerance | Filtering Name | Filter Size | Filter Passes
	
	ofile.open(outFile);
	ofile << "Program Options" << endl;
	ofile << "===============" << endl;
	ofile << "Baseline Adjustment:\t" << baseline << endl;
	ofile << "Tolerance:\t" << tolerance << endl;
	ofile << filterName << endl;
	ofile << "Filter Size:\t" << filterSize << endl;
	ofile << "Filter Passes:\t" << filterPass << endl;
	ofile << endl;

	// Integration Method -- Chosen	
	ofile << "Integration Method" << endl;
	ofile << "==================" << endl;
	ofile << integrationName << endl;
	ofile << endl;

	// Plot File Data:
	// File Name | Plot Shift for Tms
	ofile << "Plot File Data" << endl;
	ofile << "==============" << endl;
	ofile << filename << endl;
	ofile << "Plot Shifted " << tms << " ppm for TMS calibration" << endl;
	ofile << endl;
	ofile << endl;
	ofile << endl;

	// Peak Number | Peak Start X | Peak End X | Peak Manifold X | Peak Top Y | Peak Area | Whole Number Ratio
	ofile << "Peak\t" << "Begin\t" << "End\t" << "Location\t" << "Top\t" << "Area\t" << "Hydrogens" << endl;
	ofile << "===================================================================" << endl;

	ofile << endl;
	ofile << endl;
	
	// Analysis Time
	ofile << "Anaylsis took " << (double)(end - start)/CLOCKS_PER_SEC << " seconds." << endl;	
	ofile << "(From TMS Calibration through Integration)" << endl;

	ofile.close();
	return 0;
}

double findTop( double x1, double x2, spline s[], int numpoints, double baseline){

	int INDEX;

	for(int i = 0; i < numpoints; i++){
		if( x1 >= s[i].x && x1 < s[i+1].x){
			INDEX = i;
		}
	}
	double top = eqn(x1, s[INDEX], baseline);
	double tmp;
	for(double i = x1; i < x2; i += 0.000001){
		tmp = eqn(i, s[INDEX], baseline);
		if( tmp > top){
			top = tmp;
		}
	}
	return top;
}
// Guassian Quadrature Method using Pg. 233 from the texbook
double guassian( double a, double b, spline spline[], int numpoints, double baseline){

	// With 2 points, Gaussian Coefficients are 1.0000
	double x = (b+a)/2;
	double r[2] = { sqrt(3)/3, -sqrt(3)/3 };
	double w = 9.999999999999996E-001;

/*
	double sum = 0.0;
	for( int h = 0; h < 2; h++ ){
		double tmp = ((b-a)/2) * r[h] + ((b+a)/2);
		for(int i = 0; i < numpoints; i++ ){
			if( tmp >= spline[i].x && tmp < spline[+1].x ){	
				sum += w * eqn( tmp, spline[i], baseline );
			}
		}
	}
*/
	double t1, t2;
	for( int i = 0; i < numpoints; i++ ){
		if( x >= spline[i].x && x < spline[i+1].x ){
			t1 = eqn( ((b-a)*r[0] + b + a)/2, spline[i], baseline) * (b-a)/2;
			t2 = eqn( ((b-a)*r[1] + b + a)/2, spline[i], baseline) * (b-a)/2;	
		}
	}

	return t1 + t2;
}

// Adaptive Quadrature Method using Pg. 224 (Algorithm 4.3) from the text book
double adaptive( double AA, double BB, spline spline[], int numpoints, double tolerance, double baseline){

	int levsize = 50;
	double tol[levsize], a[levsize], h[levsize], fa[levsize], fc[levsize], fb[levsize], s[levsize], v[7];
	int l[levsize];
	double APP, fd, fe, s1, s2;
	int i, n, level;
	
	h[1] = 0.5 * (BB - AA);
	int INDEXA, INDEXB, INDEXAH;
	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints; i++){
		if( AA >= spline[i].x && AA < spline[i+1].x){
			INDEXA = i;
		}
		if( BB >= spline[i].x && BB < spline[i+1].x){
			INDEXB = i;
		}
		if( AA + h[1] >= spline[i].x && AA + h[1] < spline[i+1].x ){
			INDEXAH = i;
		}
	}

	// STEP ONE
	n = 500;
	APP = 0.0;
	i = 1;
	tol[i] =  tolerance;
	a[i] = AA;
	fa[i] = eqn( AA, spline[INDEXA], baseline);
	fc[i] = eqn( AA + h[i], spline[INDEXAH], baseline);
	fb[i] = eqn( BB, spline[INDEXB], baseline);	

	// APPROXIMATION FROM SIMPSONS METHOD FOR ENTIRE INTERVAL
	s[i] = h[i] * (fa[i] + 4.0 * fc[i] + fb[i]) / 3.0;
	l[i] = 1;

	bool oklevel = true;

	// STEP TWO
	while( ( i > 0 ) && ( oklevel ) ){
		// STEP THREE
		double half = a[i] + 0.5 * h[i];
		double onehalf = a[i] + 1.5 * h[i];
		for( int j = 0; j < numpoints; j++ ){
			if( half >= spline[j].x && half < spline[j+1].x ){
				fd = eqn( half, spline[j], baseline);
			}
			if( onehalf >= spline[j].x && onehalf < spline[j+1].x ){
				fe = eqn( onehalf, spline[j], baseline);
			}
		}
	
		// APPROXIMATION FROM SIMPSONS METHOD FOR HALVES OF INTERVALS
		s1 = h[i] * (fa[i] + 4.0 * fd + fc[i]) / 6.0;
		s2 = h[i] * (fc[i] + 4.0 * fe + fb[i]) / 6.0;

		// SAVE DATA AT THIS LEVEL
		v[0] = a[i];
		v[1] = fa[i];
		v[2] = fc[i];
		v[3] = fb[i];
		v[4] = h[i];
		v[5] = tol[i];
		v[6] = s[i];
		level = l[i];
	
		// STEP FOUR: DELETE THE CURRENT LEVEL
		i--;
		if( fabs( s1 + s2 - v[6]) < v[5]){
			APP += s1 + s2;
		}
		else{
			if( level >= n ){
				oklevel = false;
			}
			else{
				// ADD ONE LEVEL: DATA FOR RIGHT HAVE SUBINTERVAL
				i++;
				a[i] = v[0] + v[4];
				fa[i] = v[2];
				fc[i] = fe;
				fb[i] = v[3];
				h[i] = 0.5 * v[4];
				tol[i] = 0.5 * v[5];
				s[i] = s2;
				l[i] = level + 1;

				// DATA FOR LEFT HALF SUBINTERVAL
				i++;
				a[i] = v[0];
				fa[i] = v[1];
				fc[i] = fd;
				fb[i] = v[2];
				h[i] = h[i-1];
				tol[i] = tol[i-1];
				s[i] = s1;
				l[i] = l[i-1];
			}
		}

	}
	return APP;
}

// Romberg Integration using the algorithm outlined on Pg. 216 (Algorithm 4.2) in the textbook
double romberg (double A, double B, spline spline[], int numpoints, double baseline, double tolerance){

	int INDEXA;
	int INDEXB;

	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints; i++){
		if( A >= spline[i].x && A < spline[i+1].x){
			INDEXA = i;
		}
		if( B >= spline[i].x && B < spline[i+1].x){
			INDEXB = i;
		}
	}

	double APP;
	vector< double > R1;
	vector< double > R2;
	double h = (B - A);

	R1.push_back((h / 2) * (eqn( A, spline[INDEXA], baseline) + eqn( B, spline[INDEXB], baseline)));

	int i = 2;
	bool withinRange = false;
	while( !withinRange){

		APP = R1[0] / 2;
	
		for( int j = 1; j < pow(2,i-2)+1; j++){
			for(int k = 0; k < numpoints; k++){
				if( (A + (j - 0.5)*h) >= spline[k].x && (A + (j - 0.5)*h) < spline[k+1].x){
					APP += h * eqn(A + (j - 0.5)*h, spline[k], baseline) / 2.0;
				}
			}
		}
		if( R2.size() > 0){ R2[0] = APP; }
		else{ R2.push_back(APP);}

		for( int j = 2; j < i + 1; j++){
			APP = R2[j - 2] + (R2[j - 2] - R1[j- 2]) / (pow(4, j-1)-1);

			if( R2.size() > j){ R2[j-1] = APP; }
			else { R2.push_back(APP); }
		}

		h = h/2;

		for(int j = 0; j < i; j++){
			if( R1.size() > j){ R1[j] = R2[j]; }
			else{ R1.push_back(R2[j]); }
		}
		if( fabs(R1[R1.size() - 1] - R1[R1.size() - 2]) < tolerance && i > 2){
			withinRange = true;
		}
		i++;
	}
	return  R1[R1.size()-1];
}
// Composite Simpson Integration uses the algorithm outlined on Pg. 205 (Algorithm 4.1) in the textbook.
double compositeSimpson( double a, double b, spline spline[], int numpoints, double baseline){


	int ai, bi, ii;

	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints; i++){
		if( a >= spline[i].x && a < spline[i+1].x ){
			ai = i;
		}
		if( b >= spline[i].x && b < spline[i+1].x ){
			bi = i;
		}
	}


	double x;
	double xi = 0.0;
	double h = (b - a) / 50;
	double xi0 = eqn(a, spline[ai], baseline) + eqn( b, spline[bi], baseline);

	double xi1 = 0.0;
	double xi2 = 0.0;

	for(int j = 1; j < 50; j++){

		x = a + j*h;
		for(int i = 0; i < numpoints-1; i++){
			if( x >= spline[i].x && x < spline[i+1].x ){
				ii = i;
			}
		}
		if( j%2 == 0){
			xi2 = xi2 + eqn( x, spline[ii], baseline);
		}
		else{
			xi1 = xi1 + eqn( x, spline[ii], baseline);
		}
	}
	xi = h * (xi0 + (2*xi2) + (4*xi1)) / 3;
	return xi;
}

// Interpolates on the entire cubic spline to produce an interpolated points array necessary for plotting
vector<point> interpolateSpline ( point points[], int numpoints, spline spline[], double baseline ){

	// Interpolated points array
	vector<point> interp;

	// Interpolate on the Natural Cubic Spline
	double k;
	for(int i = 0; i < numpoints; i++){
		// Interpolate given value increase OR Decrease
		for(k = points[i].x; k < points[i+1].x; k = k + 0.001){
			point p = {k, eqn(k, spline[i], 0.0)};
			interp.push_back(p);
		}
	}
	return interp;
}


// Uses Bisection Root Finding as outlined on Pg. 49 (Algorithm 2.1) in the textbook.
// Between root endpoints r1 and r2 with spline equation s
double findRoot( point r1, point r2, spline s[], int numpoints, double baseline, double tolerance){

	double a = r1.x;
	double b = r2.x;
	double fa, p, fp;

	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints-1; i++){
		if( r1.x >= s[i].x && r1.x < s[i+1].x ){
			fa = eqn( a, s[i], baseline);
		}
	}

	for(int biIndex = 0; biIndex < 10000; biIndex++){

		p = a + (b - a) / 2;

		// Retrieves the appropriate spline index based on the 
		// x value to be calculated within the x value range
		// of each spline piece
		for(int i = 0; i < numpoints-1; i++){
			if( p >= s[i].x && p < s[i+1].x ){
				fp = eqn(p, s[i], baseline);
			}
		}

		if( (fp == 0) || ( ((b - a) / 2) < tolerance) ){
			return p;
		}

		if( fa * fp > 0){
			a = p;
			fa = fp;
		}
		else{
			b = p;
		}
	}
	return NAN;
}

// Defines the sorting process of the points array for increasing 
// x value.
bool sortIncrX( point p1, point p2){

	return (p1.x < p2.x);
}

// Returns the y value of the spline given double x, the corresponding
// spline segment, and the baseline
double eqn(double x, spline s, double baseline){

	return ( s.a + s.b * (x - s.x) + s.c * pow(x - s.x, 2)  +  s.d * pow(x - s.x, 3) ) - baseline;
}

// Returns the highest x value whose y is above the baseline
double findtmsPeak(double baseline, point points[], int numpoints){
	int tmsIndex = -1;

	for(int i = 0; i < numpoints; i++){
		if( (points[i].y > baseline) ){
			if ( ((tmsIndex != -1) && (points[i].x > points[tmsIndex].x)) || (tmsIndex == -1) ){
				tmsIndex = i;
			}
		}
	}	
	return points[tmsIndex].x;
}
