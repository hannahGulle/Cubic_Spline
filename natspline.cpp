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
double romberg( point p1, point p2, spline spline[], int numpoints, double baseline, double tolerance);
double adaptive( point p1, point p2, spline spline[], int numpoints, double tolerance, double baseline);
double compositeSimpson( point p1, point p2, spline spline[], int numpoints, double baseline);
double guassian( point p1, point p2, spline spline[], int numpoints, double tolerance, double baseline);

// CUBIC SPLINE METHODS
vector<point> interpolateSpline(point points[], int numpoints, spline spline[], double baseline);
double eqn(double x, spline s, double baseline);

// ROOT FINDING, SORTING, AND TMS PEAK METHODS
double findRoot( point r1, point r2, spline s[], int numpoints, double baseline, double tolerance);
bool sortIncrX( point p1, point p2);
double findtmsPeak(double baseline, point points[], int numpoints);

int main(int argc, char* argv[]){

	// Possible savitzky golay filter sizes with corresponding weights
	sg sgfilters = { {-3.0,12.0,17.0,12.0,-3.0}, {-2.0,3.0,6.0,7.0,6.0,3.0,-2.0}, {-36.0,9.0,44.0,69.0,84.0,89.0,84.0,69.0,44.0,9.0,-36.0} };

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

	// data file Variables
	string line;

	infile.open(dataFile.c_str());
	while(getline(infile,line)){
		numpoints++;
	}
	infile.close();
	
	// 1D Array Containing Point Structs Containing X and Y Values
	// from the data file
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
	double tms = findtmsPeak(baseline, points, numpoints);
	for(int i = 0; i < numpoints; i++){
		points[i].x = points[i].x - tms;
	}

	// *********************** APPLY FILTERS STEP **************************************
	// Three filters available: None, Boxcar, Savitzky Golay (5, 7, 11)
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
			       // Overlay a Savitzky-Golay Filter on the 2D Vector
			       if( (filterSize == 5) || (filterSize == 11) || (filterSize == 17)){

				       int m = filterSize/2;
				       double norm;
				       double sg[filterSize] = {};

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

							       ci = sg[i];
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

	int i, J;
	for(i = 0; i < numpoints; i++){
		spline[i].a = points[i].y;
		spline[i].x = points[i].x;
	}

	for (i = 1; i <= m; i++){
		h[i-1] = (points[i].x - points[i-1].x);
	}

	alpha[0] = 0.0;
	for (i = 1; i < m; i++) {
		alpha[i] =  (spline[i+1].a - spline[i].a) * 3.0 / h[i] - (spline[i].a - spline[i-1].a) * 3.0 / h[i-1];
	}

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
		spline[J].c = z[J] - u[J] * spline[J+1].c;
		spline[J].b = (spline[J+1].a - spline[J].a) / h[J] - h[J] * (spline[J+1].c + 2.0 * spline[J].c) / 3.0;
		spline[J].d = (spline[J+1].c - spline[J].c) / (3.0 * h[J]);
	} 

	for(int i = 0; i < numpoints; i++){
		cout << "a= " << spline[i].a << " b= " << spline[i].b << " c= " << spline[i].c << " d= " << spline[i].d << endl;
	}

	// ************************** Find the Roots ************************************
	// Uses Bisection method on the prev and curr root endpoints to find the exact value
	// of the root. If no root is found (NAN), no root is added to the roots point vector.
	
	vector<point> roots;
	for(int i = 0; i < numpoints-1; i++){

		point root1;
		root1.x = NAN;
		// One point must be above and the other below in either variation
		// (first below and second above OR first above and second below)
		if(points[i].y < baseline && points[i+1].y > baseline){
			point prev, curr;
			prev.x = points[i].x; curr.x = points[i+1].x;
			// Using Bisection Method
			root1.x = findRoot( prev, curr,spline, numpoints, baseline, tolerance);
		}

		if(points[i+1].y < baseline && points[i].y > baseline){
			point prev, curr;
			prev.x = points[i+1].x; curr.x = points[i].x;
			// Using Bisection Method
			root1.x = findRoot( prev, curr, spline, numpoints, baseline, tolerance);
		}

		// Push all real roots
		if( !isnan(root1.x)){
			roots.push_back(root1);
		}	
	}


	// *************************** INTEGRATE FOR AREA ***************************************
	// Four Integration Types are Offered
	// Composite Simpson, Romberg, Adaptive Quadrature, Gaussian Quadrature
	
	double area;
	// The roots point vector holds pairs of roots at either side of a peak
	// Each pair of roots is used to find the area of the peak
	// Composite Simspon
	for( int i = 0; i < roots.size(); i+=2){
		point p1 = roots[i];
		point p2 = roots[i+1];
		area = compositeSimpson(p1, p2, spline, numpoints, baseline);
		cout << "composite area= " << area << endl;
	} cout << endl;

	// Romberg
	for( int i = 0; i < roots.size(); i+=2){
		point p1 = roots[i];
		point p2 = roots[i+1];
		area = romberg(p1, p2, spline, numpoints, baseline, tolerance);
		cout << "romberg area= " << area << endl;
	} cout << endl;

	// Adaptive Quadrature
	for( int i = 0; i < roots.size(); i+=2){
		point p1 = roots[i];
		point p2 = roots[i+1];
		area = adaptive( p1, p2, spline, numpoints, tolerance, baseline);
		cout << "adaptive area= " << area << endl;
	} cout << endl;

	// Guassian Quadrature
	for( int i = 0; i < roots.size(); i+=2){
		point p1 = roots[i];
		point p2 = roots[i+1];
		area = guassian(p1, p2, spline, numpoints, tolerance, baseline);
		cout << "guassian area= " << area << endl;
	} cout << endl;

	return 0;
}

// Guassian Quadrature Method using Pg. 233 from the texbook
double guassian( point p1, point p2, spline spline[], int numpoints, double tolerance, double baseline){

	// With 2 points, Gaussian Coefficients are 1.0000

	double a = p1.x;
	double b = p2.x;
	double x = (a + b) / 2;
	double top = ((b - a) * (sqrt(3.0)/3.0) + b + a) / 2.0;
	double bot = ((b - a) * (-1 * sqrt(3.0)/3.0) + b + a) / 2.0;
	int INDEX;

	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints-1; i++){
		if( x >= spline[i].x && x < spline[i+1].x){
			INDEX = i;
	//		cout << "spline x= " << spline[i].x << " x= " << x << endl;
		}
	}

	double t1 = eqn(top, spline[INDEX], baseline) * ((b - a)/2);
	double t2 = eqn(bot, spline[INDEX], baseline) * ((b - a)/2);

//	cout << "b-a/2 = " << ((b - a)/2) << endl;
//	cout << "gauss parts: " << t1 << " " << t2 << endl;

	return t1 + t2;
}

// Adaptive Quadrature Method using Pg. 224 (Algorithm 4.3) from the text book
double adaptive( point p1, point p2, spline spline[], int numpoints, double tolerance, double baseline){

	int N = 500;
	double TOL[N], A[N], H[N], FA[N], FC[N], FB[N], S[N], V[8];
	int L[N];
	double AA, BB, APP, FD, FE, S1, S2;
	int INDEXA, INDEXB;

	AA = p1.x;
	BB = p2.x;

	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints-1; i++){
		if( AA >= spline[i].x && AA < spline[i+1].x){
			INDEXA = i;
		}
		if( BB >= spline[i].x && BB < spline[i+1].x){
			INDEXB = i;
		}
	}

	TOL[0] = tolerance;
	A[0] = AA;
	H[0] = 0.5 * (BB - AA);
	FA[0] = (eqn(AA, spline[INDEXA], baseline));
	FC[0] = (eqn((AA + H[0]), spline[INDEXA], baseline));
	FB[0] = (eqn(BB, spline[INDEXB], baseline));
	S[0] = (H[0] * (FA[0] + 4.0 * FD + FC[0]) / 3.0);
	L[0] = 1;

	int I = 1;
	while ( I > 0 && I < N){

		FD = eqn( (A[I-1] + 0.5 * H[I-1]), spline[INDEXA], baseline);
		FE = eqn( (A[I-1] + 1.5 * H[I-1]), spline[INDEXA], baseline);
		S1 = H[I-1] * (FA[I-1] + 4 * FD + FC[I-1]) / 6.0;
		S2 = H[I-1] * (FC[I-1] + 4 * FE + FB[I-1]) / 6.0;

		V[0] = A[I-1];
		V[1] = FA[I-1];
		V[2] = FC[I-1];
		V[3] = FB[I-1];
		V[4] = H[I-1];
		V[5] = TOL[I-1];
		V[6] = S[I-1];
		V[7] = L[I-1];
		I--;

		if( abs(S1 + S2 - V[6]) < V[5]){
			APP += S1 + S2;
		}
		else{
			I++;
			int z = I - 1;
			A[z] = V[0] + V[4];
			FA[z] = V[2];
			FC[z] = FE;
			FB[z] = V[3];
			H[z] = 0.5 * V[4];
			TOL[z] = 0.5 * V[5];
			S[z] = S2;
			L[z] = int(V[7] + 1);

			I++;
			z = I - 1;
			int y = I - 2;
			A[z] = V[0];
			FA[z] = V[1];
			FC[z] = FD;
			FB[z] = V[2];
			H[z] = H[y];
			TOL[z] = TOL[y];
			S[z] = S1;
			L[z] = L[y];
		}
	}
	return APP;
}

// Romberg Integration using the algorithm outlined on Pg. 216 (Algorithm 4.2) in the textbook
double romberg (point p1, point p2, spline spline[], int numpoints, double baseline, double tolerance){

	double A = p1.x;
	double B = p2.x;
	int INDEXA;
	int INDEXB;

	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints-1; i++){
		if( A >= spline[i].x && A < spline[i+1].x){
			INDEXA = i;
		}
		if( B >= spline[i].x && B < spline[i+1].x){
			INDEXB = i;
		}
	}

	double APP = 0.0;
	vector< double > R1;
	vector< double > R2;
	double h = (B - A);

	R1.push_back((h / 2.0) * (eqn( A, spline[INDEXA], baseline) + eqn( B, spline[INDEXB], baseline)));

	int i = 2;
	while( ((R1[R1.size()-1] - R1[R1.size()-2]) < tolerance) && (i > 2)){

		APP = R1[0] / 2;
		for( int j = 1; j < exp(i-2)+1; j++){
			APP += h * eqn(A + (j - 0.5)*h, spline[INDEXA], baseline) / 2.0;
		}
		if( R2.size() > 0){ R2[0] = APP; }
		else{ R2.push_back(APP);}

		for( int j = 2; j < i + 1; j++){
			APP = R2[j - 2] + (R2[j - 2] - R1[j- 2]) / (pow(4, j-1)-1);

			if( R2.size() > j){ R2[j-1] = APP; }
			else { R2.push_back(APP); }
		}

		h = h/2.0;

		for(int j = 0; j < i; j++){
			if( R1.size() > j){ R1[j] = R2[j]; }
			else{ R1.push_back(R2[j]); }
		}
		i++;
	}
	return  R1[R1.size()-1];
}
// Composite Simpson Integration uses the algorithm outlined on Pg. 205 (Algorithm 4.1) in the textbook.
double compositeSimpson( point p1, point p2, spline spline[], int numpoints, double baseline){


	double a = p1.x;
	double b = p2.x;
	int ai, bi, ii;

	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints-1; i++){
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
	for(int i = 0; i < numpoints-1; i++){
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
	for(int i = 0; i < numpoints; i++){
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

	return ( s.a + s.b * (x - s.x) + s.c * (x - s.x) * (x - s.x)  +  s.d * (x - s.x) * (x - s.x) * (x - s.x) ) - baseline;
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
