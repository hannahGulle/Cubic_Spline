#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<vector>
#include<algorithm>
#include<cmath>

struct point{
	double x;
	double y;
};

struct spline{
	double a;
	double b;
	double c;
	double d;
	double x;
};

struct sg{
	double five[5];
	double seven[7];
	double eleven[11];
};

using namespace std;

double romberg( point p1, point p2, spline spline[], int numpoints, double baseline);
double adaptive( point p1, point p2, spline spline[], int numpoints, double tolerance, double baseline);
double compositeSimpson( point p1, point p2, spline spline[], int numpoints, double baseline);
vector<point> interpolateSpline(point points[], int numpoints, spline spline[], double baseline);
double findRoot( point r1, point r2, spline s[], int numpoints, double baseline, double tolerance);
bool sortIncrX( point p1, point p2);
double findtmsPeak(double baseline, point points[], int numpoints);
double eqn(double x, spline s, double baseline);

int main(int argc, char* argv[]){

	// Possible savitzky golay filter sizes
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


	// data.dat Variables
	string line;

	infile.open(dataFile.c_str());
	while(getline(infile,line)){
		numpoints++;
	}


	// 1D Array Containing X and Y Values
	point points[numpoints];

	infile.close();
	infile.open(dataFile.c_str());

	for(int i = 0; i < numpoints; i++){
		point p;
		infile >> p.x >> p.y;
		points[i] = p;
	}
	infile.close();

	sort( points, points + numpoints, sortIncrX );



// ********************** TMS STEP ************************************************

	double tms = findtmsPeak(baseline, points, numpoints);
	for(int i = 0; i < numpoints; i++){
		points[i].x = points[i].x - tms;
	}

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
								       sum += points[i - numpoints].y;
							       }
							       else{
								       sum += points[i].y;
							       }
							       // No underflow occurs because the starting value is 1	
						       }
						       avg = sum / float(filterSize);

						       // Overflow applies to the top limit of the filter range as well
						       if(startIndex+1 > numpoints){
							       points[startIndex - numpoints].y = avg;
						       }
						       else{
							       points[startIndex+1].y = avg;
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
				       double norm;
				       double sg[filterSize] = {};

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

// ************************ Fit the Natural Cubic Spline Step **************************
	int m = numpoints - 1;

	spline spline[numpoints];
	double h[numpoints], alpha[m], l[numpoints], u[numpoints], z[numpoints];

	int i, J;
	for(i = 0; i < numpoints; i++){
		spline[i].a = points[i].y;
		spline[i].x = points[i].x;
	}

	for (i = 0; i <= m; i++){
		h[i] = (points[i+1].x - points[i].x);
	}

	alpha[0] = 0.0;
	for (i = 1; i < m; i++) {
		alpha[i] =  (spline[i+1].a - spline[i].a) * 3.0 / h[i] - (spline[i].a - spline[i-1].a) * 3.0 / h[i-1];
	}

	l[0] = 1.0;
	u[0] = 0.0;
	z[0] = 0.0;

	for (i = 1; i <= m; i++) {
		l[i] = 2.0 * (points[i+1].x - points[i-1].x) - h[i-1] * u[i-1];
		u[i] = h[i] / l[i];
		z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
	}

	l[numpoints-1] = 1.0;
	z[numpoints-1] = 0.0;
	spline[numpoints-1].c = 0.0;

	for (i = 0; i <= m; i++) {
		J = m - i;
		spline[J].c = z[J] - u[J] * spline[J+1].c;
		spline[J].b = (spline[J+1].a - spline[J].a) / h[J] - h[J] * (spline[J+1].c + 2.0 * spline[J].c) / 3.0;
		spline[J].d = (spline[J+1].c - spline[J].c) / (3.0 * h[J]);
	} 


// ***************************** INTERPOLATION STEP ************************************
	vector<point> interp;
	ofile.open("testSpline");
	// Interpolate on the Natural Cubic Spline
	interp = interpolateSpline(points, numpoints, spline, baseline);
	for(int i = 0; i < interp.size(); i++){
		ofile << interp[i].x << " " << interp[i].y << endl;
	}	
	ofile.close();

// ************************** Find the Roots ************************************
	vector<point> roots;
	for(int i = 0; i < interp.size()-1; i++){

		point root1;
		root1.x = NAN;
		if( (interp[i].y < baseline && interp[i+1].y > baseline) ||
				(interp[i+1].y < baseline && interp[i].y > baseline) ){
			point prev, curr;
			prev.x = interp[i].x; curr.x = interp[i+1].x;
			prev.y = interp[i].y; curr.y = interp[i+1].y;	
			root1.x = findRoot( prev, curr, spline, numpoints, baseline, tolerance);
		}

		if( !isnan(root1.x)){
			roots.push_back(root1);
		}	
	}

// *************************** INTEGRATE FOR AREA ***************************************

	double area;
	for( int i = 0; i < roots.size(); i+=2){
		point p1 = roots[i];
		point p2 = roots[i+1];
		area = compositeSimpson(p1, p2, spline, numpoints, baseline);
		cout << "composite area= " << area << endl;
	} cout << endl;

	for( int i = 0; i < roots.size(); i+=2){
		point p1 = roots[i];
		point p2 = roots[i+1];
		area = romberg(p1, p2, spline, numpoints, baseline);
		cout << "romberg area= " << area << endl;
	} cout << endl;

	for( int i = 0; i < roots.size(); i+=2){
		point p1 = roots[i];
		point p2 = roots[i+1];
		area = adaptive( p1, p2, spline, numpoints, tolerance, baseline);
		cout << "adaptive area= " << area << endl;
	} cout << endl;

	return 0;
}

double adaptive( point p1, point p2, spline spline[], int numpoints, double tolerance, double baseline){

	int n = 50;
	double TOL[n], A[n], H[n], FA[n], FC[n], FB[n], S[n], V[7];
	int L[n];
	double AA, BB, EPS, APP, FD, FE, S1, S2;
	int I, LEV, INDEXA, INDEXB;
	int count = 0;

	EPS = tolerance;

	AA = p1.x;
	BB = p2.x;

	for(int i = 0; i < numpoints; i++){
		if( AA >= spline[i].x && AA < spline[i+1].x){
			INDEXA = i;
		}
		if( BB >= spline[i].x && BB < spline[i+1].x){
			INDEXB = i;
		}
	}

	APP = 0.0;
	I = 1;
	TOL[I] = 10.0 * EPS;
	A[I] = AA;
	H[I] = 0.5 * (BB - AA);
	FA[I] = eqn(AA, spline[INDEXA], baseline);
	FA[I] = eqn((AA + H[I]), spline[INDEXA], baseline);
	FB[I] = eqn(BB, spline[INDEXB], baseline);
	S[I] = H[I] * (FA[I] + 4.0 * FC[I] + FB[I]) / 3.0;
	L[I] = 1;

	while ( I > 0){

		FD = eqn( (A[I] + 0.5 + H[I]), spline[INDEXA], baseline);
		FE = eqn( (A[I] + 1.5 + H[I]), spline[INDEXA], baseline);
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
	return APP;
}

double romberg (point p1, point p2, spline spline[], int numpoints, double baseline){

	double A = p1.x;
	double B = p2.x;
	int INDEXA;
	int INDEXB;

	for(int i = 0; i < numpoints; i++){
		if( A >= spline[i].x && A < spline[i+1].x){
			INDEXA = i;
		}
		if( B >= spline[i].x && B < spline[i+1].x){
			INDEXB = i;
		}
	}

	double APP = 0.0;
	int n = 21;

	double sum, part;
	int m, tmp;

	double R[2][n];
	double h = (B - A);

	R[0][0] = (h / 2.0) * (eqn( A, spline[INDEXA], baseline) + eqn( B, spline[INDEXB], baseline));

	for(int i = 2; i <= n; i++){
		for(int j = 2; j <= n; j++){
			sum = 0.0;
			m = exp((i-2) * log(2.0)) + 0.5;
			for(double k = 1; k <= m; k++){
				part = A + (h * (k - 0.5));
				sum += eqn( part, spline[INDEXA], baseline);
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
	}
	return  R[1][n-1];
}

double compositeSimpson( point p1, point p2, spline spline[], int numpoints, double baseline){


	double a = p1.x;
	double b = p2.x;
	int ai, bi, ii;

	for(int i = 0; i < numpoints-1; i++){
		if( a >= spline[i].x && a < spline[i+1].x ){
			ai = i;
		}
		if( b >= spline[i].x && b < spline[i+1].x ){
			bi = i;
		}
	}


	double x, xi;
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
	xi = h * (xi0 + (2*xi2) + (4*xi1)) / 3.0;
	return xi;
}

vector<point> interpolateSpline ( point points[], int numpoints, spline spline[], double baseline ){

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


// Uses Bisection Root Finding
// Between root endpoints r1 and r2 with spline equation s
double findRoot( point r1, point r2, spline s[], int numpoints, double baseline, double tolerance){

	double a = r1.x;
	double b = r2.x;

	int biIndex = 1;
	double fa = r1.y;
	double p, fp;

	while( biIndex <= 1000){

		p = a + (b - a) / 2.0;

		for(int i = 0; i < numpoints; i++){
			if( p >= s[i].x && p < s[i+1].x ){
				fp = eqn(p, s[i], baseline);
			}
		}

		if( (fp == baseline) || ( ((b - a) / 2.0) < tolerance) ){
			return p;
		}

		biIndex++;

		if( fa * fp > 0){
			a = p;
		}
		else{
			b = p;
		}
	}
	return NAN;
}

bool sortIncrX( point p1, point p2){

	return (p1.x < p2.x);
}

double eqn(double x, spline s, double baseline){

	return ( s.a + s.b * (x - s.x) + s.c * (x - s.x) * (x - s.x)  +  s.d * (x - s.x) * (x - s.x) * (x - s.x) ) - baseline;
}

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
