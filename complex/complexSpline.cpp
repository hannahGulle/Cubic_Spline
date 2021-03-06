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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

// Holds the double X and Y values of a point
struct point{
	gsl_complex x;
	gsl_complex y;
};

// Holds the coefficients of a cubic spline 
// segment (a, b, c, d) and the corresponding
// double x as the lower range x value
struct spline{
	gsl_complex a;
	gsl_complex b;
	gsl_complex c;
	gsl_complex d;
	gsl_complex x;
};

// Holds the 3 possible savitzky golay filter
// size arrays
struct sg{
	gsl_complex five[5];
	gsl_complex seven[7];
	gsl_complex eleven[11];
};

using namespace std;

// AREA INTEGRATION METHODS
gsl_complex romberg( gsl_complex a, gsl_complex b, spline spline[], int numpoints, gsl_complex baseline, gsl_complex tolerance);
gsl_complex adaptive( gsl_complex AA, gsl_complex BB, spline spline[], int numpoints, gsl_complex tolerance, gsl_complex baseline);
gsl_complex compositeSimpson( gsl_complex A, gsl_complex B, spline spline[], int numpoints, gsl_complex baseline);
gsl_complex guassian( gsl_complex a, gsl_complex b, spline spline[], int numpoints, gsl_complex baseline);

// CUBIC SPLINE METHODS
vector<point> interpolateSpline(point points[], int numpoints, spline spline[], gsl_complex baseline);
gsl_complex eqn(gsl_complex x, spline s, gsl_complex baseline);

// ROOT FINDING, SORTING, AND TMS PEAK METHODS
gsl_complex findRoot( point r1, point r2, spline s[], int numpoints, gsl_complex baseline, gsl_complex tolerance);
bool sortIncrX( point p1, point p2);
gsl_complex findtmsPeak(gsl_complex baseline, point points[], int numpoints);
gsl_complex findTop( gsl_complex x1, gsl_complex x2, spline s[], int numpoints, gsl_complex baseline);

int main(int argc, char* argv[]){

	// Set precision of all numeric output to 12 total digits
	cout.precision(12);

	// Possible savitzky golay filter sizes with corresponding weights
	sg sgfilters = { 
	{ gsl_complex_rect(-3.0,0), gsl_complex_rect(12.0,0), gsl_complex_rect(17.0,0), gsl_complex_rect(12.0,0), gsl_complex_rect(-3.0,0)}, 
	{ gsl_complex_rect(-2.0,0), gsl_complex_rect(3.0,0), gsl_complex_rect(6.0,0), gsl_complex_rect(7.0,0), gsl_complex_rect(6.0,0), gsl_complex_rect(3.0,0),
		gsl_complex_rect(-2.0,0)}, 
	{ gsl_complex_rect(-36.0,0), gsl_complex_rect(9.0,0), gsl_complex_rect(44.0,0), gsl_complex_rect(69.0,0), gsl_complex_rect(84.0,0),
		gsl_complex_rect(89.0,0), gsl_complex_rect(84.0,0), gsl_complex_rect(69.0,0), gsl_complex_rect(44.0,0), gsl_complex_rect(9.0,0), gsl_complex_rect(-36.0,0)} };

	// Total Number of Points in NMR File
	int numpoints = 0;

	// NMR Line Variables
	string dataFile, outFile, filterName, integrationName;
	double base, tol;
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
	infile >> base;
	infile >> tol;
	infile >> filterType;	
	infile >> filterSize;
	infile >> filterPass;
	infile >> integration;
	infile >> outFile;
	infile.close();

	// Convert Doubles to GNU Scientific Complex Numbers
	gsl_complex baseline = gsl_complex_rect( base, 0 );
	gsl_complex tolerance = gsl_complex_rect( tol, 0 );

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

	point p;
	gsl_complex cx, cy;
	double x, y;

	for(int i = 0; i < numpoints; i++){
		infile >> x >> y;
		p = { gsl_complex_rect( cx, 0 ), gsl_complex_rect( cy, 0 ) };
		points[i] = p;
	}
	infile.close();

	// Sorts the points array by increasing X value
	sort( points, points + numpoints, sortIncrX );

	// ********************** TMS STEP ************************************************
	// Finds the highest Y value above the baseline with the highest x value
	// then subtracts the tms point x value from all x values.
	ofile.open("tms");
	gsl_complex tms = findtmsPeak(baseline, points, numpoints);
	for(int i = 0; i < numpoints; i++){
		points[i].x = gsl_complex_add(points[i].x, gsl_complex_negative(tms));
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
					       gsl_complex sum = gsl_complex_rect( 0.0, 0 );	// Sum of the current points in the filter
					       gsl_complex avg;		// Average of the current points in the filter
					       bool done = false;	// Process End boolean

					       int endIndex = filterSize;
					       int startIndex = numpoints-1;

					       while(startIndex < numpoints){
						       for(int i = startIndex; i < endIndex; i++){

							       // Treats the array as a circular structure
							       // Overflow returns to the beginning of the array
							       if(i > numpoints){
								       sum = gsl_complex_add( sum, points[i - numpoints].y );
							       }
							       else{
								       sum = gsl_complex_add( sum, points[i].y );
							       }
							       // No underflow occurs because the starting value is 1	
						       }
						       avg = sum / gsl_complex(filterSize);

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
				       gsl_complex norm;			// Normalization Factor
				       gsl_complex sg[filterSize] = {};	// Filter Weights Array

					// Weight Arrays and Normalization Factors for the Corresponding Savitzky Filter Size
				       if( filterSize == 5){  
					       for(int i = 0; i < 5; i++) { sg[i] = sgfilters.five[i]; } 
					       norm = gsl_complex_rect( 35, 0 );
				       }
				       if( filterSize == 7){  
					       for(int i = 0; i < 7; i++) { sg[i] = sgfilters.seven[i]; }
					       norm = gsl_complex_rect( 21, 0 );
				       }
				       if( filterSize == 11){  
					       for(int i = 0; i < 11; i++) { sg[i] = sgfilters.eleven[i]; }
					       norm = gsl_complex_rect( 429, 0 );
				       }

					// Run the filter filterPass number of times over the entire array
					// Filtering each point as the middle index of the filter size
				       for(int z = 0; z < filterPass; z++){
					       for(int j = 0; j < numpoints; j++){
						       gsl_complex ci, next, y;
						       for(int i = 0; i < filterSize; i++){
							       if( j+i > numpoints){
								       next = points[j + i - numpoints].y;
							       }
							       else{
								       next = points[j+i].y;
							       }

							       ci = sg[i];	// Filter Weight Value
							       y = gsl_complex_add( y, gsl_complex_mul(ci , next) );
						       }
						       if(j+m > numpoints){
							       points[j+m - numpoints].y = gsl_complex_div( y , norm );
						       }
						       else{
							       points[j+m].y = gsl_complex_div( y , norm );
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
			cerr << "Invalid Filter Type; Defaulted to None. Choices:" << endl;
			cerr << "None (0), Boxcar (1), Savitzky-Golay (2)" << endl;
			break;
			}

	}

	// ************************ Fit the Natural Cubic Spline **************************
	// Coefficients of the Natural Cubic Spline stored in Spline array "Spline" with 
	// each coefficient a double property of the spline struct
	// Using the algorithm outlined in ALG034.C on page 147 of the textbook.

	int m = numpoints - 1;

	spline spline[numpoints];
	gsl_complex h[m], alpha[m], l[numpoints], u[numpoints], z[numpoints];
	gsl_complex a[numpoints], b[numpoints], c[numpoints], d[numpoints];

	int i, J;
	// STEP ZERO: SET A VALUES TO Y VALUES OF THE INITIAL DATA
	for(i = 0; i < numpoints; i++){
		a[i] = points[i].y;
	}

	// STEP ONE: SET H VALUES
	for (i = 1; i <= m; i++){
		h[i-1] = gsl_complex_add(points[i].x, gsl_complex_negative(points[i-1].x) );
	}

	// STEP TWO: SET ALPHA VALUES
	alpha[0] = 0.0;
	for (i = 1; i < m; i++) {
		alpha[i] =  gsl_complex_add(gsl_complex_div(gsl_complex_mul(gsl_complex_add(a[i+1], gsl_complex_negative(a[i])), 3.0), h[i]), gsl_complex_negative(gsl_complex_div(gsl_complex_mul(gsl_complex_add(a[i], gsl_complex_negative(a[i-1])), 3.0), h[i-1])));
	}

	// STEP THREE
	// L, U, and Z matrices for Tridiagonal Matrix Solving
	// USES METHOD DESCRIBED IN ALGORITHM 6.7
	l[0] = gsl_complex_rect(1.0, 0);
	u[0] = gsl_complex_rect(0.0, 0);
	z[0] = gsl_complex_rect(0.0, 0);

	// STEP FOUR
	for (i = 1; i < m; i++) {
		l[i] = gsl_complex_add(gsl_complex_mul(gsl_complex_rect(2.0,0), gsl_complex_add(points[i+1].x, gsl_complex_negative(points[i-1].x))), gsl_complex_negative(gsl_complex_mul(h[i-1], u[i-1])));
		u[i] = gsl_complex_div(h[i], l[i]);
		z[i] = gsl_complex_div(gsl_complex_add(alpha[i], gsl_complex_negative(gsl_complex_mul(h[i-1], z[i-1]))), l[i]);
	}

	// STEP FIVE
	l[m] = gsl_complex_rect(1.0,0);
	u[m] = gsl_complex_rect(0.0,0);
	c[m] = gsl_complex_rect(0.0,0);

	// STEP SIX: SET B, C, D COEFFICIENT VALUES
	for (J = numpoints - 2; J >= 0; J--) {
		c[J] = gsl_complex_add(z[J], gsl_complex_negative(gsl_complex_mul(u[J], c[J+1])));
		b[J] = gsl_complex_add(gsl_complex_div(gsl_complex_add(a[J+1],gsl_complex_negative( a[J])), h[J]), gsl_complex_negative(gsl_complex_div(gsl_complex_mul(h[J], gsl_complex_add(c[J+1], gsl_complex_mul(gsl_complex_rect(2.0,0), c[J]), 3.0)))));
		d[J] = gsl_complex_div(gsl_complex_add(c[J+1],gsl_complex_negative( c[J])), gsl_complex_mul(gsl_complex_rect(3.0,0), h[J]));
	} 

	// SET COEFFICIENT VALUES IN THE SPLINE ACCORDINGLY
	for(int i = 0; i < numpoints; i++){
		spline[i] = {a[i], b[i], c[i], d[i], points[i].x};
	}

	// ************************** Find the Roots ************************************
	// Uses Bisection method on the prev and curr root endpoints to find the exact value
	// of the root. If no root is found (NAN), no root is added to the roots point vector.
	
	vector<gsl_complex> roots;
	for(int i = 0; i < numpoints-1; i++){
		gsl_complex root1 = NAN;
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

	gsl_complex area[roots.size()/2];	
	switch(integration){
		// The roots point vector holds pairs of roots at either side of a peak
		// Each pair of roots is used to find the area of the peak
		case 0:{
			integrationName = "Composite Simpson Integration";
			// Composite Simspon Integration
			for( int i = 0; i < roots.size(); i+=2){
				area[i/2] = compositeSimpson(roots[i], roots[i+1], spline, numpoints, baseline);
			}
			break;
		}

		case 1:{
			integrationName = "Romberg Integration";
			// Romberg Integration 
			for( int i = 0; i < roots.size(); i+=2){
				area[i/2] = romberg(roots[i], roots[i+1], spline, numpoints, baseline, tolerance);
			}
			break;
		}

		case 2:{
			integrationName = "Adaptive Quadrature Integration";	
			// Adaptive Quadrature Integration
			for( int i = 0; i < roots.size(); i+=2){
				area[i/2] = adaptive( roots[i], roots[i+1], spline, numpoints, tolerance, baseline);
			}
			break;
		}

		case 3:{
			integrationName = "Gaussian Legendre Quadrature Integration";
			// Guassian Legendre Quadrature Integration
			for( int i = 0; i < roots.size(); i+=2){
				area[i/2] = guassian(roots[i], roots[i+1], spline, numpoints, baseline);
			}
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

	gsl_complex small = area[0];
	int hydrogens[roots.size()/2];
	// Find the smallest area
	for( int i = 0; i < roots.size()/2 ; i++ ){
		if( area[i] < small ){
			small = area[i];
		}
	}
	// Find the integer ratio of the smallest area (hydrogens) to the other areas
	for( int i = 0; i < roots.size()/2 ; i++ ){
		hydrogens[i] = int( area[i] / small );
	}

	// Find the top of each peak
	gsl_complex top[roots.size()/2];
	gsl_complex mid;
	for( int i = 0; i < roots.size(); i+=2 ){
		top[i/2] = findTop(roots[i], roots[i+1], spline, numpoints, baseline);	
	}

	// End Run Time/ Wall Time Calculation
	clock_t end = clock();

	// ************************ OUTPUT RESULTS TO ANALYSIS.DAT ************************************
	
	// Program Options as Outlined in the nmr.in file:
	// Baseline Adjustment | Tolerance | Filtering Name | Filter Size | Filter Passes
	
	ofile.open(outFile);
	ofile << "\t\t -=> NMR ANALYSIS <=-" << endl;
	ofile << endl;

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
	ofile << "File: " << filename << endl;
	ofile << "Plot Shifted " << tms << " ppm for TMS calibration" << endl;
	ofile << endl;
	ofile << endl;
	ofile << endl;

	// Peak Number | Peak Start X | Peak End X | Peak Manifold X | Peak Top Y | Peak Area | Whole Number Ratio
	ofile << "Peak\t" << "Begin\t" << "End\t" << "Location\t" << "Top\t" << "Area\t" << "Hydrogens" << endl;
	ofile << "===================================================================" << endl;
	for( int i = 0; i < roots.size(); i += 2 ){
		ofile << (i/2)+1 << "   " << roots[i] << "   " << roots[i+1] << "   " << (roots[i] + roots[i+1])/2.0 << "   " << top[i/2] << "   " << area[i/2] << "   " << hydrogens[i/2] << endl;
	}

	ofile << endl;
	ofile << endl;
	
	// Analysis Time
	ofile << "Anaylsis took " << (double)(1000 * (end - start)/CLOCKS_PER_SEC) << " milliseconds." << endl;	
	ofile << "(From TMS Calibration through Integration)" << endl;

	ofile.close();
	return 0;
}

// Finds the Peak Y Value between two roots (x values)
gsl_complex findTop( gsl_complex x1, gsl_complex x2, spline s[], int numpoints, gsl_complex baseline){

	gsl_complex top = 0.0;
	gsl_complex tmp;
	for(gsl_complex i = gsl_complex(x1,0) ; i < x2; i += 0.000001){
		for( int j = 0; j < numpoints; j++ ){
			if( i >= s[j].x && i < s[j+1].x ){
				tmp = eqn(i, s[j], baseline);
				if( tmp > top ){
					top = tmp;
				}
			}
		}
	}
	return top;
}
// Gaussian Legendre Quadrature Method using pound's code from https://github.com/drlbs/quadrature.git
// converted to c++
gsl_complex guassian( gsl_complex a, gsl_complex b, spline spline[], int numpoints, gsl_complex baseline){

	// Root Values Array
	gsl_complex x[8] = { gsl_complex_rect(-9.602898564975363E-001,0),
			gsl_complex_rect(-7.966664774136267E-001,0),
			gsl_complex_rect(-5.255324099163290E-001,0),
			gsl_complex_rect(-1.834346424956498E-001,0),
			gsl_complex_rect(1.834346424956498E-001,0),
			gsl_complex_rect(5.255324099163290E-001,0),
			gsl_complex_rect(7.966664774136267E-001,0),
			gsl_complex_rect(9.602898564975363E-001,0) };
	// Weight Values Array
	gsl_complex w[8] = {gsl_complex_rect(1.012285362903706E-001,0),
			gsl_complex_rect(2.223810344533744E-001,0),
			gsl_complex_rect(3.137066458778874E-001,0),
			gsl_complex_rect(3.626837833783621E-001,0),
			gsl_complex_rect(3.626837833783621E-001,0),
			gsl_complex_rect(3.137066458778874E-001,0),
			gsl_complex_rect(2.223810344533744E-001,0),
			gsl_complex_rect(1.012285362903706E-001,0)};	
	
	// Area Sum
	gsl_complex sum = 0.0;
	// Index of gaussian t input in spline function
	int ind;
	// Integrate to find the area
	for( int i = 0; i < 8; i++ ){
		gsl_complex tmp = ((b-a) * x[i] + a + b) / 2.0;
		for( int i = 0; i < numpoints; i++ ){
			if( tmp >= spline[i].x && tmp < spline[i+1].x ){
				ind = i;		
			}
		}
		sum += w[i] * eqn( tmp, spline[ind], baseline ) * (b-a) / 2.0;
	}
	// return the area summation
	return sum;
}

// Adaptive Quadrature Method using Pg. 224 (Algorithm 4.3) from the text book
gsl_complex adaptive( gsl_complex AA, gsl_complex BB, spline spline[], int numpoints, gsl_complex tolerance, gsl_complex baseline){

	int levsize = 50;
	gsl_complex tol[levsize], a[levsize], h[levsize], fa[levsize], fc[levsize], fb[levsize], s[levsize], v[7];
	int l[levsize];
	gsl_complex APP, fd, fe, s1, s2;
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
		gsl_complex half = a[i] + 0.5 * h[i];
		gsl_complex onehalf = a[i] + 1.5 * h[i];
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
gsl_complex romberg (gsl_complex A, gsl_complex B, spline spline[], int numpoints, gsl_complex baseline, gsl_complex tolerance){

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

	gsl_complex APP;
	vector< gsl_complex > R1;
	vector< gsl_complex > R2;
	// STEP ONE
	gsl_complex h = (B - A);

	R1.push_back((h / 2) * (eqn( A, spline[INDEXA], baseline) + eqn( B, spline[INDEXB], baseline)));

	// STEP TWO
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
		// STEP THREE
		if( R2.size() > 0){ R2[0] = APP; }
		else{ R2.push_back(APP);}

		// STEP FOUR
		for( int j = 2; j < i + 1; j++){
			APP = R2[j - 2] + (R2[j - 2] - R1[j- 2]) / (pow(4, j-1)-1);

			if( R2.size() > j){ R2[j-1] = APP; }
			else { R2.push_back(APP); }
		}

		// STEP FIVE
		h = h/2;

		// STEP SIX
		for(int j = 0; j < i; j++){
			if( R1.size() > j){ R1[j] = R2[j]; }
			else{ R1.push_back(R2[j]); }
		}
		
		// IF THE DIFFERENCE BETWEEN THE LAST TWO INDICES OF R1 IS WITHIN TOLERANCE
		// AND THIS IS NOT THE FIRST LOOP THEN THE AREA HAS BEEN ACHIEVED
		if( fabs(R1[R1.size() - 1] - R1[R1.size() - 2]) < tolerance && i > 2){
			withinRange = true;
		}
		i++;
	}
	// RETURN THE AREA
	return  R1[R1.size()-1];
}
// Composite Simpson Integration uses the algorithm outlined on Pg. 205 (Algorithm 4.1) in the textbook.
gsl_complex compositeSimpson( gsl_complex a, gsl_complex b, spline spline[], int numpoints, gsl_complex baseline){


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


	gsl_complex x;
	gsl_complex xi = gsl_complex_rect(0.0, 0);
	// STEP ONE
	gsl_complex h = gsl_complex_div(gsl_complex_add(b, gsl_complex_negative(a)), 50);
	// STEP TWO
	gsl_complex xi0 = eqn(a, spline[ai], baseline) + eqn( b, spline[bi], baseline);
	gsl_complex xi1 = 0.0;
	gsl_complex xi2 = 0.0;

	// STEP THREE
	for(int j = 1; j < 50; j++){

		// STEP FOUR
		x = a + j*h;
		for(int i = 0; i < numpoints-1; i++){
			if( x >= spline[i].x && x < spline[i+1].x ){
				ii = i;
			}
		}
		// STEP FIVE
		if( j%2 == 0){
			xi2 = gsl_complex_add(xi2, eqn( x, spline[ii], baseline));
		}
		else{
			xi1 = gsl_complex_add(xi1, eqn( x, spline[ii], baseline));
		}
	}
	// AREA HAS BEEN ACHIEVED; RETURN IT
	xi = gsl_complex_div(gsl_complex_mul(h, gsl_complex_add((xi0, gsl_complex_add(gsl_complex_mul(2,xi2), gsl_complex_mul(4,xi1))))), 3);
	return xi;
}

// Interpolates on the entire cubic spline to produce an interpolated points array necessary for plotting
vector<point> interpolateSpline ( point points[], int numpoints, spline spline[], gsl_complex baseline ){

	// Interpolated points array
	vector<point> interp;

	// Interpolate on the Natural Cubic Spline
	gsl_complex k;
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
gsl_complex findRoot( point r1, point r2, spline s[], int numpoints, gsl_complex baseline, gsl_complex tolerance){

	gsl_complex a = r1.x;
	gsl_complex b = r2.x;
	gsl_complex fa, p, fp;

	// Retrieves the appropriate spline index based on the 
	// x value to be calculated within the x value range
	// of each spline piece
	for(int i = 0; i < numpoints-1; i++){
		if( r1.x >= s[i].x && r1.x < s[i+1].x ){
			fa = eqn( a, s[i], baseline);
		}
	}

	for(int biIndex = 0; biIndex < 10000; biIndex++){

		// STEP ONE
		p = gsl_complex_add(a, gsl_complex_div(gsl_complex_add(b, gsl_complex_negative(a)), 2));

		// Retrieves the appropriate spline index based on the 
		// x value to be calculated within the x value range
		// of each spline piece
		for(int i = 0; i < numpoints-1; i++){
			if( p >= s[i].x && p < s[i+1].x ){
				fp = eqn(p, s[i], baseline);
			}
		}
		// STEP TWO
		// IF THE MID POINT BETWEEN THE ENDPOINTS IS WITHIN TOLERANCE
		// OR IF THE ROOT IS AT ZERO, THEN THE ROOT HAS BEEN FOUND
		if( (fp == 0) || ( (gsl_complex_div(gsl_complex_add(b, gsl_complex_negative(a)), 2)) < tolerance) ){
			return p;
		}
		// STEP THREE
		if( gsl_complex_mul(fa, fp) > 0){
			a = p;
			fa = fp;
		}
		else{
			b = p;
		}
	}
	// IF NO ROOT IS FOUND FOR ALL ITERATIONS OF THE BISECTION PROCESS,
	// NAN IS RETURNED.
	return NAN;
}

// Defines the sorting process of the points array for increasing 
// x value.
bool sortIncrX( point p1, point p2){

	return (p1.x < p2.x);
}

// Returns the y value of the spline given double x, the corresponding
// spline segment, and the baseline
gsl_complex eqn(gsl_complex x, spline s, gsl_complex baseline){

	return gsl_complex_add(( gsl_complex_add(s.a, gsl_complex_mul(s.b, gsl_complex_add(x, gsl_complex_negative(s.x)))) + s.c * pow(x - s.x, 2)  +  s.d * pow(gsl_complex_add(x, gsl_complex_negative(s.x)), 3) ), gsl_complex_negative(baseline));
}

// Returns the highest x value whose y is above the baseline
gsl_complex findtmsPeak(gsl_complex baseline, point points[], int numpoints){
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
