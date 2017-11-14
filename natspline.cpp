#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<vector>
#include<algorithm>

struct root{
	double x;
	int interpi;
	int initi;
};

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

bool sortIncrX( point p1, point p2);
double findtmsPeak(double baseline, point points[], int numpoints);
double eqn(double x, spline s);
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
	ofile.open("nattmsdat");
	for(int i = 0; i < numpoints; i++){
		points[i].x = points[i].x - tms;
		ofile << points[i].x << " " << points[i].y << endl;
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

	ofile.open("sgdat");
	for(int i = 0; i < numpoints; i++){
		ofile << points[i].x << " " << points[i].y << endl;
	}
	ofile.close();


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
	// Cubic Spline Made of Coefficients in Arrays a, b, c, d
	for(int i = 0; i < numpoints; i++){
		if( fabs( spline[i].a - 143.26) < 0.5){
			cout << "a= " << spline[i].a << " b= " << spline[i].b << " c= " << spline[i].c << " d= " << spline[i].d << endl;
		}
	}

	// ***************************** INTERPOLATION STEP ************************************

	vector<point> interp;

	vector<root> rootEndpts;

	ofile.open("testSpline");

	// Interpolate on the Natural Cubic Spline
	double k;
	for(i = 0; i < numpoints-1; i++){
		// Interpolate given value increase OR Decrease
		for(k = points[i].x; k < points[i+1].x; k = k + 0.001){
			point p = {k, eqn(k, spline[i])};
			interp.push_back(p);
		}
	}

	for(int i = 0; i < interp.size(); i++){
		ofile << interp[i].x << " " << interp[i].y << endl;
	}	
	ofile.close();
	
	// Find the Root Endpoints
	i = 0;
	while(i+1 < interp.size()){
		if( (interp[i].y > baseline && interp[i+1].y < baseline) ||
			(interp[i].y < baseline && interp[i+1].y > baseline)){
		
			root prev, curr;
			for(int j = 0; j < numpoints-1; j++){			

				if( (interp[i].x > points[j].x) && (interp[i].x < points[j+1].x) ){
					prev.initi = j;
					prev.interpi = i;		
				}
	
				if( (interp[i+1].x > points[j].x && interp[i+1].x < points[j+1].x)){
					curr.initi = j+1;
					curr.interpi = i+1;	
				}
			}

			prev.x = interp[i].x;
			curr.x = interp[i+1].x;
			rootEndpts.push_back(prev);
			rootEndpts.push_back(curr);
		}
		i += 2;
	}
	cout << "rootendpoints= " << rootEndpts.size() << endl;
	
	i = 0;
	int numroots = 0;
	vector<root> roots;

	while ( i < rootEndpts.size() ){

		int interpi = rootEndpts[i].interpi;
		int initi = rootEndpts[i].initi;		
		double endpta = rootEndpts[i].x;
		double endptb = rootEndpts[i+1].x;

				int biIndex = 1;
				double fa = interp[rootEndpts[i].interpi].y;
				double p;
				while( biIndex <= 1000){

					p = endpta + (endptb - endpta)/2.0;
					double fp = eqn(p, spline[initi]);

					if( (fp == baseline) || ( (endptb - endpta)/2 < tolerance) ){
						root r;
						r.x = p;
						r.interpi = interpi;
						r.initi = initi;
		//				cout << "root= " << p << " interpi= " << interpi << " initi= " << initi << endl;
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
			cout << "numroots= " << numroots << endl;
			cout << endl;

			// XXX:Composite Simpson OR Trapezoidal
			double x, xi;
			i = 0;
			while( i+1 < numroots){

				double endpointA = roots[i].x;
				double endpointB = roots[i+1].x;
				int INDEXA = roots[i].interpi;
				int INDEXB = roots[i+1].interpi;
				int initia = roots[i].initi;
				int initib = roots[i+1].initi;

				double h = (endpointB - endpointA)/500;
				double xi0 = eqn(endpointA, spline[initia]) + eqn( endpointB, spline[initib]);

				double xi1 = 0.0;
				double xi2 = 0.0;

				for(int j = 1; j < 500; j++){

					x = endpointA + j*h;
					if( j%2 == 0){
						xi2 = xi2 + eqn( x, spline[initia]);
					}
					else{
						xi1 = xi1 + eqn( x, spline[initia]);
					}
				}
				i = i + 2;
				xi = h * (xi0 + (2*xi2) + (4*xi1))/3.0;
				cout << "area= " << xi << endl;
			}
			cout << "end compsimp" << endl; 
			cout << endl;

			//XXX: Romberg
			/*
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

			   R[0][0] = (h / 2.0) * (eqn(A, interpX[INDEXA],  a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]) +
			   eqn( B, interpX[INDEXB], a[INDEXB], b[INDEXB], c[INDEXB], d[INDEXB]));

			   for(int j = 2; j <= n; j++){
			   sum = 0.0;
			   m = exp((i-2) * log(2.0)) + 0.5;
			   for(double k = 1; k <= m; k++){
			   part = A + (h * (k - 0.5));
			   sum += eqn( part, interpX[INDEXA],  a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
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
			FA[I] = eqn(AA, interpX[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
			FA[I] = eqn((AA + H[I]), interpX[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
			FB[I] = eqn(BB, interpX[INDEXB], a[INDEXB], b[INDEXB], c[INDEXB], d[INDEXB]);
			S[I] = H[I] * (FA[I] + 4.0 * FC[I] + FB[I]) / 3.0;
			L[I] = 1;

			while ( I > 0){

				FD = eqn( (A[I] + 0.5 + H[I]), interpX[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
				FE = eqn( (A[I] + 1.5 + H[I]), interpX[INDEXA], a[INDEXA], b[INDEXA], c[INDEXA], d[INDEXA]);
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
*/
return 0;
}

bool sortIncrX( point p1, point p2){

	return (p1.x < p2.x);
}

double eqn(double x, spline s){

	return( s.a + s.b * (x - s.x) + s.c * (x - s.x) * (x - s.x)  +  s.d * (x - s.x) * (x - s.x) * (x - s.x) );
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
