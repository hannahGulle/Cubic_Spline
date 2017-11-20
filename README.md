Hannah Gulle
Project: Fitting a Cubic Spline
Class: CSC 335
Language: C++

Compile Instructions:
c++ spline.cpp

Run Instructions:
./a.out nmr.in

ALGORITHM:
1.	Open and Read in NMR.in file
2.	Open and Read in data file
3. 	Find TMS peak in array
4.	Build either a Boxcar or Savitzky-Golay Filter
5.	Fit a natural cubic spline
6.	Find all roots of the spline
7.	For each pair of roots, determine the peak position using
			(Xb - Xa)/2
8.	Determine the peak manifold between A and B using one of four integration
methods defined in the NMR.in file.
9.	Find the hydrogen atom count as multiples of the least area
