/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
#include <iostream>

#include "linbox/linbox-config.h"
#include "linbox/field/PID-integer.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/det.h"
#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/algorithms/double-det.h"


using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
	// commentator.setMaxDetailLevel (-1);
	//     commentator.setMaxDepth (-1);
	//     commentator.setReportStream (std::cerr);

	PID_integer ZZ;
    
	ifstream input (argv[1]);
	if (!input) 
		{ cerr << "Error opening matrix file " << argv[1] << endl; 
			return -1; 
		}
	MatrixStream <PID_integer> ms (ZZ, input);
	BlasBlackbox <PID_integer> A (ms);
	cout << "Matrix is " << A.rowdim() << " by " << A.coldim() << endl;

	if (A.rowdim() != A.coldim() + 1){
		cerr<<"Wrong dimensions: A must be (n+1) x n"<<endl;
		exit (-1);
	}
 
	integer det1, det2;
	Timer tim; tim.clear();
	tim.start();
	bool proof=false;
	if (argc > 2)
		proof = true;

	//doubleDet (det1, det2, A, true);
	doubleDet (det1, det2, A, proof);
	tim.stop();
     
	// Check solution
	size_t n = A.coldim();
	BlasBlackbox<PID_integer> B (ZZ, n, n);
	for (size_t i=0; i<n-1; ++i)
		for (size_t j=0; j<n; ++j)
			B.setEntry (i,j,A.getEntry(i,j));
	for (size_t j=0; j<n; ++j)
		B.setEntry (n-1,j,A.getEntry(n-1,j));

	integer db, dc;
	det (db, B);
	for (size_t j=0; j<n; ++j)
		B.setEntry(n-1,j,A.getEntry(n,j));
	det(dc, B);
 
	if ((det1 != db) || (det2 != dc))
		std::cerr<<"Check FAIL"<<endl;

	cerr << "Double Det: "<<tim.usertime() << "s" << endl;
}
