
/** 
 * examples/doubledet.C
 *
 * Copyright (C) 2008, 2010 C. Pernet
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */

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
	Timer tim2;
	tim2.clear();
	tim2.start();
	det(dc, B);
	tim2.stop();

	if ((det1 != db) || (det2 != dc))
		std::cerr<<"Check FAIL"<<endl;

	cout<< "Det 1 = " << det1 << endl;
	cout<< "Det 2 = "<<det2 << endl;
	
	cerr << "Double Det: "<<tim.usertime() << "s" << endl;
	cerr << "Each single Det: "<<tim2.usertime() << "s" << endl;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
