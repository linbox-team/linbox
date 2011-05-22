
/* tests/test-solve.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
[12~ * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/dense.h"
#include "linbox/solutions/charpoly.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Charpoly of a diagonal matrix
 *
 * Constructs a random diagonal rational matrix D 
 * Computes its characteristic polynomial
 * Checks if c[0] (determinant) and c[n-1] (trace) agree 
 * 
 * Return true on success and false on failure
 */

static bool testDiagRatCharpoly (size_t n, int iterations) 
{
	commentator.start ("Testing rational charpoly of diagonal matrix ", "testNonsingularRatIntSolve", iterations);

	bool ret = true;
	int i;
	size_t j;

	GMPRationalField Q;
	SparseMatrix<GMPRationalField > A(Q,n,n);
	DenseMatrix <GMPRationalField > B(Q,n,n);
	std::vector<GMPRationalField::Element> c;

	for (i=0; i < iterations; i++) {
        	GMPRationalField::Element c0,cn;
	        Q.init(c0,1,1);
	        Q.init(cn,0,1);

		
		commentator.startIteration (i);

		for (j=0; j < n; ++j) {
			integer tmp_n, tmp_d;
			GMPRationalField::Element tmp, abstmp;
			tmp_n = (integer) rand() % (5*(i +1)) + 1;
			tmp_d = (integer) rand() % (5*(i +1)) + 1;
			if ( ( i%2) && (j % 2)) integer::negin(tmp_n);
			
			Q.init(tmp, tmp_n,tmp_d);
			
			A.setEntry(j,j,tmp);
			B.setEntry(j,j,tmp);
			

			Q.mulin(c0, tmp);
			Q.addin(cn, tmp);
		}	
		if (n%2==0) Q.negin(cn);
			
		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		charpoly (c, A);
		if ((!Q.areEqual(c[0] , c0)) || (!Q.areEqual(c[n-1] , cn) ) ) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Sparse charpoly failed" << endl;
			ret = false;
		}
		c.clear();

		charpoly (c, B);
		if ((!Q.areEqual(c[0] , c0)) || (!Q.areEqual(c[n-1] , cn) ) ) {
                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                << "ERROR: Dense charpoly failed" << endl;
                        ret = false;
                }
		
		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNonsingularRatIntSolve");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set column dimension of test matrices to N.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};
	parseArguments (argc, argv, args);
	
	commentator.start("Rational solve test suite", "solve");

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (5);
	//commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (4);

    if ( ! testDiagRatCharpoly(n,iterations) ) pass = false;

	commentator.stop("solve test suite");
    //std::cout << (pass ? "passed" : "FAILED" ) << std::endl;

	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
