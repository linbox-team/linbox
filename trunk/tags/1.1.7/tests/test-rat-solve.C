
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
#include "linbox/solutions/solve.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Solution of diagonal system
 *
 * Constructs a random nonsingular diagonal rational matrix D and a random right-hand
 * side integer b, and computes the solution to the Dx=b, checking the result
 * 
 * Return true on success and false on failure
 */

static bool testNonsingularRatIntSolve (size_t n, int iterations) 
{
	commentator.start ("Testing nonsingular solve with integer vector", "testNonsingularRatIntSolve", iterations);

	bool ret = true;
	int i;
	size_t j;

	GMPRationalField Q;
	SparseMatrix<GMPRationalField > A(Q,n,n);
	
	PID_integer Z;
	std::vector<PID_integer::Element> b(n);
	std::vector<GMPRationalField::Element> true_x(n),x(n);

	for (i=0; i < iterations; i++) {
		commentator.startIteration (i);

		for (j=0; j < n; ++j) {
			integer tmp_n, tmp_d, tmp_b;
			GMPRationalField::Element tmp;
			tmp_n = (integer) rand() % (2*(i + 1)) + 1;
			tmp_d = (integer) rand() % (2*(i + 1)) + 1;
			tmp_b = (integer) rand() % (2*(i + 1)) ;
			Q.init(tmp, tmp_n,tmp_d);
			A.setEntry(j,j,tmp);
			b[j]= tmp_b;
			if ( ( i%2) && (j % 2)) integer::negin(b[j]);
			Q.init(true_x[j] , b[j] * tmp_d, tmp_n);
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		solve (x, A, b);
		for (j=0; j < n; ++j) {
			if (!Q.areEqual(x[j] ,true_x[j])) {
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: System solution failed" << endl;
				ret = false;
			}
		}
		
		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNonsingularRatIntSolve");

	return ret;
}

/* Test 2: Solution of diagonal system
 *
 * Constructs a random nonsingular diagonal rational matrix D and a random right-hand
 * side rational b, and computes the solution to the Dx=b, checking the result
 * 
 * Return true on success and false on failure
 */

static bool testNonsingularRatRatSolve (size_t n, int iterations) 
{
	commentator.start ("Testing nonsingular solve with rational vector", "testNonsingularRatRatSolve", iterations);

	bool ret = true;
	int i;
	size_t j;

	GMPRationalField Q;
	SparseMatrix<GMPRationalField > A(Q,n,n);
	
	PID_integer Z;
	std::vector<GMPRationalField::Element> b(n);
	std::vector<GMPRationalField::Element> true_x(n),x(n);

	for (i=0; i < iterations; i++) {
		commentator.startIteration (i);

		for (j=0; j < n; ++j) {
			integer tmp_n, tmp_d, tmp_bn, tmp_bd;
			GMPRationalField::Element tmp,tmpb;
			tmp_n = (integer) rand() % (2*(i + 1)) + 1;
			tmp_d = (integer) rand() % (2*(i + 1)) + 1;
			tmp_bn = (integer) rand() % (2*(i + 1)) ;
			tmp_bd = (integer) rand() % (2*(i + 1)) + 1;
			//integer::nonzerorandom(tmp_n, 2*(i + 1) );
			//integer::nonzerorandom(tmp_d, 2*(i + 1) );
			//integer::random(tmp_bn, 2*(i + 1));
			//integer::nonzerorandom(tmp_bd, 2*(i +1) );
			if ( ( i%2) && (j % 2)) integer::negin(tmp_bn);
			Q.init(tmp, tmp_n,tmp_d);
			A.setEntry(j,j,tmp);
			Q.init(tmpb,tmp_bn,tmp_bd);
			b[j]= tmpb;
			
			Q.init(true_x[j] , tmp_bn * tmp_d, tmp_bd * tmp_n);
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		solve (x, A, b);
		
		for (j=0; j < n; ++j) {
			if (!Q.areEqual(x[j] ,true_x[j])) {
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: System solution failed" << endl;
				ret = false;
			}
		}
		
		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNonsingularRatRatSolve");

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

    if ( ! testNonsingularRatIntSolve(n,iterations) ) pass = false;
    if ( ! testNonsingularRatRatSolve(n,iterations) ) pass = false;    

	commentator.stop("solve test suite");
    //std::cout << (pass ? "passed" : "FAILED" ) << std::endl;

	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
