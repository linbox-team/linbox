/* tests/test-isposdef.C
 * Copyright (C) LinBox
 *
 * -----------------------------------------------------
 *
 * This file is part of LinBox, licensed under the GNU Lesser General
 * Public License. See COPYING for more information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/dense.h"
#include "linbox/solutions/is-positive-definite.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: positive definiteness of a random sparse matrix
 *
 * Constructs a random sparse matrix and computes its rank using Gaussian
 * elimination (direct and blas) and Wiedemann's algorithm. Checks that the results match.
 */

template <class Ring>
bool testIsPosDef(const Ring &Z, size_t n, unsigned int iterations, double sparsity = 0.05) 
{
	typedef SparseMatrix<Ring> Blackbox;

	commentator.start ("Testing isPositiveDefinite", "testIsPosDef", iterations);

	bool ret = true;
	unsigned int i;

	typename Ring::RandIter ri (Z);

	for (i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		Blackbox A (Z, n, n);
		typename Ring::Element e; Z.init(e, 1);
		for (size_t j = 0; j < n; ++j)
			A.setEntry(j, j, e);

		std::ostream & report = 
		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		Z.write( report ) << endl; 
		A.write( report ) << endl; 
		bool p;
		p = isPositiveDefinite(A); 
		report << "Positivedefiniteness on I computed by default (Hybrid) method: " << p << endl;
		if (!p) {report << "ERROR: should be pos def" << endl; ret = false;}

		Z.negin(e);
		Z. assign(A. refEntry(n/2, n/2), e);
		A.setEntry(1, 2, e);
		A.setEntry(2, 1, e);
		p = isPositiveDefinite(A); 
		report << "Matrix:\n";
		A.write( report ) << endl; 
		report << "Positivedefiniteness on indefinite example computed by default (Hybrid) method: " << p << endl;
		if (p) {report << "ERROR: should not be pos def" << endl; ret = false;}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testEliminationRank");

	return ret;
}



int main (int argc, char **argv)
{

//     commentator.setMaxDetailLevel( 100000 );
//     commentator.setMaxDepth( 100000 );
   
	bool pass = true;

	static size_t n = 80;
	static integer q = 65519U;
	//static integer q = 1000003U;
	static int iterations = 2;
        static double sparsity = 0.05;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
        { 's', "-s S", "Sparse matrices with density S.", TYPE_DOUBLE,     &sparsity },
		{ '\0' }

	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	commentator.start("IsPositiveDefinite solution test suite", "IsPositiveDefinite");
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

    PID_integer R;
        
	if (!testIsPosDef(R, n, iterations, sparsity)) pass = false;

	commentator.stop("IsPositiveDefinite solution test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
