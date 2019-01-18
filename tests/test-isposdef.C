/* tests/test-isposdef.C
 * Copyright (C) LinBox
 *
 * -----------------------------------------------------
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 * Function definitions for block Lanczos iteration
 */

/*! @file  tests/test-isposdef.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */

#include "linbox/linbox-config.h"

#include <iostream>

#include "linbox/util/commentator.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/solutions/is-positive-definite.h"

using namespace LinBox;

/* Positive definiteness check on two specific matrices
 *
 * Checks that I is found to be pos def and that a slight modification of I (to be singular and have a negative eigenvallue) is found to be not pos def.
 *
 * Sparsity is not (yet) used.
 */

template <class Ring, class Method>
bool testIsPosDef(const Ring &Z, size_t n, unsigned int iterations, Method &M, std::string methodName, double sparsity = 0.05)
{
	typedef SparseMatrix<Ring> Blackbox;

	commentator().start ("Testing isPositiveDefinite", "testIsPosDef", iterations);

	bool ret = true;
	unsigned int i;

	// iterations since there is randomness in the alg (though none in this test).
	for (i = 0; i < iterations; ++i) {
		commentator().startIteration (i);

		Blackbox A (Z, n, n);
		typename Ring::Element e; Z.assign(e, Z.one);
		for (size_t j = 0; j < n; ++j)
			A.setEntry(j, j, e);

		std::ostream & report =
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		Z.write( report ) << std::endl;
		A.write( report ) << std::endl;
		bool p;
		p = isPositiveDefinite(A,M);
		report << "Positivedefiniteness on I computed by " << methodName << ", " << p << std::endl;
		if (!p) {report << "ERROR: should be pos def" << std::endl; ret = false;}

		Z.negin(e);
		Z. assign(A. refEntry(n/2, n/2), e);
		A.setEntry(1, 2, e);
		A.setEntry(2, 1, e);
		p = isPositiveDefinite(A,M);
		report << "Matrix:\n";
		A.write( report ) << std::endl;
		report << "Positivedefiniteness on indefinite example computed by " << methodName << ", " << p << std::endl;
		if (p) {report << "ERROR: should not be pos def" << std::endl; ret = false;}

		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testEliminationRank");

	return ret;
}



int main (int argc, char **argv)
{

//     commentator().setMaxDetailLevel( 100000 );
//     commentator().setMaxDepth( 100000 );

	bool pass = true;

	static size_t n = 80;
	static unsigned int iterations = 2;
        static double sparsity = 0.05;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 's', "-s S", "Sparse matrices with density S.", TYPE_DOUBLE,     &sparsity },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	srand ((unsigned)time (NULL));

	commentator().start("IsPositiveDefinite solution test suite", "IsPositiveDefinite");
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	Givaro::ZRing<Integer> R;

	Method::Auto MH;
	pass = pass and testIsPosDef(R, n, iterations, MH, "Method::Auto", sparsity);
	Method::Elimination ME;
	pass = pass and testIsPosDef(R, n, iterations, MH, "Method::Elimination", sparsity);
	Method::Blackbox MB;
	pass = pass and testIsPosDef(R, n, iterations, MB, "Method::Blackbox", sparsity);

	commentator().stop("IsPositiveDefinite solution test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
