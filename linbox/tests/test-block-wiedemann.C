/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* tests/test-solve.C
 * evolved by -bds from test-solve.C
 *
 * --------------------------------------------------------
 *
 * Copyright (c) LinBox
 *
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
 */

/*! @file   tests/test-solve.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
//#include "linbox/algorithms/right-block-wiedemann.h"
#include "linbox/algorithms/block-wiedemann.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/vector/stream.h"

#include "test-common.h"

using namespace LinBox;

/* Solution of diagonal system
 *
 * Constructs a random nonsingular diagonal matrix D and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 *
 * F - Field over which to perform computations
 * stream1 - Vector stream for diagonal entries
 * stream2 - Vector stream for right-hand sides
 *
 * Return true on success and false on failure
 */

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 9;
	static size_t N = 16;
	static integer q = 2147483647U;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to N.", TYPE_INT,     &n },
		{ 'N', "-N N", "Set blocking factor to N.", TYPE_INT,     &N },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	typedef Modular<uint32_t> Field;
	typedef vector<Field::Element> Vector;
	typedef Diagonal <Field> Blackbox;

	parseArguments (argc, argv, args);
	Field F (q);
	VectorDomain<Field> VD (F);

	commentator.start("block wiedemann test suite", "block-wiedemann");
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	RandomDenseStream<Field> stream1 (F, n, 1), stream2 (F, n, 1);
	Vector d(n), b(n), x(n), y(n);
	stream1.next (d);
	stream2.next (b);

	VD.write (report << "Diagonal entries: ", d) << endl;
	VD.write (report << "Right-hand side:  ", b) << endl;

	Blackbox D (F, d);

#if 0
	// Yuhasz Matrix Berlekamp Massey being used
	RightBlockWiedemannSolver<Field> RBWS(F);
	RBWS.solveNonSingular(x, D, b);

	VD.write (report << "Matrix Berlekamp Massey solution:  ", x) << endl;
	D.apply (y, x);
	VD.write ( report << "Output:           ", y) << endl;

	if (!VD.areEqual (y, b)) {
		pass = false;
		report << "ERROR: Right generator based solution is incorrect" << endl;
	}
#endif

#if 1
	// Giorgi SigmaBasis Berlekamp Massey being used
	Vector z(n);
	BlockWiedemannSolver<Field> LBWS(F);
	LBWS.solveNonSingular(z, D, b);

	VD.write (report << "SigmaBasis solution:  ", z) << endl;
	D.apply (y, z);
	VD.write ( report << "Output:           ", y) << endl;

	if (!VD.areEqual (y, b)) {
		pass = false;
		report << "ERROR: Left generator based solution is incorrect" << endl;
	}

#endif

	commentator.stop("block wiedemann test suite");
    //std::cout << (pass ? "passed" : "FAILED" ) << std::endl;

	return pass ? 0 : -1;
}
