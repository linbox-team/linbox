
/* tests/test-block-wiedemann.C
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


#include "linbox-config.h"

#include <iostream>


#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#ifdef __LINBOX_HAVE_OCL
  #include "linbox/algorithms/opencl-domain.h"
#else
  #include "linbox/algorithms/blas-domain.h"
#endif
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/block-wiedemann.h"
#include "linbox/algorithms/coppersmith.h"
#include "linbox/blackbox/sparse.h"
//#include "linbox/blackbox/diagonal.h"
//#include "linbox/blackbox/scalar-matrix.h"
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
	static size_t q = 2147483647U;
	static size_t blocking = 0;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to N.", TYPE_INT,     &n },
		{ 'N', "-N N", "Set blocking factor to N.", TYPE_INT,     &N },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &q },
		{ 'b', "-b N", "Set the blocking size", TYPE_INT, &blocking },
		END_OF_ARGUMENTS
	};

	//typedef Modular<double> Field;
	typedef Modular<uint32_t> Field;
	typedef MatrixDomain<Field> MyDomain;
	typedef BlasVector<Field> Vector;

	parseArguments (argc, argv, args);
	Field F ( (uint32_t) q);
	MyDomain MD(F);
	VectorDomain<Field> VD (F);

	commentator().start("block wiedemann test suite", "block-wiedemann");
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	RandomDenseStream<Field> stream1 (F, n, 3);// stream2 (F, n, 1);
	Vector d(F,n), b(F,n), x(F,n), y(F,n);
	stream1.next (d);
	stream1.next (b);
	d[n-1] = F.zero;
	b[n-1] = F.zero;

	VD.write (report << "Right-hand side: b =  ", b) << endl;

// sparse
	typedef SparseMatrix <Field> Blackbox;
	Blackbox D (F, n, n);
	for (size_t i = 0; i < n; ++i) D.setEntry(i, i, d[i]);
	stream1.next (d);
	//for (size_t i = 0; i < n-1; ++i) D.setEntry(i, i+1, d[i]);
	//D.setEntry(0, n-1, d[n-1]);
	//... setup entries

// diagonal
	//typedef Diagonal <Field> Blackbox;
	//Blackbox D (F, d);
	//VD.write (report << "Diagonal entries: ", d) << endl;

// scalar
	//typedef ScalarMatrix<Field> Blackbox;
	//Blackbox D (F, n, F.one);
	//report << "Scalar matrix: D = " << F.one << endl;

#if 1
	// Yuhasz' Matrix Berlekamp Massey being used
	CoppersmithSolver<MyDomain> RCS(MD,blocking);
	RCS.solveNonSingular(x, D, b);

	VD.write (report << "Matrix Berlekamp Massey solution:  ", x) << endl;
	D.apply (y, x);
	VD.write ( report << "Output:           ", y) << endl;

	if (!VD.areEqual (y, b)) {
		pass = false;
		report << "ERROR: Right generator based solution is incorrect" << endl;
	}
#endif

#if 0
	// Giorgi's block method, SigmaBasis based, being used

#ifdef __LINBOX_HAVE_OCL
// using this as a test of the opencl matrix domain
	typedef OpenCLMatrixDomain<Field> Context;
// but note, shouldn't need the ifdef.  OpenCLMatrixDomain defaults to BlasMatrixDomain calls.
#else
	typedef BlasMatrixDomain<Field> Context;
#endif
	Context MD(F);
	Vector z(F,n);
	BlockWiedemannSolver<Context> LBWS(MD);
	report << "calling solver" << endl;
	LBWS.solveNonSingular(z, D, b);

	VD.write (report << "SigmaBasis solution: D^{_1}b = ", z) << endl;
	D.apply (y, z);
	VD.write ( report << "Output: D D^{-1}b = ", y) << endl;

	if (!VD.areEqual (y, b)) {
		pass = false;
		report << "ERROR: Left generator based solution is incorrect" << endl;
	}

#endif

	commentator().stop("block wiedemann test suite");
    //std::cout << (pass ? "passed" : "FAILED" ) << std::endl;

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
