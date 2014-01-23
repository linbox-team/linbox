
/* tests/test-la-block-lanczos.C
 * Copyright (C) 2004 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * --------------------------------------------------------
 *
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


/*! @file  tests/test-la-block-lanczos.C
 * @ingroup tests
 * @test NO DOC
 * @brief  no doc
 */



#include "linbox-config.h"

#include <iostream>
#include <fstream>


#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/matrix/sparse.h"
#include "linbox/vector/stream.h"
#include "linbox/algorithms/la-block-lanczos.h"

#include "test-common.h"

using namespace LinBox;
using namespace std;

/* Test 1: Test solution of random system
 */

template <class Field, class Vector1, class Vector2>
static bool testRandomSolve (const Field           &F,
			     VectorStream<Vector1> &A_stream,
			     VectorStream<Vector2> &y_stream,
			     size_t                 N)
{
	typedef LABlockLanczosSolver<Field, BlasMatrix<Field> > LABLSolver;

	commentator().start ("Testing random solve (Block Lanczos)", "testRandomSolve", y_stream.size ());

	std::ostream &report1 = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	Vector2 y, b, x1, x2;

	VectorWrapper::ensureDim (b, y_stream.dim ());
	VectorWrapper::ensureDim (y, y_stream.dim ());
	VectorWrapper::ensureDim (x1, A_stream.dim ());
	VectorWrapper::ensureDim (x2, A_stream.dim ());

	SparseMatrix2<Field> A (F, A_stream);

	report1 << "n = " << y_stream.dim () << endl;
	report1 << "N = " << N << endl;

	report << "Input matrix A:" << endl;
	A.write (report);

	typename Field::RandIter ri (F, 0, time (NULL));

	BlockLanczosTraits traits;
	traits.preconditioner (BlockLanczosTraits::NO_PRECONDITIONER);
	traits.blockingFactor (N);
	//traits.maxTries (1);

	LABLSolver lablsolver (F, traits, ri);

	while (y_stream) {
		commentator().startIteration ( (unsigned int) y_stream.pos ());

		y_stream >> y;
		A.apply (b, y);

		std::ostream &raport = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		raport << "Right-hand side b:";
		VD.write (raport, b) << endl;

		if (!lablsolver.solve (A, x2, b)) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Solve failed to solve system" << endl;
			ret = false;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	A_stream.reset ();
	y_stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandomSolve");

	return ret;
}

/* Test 2: Test sampling of nullspace of random system
 */

template <class Field, class Vector1>
static bool testSampleNullspace (const Field           &F,
				 VectorStream<Vector1> &A_stream,
				 size_t                 N,
				 unsigned int           num_iter)
{
	typedef BlasMatrix<Field> Matrix;
	typedef LABlockLanczosSolver<Field, Matrix> LABLSolver;

	commentator().start ("Testing sampling from nullspace (Block Lanczos)", "testSampleNullspace", num_iter);

	std::ostream &report1 = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	Matrix x (F,A_stream.dim (), N);

	SparseMatrix2<Field> A (F, A_stream);

	report1 << "n = " << A_stream.dim () << endl;
	report1 << "N = " << N << endl;

	report << "Input matrix A:" << endl;
	A.write (report);

	typename Field::RandIter ri (F, 0, time (NULL));

	BlockLanczosTraits traits;
	traits.preconditioner (BlockLanczosTraits::NO_PRECONDITIONER);
	traits.blockingFactor (N);
	//traits.maxTries (1);

	LABLSolver lablsolver (F, traits, ri);

	for (unsigned int i = 0; i < num_iter; ++i) {
		commentator().startIteration (i);

		unsigned int number;
		number = lablsolver.sampleNullspace (A, x);

		commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Number of nullspace vectors found: " << number << std::endl;

		commentator().stop ("done");
		commentator().progress ();
	}

	A_stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testSampleNullspace");

	return ret;
}

/* Test 3: Test rank
 */

template <class Field, class Vector1>
static bool testRank (const Field           &F,
		      VectorStream<Vector1> &A_stream,
		      size_t                 N,
		      unsigned int           num_iter)
{
	typedef BlasMatrix<Field> Matrix;
	typedef LABlockLanczosSolver<Field, Matrix> LABLSolver;

	commentator().start ("Testing rank (Block Lanczos)", "testRank", num_iter);

	std::ostream &report1 = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;

	MatrixDomain<Field> MD (F);

	SparseMatrix2<Field> A (F, A_stream);

	report1 << "n = " << A_stream.dim () << endl;
	report1 << "N = " << N << endl;

	report << "Input matrix A:" << endl;
	A.write (report);

	typename Field::RandIter ri (F, 0, time (NULL));

	BlockLanczosTraits traits;
	traits.preconditioner (BlockLanczosTraits::NO_PRECONDITIONER);
	traits.blockingFactor (N);
	//traits.maxTries (1);

	LABLSolver lablsolver (F, traits, ri);

	for (unsigned int i = 0; i < num_iter; ++i) {
		commentator().startIteration (i);

		unsigned int rank;
		rank = lablsolver.rank (A);

		commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Rank: " << rank << std::endl;

		commentator().stop ("done");
		commentator().progress ();
	}

	A_stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "rank");

	return ret;
}

int main (int argc, char **argv)
{
	static unsigned int i = 5;
	static unsigned int n = 10; // because it shows the problem
	static int k = 5;
	static unsigned int q = 2;
	static unsigned int N = 16;

	static Argument args[] = {
		{ 'i', "-i I", "Number of iterations.", TYPE_INT, &i },
		{ 'n', "-n N", "Dimension of test matrix.", TYPE_INT, &n },
		{ 'k', "-k K", "K nonzero entries per row in test matrix.", TYPE_INT, &k },
		{ 'N', "-N N", "Blocking factor.", TYPE_INT, &N },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &q },
		END_OF_ARGUMENTS
	};

	bool fail = false ;

	typedef Modular<uint8_t> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	std::cout << "Lookahead-based block Lanczos test suite" << std::endl << std::endl;

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator().getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator().getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator().getMessageClass (PROGRESS_REPORT).setMaxDepth (0);

	RandomSparseStream<Field> A_stream (F, (double) k / (double) n, n, n);
	RandomDenseStream<Field> y_stream (F, n, i);

	fail |= testRandomSolve (F, A_stream, y_stream, N);
	fail |= testSampleNullspace (F, A_stream, N, i);
	fail |= testRank (F, A_stream, N, i);

	return fail;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

