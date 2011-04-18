
/* tests/test-mg-block-lanczos.C
 * Copyright (C) 2004 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/vector/stream.h"
#include "linbox/algorithms/mg-block-lanczos.h"

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
	typedef MGBlockLanczosSolver<Field, DenseMatrixBase<typename Field::Element> > MGBLSolver;

	commentator.start ("Testing random solve (Block Lanczos)", "testRandomSolve", y_stream.size ());

	std::ostream &report1 = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	Vector2 y, b, x1, x2;

	VectorWrapper::ensureDim (b, y_stream.dim ());
	VectorWrapper::ensureDim (y, y_stream.dim ());
	VectorWrapper::ensureDim (x1, A_stream.dim ());
	VectorWrapper::ensureDim (x2, A_stream.dim ());

	SparseMatrix<Field> A (F, A_stream);

	report1 << "n = " << y_stream.dim () << endl;
	report1 << "N = " << N << endl;

	report << "Input matrix A:" << endl;
	A.write (report);

	typename Field::RandIter ri (F, 0, time (NULL));

	BlockLanczosTraits traits;
	traits.preconditioner (BlockLanczosTraits::SYMMETRIZE);
	traits.blockingFactor (N);
	//traits.maxTries (1);

	MGBLSolver mgblsolver (F, traits, ri);

	while (y_stream) {
		commentator.startIteration (y_stream.pos ());

		y_stream >> y;
		A.apply (b, y);

		std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Right-hand side b:";
		VD.write (report, b) << endl;

		if (!mgblsolver.solve (A, x2, b)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Solve failed to solve system" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	A_stream.reset ();
	y_stream.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomSolve");

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
	typedef DenseMatrixBase<typename Field::Element> Matrix;
	typedef MGBlockLanczosSolver<Field, Matrix> MGBLSolver;

	commentator.start ("Testing sampling from nullspace (Block Lanczos)", "testSampleNullspace", num_iter);

	std::ostream &report1 = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;
	unsigned int number;

	MatrixDomain<Field> MD (F);

	Matrix x (A_stream.dim (), N);

	SparseMatrix<Field> A (F, A_stream);

	report1 << "n = " << A_stream.dim () << endl;
	report1 << "N = " << N << endl;

	report << "Input matrix A:" << endl;
	A.write (report);

	typename Field::RandIter ri (F, 0, time (NULL));

	BlockLanczosTraits traits;
	traits.preconditioner (BlockLanczosTraits::SYMMETRIZE);
	traits.blockingFactor (N);
	//traits.maxTries (1);

	MGBLSolver mgblsolver (F, traits, ri);

	for (unsigned int i = 0; i < num_iter; ++i) {
		commentator.startIteration (i);

		number = mgblsolver.sampleNullspace (A, x);

		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Number of nullspace vectors found: " << number << std::endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	A_stream.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSampleNullspace");

	return ret;
}

int main (int argc, char **argv)
{
	static int i = 5; 
	static int n = 100;
	static int k = 5;
	static int q = 2;
	static int N = 16;

	bool pass = true;

	static Argument args[] = {
		{ 'i', "-i I", "Number of iterations.", TYPE_INT, &i },
		{ 'n', "-n N", "Dimension of test matrix.", TYPE_INT, &n },
		{ 'k', "-k K", "K nonzero entries per row in test matrix.", TYPE_INT, &k },
		{ 'N', "-N N", "Blocking factor.", TYPE_INT, &N },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &q },
		{ '\0', NULL, NULL, TYPE_NONE, NULL },
		{ '\0' }
	};

	typedef Modular<uint8> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	commentator.start("Montgomery block Lanczos test suite");

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (0);

	RandomSparseStream<Field> A_stream (F, (double) k / (double) n, n, n);
	RandomDenseStream<Field> y_stream (F, n, i);

	if (!testRandomSolve (F, A_stream, y_stream, N)) pass=false;;
	commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) 
		<< "	Skipping Sample Nullspace test (which has mem problems)" << std::endl;
	//if (!testSampleNullspace (F, A_stream, N, i)) pass=false;;

	commentator.stop("Montgomery block Lanczos test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
