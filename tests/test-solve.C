/* tests/test-solve.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

/*! @file   tests/test-solve.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/direct-sum.h"
//#include "linbox/blackbox/diagonal.h"
#include "linbox/vector/stream.h"
#include "linbox/solutions/solve.h"

#include "linbox/solutions/minpoly.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Solve Ix=b, check that x == b
 *
 * Constructs a black box for the inverse of an n x n identity matrix and checks
 * that that inverse is itself the identity operator by applying it to random
 * vectors.
 *
 * F - Field over which to perform computations
 * stream - Vector stream for right-hand sides
 * text - Text to appear for commentator message
 * method - Method to use for solving the system
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector, class MethodTraits>
static bool testIdentitySolve (const Field          &F,
			       VectorStream<Vector> &stream,
			       const char           *text,
			       MethodTraits          method)
{
	typedef ScalarMatrix <Field> Blackbox;

	ostringstream str;
	str << "Testing identity solve (" << text << ")";

	commentator().start (str.str ().c_str (), "testIdentitySolve", stream.m ());

	bool ret = true;

	VectorDomain<Field> VD (F);

	typename Field::Element s;
	F.assign (s, F.one);
	Blackbox I (F, stream.n (), stream.n(), s);

	size_t n = stream.n();
	Vector v(F,n), w(F,n);

	// VectorWrapper::ensureDim (v, stream.n ());
	// VectorWrapper::ensureDim (w, stream.n ());

	while (stream) {
		commentator().startIteration ((unsigned)stream.j ());

		bool iter_passed = true;

		stream.next (v);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		try {
			solve (w, I, v, method);
			//solve (I, w, v, F, traits);
		}
		catch (const InconsistentSystem<Vector>& e) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Solve reported an inconsistent system" << endl;
			ret = iter_passed = false;
		}
		catch (const LinboxError&) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Solve failed to solve system" << endl;
			ret = iter_passed = false;
		}

		if (iter_passed) {
			report << "Output vector: ";
			VD.write (report, w);
			report << endl;

			if (!VD.areEqual (w, v))
				ret = iter_passed = false;

			if (!iter_passed)
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Vectors are not equal" << endl;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testIdentitySolve");

	return ret;
}

/* Test 2: Solution of diagonal system
 *
 * Constructs a random nonsingular diagonal matrix D and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 *
 * F - Field over which to perform computations
 * stream1 - Vector stream for diagonal entries
 * stream2 - Vector stream for right-hand sides
 * text - Text to appear for commentator message
 * method - Method to use for solving the system
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector, class MethodTraits>
static bool testNonsingularSolve (const Field          &F,
				  VectorStream<Vector> &stream1,
				  VectorStream<Vector> &stream2,
				  const char           *text,
				  MethodTraits          method)
{
	typedef Diagonal <Field> Blackbox;

	ostringstream str;
	str << "Testing nonsingular solve (" << text << ")";

	commentator().start (str.str ().c_str (), "testNonsingularSolve", stream1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;

	size_t n= stream1.n ();
	Vector d(F,n), b(F,n), x(F,n), y(F,n);

	// VectorWrapper::ensureDim (d, stream1.n ());
	// VectorWrapper::ensureDim (b, stream1.n ());
	// VectorWrapper::ensureDim (x, stream1.n ());
	// VectorWrapper::ensureDim (y, stream1.n ());

	while (stream1 && stream2) {
		commentator().startIteration ((unsigned)stream1.j ());

		ActivityState state = commentator().saveActivityState ();

		bool iter_passed = true;

		stream1.next (d);
		stream2.next (b);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		report << "Right-hand side:  ";
		VD.write (report, b);
		report << endl;
		BlasVector<Field> dd(F,d.begin(),d.size());
		//BlasVector<Field> dd(d,F);
		Blackbox D (dd);

		try {
			solve (x, D, b, method);
		}
		catch (const InconsistentSystem<Vector>& e) {
			commentator().restoreActivityState (state);
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: solve reported inconsistent system" << endl;
			ret = iter_passed = false;
		}
		catch (const LinboxError&) {
			commentator().restoreActivityState (state);
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: System solution failed" << endl;
			ret = iter_passed = false;
		}

		if (iter_passed) {
			report << "System solution:  ";
			VD.write (report, x);
			report << endl;

			D.apply (y, x);

			report << "Output:           ";
			VD.write (report, y);
			report << endl;

			if (!VD.areEqual (y, b))
				ret = iter_passed = false;

			if (!iter_passed)
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Computed solution is incorrect" << endl;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	stream1.reset ();
	stream2.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testNonsingularSolve");

	return ret;
}

/* Test 3: Solution of a consistent singular diagonal system with nonsingular
 * leading principal minor
 *
 * Constructs a random diagonal matrix D of rank r and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * stream1 - Vector stream for diagonal entries
 * stream2 - Vector stream for right-hand sides
 * text - Text to appear for commentator message
 * method - Method to use for solving the system
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector, class MethodTraits>
static bool testSingularConsistentSolve (const Field          &F,
					 unsigned int          n,
					 VectorStream<Vector> &stream1,
					 VectorStream<Vector> &stream2,
					 const char           *text,
					 MethodTraits          method)
{
	typedef Diagonal <Field> Blackbox;

	ostringstream str;
	str << "Testing singular consistent solve (" << text << ")";

	commentator().start (str.str ().c_str (), "testSingularConsistentSolve", stream1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;

	Vector d1(F), b1(F), d(F), b(F), x(F), y(F);

	VectorWrapper::ensureDim (d, n);
	VectorWrapper::ensureDim (b, n);
	VectorWrapper::ensureDim (x, n);
	VectorWrapper::ensureDim (y, n);
	VectorWrapper::ensureDim (d1, n);
	VectorWrapper::ensureDim (b1, n);

	MethodTraits traits (method);
	traits.preconditioner = Preconditioner::None;

	while (stream1 && stream2) {
		commentator().startIteration ((unsigned)stream1.j ());

		ActivityState state = commentator().saveActivityState ();


		stream1.next (d1);
		stream2.next (b1);

		VD.copy (d, d1);
		VD.copy (b, b1);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		report << "Right-hand side:  ";
		VD.write (report, b);
		report << endl;

		typename Blackbox::Vector_t dd(d,F);
		Blackbox D (dd);
		//Blackbox D (d);

		try {
			bool iter_passed = true;
			solve (x, D, b, traits);

			report << "System solution:  ";
			VD.write (report, x);
			report << endl;

			D.apply (y, x);

			report << "Output:           ";
			VD.write (report, y);
			report << endl;

			if (!VD.areEqual (y, b))
				ret = iter_passed = false;

			if (!iter_passed)
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Computed solution is incorrect" << endl;
		}
		catch (const InconsistentSystem<Vector>& e) {
			commentator().restoreActivityState (state);
			ostream &Report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
			Report << "ERROR: Inconsistent system exception" << endl;

			Report << "Certificate is: ";
			VD.write (Report, e.certificate()) << endl;

			ret = false;

			commentator().restoreActivityState (state);
		}
		catch (const LinboxError&) {
			commentator().restoreActivityState (state);
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: System solution failed" << endl;
			ret = false;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	stream1.reset ();
	stream2.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testSingularConsistentSolve");

	return ret;
}

/* Test 4: Solution of an inconsistent singular diagonal system with nonsingular
 * leading principal minor
 *
 * Constructs a random diagonal matrix D of rank r and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 *
 * F - Field over which to perform computations
 * stream1 - Vector stream for diagonal entries
 * stream2 - Vector stream for right-hand sides
 * text - Text to appear for commentator message
 * method - Method to use for solving the system
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector, class MethodTraits>
static bool testSingularInconsistentSolve (const Field          &F,
					   VectorStream<Vector> &stream1,
					   VectorStream<Vector> &stream2,
					   const char           *text,
					   MethodTraits          method)
{
	typedef Diagonal <Field> Blackbox;

	ostringstream str;
	str << "Testing singular inconsistent solve (" << text << ")";

	commentator().start (str.str ().c_str (), "testSingularInconsistentSolve", stream1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;

	Vector d1(F), d(F), b(F), x(F), y(F);
	typename Field::Element uTb;

	VectorWrapper::ensureDim (d, stream2.dim ());
	VectorWrapper::ensureDim (b, stream2.dim ());
	VectorWrapper::ensureDim (x, stream2.dim ());
	VectorWrapper::ensureDim (y, stream2.dim ());
	VectorWrapper::ensureDim (d1, stream1.dim ());

	MethodTraits traits (method);
	traits.certifyInconsistency = true;
	traits.checkResult = true; // @fixme Wiedemann does know it is inconsistent if it does not check...

	while (stream1 && stream2) {
		commentator().startIteration ((unsigned)stream1.j ());

		stream1.next (d1);
		stream2.next (b);

		VD.copy (d, d1);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		report << "Right-hand side:  ";
		VD.write (report, b);
		report << endl;

		typename Blackbox::Vector_t dd(d,F);
		Blackbox D (dd);
		//Blackbox D (d);

		bool inconsistent = false;
		bool errored = false;
		try {
			// @note Certificate of inconsistency will be put in x
			solve(x, D, b, traits);
		} catch (const LinBox::LinboxMathInconsistentSystem& error) {
			inconsistent = true;
		} catch (...) {
			errored = true;
		}

		if (inconsistent) {
			Vector u(x);
			D.applyTranspose (y, u);

			report << "Certificate of inconsistency found." << endl;

			report << "Certificate is: ";
			VD.write (report, u) << endl;

			report << "u^T A = ";
			VD.write (report, y) << endl;

			VD.dot (uTb, u, b);

			report << "u^T b = ";
			F.write (report, uTb) << endl;

			if (!VD.isZero (y)) {
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: u is not in the right nullspace of D" << endl;
				ret = false;
			}

			if (F.isZero (uTb)) {
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: u^T b = 0" << endl;
				ret = false;
			}
		}
		else if (errored) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "FAILED: Solver was unable to certify inconsistency" << endl;
		}
		else {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Solver gave solution even though system is inconsistent" << endl;
			ret = false;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	stream1.reset ();
	stream2.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testSingularInconsistentSolve");

	return ret;
}

/* Test 5: Solution of an inconsistent singular diagonal system with singular
 * leading principal minor
 *
 * Constructs a random diagonal matrix D of rank r and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 *
 * F - Field over which to perform computations
 * stream1 - Vector stream for diagonal entries
 * stream2 - Vector stream for right-hand sides
 * text - Text to appear for commentator message
 * preconditioner - Preconditioner to use
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector, class SparseVector>
static bool testSingularPreconditionedSolve (const Field                  &F,
					     VectorStream<SparseVector>   &stream1,
					     VectorStream<Vector>         &stream2,
					     const char                   *text,
					     Preconditioner preconditioner)
{
	ostringstream str;
	str << "Testing singular preconditioned solve (" << text << ")";

	commentator().start (str.str ().c_str (), "testSingularPreconditionedSolve", stream1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;

	SparseVector d1;
	typename Field::Element uTb;
	typename LinBox::Vector<Field>::Dense d;
	Vector b(F), x(F), y(F);

	VectorWrapper::ensureDim (d, stream2.dim ());
	VectorWrapper::ensureDim (b, stream2.dim ());
	VectorWrapper::ensureDim (x, stream2.dim ());
	VectorWrapper::ensureDim (y, stream2.dim ());



	Method::Wiedemann traits;
	traits.preconditioner = preconditioner;

	while (stream1 && stream2) {
		commentator().startIteration ((unsigned)stream1.j ());

		stream1.next (d1);
		stream2.next (b);

		VD.copy (d, d1);

		ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		report << "Right-hand side:  ";
		VD.write (report, b);
		report << endl;

		BlasVector<Field> dd(F,d.begin(),d.size());
		Diagonal<Field> A (dd);
		//Diagonal<Field> A (F, d);

		bool inconsistent = false;
		bool errored = false;
		try {
			// @note Certificate of inconsistency will be put in x
			solve(x, A, b, traits);
		} catch (const LinBox::LinboxMathInconsistentSystem& error) {
			inconsistent = true;
		} catch (...) {
			errored = true;
		}

		if (inconsistent) {
			Vector u(x);
			A.applyTranspose (y, u);

			report << "Certificate of inconsistency found." << endl;

			report << "Certificate is: ";
			VD.write (report, u) << endl;

			report << "u^T A = ";
			VD.write (report, y) << endl;

			VD.dot (uTb, u, b);

			report << "u^T b = ";
			F.write (report, uTb) << endl;

			if (!VD.isZero (y)) {
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: u is not in the right nullspace of D" << endl;
				ret = false;
			}

			if (F.isZero (uTb)) {
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: u^T b = 0" << endl;
				ret = false;
			}
		}
		else if (errored) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "FAILED: Solver was unable to certify inconsistency" << endl;
		}
		else {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Solver gave solution even though system is inconsistent" << endl;
			ret = false;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	stream1.reset ();
	stream2.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testSingularPreconditionedSolve");

	return ret;
}

/* Test 6: Test solution of random system
 */

template <class Field, class Vector1, class Vector2, class MethodTraits>
static bool testRandomSolve (const Field                  &F,
			     VectorStream<Vector1>        &A_stream,
			     VectorStream<Vector2>        &b_stream,
			     const char                   *text,
			     MethodTraits                  method)
{
	ostringstream str;
	str << "Testing random solve (" << text << ")";

	commentator().start (str.str ().c_str (), "testRandomSolve", b_stream.size ());

	bool ret = true;

	VectorDomain<Field> VD (F);
	MatrixDomain<Field> MD (F);

	Vector2 b(F), x(F), ATAx(F), ATb(F);

	VectorWrapper::ensureDim (b, b_stream.dim ());
	VectorWrapper::ensureDim (x, A_stream.dim ());
	VectorWrapper::ensureDim (ATAx, A_stream.dim ());
	VectorWrapper::ensureDim (ATb, A_stream.dim ());

	SparseMatrix<Field> A (F, A_stream);
	SparseMatrix<Field> AT (F,A.coldim (), A.rowdim ());
	BlasMatrix<Field> ATA (A.coldim (), A.coldim ());

	A.transpose (AT);

	MD.mul (ATA, AT, A);

	ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix A:" << endl;
	A.write (report);

	report << "Matrix A^T A:" << endl;
	MD.write (report, ATA);

	MethodTraits traits (method);

	while (b_stream) {
		commentator().startIteration ((unsigned)b_stream.pos ());

		bool iter_passed = true;

		b_stream >> b;

		ostream &Report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		Report << "Right-hand side b:";
		VD.write (Report, b) << endl;

		try {
			solve (A, x, b, F, traits);
		}
		catch (const InconsistentSystem<Vector2>& e) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Solve reported an inconsistent system" << endl;
			ret = iter_passed = false;
		}
		catch (const LinboxError&) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Solve failed to solve system" << endl;
			ret = iter_passed = false;
		}

		if (iter_passed) {
			Report << "Output vector x: ";
			VD.write (Report, x) << endl;

			MD.vectorMul (ATAx, ATA, x);
			MD.vectorMul (ATb, AT, b);

			if (!VD.areEqual (ATAx, ATb)) {
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: A^T Ax != A^T b" << endl;
				ret = iter_passed = false;
			}
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	A_stream.reset ();
	b_stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandomSolve");

	return ret;
}


template <class Field>
static bool testBasicMethodsSolve (const Field &F, size_t n)
{
	// tests of Auto, Blackbox, Elimination methods
	bool ret;
	commentator().start ("Testing Basic Methods Solve");
	ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	BlasVector<Field> xd(F,n), xh(F,n), xb(F,n), xe(F,n), xs(F,n), b(F,n);
	for(size_t i = 0; i < n/2; ++i) b[i] = F.one;
	//ScalarMatrix<Field> I(F, n/2, F.one), Z(F, n/2, F.zero);
	ScalarMatrix<Field> I(F, n, n, F.one), Z(F, 0, 0, F.zero);
	DirectSum<ScalarMatrix<Field>, ScalarMatrix<Field> > A(I, Z);

	VectorDomain<Field> VD(F);
	VD.write(report<<"b ", b) << endl;
	solve(xd, A, b);
	VD.write(report<<"xd ", xd) << endl;

	solve(xh, A, b, Method::Auto());
	VD.write(report<<"xh ", xh) << endl;

	solve(xb, A, b, Method::Blackbox());
	VD.write(report<<"xb ", xb) << endl;

	solve(xe, A, b, Method::DenseElimination());
	VD.write(report<<"xe ", xe) << endl;

	solve(xs, A, b, Method::SparseElimination());
	VD.write(report<<"xs ", xs) << endl;

	ret = VD.areEqual(xd, xh) && VD.areEqual(xd, xb) && VD.areEqual(xd, xe) && VD.areEqual(xd, xs);
	commentator().stop (MSG_STATUS (ret));
	return ret;
}


int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 100;
	static size_t m = 100;
	static size_t r = 20;
	static size_t N = 16;
	static integer q = 65521U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set column dimension of test matrices to N.", TYPE_INT,     &n },
		{ 'm', "-m M", "Set row dimension of test matrices to M.", TYPE_INT,     &m },
		{ 'r', "-r R", "Set singular system rank to R.", TYPE_INT,     &r },
		{ 'N', "-N N", "Set blocking factor to N.", TYPE_INT,     &N },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	//typedef Givaro::Modular<int64_t> Field;
    typedef Givaro::Modular<double> Field;

	parseArguments (argc, argv, args);
	Field F (q); Field::RandIter gen(F); Field::NonZeroRandIter nzgen(gen);

	commentator().start("Solve test suite", "solve");

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator().getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator().getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator().getMessageClass (PROGRESS_REPORT).setMaxDepth (5);
	//commentator().getMessageClass (BRIEF_REPORT).setMaxDepth (4);

	RandomDenseStream<Field,BlasVector<Field>, Field::NonZeroRandIter> stream1 (F, nzgen, n, (unsigned int)iterations), stream2 (F, nzgen, n, (unsigned int)iterations);
	RandomDenseStream<Field> stream3 (F, gen, r, (unsigned int)iterations), stream4 (F, gen, r, (unsigned int)iterations);
	RandomSparseStream<Field> stream6 (F, gen, (double) r / (double) n, n, (unsigned int)iterations);
	RandomSparseStream<Field> A_stream (F, gen, (double) r / (double) n, n, m);
#if 1
	Method::Wiedemann WM;
	if (!testIdentitySolve               (F, stream1, "Wiedemann", WM))
		pass = false;
	if (!testNonsingularSolve            (F, stream1, stream2, "Wiedemann", WM))
		pass = false;
#endif
#if 1
	Method::BlockWiedemann BWM;
        BWM.blockingFactor = N;
	if (!testIdentitySolve               (F, stream1, "BlockWiedemann", BWM))
		pass = false;
	if (!testNonsingularSolve            (F, stream1, stream2, "BlockWiedemann", BWM))
		pass = false;
#endif

#if 1
	if (!testSingularConsistentSolve     (F, n, stream3, stream4,
					      "Wiedemann", Method::Wiedemann ()))
		pass = false;
	if (!testSingularInconsistentSolve   (F, stream3, stream2,
					      "Wiedemann", Method::Wiedemann ()))
		pass = false;
#endif

#if 1
	if (!testSingularPreconditionedSolve (F, stream6, stream2,
					      "Sparse preconditioner", Preconditioner::Butterfly))
		pass = false;
#endif

#if 0
	if (!testIdentitySolve               (F, stream1,
					      "Lanczos", Method::Lanczos ()))
		pass = false;
	if (!testNonsingularSolve            (F, stream1, stream2,
					      "Lanczos", Method::Lanczos ()))
		pass = false;
	if (!testSingularConsistentSolve     (F, n, stream3, stream4,
					      "Lanczos", Method::Lanczos ()))
		pass = false;

	Method::Lanczos traits1;
	traits1.preconditioner = Preconditioner::FullDiagonal;

	if (!testRandomSolve (F, A_stream, stream1, "Lanczos", traits1))
		pass = false;
#endif

#if 0
	Method::BlockLanczos traits2;
	traits2.preconditioner = Preconditioner::FullDiagonal;
	traits2.blockingFactor = N;

	if (!testRandomSolve (F, A_stream, stream1, "Block Lanczos", traits2))
		pass = false;
#endif
#if 1

    if ( ! testBasicMethodsSolve (F, n) )
		pass = false;
#endif
	commentator().stop("solve test suite");
    //std::cout << (pass ? "passed" : "FAILED" ) << std::endl;

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
