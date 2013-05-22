/* linbox/tests/test-blackbox.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Extracted by bds from test-generic.h written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2003-03-11 Austin Lobo <alobo2@washcoll.edu>
 *
 * Added testApply and testApplyTranspose to time critical
 *        blackbox-functions.
 * ------------------------------------
 * fix testLinearity for rectangular matrices.  2006-02-17 Hui Wang
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
 *.
 */

#ifndef __LINBOX_test_blackbox_H
#define __LINBOX_test_blackbox_H

#include <iostream>
#include <vector>
//#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/util/field-axpy.h"
#include "linbox/vector/stream.h"
//#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/matrix-domain.h"

#include "test-common.h"

using namespace std;
/*
Generic tests for black boxes

For field F, BB A, vector streams s, t:
testTranspose (F, A, s, t)
testLinearity (A, s, t)
testReadWrite(A)
testBlackbox(A, b) // calls the first 2 and, if b, calls testReadWrite.

testBB(F) has been deleted. It assumed a BB could be built from a single size param (this is never true?!).
*/

/** Generic Blackbox test 1: Random check that (u^T A) v = u^T (A v).
 *
 * Take the given black box and compute u^T A v via <A^T u, v> and <u, Av> for
 * randomly chosen u and v. Check whether the results are equal. In theory, this
 * should guarantee that tranpose is working correctly if apply and dot product
 * are also working correctly. Apply and dot product should, of course, be
 * independently checked.
 *
 * F - Field over which to perform computations
 * A - Black box of which to construct the transpose
 * iterations - Number of random vectors to which to apply matrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Blackbox, class Vector>
static bool
testTranspose (const Field                      &F,
	       Blackbox				&A,
	       LinBox::VectorStream<Vector>     &stream1,
	       LinBox::VectorStream<Vector>     &stream2)
{
	bool ret = true;

	Vector u(F), v(F), uA(F), Av(F);

	LinBox::VectorWrapper::ensureDim (u, A.rowdim ());
	LinBox::VectorWrapper::ensureDim (uA, A.coldim ());
	LinBox::VectorWrapper::ensureDim (v, A.coldim ());
	LinBox::VectorWrapper::ensureDim (Av, A.rowdim ());

	LinBox::VectorDomain <Field> VD (A.field());
	//LinBox::VectorDomain <Field> VD (F);
	typename Field::Element r1, r2;
	ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Blackbox transpose test [that u^T(Av) == (uA)^T v]" << std::endl;

	while (stream1 && stream2) {
		LinBox::commentator().startIteration ((unsigned int) stream1.j ());

		stream1.next (u);
		stream2.next (v);

		//ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		VD.write( report << "Input vector u:            ", u) << endl;
		VD.write( report << "Input vector v:            ", v) << endl;

		A.apply (Av, v);

		VD.write( report << "Result of apply:           ", Av) << endl;

		VD.dot (r1, u, Av);
		A.applyTranspose (uA, u);

		VD.write( report << "Result of transpose apply: ", uA) << endl;

		VD.dot (r2, uA, v);

		F.write( report << "<u, Av>:  ", r1) << endl;

		F.write( report << "<A^T u, v>:  ", r2) << endl;

		if (!F.areEqual (r1, r2)) {
			ret = false;
			LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Values are not equal" << endl;
		}

		LinBox::commentator().stop ("done");
		LinBox::commentator().progress ();
	}

	return ret;
}

/** Generic Blackbox test 2: Linearity of black boxes.
 *
 * Given an arbitrary black box A, compute A(x+alpha y) and Ax+alphaAy and check equality.
 *
 * F - Field over which to perform computations
 * A - Black box of which to compute the dense representation
 * stream1 - Stream for x's
 * stream2 - Stream for y's
 *
 * Return true on success and false on failure
 */

template <class BB, class Vector>
static bool
testLinearity (//const Field                             &F,
	       BB 				 &A,
	       LinBox::VectorStream<Vector>      &stream1,
	       LinBox::VectorStream<Vector>      &stream2)
{
	bool ret = true, iter_passed;

	size_t n = A.rowdim ();
	size_t m = A.coldim ();

	typedef typename BB::Field Field;
	Field F( A.field() );
	Vector x(F), y(F), xpay(F), Axpay(F), Ax(F), Ay(F), AxpaAy(F);
	LinBox::VectorDomain <Field> VD (A.field());
	typename Field::RandIter r (A.field());
	typename Field::Element alpha;

	LinBox::VectorWrapper::ensureDim (x, m);
	LinBox::VectorWrapper::ensureDim (y, m);
	LinBox::VectorWrapper::ensureDim (xpay, m);
	LinBox::VectorWrapper::ensureDim (Axpay, n);
	LinBox::VectorWrapper::ensureDim (Ax, n);
	LinBox::VectorWrapper::ensureDim (Ay, n);
	LinBox::VectorWrapper::ensureDim (AxpaAy, n);

	ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Blackbox linearity test [that A.apply to (ax + y) == a A.apply to x + A.apply to y]" << std::endl;

	while (stream1 && stream2) {
		LinBox::commentator().startIteration ((unsigned int) stream1.j ());

		iter_passed = true;

		stream1.next (x);
		stream2.next (y);

		r.random (alpha);

		VD.write( report << "Input vector x: ", x) << endl;

		VD.write( report << "Input vector y: ", y) << endl;

		A.field().write( report << "Input alpha: ", alpha) << endl;

		VD.axpy ( xpay, alpha, y, x);
		A.apply ( Axpay, xpay);

		A.apply ( Ax, x);
		A.apply ( Ay, y);
		VD.axpy ( AxpaAy, alpha, Ay, Ax);

		VD.write( report << "   x+alpha y = ", xpay) << endl;

		VD.write( report << "A(x+alpha y) = ", Axpay) << endl;

		VD.write( report << "          Ax = ", Ax) << endl;

		VD.write( report << "          Ay = ", Ay) << endl;

		VD.write( report << " Ax+alpha Ay = ", AxpaAy) << endl;

		if (!VD.areEqual (Axpay, AxpaAy))
			ret = iter_passed = false;

		if (!iter_passed)
			LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Vectors are not equal" << endl;

	//	LinBox::commentator().stop ("done");
	//	LinBox::commentator().progress ();
	}

	return ret;
}

/** Generic Blackbox test 3: black box read/write.
 *
 * write the black box A, read it back and check equality.
 *
 * F - Field over which to perform computations
 * A - Black box
 *
 * Return true on success and false on failure
 */

template <class BB>
static bool
testReadWrite(BB &A)
{ //perhaps read/write to a stringstream?
	typedef typename BB::Field Field;
	bool pass = true;
	ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Blackbox Read/Write test: write then read back" << std::endl;

	ofstream out("temp");
	if (not out) {
		pass = false;
		report << "failure to open file for writing" << std::endl;
	}
	A.write(out) << std::endl;
	BB B(A.field());
	ifstream in("temp");
	if (not in) {
		pass = false;
		report << "failure to open file for reading" << std::endl;
	}
	B.read(in);
	std::vector<typename Field::Element> x(A.coldim()), y(A.rowdim()), z(B.rowdim());
	LinBox::MatrixDomain<Field> MD(A.field());
	if (not MD.areEqual(A, B)) {
		pass = false;
		report << "failure to get same matrix back from write/read" << std::endl;
		B.write(report << "B is ") << std::endl;
	}
	if (pass) report << "PASS: successful write/read" << std::endl;
	return pass;
}
/** Generic blackbox test 3: combination of tests
 *
 * If large, time apply and applyTranspose.
 * Call testTranspose and testLinearity.
 */
template <class BB>
static bool
testBlackboxNoRW(BB &A)
{
	size_t largeThresh = 2000; // Above it do timing of apply and applyTr.
	typedef typename BB::Field Field;
	typedef LinBox::BlasVector<Field> DenseVector;
	std::ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "testBlackbox on " << A.rowdim() << " by " << A.coldim() << " matrix." << endl;

	LinBox::commentator().setMaxDepth(-1);
	bool ret = true;
	const Field& F = A.field();

	/* timing tests */
	{
		DenseVector x(F,A.coldim()), y(F,A.rowdim());
		for(size_t i = 0; i < A.coldim(); ++i) F.init(x[i], i);
		for(size_t i = 0; i < A.rowdim(); ++i) F.init(y[i], i);
		//A.apply(y, x);

		if (A.coldim() >= largeThresh)
		{
			LinBox::commentator().start ("\t--Timing Test (Av)","testApply", 1);
			A.apply(y, x);
			LinBox::commentator().stop (MSG_STATUS (true), (const char *) 0, "testApply");
		}

		if (A.rowdim() >= largeThresh)
		{
			LinBox::commentator().start ("\t--Timing Test(v^T A)",
						   "testApplyTranspose", 1);
			A.applyTranspose(x, y);
			LinBox::commentator().stop (MSG_STATUS (true), (const char *) 0, "testApplyTranspose");
		}

	} // timing test block

	size_t iterations = 1;
	typename Field::RandIter r(F);
	LinBox::RandomDenseStream<Field, DenseVector> stream1 (F, r, A.rowdim(), iterations);
	//typename Field::Element x;
	//r.random(x);
	LinBox::RandomDenseStream<Field, DenseVector> stream2 (F, r, A.coldim(), iterations);
	ret = ret && testLinearity (A, stream1, stream2);

	LinBox::RandomDenseStream<Field, DenseVector> stream3 (F, r, A.rowdim(), iterations);
	LinBox::RandomDenseStream<Field, DenseVector> stream4 (F, r, A.coldim(), iterations);
	ret = ret && testTranspose (F, A, stream3, stream4);

	return ret;
}

template <class BB>
static bool
testBlackbox(BB &A)
{
	return testBlackboxNoRW(A) and testReadWrite(A);
}
#endif // __LINBOX_test_blackbox_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
