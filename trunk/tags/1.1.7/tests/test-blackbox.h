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
 * See COPYING for license information.
 */

#ifndef __LINBOX_test_blackbox_H
#define __LINBOX_test_blackbox_H

#include <iostream>
#include <vector>
//#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/util/field-axpy.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/dense.h"

#include "test-common.h" 

using namespace std;

// Generic tests for black boxes 
/// testBlackbox combines testTranspose and testLinearity

/** Random check that (u^T A) v = u^T (A v).
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
testTranspose (Field                             &F,
			   Blackbox							 &A,
			   LinBox::VectorStream<Vector>      &stream1,
			   LinBox::VectorStream<Vector>      &stream2) 
{
	bool ret = true;

	Vector u, v, uA, Av;

	LinBox::VectorWrapper::ensureDim (u, A.rowdim ());
	LinBox::VectorWrapper::ensureDim (uA, A.coldim ());
	LinBox::VectorWrapper::ensureDim (v, A.coldim ());
	LinBox::VectorWrapper::ensureDim (Av, A.rowdim ());

	LinBox::VectorDomain <Field> VD (A.field());
	//LinBox::VectorDomain <Field> VD (F);
	typename Field::Element r1, r2;
	ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Blackbox transpose test [that u^T(Av) == (uA)^T v]" << std::endl;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		stream1.next (u);
		stream2.next (v);

		//ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

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
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Values are not equal" << endl;
		}

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
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

template <class Field, class BB, class Vector>
static bool
testLinearity (Field                             &F,
	       BB 				 &A,
	       LinBox::VectorStream<Vector>      &stream1,
	       LinBox::VectorStream<Vector>      &stream2) 
{
	bool ret = true, iter_passed;

	size_t n = A.rowdim ();
	size_t m = A.coldim ();

	Vector x, y, xpay, Axpay, Ax, Ay, AxpaAy;
	LinBox::VectorDomain <Field> VD (F);
	typename Field::RandIter r (F);
	typename Field::Element alpha;

	LinBox::VectorWrapper::ensureDim (x, m);
	LinBox::VectorWrapper::ensureDim (y, m);
	LinBox::VectorWrapper::ensureDim (xpay, m);
	LinBox::VectorWrapper::ensureDim (Axpay, n);
	LinBox::VectorWrapper::ensureDim (Ax, n);
	LinBox::VectorWrapper::ensureDim (Ay, n);
	LinBox::VectorWrapper::ensureDim (AxpaAy, n);

	ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Blackbox linearity test [that A.apply to (ax + y) == a A.apply to x + A.apply to y]" << std::endl;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (x);
		stream2.next (y);

		r.random (alpha);

		//ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		VD.write( report << "Input vector x: ", x) << endl;

		VD.write( report << "Input vector y: ", y) << endl;

		F.write( report << "Input alpha: ", alpha) << endl;

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
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	return ret;
}

/** Generic blackbox test 4: combination of tests
 * 
 * Call testTranspose and testLinearity.
 * If large, time apply and applyTranspose.
 * if small, call testSmallBlackbox.
 */
template <class BB> 
static bool 
testBlackbox(BB &A)
{
	size_t largeThresh = 2000; // Above it do timing of apply and applyTr.
	typedef typename BB::Field Field;
	typedef std::vector<typename Field::Element> DenseVector;
	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "testBlackbox on " << A.rowdim() << " by " << A.coldim() << " matrix." << endl;
	
	LinBox::commentator.setMaxDepth(-1);
	bool ret = true;
	typename BB::Field F = A.field();

	/* timing tests */   // I changed the order of all tests. Timing now is the first set of tests and then linearity and transpose
	{
		DenseVector x(A.coldim()), y(A.rowdim());
		for(size_t i = 0; i < A.coldim(); ++i) F.init(x[i], i);
		for(size_t i = 0; i < A.rowdim(); ++i) F.init(y[i], i);
		//A.apply(y, x);
		
		if (A.coldim() >= largeThresh)
			{
				LinBox::commentator.start ("\t--Timing Test (Av)","testApply", 1);
				A.apply(y, x);
				LinBox::commentator.stop (MSG_STATUS (true), (const char *) 0, "testApply");
			}
		
		if (A.rowdim() >= largeThresh)
			{
				LinBox::commentator.start ("\t--Timing Test(v^T A)", 
										   "testApplyTranspose", 1);
				A.applyTranspose(x, y);
				LinBox::commentator.stop (MSG_STATUS (true), (const char *) 0, "testApplyTranspose");
			}
		
	} // timing test block
	
#if 1 
	size_t iterations = 1; 
	typename Field::RandIter r(F);
	LinBox::RandomDenseStream<Field, DenseVector> stream1 (F, r, A.rowdim(), iterations); 
	typename Field::Element x; 
	r.random(x);
	LinBox::RandomDenseStream<Field, DenseVector> stream2 (F, r, A.coldim(), iterations); 

	LinBox::commentator.start ("\t--Testing A(ax+y) = a(Ax) + (Ay)", "testLinearity", 1);
	ret = ret && testLinearity (F, A, stream1, stream2);
	
	LinBox::commentator.stop (MSG_STATUS (ret), 
							  (const char *) 0, "testLinearity");
	
	LinBox::commentator.start ("\t--Testing u^T(Av) = (u^T A)v", 
							   "testTranspose", 1);
	
	LinBox::RandomDenseStream<Field, DenseVector> stream3 (F, r, A.rowdim(), iterations); 
	LinBox::RandomDenseStream<Field, DenseVector> stream4 (F, r, A.coldim(), iterations); 

	ret = ret && testTranspose (F, A, stream3, stream4); 
	LinBox::commentator.stop (MSG_STATUS (ret), 
							  (const char *) 0, "testTranspose");
	
#endif
	
	return ret;
}
 
/** Generic blackbox test 5: test several sizes
 * 
 * Call testTranspose and testLinearity.
 * If large, time apply and applyTranspose.
 * if small, call test
SmallBlackbox.
 */
template <class Field, class Blackbox> 
static bool 
testBB(Field& F) 
{
	bool pass = true;

	Blackbox A(10);
	if (!testBlackbox<Field, vector<typename Field::Element> >(F, A, 1))
		pass = false;
	Blackbox B(10000);
	if (!testBlackbox<Field, vector<typename Field::Element> >(F, B, 1))
		pass = false;

	return pass;
}

#endif // __LINBOX_test_blackbox_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
