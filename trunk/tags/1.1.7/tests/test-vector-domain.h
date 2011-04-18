/* linbox/tests/test-vector-domain.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 * See COPYING for license information.
 *
 * Extracted from test-generic by bds.
 * test-generic was written by Bradford Hovinen <hovinen@cis.udel.edu>
 */

#ifndef __LINBOX_test_vector_domain_H
#define __LINBOX_test_vector_domain_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/util/field-axpy.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/dense.h"
#include "linbox/integer.h"

#include "test-common.h" 


/** Test 1: Dot product of vectors
 *
 * Construct two random vectors and compute their dot product
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * stream1 - Stream for first family of vectors
 * stream2 - Stream for second family of vectors
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector1, class Vector2>
static bool testDotProduct (Field &F, const char *text, LinBox::VectorStream<Vector1> &stream1, LinBox::VectorStream<Vector2> &stream2) 
{
	std::ostringstream str;

	str << "\t--Testing " << text << " dot product" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), "testDotProduct", stream1.m ());

	bool ret = true;

	Vector1 v1;
	Vector2 v2;
	typename Field::Element sigma, rho;

	LinBox::VectorDomain<Field> VD (F);

	size_t j;

	LinBox::VectorWrapper::ensureDim (v1, stream1.n ());
	LinBox::VectorWrapper::ensureDim (v2, stream2.n ());

	LinBox::Timer timer;
	double totaltime = 0.0;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		F.init (sigma, 0);

		stream1.next (v1);
		stream2.next (v2);

		for (j = 0; j < stream1.n (); j++)
			F.axpyin (sigma,
				  LinBox::VectorWrapper::constRef<Field> (v1, j),
				  LinBox::VectorWrapper::constRef<Field> (v2, j));

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1 of size " << v1.size() << ":  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2 of size " << v2.size() << ":  ";
		VD.write (report, v2) << std::endl;

		timer.start ();
		VD.dot (rho, v1, v2);
		timer.stop ();
		totaltime += timer.realtime ();

		report << "True dot product: ";
		F.write (report, sigma) << std::endl;

		report << "Dot product from vector domain: ";
		F.write (report, rho) << std::endl;

		if (!F.areEqual (sigma, rho)) {
			ret = false;
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Dot products are not equal" << std::endl;
		}

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
		<< "Average time for dot product: " << totaltime / stream1.m () << std::endl;

	LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDotProduct");

	stream1.reset ();
	stream2.reset ();

	return ret;
}

/** Test 2: Vector-vector addition, vector-scalar multiply
 *
 * Construct two random vectors x and y and a random element a and compute (x +
 * a*y) and a*(y + a^-1*x) using vector add, sub, and mul. Check whether the
 * results are equal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testAddMul (Field &F, const char *text, LinBox::VectorStream<Vector> &stream1, LinBox::VectorStream<Vector> &stream2) 
{
	std::ostringstream str;

	str << "\t--Testing " << text << " vector add, mul" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), "testAddMul", stream1.m ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3, v4;
	typename Field::Element a;
	typename Field::Element ainv;
	typename Field::Element aneg;
	typename Field::RandIter r (F);

	LinBox::VectorWrapper::ensureDim (v1, stream1.n ());
	LinBox::VectorWrapper::ensureDim (v2, stream2.n ());
	LinBox::VectorWrapper::ensureDim (v3, stream1.n ());
	LinBox::VectorWrapper::ensureDim (v4, stream1.n ());

	LinBox::VectorDomain<Field> VD (F);

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (v1);
		stream2.next (v2);

		do r.random (a); while (F.isZero (a));

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1 of size " << v1.size() << ":  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2 of size " << v2.size() << ":  ";
		VD.write (report, v2) << std::endl;

		report << "Element a:  ";
		F.write (report, a) << std::endl;

		F.inv (ainv, a);
		F.neg (aneg, a);
		VD.mul (v3, v1, ainv);
		report << "          a^-1 * x = ";
		VD.write (report, v3) << std::endl;
		report.flush ();

		VD.addin (v3, v2);
		report << "      y + a^-1 * x = ";
		VD.write (report, v3) << std::endl;
		report.flush ();

		VD.mulin (v2, a);
		report << "             a * y = ";
		VD.write (report, v2) << std::endl;
		report.flush ();

		VD.add (v4, v1, v2);
		report << "         x + a * y = ";
		VD.write (report, v4) << std::endl;
		report.flush ();

		VD.mulin (v3, a);
		report << "a * (y + a^-1 * x) = ";
		VD.write (report, v3) << std::endl;
		report.flush ();

		if (!VD.areEqual (v3, v4))
			ret = iter_passed = false;

		if (!iter_passed)
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (x + a*y) != a*(y + a^-1*x)" << std::endl;

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAddMul");

	stream1.reset ();
	stream2.reset ();

	return ret;
}

/** Test 3: Vector-vector subtraction, vector-scalar multiply
 *
 * Construct two random vectors x and y and a random element a and compute (x -
 * a*y) and a*(a^-1*x - y) using vector add, sub, and mul. Check whether the
 * results are equal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testSubMul (Field &F, const char *text, LinBox::VectorStream<Vector> &stream1, LinBox::VectorStream<Vector> &stream2) 
{
	std::ostringstream str;

	str << "\t--Testing " << text << " vector sub, mul" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), "testSubMul", stream1.m ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3, v4;
	typename Field::Element a;
	typename Field::Element ainv;
	typename Field::Element aneg;
	typename Field::RandIter r (F);

	LinBox::VectorWrapper::ensureDim (v1, stream1.n ());
	LinBox::VectorWrapper::ensureDim (v2, stream2.n ());
	LinBox::VectorWrapper::ensureDim (v3, stream1.n ());
	LinBox::VectorWrapper::ensureDim (v4, stream1.n ());

	LinBox::VectorDomain<Field> VD (F);

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (v1);
		stream2.next (v2);

		do r.random (a); while (F.isZero (a));

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1 of size " << v1.size() << ":  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2 of size " << v2.size() << ":  ";
		VD.write (report, v2) << std::endl;

		report << "Element a:  ";
		F.write (report, a) << std::endl;

		F.inv (ainv, a);
		F.neg (aneg, a);
		VD.mul (v3, v1, ainv);
		report << "          a^-1 * x = ";
		VD.write (report, v3) << std::endl;
		report.flush ();

		VD.subin (v3, v2);
		report << "      a^-1 * x - y = ";
		VD.write (report, v3) << std::endl;
		report.flush ();

		VD.mulin (v2, a);
		report << "             a * y = ";
		VD.write (report, v2) << std::endl;
		report.flush ();

		VD.sub (v4, v1, v2);
		report << "         x - a * y = ";
		VD.write (report, v4) << std::endl;
		report.flush ();

		VD.mulin (v3, a);
		report << "a * (y - a^-1 * x) = ";
		VD.write (report, v4) << std::endl;
		report.flush ();

		if (!VD.areEqual (v3, v4))
			ret = iter_passed = false;

		if (!iter_passed)
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (x - a*y) != a*(a^-1*x - y)" << std::endl;

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSubMul");

	stream1.reset ();
	stream2.reset ();

	return ret;
}

/** Test 4: Vector-vector axpy
 *
 * Construct two random vectors x and y and a random element a and compute (x +
 * a*y) - a*(y + a^-1*x) using vector axpy. Check whether the result is 0.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testAXPY (Field &F, const char *text, LinBox::VectorStream<Vector> &stream1, LinBox::VectorStream<Vector> &stream2) 
{
	std::ostringstream str;
	str << "\t--Testing " << text << " vector axpy" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), "testAXPY", stream1.m ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3, v4;
	typename Field::Element a;
	typename Field::Element ainv;
	typename Field::Element aneg;
	typename Field::RandIter r (F);

	LinBox::VectorWrapper::ensureDim (v1, stream1.n ());
	LinBox::VectorWrapper::ensureDim (v2, stream2.n ());
	LinBox::VectorWrapper::ensureDim (v3, stream1.n ());
	LinBox::VectorWrapper::ensureDim (v4, stream1.n ());

	LinBox::VectorDomain<Field> VD (F);

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (v1);
		stream2.next (v2);

		do r.random (a); while (F.isZero (a));

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1 of size " << v1.size() << ":  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2 of size " << v2.size() << ":  ";
		VD.write (report, v2) << std::endl;

		report << "Element a:  ";
		F.write (report, a) << std::endl;

		F.inv (ainv, a);
		F.neg (aneg, a);
		VD.axpy (v3, a, v2, v1);
		VD.axpy (v4, ainv, v1, v2);
		VD.axpyin (v3, aneg, v4);

		report << "Output vector:  ";
		VD.write (report, v3) << std::endl;

		if (!VD.isZero (v3))
			ret = iter_passed = false;

		if (!iter_passed)
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (x + a*y) - a*(y + a^-1*x) != 0" << std::endl;

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAXPY");

	stream1.reset ();
	stream2.reset ();

	return ret;
}

/** Test 5: Copy and areEqual
 *
 * Constructs a random vector and copies it to another vector. Then checks equality.
 *
 * F - Field over which to perform computations
 * text - Text to use for test
 * stream - Stream generating vectors
 * stream2 - Dummy stream of second vector type to trick the compiler
 *
 * Return true on success and false on failure
 */
template <class Field, class Vector1, class Vector2>
static bool testCopyEqual (Field &F, const char *text, LinBox::VectorStream<Vector1> &stream, LinBox::VectorStream<Vector2> &stream2) 
{
	std::ostringstream str;

	str << "\t--Testing " << text << " vector copy, areEqual" << std::ends;
	LinBox::commentator.start (str.str ().c_str (), "testCopyEqual", stream.m ());

	bool ret = true;
	bool iter_passed;

	Vector1 v;
	Vector2 w;

	LinBox::VectorWrapper::ensureDim (v, stream.n ());
	LinBox::VectorWrapper::ensureDim (w, stream.n ());

	LinBox::VectorDomain<Field> VD (F);

	while (stream) {
		LinBox::commentator.startIteration (stream.j ());

		iter_passed = true;

		stream.next (v);
		VD.copy (w, v);

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:   ";
		VD.write (report, v) << std::endl;

		report << "Output vector:  ";
		VD.write (report, w) << std::endl;

		if (!VD.areEqual (v, w))
			ret = iter_passed = false;

		if (!iter_passed)
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << std::endl;

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testCopyEqual");

	stream.reset ();
	stream2.reset ();

	return ret;
}

#endif // __LINBOX_test_vector_domain_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
