/* -*- mode: C++; style: linux -*- */

/* tests/test-vector-domain.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <strstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/vector-domain.h"
#include "linbox/util/vector-factory.h"
#include "linbox/util/field-axpy.h"

#include "test-common.h"

using namespace std;
using namespace LinBox;

/* Test 1: Dot product of vectors
 *
 * Construct two random vectors and compute their dot product
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * factory1 - Factory for first family of vectors
 * factory2 - Factory for second family of vectors
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector1, class Vector2>
static bool testDotProduct (Field &F, const char *text, VectorFactory<Vector1> &factory1, VectorFactory<Vector2> &factory2) 
{
	char buf[128];
	ostrstream str (buf, 128);
	str << "Testing " << text << " dot product" << ends;
	commentator.start (buf, "testDotProduct", factory1.m ());

	bool ret = true;

	Vector1 v1;
	Vector2 v2;
	typename Field::Element sigma, rho;

	VectorDomain<Field> VD (F);

	int j;

	VectorWrapper::ensureDim (v1, factory1.n ());
	VectorWrapper::ensureDim (v2, factory2.n ());

	Timer timer;
	double totaltime = 0.0;

	while (factory1 && factory2) {
		commentator.startIteration (factory1.j ());

		F.init (sigma, 0);

		factory1.next (v1);
		factory2.next (v2);

		for (j = 0; j < factory1.n (); j++)
			F.axpyin (sigma,
				  VectorWrapper::constRef<Field, Vector1> (v1, j),
				  VectorWrapper::constRef<Field, Vector2> (v2, j));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		timer.start ();
		VD.dot (rho, v1, v2);
		timer.stop ();
		totaltime += timer.realtime ();

		commentator.indent (report);
		report << "True dot product: ";
		F.write (report, sigma);
		report << endl;

		commentator.indent (report);
		report << "Dot product from vector domain: ";
		F.write (report, rho);
		report << endl;

		if (!F.areEqual (sigma, rho)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Dot products are not equal" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
		<< "Average time for dot product: " << totaltime / factory1.m () << endl;

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDotProduct");

	return ret;
}

/* Test 2: Vector-vector addition, vector-scalar multiply
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
static bool testAddMul (Field &F, const char *text, VectorFactory<Vector> &factory1, VectorFactory<Vector> &factory2) 
{
	char buf[128];
	ostrstream str (buf, 128);
	str << "Testing " << text << " vector add, mul" << ends;
	commentator.start (buf, "testAddMul", factory1.m ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3, v4;
	typename Field::Element a;
	typename Field::Element ainv;
	typename Field::Element aneg;
	typename Field::RandIter r (F);

	VectorWrapper::ensureDim (v1, factory1.n ());
	VectorWrapper::ensureDim (v2, factory2.n ());
	VectorWrapper::ensureDim (v3, factory1.n ());
	VectorWrapper::ensureDim (v4, factory1.n ());

	VectorDomain<Field> VD (F);

	int i, j;

	while (factory1 && factory2) {
		commentator.startIteration (factory1.j ());

		iter_passed = true;

		factory1.next (v1);
		factory2.next (v2);

		do r.random (a); while (F.isZero (a));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		commentator.indent (report);
		report << "Element a:  ";
		F.write (report, a);
		report << endl;

		F.inv (ainv, a);
		F.neg (aneg, a);
		VD.mul (v3, v1, ainv);
		commentator.indent (report);
		report << "          a^-1 * x = ";
		printVector<Field> (F, report, v3);
		report.flush ();

		VD.addin (v3, v2);
		commentator.indent (report);
		report << "      y + a^-1 * x = ";
		printVector<Field> (F, report, v3);
		report.flush ();

		VD.mulin (v2, a);
		commentator.indent (report);
		report << "             a * y = ";
		printVector<Field> (F, report, v2);
		report.flush ();

		VD.add (v4, v1, v2);
		commentator.indent (report);
		report << "         x + a * y = ";
		printVector<Field> (F, report, v4);
		report.flush ();

		VD.mulin (v3, a);
		commentator.indent (report);
		report << "a * (y + a^-1 * x) = ";
		printVector<Field> (F, report, v3);
		report.flush ();

		if (!VD.areEqual (v3, v4))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (x + a*y) != a*(y + a^-1*x)" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAddMul");

	return ret;
}

/* Test 3: Vector-vector subtraction, vector-scalar multiply
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
static bool testSubMul (Field &F, const char *text, VectorFactory<Vector> &factory1, VectorFactory<Vector> &factory2) 
{
	char buf[128];
	ostrstream str (buf, 128);
	str << "Testing " << text << " vector sub, mul" << ends;
	commentator.start (buf, "testSubMul", factory1.m ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3, v4;
	typename Field::Element a;
	typename Field::Element ainv;
	typename Field::Element aneg;
	typename Field::RandIter r (F);

	VectorWrapper::ensureDim (v1, factory1.n ());
	VectorWrapper::ensureDim (v2, factory2.n ());
	VectorWrapper::ensureDim (v3, factory1.n ());
	VectorWrapper::ensureDim (v4, factory1.n ());

	VectorDomain<Field> VD (F);

	int i, j;

	while (factory1 && factory2) {
		commentator.startIteration (factory1.j ());

		iter_passed = true;

		factory1.next (v1);
		factory2.next (v2);

		do r.random (a); while (F.isZero (a));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		commentator.indent (report);
		report << "Element a:  ";
		F.write (report, a);
		report << endl;

		F.inv (ainv, a);
		F.neg (aneg, a);
		VD.mul (v3, v1, ainv);
		commentator.indent (report);
		report << "          a^-1 * x = ";
		printVector<Field> (F, report, v3);
		report.flush ();

		VD.subin (v3, v2);
		commentator.indent (report);
		report << "      a^-1 * x - y = ";
		printVector<Field> (F, report, v3);
		report.flush ();

		VD.mulin (v2, a);
		commentator.indent (report);
		report << "             a * y = ";
		printVector<Field> (F, report, v2);
		report.flush ();

		VD.sub (v4, v1, v2);
		commentator.indent (report);
		report << "         x - a * y = ";
		printVector<Field> (F, report, v4);
		report.flush ();

		VD.mulin (v3, a);
		commentator.indent (report);
		report << "a * (y - a^-1 * x) = ";
		printVector<Field> (F, report, v4);
		report.flush ();

		if (!VD.areEqual (v3, v4))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (x - a*y) != a*(a^-1*x - y)" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSubMul");

	return ret;
}

/* Test 4: Vector-vector axpy
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
static bool testAXPY (Field &F, const char *text, VectorFactory<Vector> &factory1, VectorFactory<Vector> &factory2) 
{
	char buf[128];
	ostrstream str (buf, 128);
	str << "Testing " << text << " vector axpy" << ends;
	commentator.start (buf, "testAXPY", factory1.m ());

	bool ret = true;
	bool iter_passed;

	Vector v1, v2, v3, v4;
	typename Field::Element a;
	typename Field::Element ainv;
	typename Field::Element aneg;
	typename Field::RandIter r (F);

	VectorWrapper::ensureDim (v1, factory1.n ());
	VectorWrapper::ensureDim (v2, factory2.n ());
	VectorWrapper::ensureDim (v3, factory1.n ());
	VectorWrapper::ensureDim (v4, factory1.n ());

	VectorDomain<Field> VD (F);

	int i, j;

	while (factory1 && factory2) {
		commentator.startIteration (factory1.j ());

		iter_passed = true;

		factory1.next (v1);
		factory2.next (v2);

		do r.random (a); while (F.isZero (a));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		commentator.indent (report);
		report << "Element a:  ";
		F.write (report, a);
		report << endl;

		F.inv (ainv, a);
		F.neg (aneg, a);
		VD.axpy (v3, v1, a, v2);
		VD.axpy (v4, v2, ainv, v1);
		VD.axpyin (v3, aneg, v4);

		commentator.indent (report);
		report << "Output vector:  ";
		printVector<Field> (F, report, v3);

		if (!VD.isZero (v3))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (x + a*y) - a*(y + a^-1*x) != 0" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAXPY");

	return ret;
}

/* Test 5: Copy and areEqual
 *
 * Constructs a random vector and copies it to another vector. Then checks equality.
 *
 * F - Field over which to perform computations
 * text - Text to use for test
 * factory - Factory generating vectors
 * factory2 - Dummy factory of second vector type to trick the compiler
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector1, class Vector2>
static bool testCopyEqual (Field &F, const char *text, VectorFactory<Vector1> &factory, VectorFactory<Vector2> &factory2) 
{
	char buf[128];
	ostrstream str (buf, 128);
	str << "Testing " << text << " vector copy, areEqual" << ends;
	commentator.start (buf, "testCopyEqual", factory.m ());

	bool ret = true;
	bool iter_passed;

	Vector1 v;
	Vector2 w;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, factory.n ());

	VectorDomain<Field> VD (F);

	int i, j;

	while (factory) {
		commentator.startIteration (factory.j ());

		iter_passed = true;

		factory.next (v);
		VD.copy (w, v);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:   ";
		printVector<Field> (F, report, v);

		commentator.indent (report);
		report << "Output vector:  ";
		printVector<Field> (F, report, w);

		if (!VD.areEqual (v, w))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testCopyEqual");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long n = 100;
	static integer q = 101;
	static int iterations = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to N (default 100)",   TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",   TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Vector domain test suite" << endl << endl;
	cout.flush ();

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (2);

	RandomDenseVectorFactory<Modular<long> > factory1 (F, n, iterations), factory2 (F, n, iterations);
	RandomSparseSeqVectorFactory<Modular<long> > factory3 (F, n, n / 10, iterations), factory4 (F, n, n / 10, iterations);
	RandomSparseMapVectorFactory<Modular<long> > factory5 (F, n, n / 10, iterations), factory6 (F, n, n / 10, iterations);

	if (!testDotProduct (F, "dense/dense", factory1, factory2)) pass = false;

	factory1.reset ();
	if (!testDotProduct (F, "sparse sequence/dense", factory3, factory1)) pass = false;

	factory1.reset ();
	if (!testDotProduct (F, "sparse associative/dense", factory5, factory1)) pass = false;

	factory3.reset ();
	if (!testDotProduct (F, "sparse sequence/sparse sequence", factory3, factory4)) pass = false;

	factory3.reset ();
	factory5.reset ();
	if (!testDotProduct (F, "sparse associative/sparse sequence", factory5, factory3)) pass = false;

	factory5.reset ();
	if (!testDotProduct (F, "sparse associative/sparse associative", factory5, factory6)) pass = false;

	factory1.reset ();
	factory2.reset ();
	if (!testAddMul (F, "dense", factory1, factory2)) pass = false;

	factory3.reset ();
	factory4.reset ();
	if (!testAddMul (F, "sparse sequence", factory3, factory4)) pass = false;

	factory5.reset ();
	factory6.reset ();
	if (!testAddMul (F, "sparse associative", factory5, factory6)) pass = false;

	factory1.reset ();
	factory2.reset ();
	if (!testSubMul (F, "dense", factory1, factory2)) pass = false;

	factory3.reset ();
	factory4.reset ();
	if (!testSubMul (F, "sparse sequence", factory3, factory4)) pass = false;

	factory5.reset ();
	factory6.reset ();
	if (!testSubMul (F, "sparse associative", factory5, factory6)) pass = false;

	factory1.reset ();
	factory2.reset ();
	if (!testAXPY (F, "dense", factory1, factory2)) pass = false;

	factory3.reset ();
	factory4.reset ();
	if (!testAXPY (F, "sparse sequence", factory3, factory4)) pass = false;

	factory5.reset ();
	factory6.reset ();
	if (!testAXPY (F, "sparse associative", factory5, factory6)) pass = false;

	factory1.reset ();
	if (!testCopyEqual (F, "dense/dense", factory1, factory1)) pass = false;

	factory1.reset ();
	if (!testCopyEqual (F, "dense/sparse sequence", factory1, factory3)) pass = false;

	factory1.reset ();
	if (!testCopyEqual (F, "dense/sparse associative", factory1, factory5)) pass = false;

	factory3.reset ();
	if (!testCopyEqual (F, "sparse sequence/dense", factory3, factory1)) pass = false;

	factory3.reset ();
	if (!testCopyEqual (F, "sparse sequence/sparse sequence", factory3, factory3)) pass = false;

	factory3.reset ();
	if (!testCopyEqual (F, "sparse sequence/sparse associative", factory3, factory5)) pass = false;

	factory5.reset ();
	if (!testCopyEqual (F, "sparse associative/dense", factory5, factory1)) pass = false;

	factory5.reset ();
	if (!testCopyEqual (F, "sparse associative/sparse sequence", factory5, factory3)) pass = false;

	factory5.reset ();
	if (!testCopyEqual (F, "sparse associative/sparse associative", factory5, factory5)) pass = false;

	return pass ? 0 : -1;
}
