/* -*- mode: C++; style: linux -*- */

/* linbox/tests/test-generic.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-11 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added test testFieldAXPY
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __TEST_GENERIC_H
#define __TEST_GENERIC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/vector-factory.h"
#include "linbox/field/vector-domain.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/dense-matrix.h"
#include "linbox/integer.h"

#include "test-common.h"

/* Generic test 1: Test of field operations
 *
 * Test various field oeprations
 *
 * F - Field over which to perform computations
 * title - String to use as the descriptive title of this test
 *
 * Return true on success and false on failure
 */

template<class Field>
bool testField (Field &F, const char *title) 
{
	typename Field::Element zero, one, two, three;
	typename Field::Element a, b, c, d, e, f;

	commentator.start (title, "testField", 5);

	// there is an extra char in the output - bds 3/02
	ostream &report = commentator.report (LinBox::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Field self description: " << F.write (report) << endl;
	commentator.indent (report);
//	report << "field Element 2: " << F.write (report, two) << endl;

	LinBox::integer n, m;
	bool pass = true, part_pass;

	commentator.start ("Testing characteristic/cardinality match");
	part_pass = true;

	F.characteristic (n); 
	F.cardinality (m);

	if (n > 0 && !isPower (m, n)) {
		pass = part_pass = false; 
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Characteristic, cardinality mismatch" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	/* tests for presence of members with minimal check of semantics */
	// these checks need improvement 

	commentator.start ("Testing correctness of 0 and 1");
	part_pass = true;

	F.init (zero, 0);
	F.init (one, 1);

	if (!F.isZero (zero)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: isZero (0) is false" << endl;
	}

	if (F.isZero (one)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: isZero (1) is true" << endl;
	}

	if (F.isOne (zero)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: isOne (0) is true" << endl;
	}

	if (!F.isOne (one)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: isOne (1) is false" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing init/convert");
	part_pass = true;

	if (F.cardinality (m) <= 0)
	  n = 49193295;   // Just using some odd value
	else
	  n = m-1;

	F.init (a, n);  // !!!!!!! non-generic with a finite field ...
	
	{
		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Initial value: " << n << endl;
		commentator.indent (report);
		report << "Result of conversion: ";
		F.write (report, a);
		report << endl;
	}

	F.convert (m, a);

	if (m != n) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: F.convert (m, F.init (a, n)) != n" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing field arithmetic");
	part_pass = true;

	F.init (two, 2);
	F.init (three, 3);

	F.init (b, n-2);
	F.init (d, n-2);
	F.init (e, 3);

	F.add (a, three, two); 
	F.assign (d, three);
	F.addin (d, two);

	{
		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Result of 2 + 3: ";
		F.write (report, a);
		report << endl;

		commentator.indent (report);
		report << "Result of 2 + 3 (inplace): ";
		F.write (report, d);
		report << endl;
	}

	if (!F.areEqual (a, F.init (f, 5)) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of add are incorrect" << endl;
	}

	F.neg (a, two); 
	F.assign (d, two);
	F.negin (d);

	{
		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Result of -2: ";
		F.write (report, a);
		report << endl;

		commentator.indent (report);
		report << "Result of -2 (inplace): ";
		F.write (report, d);
		report << endl;
	}

	if (!F.areEqual (a, F.init (f, -2)) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of neg are incorrect" << endl;
	}

	F.sub (a, three, two);
	F.init (d, 3);
	F.subin (d, two);

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of sub are incorrect" << endl;
	}

	F.mul (a, two, three);
	F.assign (d, two);
	F.mulin (d, three);
	F.init (f, 6);

	if (!F.areEqual (a, f) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of mul are incorrect" << endl;
	}

	F.inv (a, one);
	F.assign (d, one);
	F.invin (d);

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of inv are incorrect" << endl;
	}

	F.div (a, two, two);
	F.assign (d, three);
	F.divin (d, three);
	
	{
		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Result of 2/2: ";
		F.write (report, a);
		report << endl;

		commentator.indent (report);
		report << "Result of 3/3: ";
		F.write (report, d);
		report << endl;
	}

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results of div are incorrect" << endl;
	}

	F.axpy (a, two, three, two); 
	F.assign (d, two);
	F.axpyin (d, two, three);

	if ( !F.areEqual (a, F.init (f, 8)) || !F.areEqual (d, a) ) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Resutls of axpy are incorrect" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("Testing summation of powers of 2");

	//,..
	// 2^101 - 1 vs 1 + 2 + 4 + ... + 2^100

	F.init (a, 1);
	F.init (b, 2);
	F.init (c, 0);

	for (int i = 1; i <= 101; ++i) {
		F.addin (c, a);
		F.mulin (a, b);
	}

	F.subin (a, F.init (f, 1));

	if (!F.areEqual (a, c)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results are incorrect using mul/add" << endl;
	}

	F.assign (d, one);
	for (int i = 1; i < 101; ++i)
		F.axpy (d, two, d, one);

	if (!F.areEqual (a, d)) {
		pass = part_pass = false;
		commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Results are incorrect using axpy" << endl;
	}

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	/* untested so far
	   ostream &write (ostream &os) const 
	   istream &read (istream &is)
	   ostream &write (ostream &os, const Element &x) const 
	   istream &read (istream &is, Element &x) const
	   FieldArchetype (FieldAbstract*, ElementAbstract*, RandIterAbstract* = 0)
	*/

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "testField");

	return pass;
}

/* Generic test 2: Consistency of FieldAXPY object with Field::axpy
 *
 * Constructs two random vectors of n elements each and computes their dot
 * product manually (i.e., not using VectorDomain), using both FieldAXPY and
 * Field::axpy. Checks whether the results are the same. Useful for performance
 * comparisons in regards to the delayed modding out feature.
 *
 * F - Field over which to perform computations
 * n - Number of elements to make vectors
 * iterations - Number of random vectors to check
 * title - String to use as the descriptive title of this test
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testFieldAXPY (Field &F, long n, int iterations, const char *title) 
{
	typedef vector <typename Field::Element> Vector;

	commentator.start (title, "testFieldAXPY", iterations);

	bool ret = true;

	int i, j;

	Vector u(n), v(n);
	typename Field::RandIter r (F);
	typename Field::Element r1, r2;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		for (j = 0; j < n; j++) {
			r.random (u[j]);
			r.random (v[j]);
		}

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector u:  ";
		printVector<Field> (F, report, u);

		commentator.indent (report);
		report << "Input vector v:  ";
		printVector<Field> (F, report, v);

		F.init (r1, 0);

		commentator.start ("Field::axpy");

		for (j = 0; j < n; j++)
			F.axpyin (r1, u[j], v[j]);

		commentator.stop ("Done");

		LinBox::FieldAXPY<Field> field_axpy (F);

		commentator.start ("FieldAXPY");

		for (j = 0; j < n; j++)
			field_axpy.accumulate (u[j], v[j]);

		field_axpy.get (r2);

		commentator.stop ("Done");

		commentator.indent (report);
		report << "Result of Field::axpy: ";
		F.write (report, r1);
		report << endl;

		commentator.indent (report);
		report << "Result of FieldAXPY: ";
		F.write (report, r2);
		report << endl;

		if (!F.areEqual (r1, r2)) {
			ret = false;
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Inner products are not equal" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldAXPY");

	return ret;
}

/* Generic test 3: Application of transpose of a matrix
 *
 * Take the given black box and compute u^T A v via <A^T u, v> and <u, Av> for
 * randomly chosen u and v; check whether the results are equal. In theory, this
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

template <class Field, class Vector>
static bool
testTranspose (Field                             &F,
	       LinBox::BlackboxArchetype<Vector> &A,
	       LinBox::VectorFactory<Vector>     &factory1,
	       LinBox::VectorFactory<Vector>     &factory2) 
{
	bool ret = true;

	int i, j;

	Vector u, v, w;

	LinBox::VectorWrapper::ensureDim (u, A.rowdim ());
	LinBox::VectorWrapper::ensureDim (v, A.coldim ());
	LinBox::VectorWrapper::ensureDim (w, A.coldim ());

	LinBox::VectorDomain <Field> VD (F);
	typename Field::Element r1, r2;

	while (factory1 && factory2) {
		commentator.startIteration (factory1.j ());

		factory1.next (u);
		factory2.next (v);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector u:            ";
		printVector<Field> (F, report, u);

		commentator.indent (report);
		report << "Input vector v:            ";
		printVector<Field> (F, report, v);

		A.apply (w, v);

		commentator.indent (report);
		report << "Result of apply:           ";
		printVector<Field> (F, report, w);

		VD.dot (r1, u, w);

		A.applyTranspose (w, u);

		commentator.indent (report);
		report << "Result of transpose apply: ";
		printVector<Field> (F, report, w);

		VD.dot (r2, w, v);

		commentator.indent (report);
		report << "<u, Av>:  ";
		F.write (report, r1);
		report << endl;

		commentator.indent (report);
		report << "<A^T u, v>:  ";
		F.write (report, r2);
		report << endl;

		if (!F.areEqual (r1, r2)) {
			ret = false;
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	return ret;
}

/* Generic test 4: Linearity of black boxes
 *
 * Given an arbitrary black box A, compute A(x+alpha y) and Ax+alphaAy and check equality.
 *
 * F - Field over which to perform computations
 * A - Black box of which to compute the dense representation
 * factory1 - Factory for x's
 * factory2 - Factory for y's
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool
testLinearity (Field                              &F,
	       LinBox::BlackboxArchetype <Vector> &A,
	       LinBox::VectorFactory<Vector>      &factory1,
	       LinBox::VectorFactory<Vector>      &factory2) 
{
	bool ret = true, iter_passed;

	int i, j;

	size_t n = A.rowdim ();
	size_t m = A.coldim ();

	Vector x, y, xpay, Axpay, Ax, Ay, AxpaAy;
	LinBox::VectorDomain <Field> VD (F);
	typename Field::RandIter r (F);
	typename Field::Element alpha;

	LinBox::VectorWrapper::ensureDim (x, n);
	LinBox::VectorWrapper::ensureDim (y, n);
	LinBox::VectorWrapper::ensureDim (xpay, n);
	LinBox::VectorWrapper::ensureDim (Axpay, m);
	LinBox::VectorWrapper::ensureDim (Ax, m);
	LinBox::VectorWrapper::ensureDim (Ay, m);
	LinBox::VectorWrapper::ensureDim (AxpaAy, m);

	while (factory1 && factory2) {
		commentator.startIteration (factory1.j ());

		iter_passed = true;

		factory1.next (x);
		factory2.next (y);

		r.random (alpha);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector x: ";
		printVector<Field> (F, report, x);

		commentator.indent (report);
		report << "Input vector y: ";
		printVector<Field> (F, report, y);

		commentator.indent (report);
		report << "Input alpha: ";
		F.write (report, alpha) << endl;

		VD.axpy (xpay, x, alpha, y);
		A.apply (Axpay, xpay);

		A.apply (Ax, x);
		A.apply (Ay, y);
		VD.axpy (AxpaAy, Ax, alpha, Ay);

		commentator.indent (report);
		report << "   x+alpha y = ";
		printVector<Field> (F, report, xpay);

		commentator.indent (report);
		report << "A(x+alpha y) = ";
		printVector<Field> (F, report, Axpay);

		commentator.indent (report);
		report << "          Ax = ";
		printVector<Field> (F, report, Ax);

		commentator.indent (report);
		report << "          Ay = ";
		printVector<Field> (F, report, Ay);

		commentator.indent (report);
		report << " Ax+alpha Ay = ";
		printVector<Field> (F, report, AxpaAy);

		if (!VD.areEqual (Axpay, AxpaAy))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	return ret;
}

#endif // __TEST_GENERIC_H
