/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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
#include <sstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/util/field-axpy.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/dense.h"
#include "linbox/integer.h"

#include "test-common.h"

/* Modular exponentiation */
using namespace std;

template <class Field>
typename Field::Element expt (const Field &F, typename Field::Element &res, const typename Field::Element &a, LinBox::integer &n) 
{
	if (n == 0) {
		F.init (res, 1);
	}
	else if (n == 1) {
		F.assign (res, a);
	}
	else if (n[0] & 1) {
		n -= 1;
		expt (F, res, a, n);
		F.mulin (res, a);
	} else {
		n /= 2;
		expt (F, res, a, n);
		F.mulin (res, res);
	}

	return res;
}

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

		report << "Result of -2 (inplace): ";
		F.write (report, d);
		report << endl;
	}

	if (!F.areEqual (a, F.init (f, -2)) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
		report << "ERROR: Results of neg are incorrect (";
		F.write (report, f);
		report << ")" << endl;
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

		report << "Result of 3/3: ";
		F.write (report, d);
		report << endl;
	}

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) {
		pass = part_pass = false;
		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
		report << "ERROR: Results of div are incorrect" << endl;
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
	F.init (f, 1);
	F.subin (a, f);

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


/* Tests of algebraic properties of fields */

/* Generic test 6: Negation of elements
 *
 * Negates random elements and checks that they are true negatives
 */

template <class Field>
bool testFieldNegation (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " negation" << ends;
	commentator.start (str.str ().c_str (), "testFieldNegation", iterations);

	typename Field::Element a, neg_a, neg_a_a, zero;
	typename Field::RandIter r (F);

	bool ret = true;

	F.init (zero, 0);

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);
		
		r.random (a);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random element a: ";
		F.write (report, a) << endl;

		F.neg (neg_a, a);

		report << "-a = ";
		F.write (report, neg_a) << endl;

		F.add (neg_a_a, neg_a, a);

		report << "a + -a = ";
		F.write (report, neg_a_a) << endl;

		if (!F.areEqual (neg_a_a, zero)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: a + -a != 0" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldNegation");

	return ret;
}

/* Generic test 5: Inversion of elements
 *
 * Inverts random elements and checks that they are true inverses
 */

template <class Field>
bool testFieldInversion (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " inversion" << ends;
	commentator.start (str.str ().c_str (), "testFieldInversion", iterations);

	typename Field::Element a, ainv, aainv, one;
	typename Field::RandIter r (F);

	bool ret = true;

	F.init (one, 1);

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		do r.random (a); while (F.isZero (a));

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random element a: ";
		F.write (report, a) << endl;

		F.inv (ainv, a);

		report << "a^-1 = ";
		F.write (report, ainv) << endl;

		F.mul (aainv, ainv, a);

		report << "a a^-1 = ";
		F.write (report, aainv) << endl;

		if (!F.areEqual (aainv, one)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: a a^-1 != 1" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldInversion");

	return ret;
}

/* Generic test 7: Commutativity and distributivity of addition
 * and multiplication
 *
 * Given random field elements 'a', 'b', and 'c', checks that
 * (a + b) * c = a * c + b * c = c * (a + b) = b * c + a * c
 */

template <class Field>
bool testFieldAxioms (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " commutativity, distributivity" << ends;
	commentator.start (str.str ().c_str (), "testFieldAxioms", iterations);

	typename Field::Element a, b, c, a_b, a_bc, ac, bc, ac_bc, ca_b, bc_ac;
	typename Field::RandIter r (F);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);
		r.random (c);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << ", c = ";
		F.write (report, c) << endl;

		F.add (a_b, a, b);
		F.mul (a_bc, a_b, c);
		F.mul (ca_b, c, a_b);
		F.mul (ac, a, c);
		F.mul (bc, b, c);
		F.add (ac_bc, ac, bc);
		F.add (bc_ac, bc, ac);

		report << "(a + b) * c = ";
		F.write (report, a_bc) << endl;

		report << "a * c + b * c = ";
		F.write (report, ac_bc) << endl;

		report << "c * (a + b) = ";
		F.write (report, ca_b) << endl;

		report << "b * c + a * c = ";
		F.write (report, bc_ac) << endl;

		if (!F.areEqual (a_bc, ac_bc) || !F.areEqual (a_bc, ca_b) || !F.areEqual (a_bc, bc_ac)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Results are not equal" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldAxioms");

	return ret;
}

/* Generic test 7: Associativity of addition and multiplication
 *
 * Given random field elements 'a', 'b', and 'c', checks that
 * (a * b) * c = a * (b * c) and (a + b) + c = a + (b + c)
 */

template <class Field>
bool testFieldAssociativity (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " associativity" << ends;
	commentator.start (str.str ().c_str (), "testFieldAssociativity", iterations);

	typename Field::Element a, b, c, a_b, b_c, a_bc, ab_c;
	typename Field::RandIter r (F);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);
		r.random (c);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << ", c = ";
		F.write (report, c) << endl;

		F.add (a_b, a, b);
		F.add (ab_c, a_b, c);
		F.add (b_c, b, c);
		F.add (a_bc, a, b_c);

		report << "(a + b) + c = ";
		F.write (report, ab_c) << endl;

		report << "a + (b + c) = ";
		F.write (report, a_bc) << endl;

		if (!F.areEqual (ab_c, a_bc)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Results are not equal" << endl;
			ret = false;
		}

		F.mul (a_b, a, b);
		F.mul (ab_c, a_b, c);
		F.mul (b_c, b, c);
		F.mul (a_bc, a, b_c);

		report << "(a * b) * c = ";
		F.write (report, ab_c) << endl;

		report << "a * (b * c) = ";
		F.write (report, a_bc) << endl;

		if (!F.areEqual (ab_c, a_bc)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Results are not equal" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldAssociativity");

	return ret;
}

/* Generic test 2: Geometric summation
 *
 * Generates a random field element 'a' and raises it through repeated
 * exponentiation to the power n. Takes the sum k of all intermediate values and
 * checks that a^n = (k-1)/(a-1).
 */

template <class Field>
bool testGeometricSummation (const Field &F, const char *name, unsigned int iterations, unsigned int n) 
{
	std::ostringstream str;
	str << "Testing " << name << " geometric summation" << ends;
	commentator.start (str.str ().c_str (), "testGeometricSummation", iterations);

	typename Field::Element a, a_n, k, zero, one;
	typename Field::RandIter r (F);

	F.init (zero, 0);
	F.init (one, 1);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		do r.random (a); while (F.areEqual (a, one));

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random element a: ";
		F.write (report, a) << endl;

		F.assign (k, one);
		F.assign (a_n, a);

		for (unsigned int j = 0; j < n; ++j) {
			F.addin (k, a_n);
			F.mulin (a_n, a);
		}

		report << "a^n = ";
		F.write (report, a_n) << endl;

		report << "sum(a^i, i = 0..n-1) = ";
		F.write (report, k) << endl;

		F.subin (a_n, one);
		F.subin (a, one);
		F.divin (a_n, a);

		report << "(a^n - 1) / (a - 1) = ";
		F.write (report, a_n) << endl;

		if (!F.areEqual (k, a_n)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Field elements are not equal" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testGeometricSummation");

	return ret;
}

/* Generic test 3: Test of field characteristic
 *
 * Take random field elements and add them p times, where p is the
 * characteristic of the field. Checks that the sum is 0. The test is not too
 * useful when the characteristic of the field is 0, but it should still work
 * correctly.
 */

template <class Field>
bool testFieldCharacteristic (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " characteristic" << ends;
	commentator.start (str.str ().c_str (), "testFieldCharacteristic", iterations);

	LinBox::integer p, j;
	typename Field::Element a, sigma, zero;
	typename Field::RandIter r (F);

	F.characteristic (p);
	F.init (zero, 0);

	bool ret = true;

	ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Field characteristic: " << p << endl;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random element a: ";
		F.write (report, a) << endl;

		F.assign (sigma, zero);

		for (j = 0; j < p; j += 1)
			F.addin (sigma, a);

		report << "p a = ";
		F.write (report, sigma) << endl;

		if (!F.isZero (sigma)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: p a != 0" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldCharacteristic");

	return ret;
}

/* Generic test 4: The Freshman's Dream
 *
 * Generates two random field elements 'a' and 'b', and checks whether
 * (a + b)^p = a^p + b^p, where p is the characteristic of the field. Bails out
 * with an error if the field is of characteristic 0.
 */

template <class Field>
bool testFreshmansDream (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "Testing " << name << " Freshman's Dream" << ends;
	commentator.start (str.str ().c_str (), "testFreshmansDream", iterations);

	LinBox::integer c, j;

	F.characteristic (c);

	if (iszero (c)) {
		commentator.stop ("skipping", "Field characteristic is 0, so this test makes no sense", "testFreshmansDream");
		return true;
	}

	bool ret = true;

	typename Field::RandIter r (F);
	typename Field::Element a, b, a_b, a_b_p, a_p, b_p, a_p_b_p;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << endl;

		F.add (a_b, a, b);

		j = c; expt (F, a_b_p, a_b, j);
		j = c; expt (F, a_p, a, j);
		j = c; expt (F, b_p, b, j);

		F.add (a_p_b_p, a_p, b_p);

		report << "(a + b)^p = ";
		F.write (report, a_b_p);
		report << endl;

		report << "      a^p = ";
		F.write (report, a_p);
		report << endl;

		report << "      b^p = ";
		F.write (report, b_p);
		report << endl;

		report << "a^p + b^p = ";
		F.write (report, a_p_b_p);
		report << endl;

		if (!F.areEqual (a_b_p, a_p_b_p)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (a + b)^p != a^p + b^p" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFreshmansDream");

	return ret;
}


/* Tests of field features */ 

/* Generic test 7: Consistency of in-place and out-of-place arithmetic
 *
 * Generates random elements 'a' and 'b' and performs all basic arithmetic
 * operations in-place and out-of-place, checking for consistency
 */

template <class Field>
bool testArithmeticConsistency (const Field &F, const char *name, unsigned int iterations)
{
	std::ostringstream str;
	str << "Testing " << name << " in-place/out-of-place arithmetic consistency" << ends;
	commentator.start (str.str ().c_str (), "testArithmeticConsistency", iterations);

	bool ret = true;

	typename Field::RandIter r (F);
	typename Field::Element a, b, c1, c2;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);

		// This should be immaterial, since we have already "proven" commutativity
		if (F.isZero (a) && !F.isZero (b))
			std::swap (a, b);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << endl;

		F.add (c1, a, b);
		F.assign (c2, a);
		F.addin (c2, b);

		report << "a + b = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Consistency failure for addition" << endl;
			ret = false;
		}

		F.sub (c1, a, b);
		F.assign (c2, a);
		F.subin (c2, b);

		report << "a - b = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Consistency failure for subtraction" << endl;
			ret = false;
		}

		F.neg (c1, a);
		F.assign (c2, a);
		F.negin (c2);

		report << "-a = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Consistency failure for negation" << endl;
			ret = false;
		}

		F.mul (c1, a, b);
		F.assign (c2, a);
		F.mulin (c2, b);

		report << "a * b = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Consistency failure for multiplication" << endl;
			ret = false;
		}

		if (!F.isZero (a)) {
			F.div (c1, b, a);
			F.assign (c2, b);
			F.divin (c2, a);

			report << "a * b = (out-of-place) ";
			F.write (report, c1) << ", (in-place) ";
			F.write (report, c2) << endl;

			if (!F.areEqual (c1, c2)) {
				commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Consistency failure for multiplication" << endl;
				ret = false;
			}

			F.inv (c1, a);
			F.assign (c2, a);
			F.invin (c2);

			report << "a^-1 = (out-of-place) ";
			F.write (report, c1) << ", (in-place) ";
			F.write (report, c2) << endl;

			if (!F.areEqual (c1, c2)) {
				commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Consistency failure for inversion" << endl;
				ret = false;
			}
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testArithmeticConsistency");

	return ret;
}

/* Generic test 8: Consistency of axpy
 *
 * Generates random elements 'a', 'x', and 'y' and checks that a * x + y is the
 * same for axpy, axpyin, add/mul
 */

template <class Field>
bool testAxpyConsistency (const Field &F, const char *name, unsigned int iterations)
{
	std::ostringstream str;
	str << "Testing " << name << " axpy/add-mul consistency" << ends;
	commentator.start (str.str ().c_str (), "testAxpyConsistency", iterations);

	bool ret = true;

	typename Field::RandIter r (F);
	typename Field::Element a, x, y, c1, c2, c3;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (x);
		r.random (y);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", x = ";
		F.write (report, x) << ", y = ";
		F.write (report, y) << endl;

		F.mul (c1, a, x);
		F.addin (c1, y);
		F.axpy (c2, a, x, y);
		F.assign (c3, y);
		F.axpyin (c3, a, x);

		report << "a * x + y = (add-mul) ";
		F.write (report, c1) << ", (out-of-place) ";
		F.write (report, c2) << ", (in-place) ";
		F.write (report, c3) << endl;

		if (!F.areEqual (c1, c2) || !F.areEqual (c1, c3)) {
			commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Consistency failure for axpy" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAxpyConsistency");

	return ret;
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
	LinBox::VectorDomain <Field> VD (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		for (j = 0; j < n; j++) {
			r.random (u[j]);
			r.random (v[j]);
		}

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector u:  ";
		VD.write (report, u);
		report << endl;

		report << "Input vector v:  ";
		VD.write (report, v);
		report << endl;

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

		report << "Result of Field::axpy: ";
		F.write (report, r1);
		report << endl;

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


/* Convenience function to run all of the field tests on a given field */

template <class Field>
bool runFieldTests (const Field &F, const char *desc, unsigned int iterations, size_t n, bool runCharacteristicTest) 
{
	bool pass = true;
	ostringstream str1, str2;

	str1 << "Testing " << desc << " field" << ends;
	str2 << "Testing " << desc << " FieldAXPY" << ends;

	commentator.start (str1.str ().c_str (), "runFieldTests", runCharacteristicTest ? 11 : 10);
	
	if (!testField                 (F, str1.str ().c_str ()))                pass = false; commentator.progress ();
	if (!testFieldNegation         (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldInversion        (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldAxioms           (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldAssociativity    (F, desc, iterations))                    pass = false; commentator.progress ();

	if (runCharacteristicTest) {
		if (!testFieldCharacteristic (F, desc, iterations))
			pass = false;

		commentator.progress ();
	}

	if (!testGeometricSummation    (F, desc, iterations, 100))               pass = false; commentator.progress ();
	if (!testFreshmansDream        (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testArithmeticConsistency (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testAxpyConsistency       (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldAXPY             (F, n, iterations, str2.str ().c_str ())) pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runFieldTests");

	return pass;
}


/* Generic tests for black boxes */

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
	       LinBox::VectorStream<Vector>     &stream1,
	       LinBox::VectorStream<Vector>     &stream2) 
{
	bool ret = true;

	Vector u, v, w;

	LinBox::VectorWrapper::ensureDim (u, A.rowdim ());
	LinBox::VectorWrapper::ensureDim (v, A.coldim ());
	LinBox::VectorWrapper::ensureDim (w, A.coldim ());

	LinBox::VectorDomain <Field> VD (F);
	typename Field::Element r1, r2;
	ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Blackbox transpose test [that u^T(Av) == (uA)^T v]" << std::endl;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		stream1.next (u);
		stream2.next (v);

		ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector u:            ";
		VD.write (report, u);
		report << endl;

		report << "Input vector v:            ";
		VD.write (report, v);
		report << endl;

		A.apply (w, v);

		report << "Result of apply:           ";
		VD.write (report, w);
		report << endl;

		VD.dot (r1, u, w);

		A.applyTranspose (w, u);

		report << "Result of transpose apply: ";
		VD.write (report, w);
		report << endl;

		VD.dot (r2, w, v);

		report << "<u, Av>:  ";
		F.write (report, r1);
		report << endl;

		report << "<A^T u, v>:  ";
		F.write (report, r2);
		report << endl;

		if (!F.areEqual (r1, r2)) {
			ret = false;
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;
		}

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	return ret;
}

/* Generic test 4: Linearity of black boxes
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

template <class Field, class Vector>
static bool
testLinearity (Field                              &F,
	       LinBox::BlackboxArchetype <Vector> &A,
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

	LinBox::VectorWrapper::ensureDim (x, n);
	LinBox::VectorWrapper::ensureDim (y, n);
	LinBox::VectorWrapper::ensureDim (xpay, n);
	LinBox::VectorWrapper::ensureDim (Axpay, m);
	LinBox::VectorWrapper::ensureDim (Ax, m);
	LinBox::VectorWrapper::ensureDim (Ay, m);
	LinBox::VectorWrapper::ensureDim (AxpaAy, m);

	ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Blackbox linearity test [that A.apply to (ax + y) == a A.apply to x + A.apply to y]" << std::endl;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (x);
		stream2.next (y);

		r.random (alpha);

		ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector x: ";
		VD.write (report, x);
		report << endl;

		report << "Input vector y: ";
		VD.write (report, y);
		report << endl;

		report << "Input alpha: ";
		F.write (report, alpha) << endl;

		VD.axpy (xpay, alpha, y, x);
		A.apply (Axpay, xpay);

		A.apply (Ax, x);
		A.apply (Ay, y);
		VD.axpy (AxpaAy, alpha, Ay, Ax);

		report << "   x+alpha y = ";
		VD.write (report, xpay);
		report << endl;

		report << "A(x+alpha y) = ";
		VD.write (report, Axpay);
		report << endl;

		report << "          Ax = ";
		VD.write (report, Ax);
		report << endl;

		report << "          Ay = ";
		VD.write (report, Ay);
		report << endl;

		report << " Ax+alpha Ay = ";
		VD.write (report, AxpaAy);
		report << endl;

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

/* 
/// test 5 testSmallBlackbox - equivalence to dense matrix obtained with n applies
template <class Field, class Vector, class Blackbox = LinBox::BlackboxArchetype<Vector> >
bool testSmallBlackbox(Field& F, Blackbox& A)
{
typedef std::vector<typename Field::Element> DenseVector;

DenseMatrix<?> B(F, A.rowdim(), A.coldim());
B::RowofColsIterator i = B.begin();
BasisVectorStream<?> j(A.coldim());
for(i = B.begin(); i < B.end(); ++i, ++j) A.apply(*i, j.next());
// display B in report

int iterations = 1; // could be higher if cardinality is small.
RandomDenseStream<Field, DenseVector> stream1 (F, A.rowdim(), iterations), stream2 (F, A.coldim(), iterations);
Vector y(A.rowdim()), z(A.rowdim()), x(stream1.next());
A.apply(y, x); B.apply(z, x);
// display x, y, z in report

return y==z;
//&& testLinearity(F, A, ..) && testTranspose(F, A, ..);
}
*/

/// test 6 testBlackbox - call testTranspose and testLinearity
template <class Field, class Vector>
bool testBlackbox(Field& F, LinBox::BlackboxArchetype <Vector> &A)
{
	typedef std::vector<typename Field::Element> DenseVector;

	int iterations = 1; 
	
	LinBox::commentator.start ("Testing A(ax+y) = a(Ax) + (Ay)", "testLinearity", 1);
	LinBox::RandomDenseStream<Field, DenseVector>
		stream1 (F, A.rowdim(), iterations), stream2 (F, A.coldim(), iterations);
	bool ret = testLinearity (F, A, stream1, stream2);
	LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testLinearity");

	LinBox::commentator.start ("Testing u(Av) = (uA)v", "testTranspose", 1);
	LinBox::RandomDenseStream<Field, DenseVector>
		stream3 (F, A.rowdim(), iterations), stream4 (F, A.coldim(), iterations);
	ret = ret && testTranspose (F, A, stream3, stream4); 
	LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testTranspose");

	return ret;
}

#endif // __TEST_GENERIC_H
