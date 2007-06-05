/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */

/* linbox/tests/test-field.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 * See COPYING for license information.
 *
 * Extracted by bds from test-generic.h, written by Bradford Hovinen <hovinen@cis.udel.edu>
 */

#ifndef __TEST_FIELD_H
#define __TEST_FIELD_H

#include <iostream>
//#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/util/field-axpy.h"
//#include "linbox/vector/stream.h"
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

/// @name Generic field tests
//@{
/** Generic test 1: Test of field operations
 *
 * Test various field operations
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

	ostream &report = commentator.report (LinBox::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Field self description: " << F.write (report) << endl;
	//	report << "field Element 2: " << F.write (report, two) << endl;

	LinBox::integer n, m;
	bool pass = true, part_pass;

	commentator.start ("\t--Testing characteristic/cardinality match");
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

	commentator.start ("\t--Testing correctness of 0 and 1");
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

	commentator.start ("\t--Testing init/convert");
	part_pass = true;

	if (F.cardinality (m) <= 0)
		n = 49193295;   // Just using some odd value
	else
		n -= 1;

	F.init (a, n);  
	
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

	commentator.start ("\t--Testing field arithmetic");
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
		F.write( F.write( F.write( 
                    commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                         << "ERROR: Results of inv are incorrect : 1/", one)
                         << " != ", a) 
                         << " or != ", d) << endl;
	}

        F.write( commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << " two : ", two) << ", is zero ? : " << F.isZero(two) << std::endl;
        
	if ( ! F.isZero(two) )
	{
	        F.div (a, two, two);
		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Result of 2/2: ";
		F.write (report, a);
		report << endl;
	        if (!F.areEqual (a, one) ) {
		    pass = part_pass = false;
		    report << "ERROR: Result of div is incorrect" << endl;
	        }
	}

        F.write( commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << " three : ", three) << ", is zero ? : " << F.isZero(three) << std::endl;
	if ( ! F.isZero(three) ) {
		F.assign (d, three);
		F.divin (d, three);
		ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Result of 3/3: ";
		F.write (report, d);
		report << endl;
		if (!F.areEqual (d, one)) {
		    pass = part_pass = false;
		    report << "ERROR: Result of divin is incorrect" << endl;
		}
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

	commentator.start ("\t--Testing summation of powers of 2");

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


/** Tests of algebraic properties of fields */

/* Generic test 6: Negation of elements
 *
 * Negates random elements and checks that they are true negatives
 */

template <class Field>
bool testFieldNegation (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " negation" << ends;
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

/** Generic test 5: Inversion of elements
 *
 * Inverts random elements and checks that they are true inverses
 */

template <class Field>
bool testFieldInversion (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " inversion" << ends;
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

		report << "a^{-1} = ";  F.write (report, ainv) << endl;

		F.mul (aainv, ainv, a);

		report << "a a^{-1} = ";  F.write (report, aainv) << endl;

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

/** @brief Generic test 7: Commutativity and distributivity of addition and multiplication

 * Given random field elements 'a', 'b', and 'c', checks that
 * (a + b) * c = a * c + b * c = c * (a + b) = b * c + a * c
 */

template <class Field>
bool testFieldAxioms (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " commutativity, distributivity" << ends;
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

/** Generic test 7: Associativity of addition and multiplication
 *
 * Given random field elements 'a', 'b', and 'c', checks that
 * (a * b) * c = a * (b * c) and (a + b) + c = a + (b + c)
 */

template <class Field>
bool testFieldAssociativity (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " associativity" << ends;
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

/** Generic test 2: Geometric summation
 *
 * Generates a random field element 'a' and raises it through repeated
 * exponentiation to the power n. Takes the sum k of all intermediate values and
 * checks that a^n = (k-1)/(a-1).
 */

template <class Field>
bool testGeometricSummation (const Field &F, const char *name, unsigned int iterations, unsigned int n) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " geometric summation" << ends;
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

		report << "n = " << n << " a^n = ";
		F. write (report, a_n) << endl;
		F. write(report);
		report<<std::endl;

		report << "sum(a^i, i = 0..n-1) = ";
		F.write (report, k) << endl;

		F.subin (a_n, one);
		report << "(a^n - 1) = ";
		F.write (report, a_n) << endl;

		F.subin (a, one);
		report << "(a - 1) = ";
		F.write (report, a) << endl;

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

/** Generic test 3: Test of field characteristic
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
	str << "\t--Testing " << name << " characteristic" << ends;
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

/** Generic test 4: The Freshman's Dream
 *
 * Generates two random field elements 'a' and 'b', and checks whether
 * (a + b)^p = a^p + b^p, where p is the characteristic of the field. Bails out
 * (returning true) if the field is of characteristic 0.
 */

template <class Field>
bool testFreshmansDream (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " Freshman's Dream" << ends;
	commentator.start (str.str ().c_str (), "testFreshmansDream", iterations);

	LinBox::integer c, j;
	typename Field::Element cp; F.init(cp, c);

	F.characteristic (c);

	if (F.isZero (cp)) {
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

/** Generic test 7: Consistency of in-place and out-of-place arithmetic
 *
 * Generates random elements 'a' and 'b' and performs all basic arithmetic
 * operations in-place and out-of-place, checking for consistency
 */

template <class Field>
bool testArithmeticConsistency (const Field &F, const char *name, unsigned int iterations)
{
	std::ostringstream str;
	str << "\t--Testing " << name << " in-place/out-of-place arithmetic consistency" << ends;
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

			report << "b / a = (out-of-place) ";
			F.write (report, c1) << ", (in-place) ";
			F.write (report, c2) << endl;

			if (!F.areEqual (c1, c2)) {
				commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Consistency failure for division" << endl;
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

/** Generic test 8: Consistency of axpy
 *
 * Generates random elements 'a', 'x', and 'y' and checks that a * x + y is the
 * same for axpy, axpyin, add/mul
 */

template <class Field>
bool testAxpyConsistency (const Field &F, const char *name, unsigned int iterations)
{
	std::ostringstream str;
	str << "\t--Testing " << name << " axpy/add-mul consistency" << ends;
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

/** Generic test 9: Basic concept check of RandIter
 *
 * In a loop, generates random element 'a', and fails
 * if it is always zero.
 */
template <class Field>
bool testRanditerBasic(const Field &F, const char *name, unsigned int iterations)
{
	bool ret=false;
	std::ostringstream str;
	str << "\t--Testing " << name << " randiter basic operation " << ends;
	commentator.start (str.str ().c_str (), "testAxpyConsistency", iterations);

	typename Field::RandIter r (F);
	typename Field::Element a;

	for (unsigned int i = 0; i < iterations; i++) {
		r.random (a);
		if ( ! F.isZero(a) ) {ret = true; break;}

	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRanditerBasic");

	return ret;
}


/* Convenience function to run all of the field tests on a given field */

template <class Field>
bool runFieldTests (const Field &F, const char *desc, unsigned int iterations, size_t n, bool runCharacteristicTest) 
{
	bool pass = true;
	ostringstream str1, str2;

	str1 << "\t--Testing " << desc << " field" << ends;
	str2 << "\t--Testing " << desc << " FieldAXPY" << ends;

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
	if (!testRanditerBasic       (F, desc, iterations))                    pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runFieldTests");

	return pass;
}
//@}

/// @name Generic field tests
//@{
/** Random number test
 *
 * Test that the random iterator over the given field works.
 *
 * Test up to five times, accepting any one, to increase probability of 
 * passing statistical tests.
 */
template <class Field>
bool testRandomIterator (const Field &F, const char *text,
			 unsigned int num_trials,
			 unsigned int num_categories,
			 unsigned int hist_len) 
{
	std::ostringstream str;

	str << "\t--Testing " << text << "::RandIter" << std::ends;

	LinBox::commentator.start (str.str ().c_str (), "testRandomIterator");

	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	/* This test either passes or runs forever */
	for (int i = 1; 
	     !  testRandomIteratorStep (F, text, num_trials, num_categories, hist_len) ;
	     ++i ){
		if (0 == i % 10)  
			report << "Warning! Probable failure of uniformity" << std::endl;
		};

	LinBox::commentator.stop (MSG_STATUS (true), (const char *) 0, "testRandomIterator");
	return true;

}

/* Random number test
 *
 * Test that the random iterator over the given field works
 */

template <class Field>
bool testRandomIteratorStep (const Field &F,
			 const char *text,
			 unsigned int num_trials,
			 unsigned int num_categories,
			 unsigned int hist_len) 
{
	//std::ostringstream str;

	//str << "\t--Testing " << text << "::RandIter" << std::ends;

	//LinBox::commentator.start (str.str ().c_str (), "testRandomIteratorStep");
	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	bool ret = true;

	LinBox::integer card;
	unsigned int i;
	std::vector<int> categories1 (num_categories, 0);
	std::vector<int> categories2 (num_categories, 0);
	std::list<std::vector<int> > diff_categories;
	std::list<typename Field::Element> x_queue;

	F.cardinality (card);

	typename Field::RandIter iter (F);
	typename Field::Element x, x_prev, x_prev2, d;

	std::list<std::vector<int> >::iterator diff_cat_iter;

	for (i = 0; i < hist_len; ++i)
		diff_categories.push_back (std::vector<int> (num_categories, 0));

	// I make the simplifying assumption that field elements are actually
	// C++ ints. Otherwise, I don't know how to place the numbers into
	// categories in any well-defined manner.
	for (i = 0; i < num_trials; ++i) {
		F.assign (x_prev2, x_prev);
		F.assign (x_prev, x);

                integer ix, id;
		F.convert(ix, iter.random (x));

   
		categories1[ix % num_categories]++;
		categories2[(unsigned int) (double (ix) / double (card) * num_categories)]++;

		typename std::list<typename Field::Element>::iterator x_queue_iter = x_queue.begin ();
		diff_cat_iter = diff_categories.begin ();

		for (; x_queue_iter != x_queue.end (); ++x_queue_iter, ++diff_cat_iter) {
			F.convert(id, F.sub (d, *x_queue_iter, x));
			(*diff_cat_iter)[id % num_categories]++;
		}

		x_queue.push_front (x);

		while (x_queue.size () > hist_len)
			x_queue.pop_back ();
	}

	double p, chi_squared = 0.0;

	for (i = 0; i < num_categories; ++i)
		chi_squared += pow (double (categories1[i]) -
				    double (num_trials) / double (num_categories), 2);

	p = chiSquaredCDF (chi_squared * (double)num_categories / (double)num_trials, (double)num_categories - 1.0);

	report << "Test of distribution uniformity (low-order): chi^2 = "
	       << chi_squared * num_categories / num_trials << std::endl;
	report << "Test of distribution uniformity (low-order):     p = " << p << std::endl;

	if (p < 0.05 || p > 0.95) {
		LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Random iterator's values do not appear to be uniformly distributed"
			<< std::endl;
		ret = false;
	}

	chi_squared = 0.0;

	for (i = 0; i < num_categories; ++i)
		chi_squared += pow (double (categories2[i]) -
				    double (num_trials) / double (num_categories), 2);

	p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

	report << "Test of distribution uniformity (high-order): chi^2 = "
	       << chi_squared * num_categories / num_trials << std::endl;
	report << "Test of distribution uniformity (high-order):     p = " << p << std::endl;

	if (p < 0.05 || p > 0.95) {
		LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Random iterator's values do not appear to be uniformly distributed"
			<< std::endl;
		ret = false;
	}

	diff_cat_iter = diff_categories.begin ();

	int idx = 0;

	for (; diff_cat_iter != diff_categories.end (); ++diff_cat_iter, ++idx) {
		chi_squared = 0.0;

		for (i = 0; i < num_categories; ++i)
			chi_squared += pow (double ((*diff_cat_iter)[i]) -
					    double (num_trials) / double (num_categories), 2);

		p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

		report << "Test of " << idx + 1 << " spacing: chi^2 = "
		       << chi_squared * num_categories / num_trials << std::endl;
		report << "Test of " << idx + 1 << " spacing:     p = " << p << std::endl;

		if (p < 0.05 || p > 0.95) {
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Difference values do not appear to be uniformly distributed"
				<< std::endl;
			ret = false;
		}
	}

	//LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomIteratorStep");
	return ret;
}
//@}
#endif // __TEST_FIELD_H
