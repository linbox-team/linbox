/* linbox/tests/test-generic-for-quad.h
 * Copyright (C) LinBox
 *
 * Evolved from test-generic.h by Hui Wang, spring/summer 2006
 * The test-blackbox is substantially different.  
 * It is setup to compare blackbox/zo.h with blackbox/quad-matrix.h, 
 * two representations of zero-one matrices. (quad wins)
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_test_generic_quad_H
#define __LINBOX_test_generic_quad_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>

#include "linbox/util/timer.h"
#include "linbox/util/commentator.h"
#include "linbox/util/field-axpy.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/vector-domain.h"
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

/// @name Generic field tests
//@{
/** Generic test 1: Test of field operations
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
 * with an error if the field is of characteristic 0.
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


///@name Generic tests for black boxes 
//@{
/** Generic Blackbox test 1: (u^T A) v = u^T (A v).
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

	LinBox::VectorDomain <Field> VD (F);
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

/** Generic blackbox test 3: compare to a dense matrix.
 *
 * An equivalent dense matrix B is obtained with n applies.
 * Then it's behaviour is compared to A's.
 *
 * F - Field over which to perform computations
 * A - Black box of which to compute the dense representation
 */
template <class Field, class Blackbox> 
static bool 
testSmallBlackbox(Field& F, Blackbox& A)
{
	size_t m = A.rowdim(), n = A.coldim();
	typedef std::vector<typename Field::Element> Vector;

	// e for cols of identity
	typename Field::Element zero, one; 
	F.init(zero, 0); F.init(one, 1);
	Vector e(n, zero);
	Vector v(m);

	// construct dense matrix
	LinBox::DenseMatrix<Field> B(F, m, n);
	for(size_t j = 0; j < n; ++j)
	{	e[(n + j-1)%n] = zero;
		e[j] = one;
		A.apply(v, e);
		for (size_t i = 0; i < m; ++i) B.setEntry(i, j, v[i]);
	}

// really, you have to look at B to see if it is what is intended, else this is just another linearity test.

	// display B in report
	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	B.write(report);

	// compare blackbox A and dense B on random vector
	int iterations = 1; // could be higher if cardinality is small.
	LinBox::RandomDenseStream<Field, Vector> stream (F, n, iterations);
	Vector y(m), z(m), x(n);
	stream.next(x);
	A.apply(y, x); B.apply(z, x);

	// display x, y, z in report
	LinBox::VectorDomain<Field> VD(F); 
	VD.write(report, x); report << " is x" << endl;
	VD.write(report, y); report << " is y" << endl;
	VD.write(report, z); report << " is z" << endl;

	return VD.areEqual(y, z);
}

/** Generic blackbox test 4: combination of tests
 * 
 * Call testTranspose and testLinearity.
 * If large, time apply and applyTranspose.
 * if small, call testSmallBlackbox.
 */
template <class Field, class BB> 
//template <class Field, class AA, class BB> 
static bool 
testBlackbox(Field& F, BB &A)
//testBlackbox(Field& F, AA &A, BB &B)
{
	size_t smallThresh = 20; // Below it do dense matrix comparison.
	size_t largeThresh = 2000; // Above it do timing of apply and applyTr.
	typedef std::vector<typename Field::Element> DenseVector;
	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "testBlackbox on " << A.rowdim() << " by " << A.coldim() << " matrix." << endl;
	if (largeThresh < smallThresh) 
		std::cout << "telling compiler not to issue warning";

	//LinBox::VectorDomain<Field> VD (F);
	
	LinBox::commentator.setMaxDepth(-1);
	bool ret = true;
	UserTimer t;

#if 0
	size_t iterations = 1; 
	LinBox::commentator.start ("\t--Testing A(ax+y) = a(Ax) + (Ay)", 
							   "testLinearity", 1);
	typename Field::RandIter r(F);
	LinBox::RandomDenseStream<Field, DenseVector> stream1 (F, r, A.rowdim(), iterations); 
	typename Field::Element x; 
	r.random(x);
	LinBox::RandomDenseStream<Field, DenseVector> stream2 (F, r, A.coldim(), iterations); 

	ret = ret && testLinearity (F, A, stream1, stream2);
	
	LinBox::commentator.stop (MSG_STATUS (ret), 
							  (const char *) 0, "testLinearity");
	
	/*
	LinBox::commentator.start ("\t--Testing u^T(Av) = (u^T A)v", 
							   "testTranspose", 1);
	
	r.random(x);
	r.random(x);
	LinBox::RandomDenseStream<Field, DenseVector> stream3 (F, r, A.rowdim(), iterations); 
	r.random(x);
	r.random(x);
	r.random(x);
	r.random(x);
	LinBox::RandomDenseStream<Field, DenseVector> stream4 (F, r, A.coldim(), iterations); 

	ret = ret && testTranspose (F, A, stream3, stream4); 
	LinBox::commentator.stop (MSG_STATUS (ret), 
							  (const char *) 0, "testTranspose");
	*/
	
#endif
#if 1
	/* timing tests */   // I changed the order of all tests. Timing now is the first set of tests and then linearity and transpose
	{
		//{
			DenseVector x(A.coldim()), ya(A.rowdim());
			//DenseVector yb(B.rowdim());
			for(size_t i = 0; i < A.coldim(); ++i) F.init(x[i], i);
			for(size_t i = 0; i < A.rowdim(); ++i) F.init(ya[i], i);
			//for(size_t i = 0; i < A.rowdim(); ++i) F.init(yb[i], i);
			//std::cout << " -- testblackbox: correctness test --> ";
			A.apply(ya, x);
			//B.apply(yb, x);
	    //}
			//if (!VD.areEqual(ya, yb)) {std::cout << " wrong computation -- " << endl;  return ret = false;}
			//else std::cout << " correct computation -- " << endl;

		//if (A.coldim() >= largeThresh)
			{
				//LinBox::commentator.start ("\t--Timing Test (Av)","testApply", 1);
				t.clear();
				t.start();
				//std::cout << " -- before call to apply -- " << std::endl;
				for(size_t i = 0; i < 10; ++i )
					{
				        //DenseVector x(A.coldim()), y(A.rowdim());
						//for(size_t i = 0; i < A.coldim(); ++i) F.init(x[i], i);
						//for(size_t i = 0; i < A.rowdim(); ++i) F.init(y[i], i);
						A.apply(ya, x);
						//for(int j = 0; j < A.rowdim(); ++j)
						//	report << ya[j] << " ";
						//report << std::endl;
					}
				//LinBox::commentator.stop (MSG_STATUS (true), (const char *) 0, "testApply");
				t.stop();
				std::cout << t.time() << ", ";
				//std::cout << " -- testblackbox: time of quad -- " << t.time() << std::endl;
				//t.clear();
				//t.start();
				//for(size_t i = 0; i < 3; ++i )
				//B.apply(yb, x);
				//t.stop();
				//std::cout << t.time() << ")" << endl;
				//std::cout << " -- testblackbox: time of zo  -- " << t.time() << std::endl;
			}

		//if (A.rowdim() >= largeThresh)
			{
				//LinBox::commentator.start ("\t--Timing Test(v^T A)", "testApplyTranspose", 1);
				//A.applyTranspose(x, y);
	  			//LinBox::commentator.stop (MSG_STATUS (true), (const char *) 0, "testApplyTranspose");
	   		}
	   	
	} // timing test block
#endif
#if 0	
	/*  Testing against constructed dense matrix doesn't really add much.  
	    May be useful to see the matrix in the report.
	*/
	if (A.rowdim() <= smallThresh && A.coldim() <= smallThresh)
		{
			LinBox::commentator.start ("\t--Testing A behaves like Dense A", 
									   "testSmallBlackbox", 1);
			ret = ret && testSmallBlackbox(F, A);
			LinBox::commentator.stop (MSG_STATUS (ret), 
									  (const char *) 0, "testSmallBlackbox");
		}
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
			 const char *, // text
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

		iter.random (x);

		categories1[x % num_categories]++;
		categories2[(unsigned int) (double (x) / double (card) * num_categories)]++;

		typename std::list<typename Field::Element>::iterator x_queue_iter = x_queue.begin ();
		diff_cat_iter = diff_categories.begin ();

		for (; x_queue_iter != x_queue.end (); ++x_queue_iter, ++diff_cat_iter) {
			F.sub (d, *x_queue_iter, x);
			(*diff_cat_iter)[d % num_categories]++;
		}

		x_queue.push_front (x);

		while (x_queue.size () > hist_len)
			x_queue.pop_back ();
	}

	double p, chi_squared = 0.0;

	for (i = 0; i < num_categories; ++i)
		chi_squared += pow (double (categories1[i]) -
				    double (num_trials) / double (num_categories), 2);

	p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

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
/// @name Vector operation tests
//@{

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
		report << "Input vector 1:  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2:  ";
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
		report << "Input vector 1:  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2:  ";
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
		report << "Input vector 1:  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2:  ";
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
		report << "Input vector 1:  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2:  ";
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
//@}

#endif // __LINBOX_test_generic_quad_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
