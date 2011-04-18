/* linbox/tests/test-field.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 * See COPYING for license information.
 *
 * Extracted by bds from test-generic.h, written by Bradford Hovinen <hovinen@cis.udel.edu>
 */

#ifndef __LINBOX_test_field_H
#define __LINBOX_test_field_H

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
typename Field::Element& expt (const Field &F, typename Field::Element &res, const typename Field::Element &a, LinBox::integer &n) 
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
                typename Field::Element tmp;
                F.init(tmp,0);
		res = F.mul (tmp, res, res);
	}

	return res;
}

bool reportError(string rep, bool& flag)
{
	ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
	report << "ERROR: " << rep << endl;
	return flag = false;
}

/** Check each field or ring operation.
 *
 * Test various field operations
 *
 * F - Field over which to perform computations
 * title - String to use as the descriptive title of this test
 * fieldp - use true if inv and div must work for all nonzero denominators
 *
 * Return true on success and false on failure
 */

template<class Field>
bool testField (Field &F, const char *title, bool fieldp = true) 
{
	commentator.start (title, "testField", 5);
	ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	typename Field::Element zero, one, two, three;
	F.init(zero, 0); F.init(one, 1); F.init(two, 2); F.init(three, 3);

	typename Field::Element a, b, c, d, e, f;
	F.init(a,0); F.init(b,0); F.init(c,0); F.init(d,0); F.init(e,0); F.init(f,0);


	report << "Field self description: " << F.write (report) << endl;
	//	report << "field Element 2: " << F.write (report, two) << endl;

	LinBox::integer n, m;
	bool pass = true, part_pass = true;

	commentator.start ("\t--Testing characteristic/cardinality match");

	F.characteristic (n); 
	F.cardinality (m);

	if (n > 0 && !isPower (m, n)) part_pass = reportError("Characteristic, cardinality mismatch", pass);

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	/* tests for presence of members with minimal check of semantics */
	// these checks need improvement 

	commentator.start ("\t--Testing correctness of 0 and 1");
	part_pass = true;

	if (!F.isZero (zero)) {
		part_pass = reportError( "isZero (0) is false", pass);
	}
	if (F.isZero (one)) part_pass = reportError( "isZero (1) is true", pass);
	if (F.isOne (zero)) part_pass = reportError( "isOne (0) is true", pass);
	if (!F.isOne (one)) {
		part_pass = reportError( "isOne (1) is false", pass);
	}


	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("\t--Testing init/convert");
	part_pass = true;

	if (F.cardinality (m) <= 0)
		n = 49193295;   // Just using some odd value
	else
		n -= 1;

	
	report << "Initial integer: " << n << endl;
	F.init (a, n);  F.write ( report << "Result of init: ", a) << endl;

	F.convert (m, a);
	report << "Result of convert: " << m << endl;

	if (m != n) part_pass = reportError( "F.convert (m, F.init (a, n)) != n", pass);

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("\t--Testing field arithmetic");
	part_pass = true;

	F.init (b, n-2);
	F.init (d, n-2);
	F.init (e, 3);

	F.add (a, three, two); F.write (report << "Result of 2 + 3: ", a) << endl;
	F.assign (d, three);
	F.addin (d, two); F.write (report << "Result of 2 + 3 (inplace): ", d) << endl;

	if (!F.areEqual (a, F.init (f, 5)) || !F.areEqual (d, a)) 
		part_pass = reportError( "Results of add are incorrect", pass);

	F.neg (a, two); F.write (report << "Result of -2: ", a) << endl;
	F.assign (d, two);
	F.negin (d); F.write (report << "Result of -2 (inplace): ", d) << endl;
	F.init (f, -2); F.write( report << "-2 via init: ", f) << endl;

	if (!F.areEqual (a, f) || !F.areEqual (d, a)) 
		part_pass = reportError( "Results of neg are incorrect", pass);

	F.sub (a, three, two); F.write (report << "Result of 3 - 2: ", a) << endl;
	F.init (d, 3);
	F.subin (d, two); F.write (report << "Result of 3 - 2 (inplace): ", d) << endl;

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) 
		part_pass = reportError( "Results of neg sub incorrect", pass);

	F.mul (a, two, three); F.write (report << "Result of 2 * 3: ", a) << endl;
	F.assign (d, two);
	F.mulin (d, three); F.write (report << "Result of 2 * 3 (inplace): ", d) << endl;
	F.init (f, 6);

	if (!F.areEqual (a, f) || !F.areEqual (d, a)) 
		part_pass = reportError( "Results of mul incorrect", pass);

	F.inv (a, one); F.write (report << "Result of inverting 1: ", a) << endl;
	F.assign (d, one);
	F.invin (d); F.write (report << "Result of inverting 1 (inplace): ", d) << endl;

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) 
		part_pass = reportError( "Results of inv incorrect", pass);

	if ( F.isZero(two) ) F.assign(f, three); else F.assign(f, two);
	F.div (a, f, f); F.write ( report << "Result of f/f: ", a) << endl;
	F.assign(d, f);
	F.divin (d, f); F.write ( report << "Result of f/f (inplace): ", d) << endl;

	if (!F.areEqual (a, one) || !F.areEqual (d, a)) 
		part_pass = reportError( "Results of div incorrect", pass);

	F.axpy (a, two, three, two); F.write ( report << "Result of axpy 2*3 + 2: ", a) << endl;
	F.assign (d, two);
	F.axpyin (d, two, three); F.write ( report << "Result of axpy 2*3 + 2 (inplace): ", d) << endl;
	F.init (f, 8);

	if ( !F.areEqual (a, f) || !F.areEqual (d, a) ) 
		part_pass = reportError( "Results of axpy incorrect", pass);

	commentator.stop (MSG_STATUS (part_pass));
	commentator.progress ();

	commentator.start ("\t--Testing summation of powers of 2");

	//,..
	// 2^101 - 1 vs 1 + 2 + 4 + ... + 2^100

	F.init (a, 1);
	F.init (b, 1);
	F.init (c, 0);
	
	n = 101;
	expt(F, a, two, n);
	F.subin (a, one); F.write( report << "using expt, 2^101 - 1: ", a) << endl;

	for (int i = 1; i <= 101; ++i) {
		F.addin (c, b);
		F.mulin (b, two);
	}

	if (!F.areEqual (a, c)) 
		part_pass = reportError( "2^101 - 1 != 1 + 2 + .. + 2^100", pass);

	// 1 + 2*(1 + 2*( ... (1) ... )), 100 times.
	F.assign (d, one);
	for (int i = 1; i < 101; ++i)
	{	F.axpy (b, two, d, one); F.assign(d, b); }
    F.write( report << "using axpy, 1 + 2(1 + ... + 2(1)...), with 100 '+'s: ", d) << endl;

	if (!F.areEqual (a, d)) 
	part_pass = reportError( "2^101 - 1 != 1 + 2(1 + ... + 2(1)...), with 100 '+'s: ", pass);

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
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testFieldNegation", iterations);

	typename Field::Element a, neg_a, neg_a_a, zero;
	F.init(a,0); F.init(neg_a,0); F.init(neg_a_a,0); F.init (zero, 0);
	typename Field::RandIter r (F);

	bool ret = true;

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

		if (!F.areEqual (neg_a_a, zero)) reportError("a + -a != 0" , ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldNegation");
	delete[] st;
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
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testFieldInversion", iterations);

	typename Field::Element a, ainv, aainv, one;
        F.init (a,0); F.init (ainv,0); F.init (aainv,0);
	F.init (one, 1);
	typename Field::RandIter r (F);

	bool ret = true;


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

		if (!F.areEqual (aainv, one)) reportError("a a^-1 != 1", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldInversion");
	delete[] st;
	return ret;
}

/** @brief Generic test 7a: Distributivity of multiplication over addition

 * Given random field elements 'a', 'b', and 'c', checks that
 * (a + b) * c = a * c + b * c  and  c * (a + b) = c * a + c * b
 */


template <class Field>
bool testFieldDistributivity(const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " distributivity" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testFieldDistributivity", iterations);

	typename Field::Element a, b, c, a_b, a_bc, ac, bc, ac_bc, ca_b, ca, cb, ca_cb;
        F.init (a,0); F.init (b,0); F.init (c,0); 
        F.init (a_b,0); F.init (a_bc,0); F.init (ac,0); F.init (bc,0);
        F.init (ac_bc,0); 
		F.init (ca_b,0); F.init (ca,0); F.init (cb,0); F.init (ca_cb,0); 
        
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

		F.add (a_b, a, b);//a + b
		F.mul (a_bc, a_b, c);//(a+b)*c
		F.mul (ca_b, c, a_b);//c*(a+b)
		F.mul (ac, a, c); //a*c
		F.mul (bc, b, c); //b*c
		F.mul (ca, c, a); //c*a
		F.mul (cb, c, b); //c*b
		F.add (ac_bc, ac, bc); //a*c + b*c
		F.add (ca_cb, ca, cb); //c*a + c*b

		report << "(a + b) * c = ";
		F.write (report, a_bc) << endl;

		report << "a * c + b * c = ";
		F.write (report, ac_bc) << endl;

		report << "c * (a + b) = ";
		F.write (report, ca_b) << endl;

		report << "c * a + c * b = ";
		F.write (report, ca_cb) << endl;
		if (!F.areEqual (a_bc, ac_bc) || !F.areEqual (ca_b, ca_cb)) 
			reportError("Operations were not distributative", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldDistributivity");
	delete[] st;
	return ret;
}


/** @brief Generic test 7b: Commutativity of multiplication and addition

 * Given random field elements 'a', 'b', checks that
 * a*b = b*a
 * a+b = b+a
 */


template <class Field>
bool testFieldCommutativity (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " commutativity," << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testFieldCommutativity", iterations);

	typename Field::Element a, b, ab, ba, a_b, b_a;
        F.init (a,0); F.init (b,0);
        F.init (ab,0); F.init (ba,0);
        F.init (a_b,0); F.init (b_a,0);

        
	typename Field::RandIter r (F);

	bool ret = true;

	for (unsigned int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		r.random (a);
		r.random (b);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << endl;

		F.mul (ab, a, b);
		F.mul (ba, b, a);

		report << "a*b = ";
		F.write (report, ab) << endl;

		report << "b*a = ";
		F.write (report, ba) << endl;

		if (!F.areEqual (ab, ba)) 
			reportError("Multiplication was not commutative", ret);


		F.add (a_b, a, b);
		F.add (b_a, b, a);

		report << "a+b = ";
		F.write (report, a_b) << endl;

		report << "b+a = ";
		F.write (report, b_a) << endl;

		if (!F.areEqual (a_b, b_a)) 
			reportError("Addition was not commutative", ret);



		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldCommutativity");
	delete[] st;
	return ret;
}


/** Generic test 7c: Associativity of addition and multiplication
 *
 * Given random field elements 'a', 'b', and 'c', checks that
 * (a * b) * c = a * (b * c) and (a + b) + c = a + (b + c)
 */

template <class Field>
bool testFieldAssociativity (const Field &F, const char *name, unsigned int iterations) 
{
	std::ostringstream str;
	str << "\t--Testing " << name << " associativity" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testFieldAssociativity", iterations);

	typename Field::Element a, b, c, a_b, b_c, a_bc, ab_c;
        F.init (a,0); F.init (b,0); F.init (c,0);
        F.init (a_b,0); F.init (b_c,0); F.init (a_bc,0); F.init (ab_c,0);
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

		if (!F.areEqual (ab_c, a_bc)) reportError( "Results are not equal", ret);

		F.mul (a_b, a, b);
		F.mul (ab_c, a_b, c);
		F.mul (b_c, b, c);
		F.mul (a_bc, a, b_c);

		report << "(a * b) * c = ";
		F.write (report, ab_c) << endl;

		report << "a * (b * c) = ";
		F.write (report, a_bc) << endl;

		if (!F.areEqual (ab_c, a_bc)) reportError( "Results are not equal", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldAssociativity");
	delete[] st;
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
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testGeometricSummation", iterations);

	typename Field::Element a, a_n, k, zero, one;
	typename Field::RandIter r (F);

	F.init (zero, 0);
	F.init (one, 1);
        F.init (a,0); F.init (a_n,0); F.init (k,0);

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

		if (!F.areEqual (k, a_n)) reportError("Field elements are not equal", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testGeometricSummation");
	delete[] st;
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
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (string(str.str()).c_str(), "testFieldCharacteristic", iterations);

	LinBox::integer p, j;
	typename Field::Element a, sigma, zero;
	typename Field::RandIter r (F);

	F.characteristic (p);
	F.init (zero, 0);
        F.init (a,0); F.init (sigma,0);
        
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

// suggestion: make this run in time O(lg(p)), then take the conditional out of the run...Tests
		for (j = 0; j < p; j += 1)
			F.addin (sigma, a);

		report << "p a = ";
		F.write (report, sigma) << endl;

		if (!F.isZero (sigma)) reportError("p a != 0", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFieldCharacteristic");
	delete[] st;
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
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testFreshmansDream", iterations);

	LinBox::integer c, j;

	F.characteristic (c);

	if (c == 0) {
		commentator.stop ("skipping", "Field characteristic is 0, so this test makes no sense", "testFreshmansDream");
		delete[] st;
		return true;
	}

	bool ret = true;

	typename Field::RandIter r (F);
	typename Field::Element a, b, a_b, a_b_p, a_p, b_p, a_p_b_p;

        F.init (a,0); F.init (b,0); F.init (a_b,0);
        F.init (a_b_p,0); F.init (a_p,0); F.init (b_p,0); F.init (a_p_b_p,0);

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

		if (!F.areEqual (a_b_p, a_p_b_p)) reportError("(a + b)^p != a^p + b^p", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testFreshmansDream");
	delete[]  st;
	return ret;
}


/* Tests of field features */ 

/** Generic test 7: Consistency of in-place and out-of-place arithmetic
 *
 * Generates random elements 'a' and 'b' and performs all basic arithmetic
 * operations in-place and out-of-place, checking for consistency.
 *
 * Div and inv are checked in a separate function.
 */

template <class Field>
bool testArithmeticConsistency (const Field &F, const char *name, unsigned int iterations)
{  return testRingArithmeticConsistency(F, name, iterations) 
	   && testInvDivConsistency(F, name, iterations); 
}
template <class Field>
bool testRingArithmeticConsistency (const Field &F, const char *name, unsigned int iterations)
{
	std::ostringstream str;
	str << "\t--Testing " << name << " in-place/out-of-place arithmetic consistency" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRingArithmeticConsistency", iterations);

	bool ret = true;

	typename Field::RandIter r (F);
	typename Field::Element a, b, c1, c2;
        F.init (a,0); F.init (b,0); F.init (c1,0); F.init (c2,0);

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

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for addition", ret);

		F.sub (c1, a, b);
		F.assign (c2, a);
		F.subin (c2, b);

		report << "a - b = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for subtraction", ret);

		F.neg (c1, a);
		F.assign (c2, a);
		F.negin (c2);

		report << "-a = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for negation", ret);

		F.mul (c1, a, b);
		F.assign (c2, a);
		F.mulin (c2, b);

		report << "a * b = (out-of-place) ";
		F.write (report, c1) << ", (in-place) ";
		F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for multiplication", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRingArithmeticConsistency");
	delete[] st;
	return ret;
}

template <class Field>
bool testInvDivConsistency (const Field &F, const char *name, unsigned int iterations)
{
    std::ostringstream str;
    str << "\t--Testing " << name << " in-place/out-of-place inv and div consistency" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testInvDivConsistency", iterations);

    bool ret = true;

    typename Field::RandIter r (F);
    typename Field::Element a, b, c1, c2;
    F.init (a,0); F.init (b,0); F.init (c1,0); F.init (c2,0);

    for (unsigned int i = 0; i < iterations; i++) {
        commentator.startIteration (i);

        r.random (a); r.random (b);

		// This should be immaterial, since we have already "proven" commutativity
		if (F.isZero (a) && !F.isZero (b))
		std::swap (a, b);

		ostream &report = commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Random elements a = ";
		F.write (report, a) << ", b = ";
		F.write (report, b) << endl;


		if (!F.isZero (a)) {
			F.div (c1, b, a);
			F.assign (c2, b);
			F.divin (c2, a);

			report << "b / a = (out-of-place) ";
			F.write (report, c1) << ", (in-place) ";
			F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for division", ret);

			F.inv (c1, a);
			F.assign (c2, a);
			F.invin (c2);

			report << "a^-1 = (out-of-place) ";
			F.write (report, c1) << ", (in-place) ";
			F.write (report, c2) << endl;

		if (!F.areEqual (c1, c2)) reportError("Consistency failure for inversion", ret);
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testInvDivConsistency");
	delete[] st;
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
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testAxpyConsistency", iterations);

	bool ret = true;

	typename Field::RandIter r (F);
	typename Field::Element a, x, y, c1, c2, c3;
        F.init (a,0); F.init (x,0); F.init (y,0); 
        F.init (c1,0); F.init (c2,0); F.init (c3,0);

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

		if (!F.areEqual (c1, c2) || !F.areEqual (c1, c3)) reportError("Consistency failure for axpy", ret);

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testAxpyConsistency");
	delete[] st;
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
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());
	commentator.start (st, "testRanditerBasic", iterations);

	 typename Field::RandIter r (F);
	 typename Field::Element a;
		 F.init (a,0);

	 if (iterations < 20) iterations = 20;
	 for (unsigned int i = 0; i < iterations; i++) {
		 r.random (a);
		 if ( ! F.isZero(a) ) {ret = true; break;}

	 }

	 commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRanditerBasic");
	 delete[] st;
	 return ret;
 }

 
 /* Convenience function to run all of the field tests on a given field */

 template <class Field>
 bool runFieldTests (const Field &F, const char *desc, unsigned int iterations, size_t n, bool runCharacteristicTest = true) 
 // n is not used.
 {	ostringstream str;

	 str << "\t--Testing " << desc << " field" << ends;
	 char * st = new char[str.str().size()];
	 strcpy (st, str.str().c_str());
	 commentator.start (st, "runFieldTests");
	 bool ret =  runBasicRingTests(F, desc, iterations, runCharacteristicTest)
		 && testInvDivConsistency(F, desc, iterations) 
		 && testFieldInversion (F, desc, iterations)
		 && testFieldCommutativity (F, desc, iterations)
		 && testFreshmansDream(F, desc, iterations);

	 commentator.stop (MSG_STATUS (ret));
	 delete[] st;
	 return ret;
}

template <class Field>
bool runBasicRingTests (const Field &F, const char *desc, unsigned int iterations, bool runCharacteristicTest = true) 
{
	bool pass = true;
	ostringstream str;

	str << "\t--Testing " << desc << " ring" << ends;
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());

	commentator.start (st, "runBasicRingTests", runCharacteristicTest ? 11 : 10);
	
	if (!testField                 (F, string(str.str()).c_str()))                pass = false; commentator.progress ();
	if (!testFieldNegation         (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldDistributivity           (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldAssociativity    (F, desc, iterations))                    pass = false; commentator.progress ();

	if (runCharacteristicTest) {
		if (!testFieldCharacteristic (F, desc, iterations))
			pass = false;

		commentator.progress ();
	}

	if (!testGeometricSummation    (F, desc, iterations, 100))               pass = false; commentator.progress ();
	if (!testRingArithmeticConsistency (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testAxpyConsistency       (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testRanditerBasic       (F, desc, iterations))                    pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runBasicRingTests");
	delete[] st;
	return pass;
}

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
	char * st = new char[str.str().size()];
	strcpy (st, str.str().c_str());

	LinBox::commentator.start (st, "testRandomIterator");

	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	/* This test either passes or runs a lot of times */
	for (int i = 1; 
	     (!  testRandomIteratorStep (F, text, num_trials, num_categories, hist_len)) && (i < 20) ;
	     ++i ){
		if (0 == i % 10)  
			report << "Warning! Probable failure of uniformity" << std::endl;
		};

	LinBox::commentator.stop (MSG_STATUS (true), (const char *) 0, "testRandomIterator");

	delete[] st;
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
	typename Field::Element x,  d;

	std::list<std::vector<int> >::iterator diff_cat_iter;

	for (i = 0; i < hist_len; ++i)
		diff_categories.push_back (std::vector<int> (num_categories, 0));

	// I make the simplifying assumption that field elements are actually
	// C++ ints. Otherwise, I don't know how to place the numbers into
	// categories in any well-defined manner.
	for (i = 0; i < num_trials; ++i) {
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

	if (p < 0.05 || p > 0.95) 
		reportError("Random iterator's values do not appear to be uniformly distributed", ret);

	chi_squared = 0.0;

	for (i = 0; i < num_categories; ++i)
		chi_squared += pow (double (categories2[i]) -
				    double (num_trials) / double (num_categories), 2);

	p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

	report << "Test of distribution uniformity (high-order): chi^2 = "
	       << chi_squared * num_categories / num_trials << std::endl;
	report << "Test of distribution uniformity (high-order):     p = " << p << std::endl;

	if (p < 0.05 || p > 0.95) 
		reportError("Consistency failure for addition", ret);

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

		if (p < 0.05 || p > 0.95) 
			reportError("Difference values do not appear to be uniformly distributed", ret);
	}

	//LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomIteratorStep");
	return ret;
}
//@}
#endif // __LINBOX_test_field_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
