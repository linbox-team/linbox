/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-scalar-matrix.C
 * Evolved from test-diagonal.C  -bds
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

//#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/solutions/minpoly.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Application of identity matrix onto random vectors
 *
 * Construct the identity matrix and a series of randomly-generated
 * vectors. Apply the identity to each vector and test whether the input and
 * output are equal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testScalarApply (Field &F, size_t n, ostream &report) 
{
	typedef typename Field::Element Element;
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef ScalarMatrix <Field, Vector> Blackbox;

	cout << "Testing scalar matrix apply with 1, 0, random ...";
	cout.flush ();
	report << "Testing identity apply:" << endl;

	bool ret = true;

	Vector u(n), v(n), w(n);
	typename Field::RandIter r (F);
	for (size_t j = 0; j < n; j++) r.random (v[j]);

	Element d;
	/*  indentity matrix */
	F.init (d, 1);
	Blackbox I (F, n, d);

	if (I.rowdim() != n || n != I.coldim())
	  { ret = false; report << "    Bad dimension(s), should be " << n << endl; }

	report << "    Input vector:  ";
	printVector<Field> (F, report, v);

	I.apply (w, v);

	report << "    Output vector (should equal input): ";
	printVector<Field> (F, report, w);

	for (size_t j = 0; j < n; j++)
		if (!F.areEqual (w[j], v[j])) {ret = false; break;}

	/*  square zero matrix */
	F.init (d, 0);
	Blackbox Z (F, n, d);

	report << "    Input vector:  ";
	printVector<Field> (F, report, v);

	Z.apply (w, v);

	report << "    Output vector (should be zero): ";
	printVector<Field> (F, report, w);

	for (size_t j = 0; j < n; j++)
		if (!F.isZero(w[j]) ) {ret = false; break;}

	/* random scalar matrix */
	do {r.random (d);} while ( F.isZero(d) );
	Blackbox R (F, n, d);

	report << "    Input vector:  ";
	printVector<Field> (F, report, v);

	R.apply (w, v);

	report << "    apply Output vector (random): ";
	printVector<Field> (F, report, w);

	R.applyTranspose (u, v);

	report << "    applyTranspose output vector (should equal apply output): ";
	printVector<Field> (F, report, u);

	for (size_t j = 0; j < n; j++)
		if (!F.areEqual (u[j], w[j])) {ret = false; break;}

	F.invin(d);
	Blackbox Rinv (F, n, d);

	Rinv.apply (u, w);

	report << "    Inverse output vector (should equal input): ";
	printVector<Field> (F, report, u);

	for (size_t j = 0; j < n; j++)
		if (!F.areEqual (u[j], v[j])) {ret = false; break;}

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();
	return ret;
}

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Scalar matrix black box test suite" << endl << endl;

	pass = pass && testScalarApply<Modular<long> > (F, n, report);

	//if (!testRandomMinpoly<Modular<long> >    (F, n, report, iterations)) pass = false;
	//if (!testRandomTranspose<Modular<long> >  (F, n, report, iterations)) pass = false;

	return pass ? 0 : -1;
}
