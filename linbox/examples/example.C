/* -*- mode: c; style: linux -*- */

/* examles/example.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 */

/* Simple example on Linbox use. Demonstrates loading and application of
 * blackbox matrix to a vector and computation of the minimal polynomial.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>

// Field we are working with
#include "linbox/field/modular.h"

// Black box classes we are going to work with
#include "linbox/blackbox/sparse0.h"

// Minimal polynomial algorithm
#include "linbox/solutions/minpoly.h"

#include "linbox/field/vector-domain.h"

using namespace LinBox;

// This is the field we are going to be working with - integers mod q
typedef Modular<long> Field;

// Some typedefs to make the type names less daunting
typedef vector <Field::element> Vector;
typedef vector <Field::element> Polynomial;
typedef vector <pair <size_t, Field::element> > Row;
typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

// Constants: we are working with an n x n matrix over GF(q)
const int n = 10;
const int q = 101;

void printVector (const Field &F, const Vector &v) 
{
	int i;

	cout << '(';
	for (i = 0; i < v.size (); i++) {
		F.write (cout, v[i]);
		if (i < v.size () - 1)
			cout << ", ";
	}
	cout << ')' << endl;
}

void printPolynomial (const Field &F, const Polynomial &v) 
{
	int i;

	for (i = v.size () - 1; i >= 0; i--) {
		F.write (cout, v[i]);

		if (i > 0)
			cout << " x^" << i << " + ";
	}

	cout << endl;
}

void makeTestVector (const Field &F, const Blackbox &A, Vector &v) 
{
	int i;

	v.resize (A.coldim ());

	for (i = 0; i < v.size (); i++)
		F.init (v[i], i);
}

void testApply (const Field &F, const Blackbox &A, const Vector &v)
{
	Vector w (A.coldim ());

	cout << "v = ";
	printVector (F, v);

	A.apply (w, v);

	cout << "Av = ";
	printVector (F, w);
	cout << endl;
}

void testMinpoly (const Field &F, const Blackbox &A) 
{
	Polynomial m_A;

	minpoly<Field, Polynomial, Vector> (m_A, A, F);

	cout << "Minimal polynomial m_A of A is: ";
	printPolynomial (F, m_A);
	cout << endl;
}

int main (int argc, char **argv)
{
	srand (time (NULL));

	// Construct the field GF(q) and a vector over GF(q)^n to
	// which to apply the matrix
	Field F (q);
	Vector v (n);

	// Construct and load the sparse test matrix
	Blackbox A (F, n, n);
	ifstream input (EXAMPLE_DATADIR "/test.matrix");
	A.read (input);

	// Run tests
	makeTestVector (F, A, v);
	testApply (F, A, v);
	testMinpoly (F, A);

	return 0;
}
