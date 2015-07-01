/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// Copyright (C) 2001, 2002 Bradford Hovinen
/** @name examples/blackbox/example.C
 *
 * @author Bradford Hovinen <hovinen@cis.udel.edu>
 * @memo 
 * Simple example on Linbox use. Demonstrates loading and application of
 * blackbox matrix to a vector and computation of the minimal polynomial.
 * @doc 
 * FIXME what is shown different that in other minpoly example?
 *
 */
//@{
#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

// Field we are working with
//#include "linbox/field/modular.h"
#include "linbox/field/givaro-gfq.h"

// Black box classes we are going to work with
#include "linbox/blackbox/sparse.h"

// Minimal polynomial algorithm
#include "linbox/solutions/minpoly.h"

#include "linbox/vector/vector-domain.h"

//using namespace LinBox;
using namespace std;

// This is the field we are going to be working with - integers mod q
//typedef LinBox::Modular<LinBox::uint32> Field;
typedef LinBox::GivaroGfq Field;

// Some typedefs to make the type names less daunting
typedef vector <Field::Element> Vector;
typedef vector <Field::Element> Polynomial;
typedef vector <pair <size_t, Field::Element> > Row;
typedef LinBox::SparseMatrix <Field, Row> Blackbox;

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

	LinBox::minpoly (m_A, A, F);

	cout << "Minimal polynomial m_A of A is: ";
	printPolynomial (F, m_A);
	cout << endl;
}

/// no command line args
int main (int argc, char **argv)
{
	srand (time (NULL));

	// Construct the field GF(q) and a vector over GF(q)^n to
	// which to apply the matrix
	//Field F (q);
	Field F (2, 15);
	Vector v (n);

	// Construct and load the sparse test matrix
	Blackbox A (F, n, n);
	//ifstream input (EXAMPLE_DATADIR "/test.matrix");
	ifstream input ("mat.txt");
	A.read (input);

	// Run tests
	makeTestVector (F, A, v);
	testApply (F, A, v);
	testMinpoly (F, A);

	return 0;
}
//@}
