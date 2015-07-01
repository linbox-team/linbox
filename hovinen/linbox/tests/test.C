/* -*- mode: c; style: linux -*- */

/* linbox/tests/test.C
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>

// Field we are working with
#include "linbox/field/param-modular.h"

// Black box classes we are going to work with
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/blackbox/hilbert.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/blackbox/sparse-matrix.h"

#include "linbox/solutions/minpoly.h"

using namespace LinBox;

template <class Blackbox, class Field>
class TestBlackbox {
	const Field     &F;
	const Blackbox  &blackbox;
	ostream         &output;

public:
	typedef vector<typename Field::Element> Polynomial;

	TestBlackbox (const Field &_F, const Blackbox &_blackbox, ostream &_output) 
		: F (_F), blackbox (_blackbox), output (_output)
	{
		output << "Testing black box functionality" << endl;
		output << "Field: ";
		F.write (output);
		output << endl;
	}

	void testApply (const vector <typename Field::Element> &v)
	{
		vector <typename Field::Element> w (blackbox.coldim ());

		output << "Testing apply ()..." << endl;

		output << "v = ";
		printVector (v);

		blackbox.apply (w, v);

		output << "Av = ";
		printVector (w);
		output << endl;
	}

	void testMinpoly (void) 
	{
		Polynomial m_A;

		minpoly<Field, Polynomial, vector<typename Field::Element>, MethodTrait::Wiedemann> (m_A, blackbox, F);

		output << "Minimal polynomial m_A of A is: ";
		printPolynomial (m_A);
		output << endl;
	}

	void makeTestVector (vector <typename Field::Element> &v) 
	{
		int i;

		v.resize (blackbox.coldim ());

		for (i = 0; i < v.size (); i++)
			F.init (v[i], i);
	}

	void printVector (const vector <typename Field::Element> &v) 
	{
		int i;

		output << '(';
		for (i = 0; i < v.size (); i++) {
			F.write (output, v[i]);
			if (i < v.size () - 1)
				output << ", ";
		}
		output << ')' << endl;
	}

	void printPolynomial (const Polynomial &v) 
	{
		int i;

		for (i = v.size () - 1; i >= 0; i--) {
			F.write (output, v[i]);

			if (i > 0)
				output << " x^" << i << " + ";
		}

		output << endl;
	}
};

int main (int argc, char **argv)
{
	// This is the field we are going to be working with - integers mod q
	typedef ParamModular Field;

	// Some typedefs to make the type names less daunting
	typedef vector <Field::Element> Vector;
	typedef vector <pair <size_t, Field::Element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	// Constants: we are working with an n x n matrix over GF(q)
	const int n = 10;
	const int q = 23;

	// Construct the field GF(q) and a vector over GF(q)^n to which to apply the matrix
	Field F (q);
	Vector v (n);

	// Construct and load the sparse test matrix
	Blackbox A (F, n, n);
	ifstream input ("test.matrix");
	A.read (input);

	// Construct the test object defined above
	TestBlackbox <Blackbox, Field> test (F, A, cout);

	// Run tests
	test.makeTestVector (v);
	test.testApply (v);
	test.testMinpoly ();

	return 0;
}
