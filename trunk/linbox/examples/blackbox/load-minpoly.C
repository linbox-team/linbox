/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// Copyright (C) 2001, 2002 Bradford Hovinen
// See COPYING for license information.
/** @name examples/blackbox/load-minpoly.C
 *
 * @author Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * @memo 
 * Small program that loads and computes the minimal polynomial of a matrix
 * whose filename is given on the command line.
 */
//@{

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/minpoly.h"

using namespace LinBox;
using namespace std;

// Select our field: integers modulo a word-size (max. 31-bit) modulus
typedef Modular<uint32> Field;

// Select our black box: a sparse matrix over the above-mentioned field with
// default application vector and row representation types
typedef SparseMatrix<Field> Blackbox;

// We are using dense vectors to represent polynomials
typedef Vector<Field>::Dense Polynomial;

// Constants: we are working with an n x n matrix over GF(q)
const int n = 1000;
const int q = 65521U;

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

/// load-minpoly matrix-file
int main (int argc, char **argv)
{
	Polynomial m_A;

	if (argc < 2) {
		cerr << "Usage: load-minpoly <matrix>" << endl;
		return -1;
	}

	commentator.setMaxDepth (2);
	commentator.setReportStream (cout);

	Field F (q);

	Blackbox A (F, n, n);
	ifstream input (argv[1]);

	if (!input) {
		cerr << "Error: Cannot load matrix " << argv[1] << endl;
		return -1;
	}

	A.read (input);

	minpoly (m_A, A, F);

	cout << "Minimal polynomial is ";
	printPolynomial (F, m_A);

	return 0;
}
//@}
