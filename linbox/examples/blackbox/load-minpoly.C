
/*
 * examples/blackbox/load-minpoly.C
 *
 * Copyright (C) 2001, 2002 Bradford Hovinen <hovinen@cis.udel.edu>
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/** @name examples/blackbox/load-minpoly.C
 *
 * @author Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * @memo
 * Small program that loads and computes the minimal polynomial of a matrix
 * whose filename is given on the command line.
 */
//@{

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/minpoly.h"

using namespace LinBox;
using namespace std;

// Select our field: integers modulo a word-size (max. 31-bit) modulus
typedef Modular<uint32_t> Field;

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

	commentator().setMaxDepth (2);
	commentator().setReportStream (cout);

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

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

