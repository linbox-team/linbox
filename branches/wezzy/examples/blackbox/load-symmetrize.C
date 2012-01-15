/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* -*- mode: C++; style: linux -*- */

/*
 * examples/blackbox/load-sylletrize.C
 *
 * Copyright (C) 2001, 2002, 2010 Bradford Hovinen <hovinen@cis.udel.edu>
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

/** @name examples/blackbox/load-symmetrize.C
 *
 * @author Bradford Hovinen <hovinen@cis.udel.edu>
 * @memo
 * Small program that loads and computes the minimal polynomial of A^T A, where
 * A is a matrix whose filename is given on the command line
 */
//@{

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/algorithms/blackbox-container-symmetrize.h"
#include "linbox/algorithms/massey-domain.h"

using namespace LinBox;
using namespace std;

typedef Modular<uint32_t> Field;

typedef vector <Field::Element> Vector;
typedef vector <Field::Element> Polynomial;
typedef SparseMatrix <Field, Vector> Blackbox;

// Constants: we are working with an n x n matrix over GF(q)
const int n = 1000;
const int q = 65521U;

/// load-symmetrize matrix-file
int main (int argc, char **argv)
{
	if (argc < 2) {
		cerr << "Usage: load-symmetrize <matrix>" << endl;
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

	unsigned long               deg;
	std::vector<Field::Element> P;

	BlackboxContainerSymmetrize<Field, Vector> TF (&A, F);
	MasseyDomain< Field, BlackboxContainerSymmetrize<Field, Vector> > WD (&TF, 20);

	WD.pseudo_minpoly (P, deg);

	cout << "Determinant is " << P[0] << endl;
	cout << "Degree is " << P.size () - 1 << endl;

	return 0;
}
//@}
