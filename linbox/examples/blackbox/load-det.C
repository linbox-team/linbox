/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* examples/blackbox/load-det.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------
 *
 * See COPYING for license information.
 *
 * --------------------------------------------------
 * Small program that loads and computes the determinant of a matrix whose
 * filename is given on the command line.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/det.h"

using namespace LinBox;
using namespace std;

typedef Modular<long> Field;

typedef vector <Field::Element> Vector;
typedef vector <Field::Element> Polynomial;
typedef vector <pair <size_t, Field::Element> > Row;
typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

// Constants: we are working with an n x n matrix over GF(q)
const int n = 1000;
const int q = 2147483647U;

int main (int argc, char **argv)
{
	long det_A;

	if (argc < 2) {
		cerr << "Usage: load-det <matrix>" << endl;
		return -1;
	}

	srand (time (NULL));

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

	det (det_A, A, F);

	cout << "Determinant is " << det_A << endl;

	return 0;
}
