/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// Copyright (C) 2001, 2002 Bradford Hovinen
// See COPYING for license information.
/** @name examples/blackbox/load-det.C
 *
 * @author Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * @memo 
 * Small program that loads and computes the determinant of a matrix whose
 * filename is given on the command line.
 */
//@{
#include "linbox-config.h"

#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/dense.h"
#include "linbox/solutions/det.h"

using namespace LinBox;
using namespace std;

// Select our field: integers modulo a word-size (max. 31-bit) modulus
typedef Modular<uint32> Field;

// Select our black box: a sparse matrix over the above-mentioned field with
// default application vector and row representation types
typedef DenseMatrix<Field> Blackbox;

// Constants: we are working with a matrix over GF(q)

/// load-det matrix-file
int main (int argc, char **argv)
{
	Field::Element det_A;
    int q = 65521U;

	if (argc < 2 || argc > 3) {
		cerr << "Usage: load-det <matrix> [<p>]" << endl;
		return -1;
	}

	commentator.setMaxDepth (2);
	commentator.setReportStream (cout);

	if (argc == 3) q = atoi(argv[2]);
	Field F (q);

	Blackbox A (F);
	ifstream input (argv[1]);

	if (!input) {
		cerr << "Error: Cannot load matrix " << argv[1] << endl;
		return -1;
	}

	A.read (input); // size is determined by the input.
	cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

	det (det_A, A, F);

	cout << "Determinant is " << det_A << " mod " << q << endl;

	return 0;
}
//@}
