/* -*- mode: C++; style: linux -*- */

// Copyright (C) 2001, 2002 Bradford Hovinen
// See COPYING for license information.
/** @name examples/blackbox/load-symmetrize.C
 *
 * @author Bradford Hovinen <hovinen@cis.udel.edu>
 * @memo 
 * Small program that loads and computes the minimal polynomial of A^T A, where
 * A is a matrix whose filename is given on the command line
 */
//@{

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/algorithms/blackbox-container-symmetrize.h"
#include "linbox/algorithms/massey-domain.h"

using namespace LinBox;
using namespace std;

typedef Modular<uint32> Field;

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
