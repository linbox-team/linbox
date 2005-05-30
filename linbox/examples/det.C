/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
#include "linbox-config.h"

#include "linbox/field/modular-double.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/det.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{

	if (argc < 2 || argc > 3) {
		cerr << "Usage: load-det <matrix-file-in-SMS-format> [<p>]" << endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) {
		cerr << "Error: Cannot load matrix " << argv[1] << endl;
		return -1;
	}

	if (argc == 2) {

		SparseMatrix<GMP_Integers> A (F);
		A.read (input);
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		GMP_Integers::Element det_A;
		det (det_A, A)

		cout << "Determinant is " << det_A << endl;
	}
	if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		SparseMatrix<Field> B (F);
		B.read (input);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;

		Field::Element det_B;
		det (det_B, B)

		cout << "Determinant is ";
		F.write(cout, det_B) << " mod " << q << endl;
	}

	return 0;
}
