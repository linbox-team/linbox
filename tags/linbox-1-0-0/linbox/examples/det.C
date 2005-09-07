/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/**\file examples/det.C examples/det.C
\brief Determinant of sparse matrix over Z or Zp.
\ingroup examples
*/
//#include "linbox-config.h"
#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/det.h"
#include "linbox/util/matrix-stream.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{

	if (argc < 2 || argc > 3) {
		cerr << "Usage: det <matrix-file-in-supported-format> [<p>]" << endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	if (argc == 2) {

		GMP_Integers ZZ;
		MatrixStream< GMP_Integers > ms( ZZ, input );
		SparseMatrix<GMP_Integers> A (ms);
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		GMP_Integers::Element det_A;
		det (det_A, A);

		cout << "Determinant is ";
		ZZ.write(cout, det_A) << endl;
	}
	if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		MatrixStream< Field > ms ( F, input );
		SparseMatrix<Field> B (ms);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;

		Field::Element det_B;
		det (det_B, B);

		cout << "Determinant is ";
		F.write(cout, det_B) << " mod " << q << endl;
	}

	return 0;
}
