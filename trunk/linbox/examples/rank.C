/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
#include "linbox-config.h"

#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/rank.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{

	if (argc < 2 || argc > 3) 
	{	cerr << "Usage: rank <matrix-file-in-SMS-format> [<p>]" << endl; return -1; }

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

	long unsigned int r;

	if (argc == 2) { // rank over the integers.

	   /* We could pick a random prime and work mod that prime, But the point here 
	   is that the rank function in solutions/ handles that issue.  Our matrix here 
	   is an integer matrix and our concept is that we are getting the rank of that 
	   matrix by some blackbox magic inside linbox.
	   */

		GMP_Integers ZZ;
		SparseMatrix<GMP_Integers> A (ZZ);
		A.read (input);
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		rank (r, A);

		cout << "Rank is " << r << endl;
	}
	if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		SparseMatrix<Field> B (F);
		B.read (input);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;

		rank (r, B);

		cout << "Rank is " << r << endl;
	}

	return 0;
}
