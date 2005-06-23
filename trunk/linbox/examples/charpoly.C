/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */


// Warning : example under development.
// integer computation must be completed
#include <iostream>
#include <iomanip>
#include "linbox/field/modular-double.h"
#include "linbox/field/PID-integer.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/solutions/charpoly.h"
#include "Matio.h"
#include "linbox/ring/givaro-polynomial.h"
using namespace LinBox;
using namespace std;

template <class Field, class Polynomial>
void printPolynomial (const Field &F, const Polynomial &v) 
{
	for (int i = v.size () - 1; i >= 0; i--) {
		F.write (cout, v[i]);
		if (i > 0)
			cout << " x^" << i << " + ";
	}
	cout << endl;
}

int main (int argc, char **argv)
{
	cout<<setprecision(8);
	cerr<<setprecision(8);
	if (argc < 2 || argc > 3) {
		cerr << "Usage: charpoly <matrix-file-in-SMS-format> [<p>]" << endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	if (argc == 2) {

		PID_integer ZZ;
		SparseMatrix<PID_integer> A (ZZ);
		A.read (input);
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
		typedef GivPolynomialRing<PID_integer,Dense> IntPolRing;
		IntPolRing::Element c_A;
		charpoly (c_A, A);


		cout << "Characteristic Polynomial is ";
		printPolynomial (ZZ, c_A);
	}
	if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		SparseMatrix<Field> B (F);
		B.read (input);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
		vector<Field::Element> c_B;
		charpoly (c_B, B);
		cout << "Characteristic Polynomial is ";
		printPolynomial (F, c_B);
	}

	return 0;
}
