/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/** \file examples/minpoly.C examples/minpoly.C 
\brief Minimal polynomial of a sparse matrix.
\ingroup examples
*/
#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/minpoly.h"

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

	if (argc < 2 || argc > 3) {
		cerr << "Usage: minpoly <matrix-file-in-SMS-format> [<p>]" << endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	if (argc == 2) {
	   int process = 0;
	   
           Method::Blackbox M;
#ifdef __LINBOX_HAVE_MPI
           Communicator C(&argc, &argv);
           process = C.rank();
           M.communicatorp(&C);
#endif
           
           PID_integer ZZ;
           SparseMatrix<PID_integer> A (ZZ);
           A.read (input);
           
           if(process == 0)
               cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
           
           vector<PID_integer::Element> m_A;
           minpoly (m_A, A, M); 
           
           
           
           if(process == 0){
               cout << "Minimal Polynomial is ";
               printPolynomial (ZZ, m_A);
           }
	}
	if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		SparseMatrix<Field> B (F);
		B.read (input);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;

		vector<Field::Element> m_B;
		minpoly (m_B, B);

		cout << "Minimal Polynomial is ";
		printPolynomial (F, m_B);
	}

	return 0;
}
