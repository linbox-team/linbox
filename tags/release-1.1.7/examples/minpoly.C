
/** 
 * examples/minpoly.C
 *
 * Copyright (C) 2005, 2010 D Saunders, J-G Dumas
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */

/** \file examples/minpoly.C examples/minpoly.C 
\brief Minimal polynomial of a sparse matrix.
\ingroup examples
*/
#include <iostream>
template <class Field, class Polynomial>
void printPolynomial (const Field &F, const Polynomial &v) 
{
	for (int i = v.size () - 1; i >= 0; i--) {
		F.write (std::cout, v[i]);
		if (i > 0)
			std::cout << " x^" << i << " + ";
	}
	std::cout << std::endl;
}

#include "linbox/field/modular-double.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/minpoly.h"

using namespace LinBox;
using namespace std;


int main (int argc, char **argv)
{
        commentator.setMaxDetailLevel (-1);
        commentator.setMaxDepth (-1);
        commentator.setReportStream (std::cerr);
		  
		  int a = argc;
		  while(a--){
			  cerr << "argv[" << a << "] = " << argv[a] << endl;
		  }
	if (argc < 2) {
		cerr << "Usage: minpoly <matrix-file-in-SMS-format> [<p>]" << endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	if (argc != 3) {
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

//            minpoly (m_A, A, Method::BlasElimination() );
           
//            if(process == 0){
//                cout << "Minimal Polynomial is ";
//                printPolynomial (ZZ, m_A);
//            }
           
	}
	else{
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
                
//                 minpoly (m_A, A, Method::BlasElimination() );
//  		cout << "Minimal Polynomial is ";
// 		printPolynomial (F, m_B);
                
	}

	return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
