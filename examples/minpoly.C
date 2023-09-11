/*
 * examples/minpoly.C
 *
 * Copyright (C) 2005, 2010 D Saunders, J-G Dumas
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

/** \file examples/minpoly.C
 * @example  examples/minpoly.C
  \brief Minimal polynomial of a sparse matrix.
  \ingroup examples
  */

#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif

#define __LB_VALENCE_REPORTING__ 1
#define __LB_CRA_REPORTING__ 1
#define __LB_CRA_TIMING__ 1
#define __LINBOX_HEURISTIC_CRA 1



#include <linbox/linbox-config.h>
#include <iostream>
//! @bug this should be elsewhere
template <class Field, class Polynomial>
std::ostream& printPolynomial (std::ostream& out, const Field &F, const Polynomial &v)
{
	for (int i = (int)v.size () ; i-- ; ) {
		F.write (out, v[(size_t)i]);
		if (i > 0)
			out << " x^" << i << " + ";
	}
    return out;
}

#include <linbox/ring/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/blackbox/compose.h>
#include <linbox/polynomial/dense-polynomial.h>
#include <linbox/solutions/minpoly.h>

using namespace LinBox;
using namespace std;


int main (int argc, char **argv)
{
	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::clog);

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
        LinBox::Timer chrono;
        Method::Auto M;

		Givaro::ZRing<Integer> ZZ;
		DensePolynomial<Givaro::ZRing<Integer> > m_A(ZZ);
		typedef SparseMatrix<Givaro::ZRing<Integer>> SpMat;
		SpMat B (ZZ);
		B.read (input);
        if (B.rowdim() == B.coldim()) {

			std::clog << "B is " << B.rowdim() << " by " << B.coldim() << endl;

            chrono.start();
            PAR_BLOCK {
                minpoly (m_A, B, M);
            }

            chrono.stop();


        } else {
                // Minpoly of B . Transpose(B) or Transpose(B) . B
            typedef Transpose<SpMat> BB2;
            BB2 BT(B);
            if (B.rowdim() > B.coldim()) {
                Compose<BB2, SpMat> A(BT,B);
                std::clog << "(B^T . B) is " << A.rowdim() << " by " << A.coldim() << endl;

                chrono.start();
                PAR_BLOCK {
                    minpoly (m_A, A, M);
                }
                chrono.stop();
            } else {
                Compose<SpMat, BB2 > A(B,BT);
                std::clog << "(B . B^T) is " << A.rowdim() << " by " << A.coldim() << endl;

                chrono.start();
                PAR_BLOCK {
                    minpoly (m_A, A, M);
                }
                chrono.stop();
            }
        }


        std::clog << "Minimal Polynomial, " << chrono << ", is: ";
        printPolynomial (std::cout, ZZ, m_A) << std::endl;

	}
	else{
		typedef Givaro::Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		SparseMatrix<Field> B (F);
        Method::Auto M;
		B.read (input);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;

		DensePolynomial<Field> m_B(F);
		minpoly (m_B, B, M);

		cout << "Minimal Polynomial is ";
		printPolynomial (std::cout, F, m_B) << std::endl;

	}

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
