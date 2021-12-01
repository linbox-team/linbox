#define _LB_CRATIMING 1
//#define __LINBOX_HEURISTIC_CRA
/*
 * examples/charpoly.C
 *
 * Copyright (C) 2005, 2010 D. Saunders, C. Pernet, J-G. Dumas
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

/** \file examples/charpoly.C
 * @example  examples/charpoly.C
 \brief Characteristic polynomial of matrix over Z or Zp.
 \ingroup examples
*/
#include <linbox/linbox-config.h>

#include <iostream>
#include <iomanip>

#include <linbox/util/timer.h>
#include <linbox/ring/modular.h>
#include <givaro/zring.h>
#include <linbox/matrix/sparse-matrix.h>
using namespace std;

#include <linbox/solutions/charpoly.h>
#include <linbox/ring/polynomial-ring.h>
using namespace LinBox;

template <class Field, class Polynomial>
std::ostream& printPolynomial (std::ostream& out, const Field &F, const Polynomial &v)
{
	for (int i = (int)(v.size () - 1); i >= 0; i--) {
		F.write (out, v[(size_t)i]);
		if (i > 0)
			out << " X^" << i << " + ";
	}
	return out;
}
template <class Field, class Polynomial>
std::ostream& prettyprintIntegerPolynomial (std::ostream& out, const Field &F, const Polynomial &v)
{
	size_t n = v.size()-1;
	if (n == 0) {
		F.write(out, v[0]);
	}
	else {
		if(v[n] != 0) {
			if (v[n] != 1) F.write(out, v[n]) << '*';
			out << 'X';
			if (n > 1) out << '^' << n;
			for (int i = (int)n - 1; i > 0; i--) {
				if (v[(size_t)i] != 0) {
					if (v[(size_t)i] >0) out << " + ";
					if (v[(size_t)i] != 1) F.write (out, v[(size_t)i]) << '*';
					out << 'X';
					if (i > 1) out << '^' << i;
				}
			}
			if (v[0] != 0) {
				if (v[0] >0) out << " + ";
				F.write(out, v[0]);
			}
		}
	}
	return out;
}
template <class Field, class Factors, class Exponents>
std::ostream& printFactorization (std::ostream& out, const Field &F, const Factors &f, const Exponents& exp)
{
	typename Factors::const_iterator itf = f.begin();
	typename Exponents::const_iterator ite = exp.begin();
	for ( ; itf != f.end(); ++itf, ++ite) {
		prettyprintIntegerPolynomial(out << '(', F, *itf) << ')';
		if (*ite > 1) out << '^' << *ite;
		out << endl;
	}
	return out;
}

typedef Givaro::ZRing<Givaro::Integer> IntDom;

int main (int argc, char **argv)
{
	commentator().setMaxDetailLevel (2);
	commentator().setMaxDepth (2);
	commentator().setReportStream (std::cerr);

	cout<<setprecision(8);
	cerr<<setprecision(8);
	if (argc < 2 || argc > 3) {
		cerr << "Usage: charpoly <matrix-file-in-SMS-format> [<p>]" << endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	if (argc == 2) {

		IntDom ZZ;
		DenseMatrix<IntDom > A (ZZ);
		A.read (input);
        DensePolynomial<IntDom> c_A(ZZ);

		Timer tim; tim.clear();tim.start();
		charpoly (c_A, A);
		tim.stop();

		clog << "Characteristic Polynomial is ";
		printPolynomial (clog, ZZ, c_A) << endl;
		cout << tim << endl;

#ifdef __LINBOX_HAVE_NTL
		clog << "Do you want a factorization (y/n) ? ";
		char tmp;
		cin >> tmp;
		if (tmp == 'y' || tmp == 'Y') {
			commentator().start("Integer Polynomial factorization by NTL", "NTLfac");
			vector<DensePolynomial<IntDom> > intFactors;
			vector<uint64_t> exp;
			PolynomialRing<IntDom> IPD(ZZ);
			tim.start();
			IPD.factor (intFactors, exp, c_A);
			tim.stop();
			commentator().stop("done", NULL, "NTLfac");
			printFactorization(clog << intFactors.size() << " integer polynomial factors:" << endl, ZZ, intFactors, exp) << endl;

			cout << tim << endl;

		}
#endif
	}
	if (argc == 3) {

		typedef Givaro::Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		DenseMatrix<Field> B (F);
		B.read (input);
		clog << "B is " << B.rowdim() << " by " << B.coldim() << endl;
        DensePolynomial<Field> c_B(F);
		Timer tim; tim.clear();tim.start();
		charpoly (c_B, B);
		tim.stop();

		clog << "Characteristic Polynomial is ";
		printPolynomial (clog, F, c_B) << endl;
		cout << tim << endl;
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
