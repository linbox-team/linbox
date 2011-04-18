/** 
 * examples/charpoly.C
 *
 * Copyright (C) 2005, 2010 D. Saunders, C. Pernet, J-G. Dumas 
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

/** \file examples/charpoly.C
    \brief Characteristic polynomial of matrix over Z or Zp.
    \ingroup examples
*/
#include <iostream>
#include <iomanip>

#include "linbox/util/timer.h"
#include "linbox/field/modular-double.h"
#include "linbox/field/unparametric.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/blas-blackbox.h"
using namespace std;

#include "linbox/solutions/charpoly.h"
#include "linbox/ring/givaro-polynomial.h"
using namespace LinBox;

template <class Field, class Polynomial>
std::ostream& printPolynomial (std::ostream& out, const Field &F, const Polynomial &v) 
{
	for (int i = v.size () - 1; i >= 0; i--) {
		F.write (out, v[i]);
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
  } else {
    if(v[n] != 0) {
      if (v[n] != 1) F.write(out, v[n]) << '*';
      out << 'X';
      if (n > 1) out << '^' << n;
      for (int i = n - 1; i > 0; i--) {
	  if (v[i] != 0) {
	    if (v[i] >0) out << " + ";
	    if (v[i] != 1) F.write (out, v[i]) << '*';
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
    prettyprintIntegerPolynomial(out << '(', F, *(*itf)) << ')';
    if (*ite > 1) out << '^' << *ite;
    out << endl;
  }
  return out;
}

int main (int argc, char **argv)
{
	commentator.setMaxDetailLevel (2);
	commentator.setMaxDepth (2);
	commentator.setReportStream (std::cerr);

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
		DenseMatrix<PID_integer > A (ZZ);
		A.read (input);
		typedef GivPolynomialRing<PID_integer,Dense> IntPolRing;
		IntPolRing::Element c_A;

		Timer tim; tim.clear();tim.start();charpoly (c_A, A, Method::Blackbox());
                tim.stop();

		cout << "Characteristic Polynomial is ";
		    //printPolynomial (cout, ZZ, c_A) << endl;
                cout << tim << endl;

#ifdef __LINBOX_HAVE_NTL
		cout << "Do you want a factorization (y/n) ? ";
		char tmp;
		cin >> tmp;
		if (tmp == 'y' || tmp == 'Y') {
                  commentator.start("Integer Polynomial factorization by NTL", "NTLfac");
		  vector<IntPolRing::Element*> intFactors;    
		  vector<unsigned long> exp;
		  IntPolRing IPD(ZZ);
                  tim.start();
		  IPD.factor (intFactors, exp, c_A);
                  tim.stop();
                  commentator.stop("done", NULL, "NTLfac");
		  printFactorization(cout << intFactors.size() << " integer polynomial factors:" << endl, ZZ, intFactors, exp) << endl;
		  vector<IntPolRing::Element*>::const_iterator itf = intFactors.begin();
		  for ( ; itf != intFactors.end(); ++itf) 
			  delete *itf;
			  
cout << tim << endl;
		  
		}
#endif
	}
	if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		DenseMatrix<Field> B (F);
		B.read (input);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
		GivPolynomialRing<Field,Dense>::Element c_B;
                Timer tim; tim.clear();tim.start();
		charpoly (c_B, B);
		tim.stop();
		
		cout << "Characteristic Polynomial is ";
		    //printPolynomial (cout, F, c_B) << endl;
                cout << tim << endl;
	}

	return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
