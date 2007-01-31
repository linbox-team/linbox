/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/** \file examples/charpoly.C
\brief Characteristic polynomial of matrix over Z or Zp.
 \ingroup examples
*/

// Warning : example under development.
// integer computation must be completed
#include <iostream>
#include <iomanip>
#include "Matio.h"

#include "linbox/field/modular-double.h"
#include "linbox/field/unparametric.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/blas-blackbox.h"
using namespace std;

#include "linbox/solutions/charpoly.h"
#include "linbox/ring/givaro-polynomial.h"
#include "linbox/element/givaro-polynomial.h"
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
	cout<<setprecision(8);
	cerr<<setprecision(8);
	if (argc < 2 || argc > 3) {
		cerr << "Usage: charpoly <matrix-file-in-SMS-format> [<p>]" << endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	if (argc == 2) {

		UnparametricField<integer> ZZ;
		SparseMatrix<UnparametricField<integer> > A (ZZ);
		A.read (input);
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
		typedef GivPolynomialRing<UnparametricField<integer>,Dense> IntPolRing;
		IntPolRing::Element c_A;
		charpoly (c_A, A);


		cout << "Characteristic Polynomial is ";
		printPolynomial (cout, ZZ, c_A) << endl;

#ifdef __LINBOX_HAVE_NTL
		cout << "Do you want a factorization (y/n) ? ";
		char tmp;
		cin >> tmp;
		if (tmp == 'y' || tmp == 'Y') {
		  vector<IntPolRing::Element*> intFactors;    
		  vector<unsigned long> exp;
		  IntPolRing IPD(ZZ);
		  IPD.factor (intFactors, exp, c_A);
		  printFactorization(cout << intFactors.size() << " integer polynomial factors:" << endl, ZZ, intFactors, exp) << endl;
		}
#endif
	}
	if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		SparseMatrix<Field> B (F);
		B.read (input);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
		GivPolynomial<Field::Element> c_B;
		charpoly (c_B, B);
		cout << "Characteristic Polynomial is ";
		printPolynomial (cout, F, c_B) << endl;
	}

	return 0;
}
