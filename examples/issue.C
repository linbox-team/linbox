/** compilation: 
* g++ issue_sage_35846_charpoly.C -c -o issue_charpoly.o
* g++ -o issue_charpoly issue_charpoly.o -fopenmp -llinbox -lntl -lblis -lgivaro -lgmp
*/

#include <givaro/modular.h>
#include <givaro/modular-inttype.h>

#include <fflas-ffpack/field/field-traits.h>
#include <iostream>

#include <linbox/linbox-config.h>
#include <linbox/solutions/charpoly.h>
#include <linbox/matrix/dense-matrix.h>

using namespace LinBox;
using namespace std;

template <class Field, class Polynomial>
std::ostream& prettyprintIntPoly (std::ostream& out, const Field &F, const Polynomial &v)
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

typedef Givaro::ZRing<Givaro::Integer> CoeffDomain;
//typedef Givaro::Modular<double> CoeffDomain;

int main()
{
	clog << "***********   ICI " << endl;

    ifstream input("issue_sage_35846_matrix.sms");

	clog << "***********   ICI " << endl;

    CoeffDomain R;   // <-- integer
    //CoeffDomain R(CoeffDomain::maxCardinality()); // <-- modular
    BlasMatrix<CoeffDomain> mat(R);
    mat.read(input);

    mat.write(cerr<<"LinBox::Echelon = "<<endl)<<endl;

	clog << "***********   ICI " << endl;

    DensePolynomial<CoeffDomain> cp(R);
    charpoly(cp, mat);

	cerr << "***********   ICI " << endl;

    clog << "Characteristic Polynomial is ";
    prettyprintIntPoly(clog, R, cp) << endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
