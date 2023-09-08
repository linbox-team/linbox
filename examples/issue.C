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

typedef Givaro::ZRing<Givaro::Integer> Z_Domain;
//typedef Givaro::Modular<double> CoeffDomain;

int main()
{

    //ifstream input("issue_sage_35846_matrix.sms");
    ifstream input("mat.sms");


    Z_Domain Z;   // <-- integer
    //CoeffDomain R(CoeffDomain::maxCardinality()); // <-- modular
    BlasMatrix<Z_Domain> mat(Z);
    mat.read(input);

  


    DensePolynomial<Z_Domain> cp(Z);
    charpoly(cp, mat);


    clog << "Characteristic Polynomial is ";
    prettyprintIntPoly(clog, Z, cp) << endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
