/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/** @name mmZZ_pE.C
 * @memo The wrappered NTL::ZZ_pE
 * @doc
 * showing a use of NTL's finite extension fields
 *
 * FIXME - What does it do? Something more should be written here
 *
 * @author Zhendong Wan <wan@mail.eecis.udel.edu>
 */
// See COPYING for license information
//@{ 

// FIXME  it appears DenseMatrix has changed. -bds

#include <iostream>
#include <linbox/blackbox/dense.h>
#include <linbox/field/FieldBLAS.h>
#include <linbox/field/ntl-ZZ_pE.h>
#include <iterator>
#include <algorithm>

using namespace LinBox;
using namespace std;

/// usage: mmZZ_pE n p e, for using n by n matrix over GF(p^e) 
int main (int argc, const char* argv[])
{
        // argument parsing.
	int n, p, e; // note: larger p is possible with code change.

	if (argc != 4) {
		std::cout << "Usage: " << argv[0] << " n p e\n";
		std::cout << "returns k, after computing over GF(p^e) a product of k n by n random matrices\n"; 
		std::cout << "over GF(p^e)\n which is zero, with no leading subproducts zero.\n";
		return 0;
	}
	else {
		n = atoi(argv[1]);
		p = atoi(argv[2]);
		e = atoi(argv[3]);
	}

	// types
	typedef UnparametricField<NTL::ZZ_pE> Field;
	typedef DenseMatrix<Field> Matrix;

	// generate the field
	srand(time(0));

	NTL::ZZ_p::init(NTL::to_ZZ(p)); // prime field
	NTL::ZZ_pX Poly;
        NTL::BuildIrred(Poly, e); // generate an irreducible polynomial P
	// of degree e over GF(p)

        NTL::ZZ_pE::init(Poly);
	Field F;  // F represents GF(p^e).

	// generate random matrix P as initial product.

	Field::RandIter r(F,16,0); // what??
	Matrix P(F, n, n, r); // n x n matrix, entries uniformly random from F.
	int counter = 1;
	std::cout<<"matrix " << counter << ": ";
	std::copy(P.rawBegin(),P.rawEnd(),std::ostream_iterator<NTL::ZZ_pE>(std::cout, " "));
        std::cout<<"\n";



	// generate matrix domain which provides the matrix functions isZero() and mul(). 

	FieldBLAS<Field> MD(F);

	// iterate until product is zero.

	cout << "\nproduct length: " << counter;
	cout.flush();

	for ( ;  ! MD.isZero(P) ; ++counter) {
		Matrix M(F, n, n, r);
		MD.mulin(P, M);
		if (! (counter % 100)) {
			cout << "\rproduct length: " << counter;
			cout.flush();
		}

		if (n < 4) {
			std::cout<<"matrix " << counter << ": ";
			std::copy(M.rawBegin(),M.rawEnd(),std::ostream_iterator<NTL::ZZ_pE>(std::cout, " "));
			std::cout<<",        product:";
			std::copy(P.rawBegin(),P.rawEnd(),std::ostream_iterator<NTL::ZZ_pE>(std::cout, " "));
			std::cout<<"\n";
		}
	}
	std::cout << "steps to zero matrix: " << counter << "\n";
	return 0;
}
//@}
