/*
 * Copyright (C) 2012 the LinBox group
 *
 * written by BB <bboyer@imag.fr>
 *
 *
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
 *.
 */

/*! @file  tests/test-toom-cook.C
 * @ingroup tests
 * @brief toom-cook multiplication routine
 * @test toom-cook multiplication routine
 */


#include <iostream>
#include <linbox-config.h>
#include <linbox/field/modular.h>
#include <linbox/matrix/dense-matrix.h>
#include <linbox/matrix/random-matrix.h>
// #include <fflas-ffpack/fflas/fflas.h>
#include <linbox/field/givaro.h>
#include <linbox/util/timer.h>

#include "linbox/algorithms/matrix-blas3/mul.h"

#include "test-common.h"
int main(int ac, char ** av) {
	static int p = 1009;
	static int e = 3 ;
	static size_t m = 10, n = 10 , k = 10;
	static size_t b = 10 ;

	static Argument as[] = {
		{ 'n', "-n N", "Set cols of C .",                 TYPE_INT,     &n },
		{ 'm', "-m N", "Set rows of C .",                 TYPE_INT,     &m },
		{ 'k', "-k N", "Set rows of B .",                 TYPE_INT,     &k },
		{ 'p', "-p N", "Set characteristic.",                 TYPE_INT,     &p },
		{ 'e', "-e N", "Set degree.",                 TYPE_INT,     &e },
		{ 'b', "-b N", "Set length of integers.",                 TYPE_INT,     &b },
		END_OF_ARGUMENTS
	};


	LinBox::parseArguments (ac, av, as);

	LinBox::Timer Tim ;
	{ /* Toom Cook over GivarorExtension */
		typedef LinBox::Modular<int64_t> Zpz;
		// typedef LinBox::Modular<double> Zpz;
		typedef LinBox::GivaroExtension< Zpz > GFpe ;

		// Z/pZ
		Zpz F(p);
		// GF(p^e) ;
		GFpe GF (F, 2);
		LinBox::MatrixDomain<GFpe> MD(GF);

		GF.write(std::cout << "This is the field with " << (LinBox::Integer)pow((LinBox::Integer)p,e) << " elements: ") << ", using: "   << GF.irreducible() << " as irreducible polynomial" << std::endl;
		std::cout << "matrices are " << m << 'x' << k << " and " << k << 'x' << n <<  std::endl;

		LinBox::BlasMatrix<GFpe> A(GF,m,k);
		LinBox::BlasMatrix<GFpe> B(GF,k,n);
		LinBox::BlasMatrix<GFpe> C(GF,m,n);

		typedef GFpe::RandIter Randiter;
		Randiter R(GF);
		LinBox::RandomDenseMatrix<Randiter,GFpe> randomizer(GF,R) ;
		randomizer.random(A);
		// std::cout << "A[0,0] = " << A.getEntry(0,0) << std::endl;
		randomizer.random(B);

		std::cout << "naive over GFq" << std::endl;
		Tim.clear(); Tim.start();
		LinBox::BLAS3::mul(C,A,B,LinBox::BLAS3::mulMethod::naive());
		Tim.stop();
		std::cout << Tim << '(' << C.getEntry(0,0) << ')' << std::endl;

		// for (size_t i = 0 ; i< m ; ++i)
		// for (size_t j = 0 ; j< n ; ++j)
		// C.setEntry(i,j,GF.zero);

		// std::cout << "A[0,0] = " << A.getEntry(0,0) << std::endl;
		{
			std::cout << "ToomCook low mem" << std::endl;
			LinBox::BlasMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			LinBox::BLAS3::mul(D,A,B,LinBox::BLAS3::mulMethod::ToomCook<GFpe>(GF,false));
			Tim.stop();
			std::cout << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				std::cout << "error" << std::endl;
				return 1;
			}
		}

		{
			std::cout << "ToomCook high mem" << std::endl;
			LinBox::BlasMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			LinBox::BLAS3::mul(D,A,B,LinBox::BLAS3::mulMethod::ToomCook<GFpe>(GF,true));
			Tim.stop();
			std::cout << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				std::cout << "error" << std::endl;
				return 1;
			}
		}

		{
			std::cout << "Matrix Domain" << std::endl;
			LinBox::BlasMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			MD.mul(D,A,B);
			Tim.stop();
			std::cout << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				std::cout << "error" << std::endl;
				return 1;
			}
		}
	}

	{ /* ZZ mat mul */

		LinBox::PID_integer ZZ ;
		LinBox::MatrixDomain<LinBox::PID_integer> MD(ZZ);
		LinBox::BlasMatrix<LinBox::PID_integer> A(ZZ,m,k) ;
		LinBox::BlasMatrix<LinBox::PID_integer> B(ZZ,k,n) ;
		LinBox::BlasMatrix<LinBox::PID_integer> C(ZZ,m,n) ;

		A.random((unsigned)b);
		B.random((unsigned)b);

		std::cout << "Naïve " << std::endl ;
		Tim.clear() ; Tim.start() ;
		LinBox::BLAS3::mul(C,A,B,LinBox::BLAS3::mulMethod::naive());
		Tim.stop();
		std::cout << Tim << '(' << C.getEntry(0,0) << ')' << std::endl;

#ifdef __LINBOX_HAVE_FLINT
		{
			std::cout << "FLINT " << std::endl;
			LinBox::BlasMatrix<LinBox::PID_integer> D(ZZ,m,n);
			Tim.clear(); Tim.start();
			LinBox::BLAS3::mul(D,A,B,LinBox::BLAS3::mulMethod::FLINT());
			Tim.stop();
			std::cout << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				std::cout << "error" << std::endl;
				return 1;
			}
		}
#endif // __LINBOX_HAVE_FLINT

		{
			std::cout << "Matrix Domain" << std::endl;
			LinBox::BlasMatrix<LinBox::PID_integer> D(ZZ,m,n);
			Tim.clear(); Tim.start();
			MD.mul(D,A,B);
			Tim.stop();
			std::cout << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				std::cout << "error" << std::endl;
				return 1;
			}
		}
	}


	return 0;
}
