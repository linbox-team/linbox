
#include <iostream>
#include <linbox/linbox-config.h>
#include <linbox/field/modular.h>
#include <linbox/matrix/blas-matrix.h>
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

	static Argument as[] = {
		{ 'n', "-n N", "Set cols of C .",                 TYPE_INT,     &n },
		{ 'm', "-m N", "Set rows of C .",                 TYPE_INT,     &m },
		{ 'k', "-k N", "Set rows of B .",                 TYPE_INT,     &k },
		{ 'p', "-p N", "Set characteristic.",                 TYPE_INT,     &p },
		{ 'e', "-e N", "Set degree.",                 TYPE_INT,     &e },
		END_OF_ARGUMENTS
	};


	LinBox::parseArguments (ac, av, as);

	if (e < 2) {
		std::cout << " e > 1 please ! " << std::endl;
		return 0 ;
	}


	typedef LinBox::Modular<double> Zpz;
	typedef LinBox::GivaroExtension< Zpz > GFpe ;

	// Z/pZ
	Zpz F(p);
	// GF(p^e) ;
	GFpe GF (F, 2);

	GF.write(std::cout << "This is the field with " << (LinBox::Integer)pow((LinBox::Integer)p,e) << " elements: ") << ", using: "   << GF.irreducible() << " as irreducible polynomial" << std::endl;

	LinBox::BlasMatrix<GFpe> A(GF,m,k);
	LinBox::BlasMatrix<GFpe> B(GF,k,n);
	LinBox::BlasMatrix<GFpe> C(GF,m,n);

	typedef typename GFpe::RandIter Randiter;
	Randiter R(GF);
	LinBox::RandomDenseMatrix<Randiter,GFpe> randomizer(GF,R) ;
	randomizer.random(A);
	// std::cout << "A[0,0] = " << A.getEntry(0,0) << std::endl;
	randomizer.random(B);

	LinBox::Timer Tim ;
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
		if (m*n < 100) {
			for (size_t i = 0 ; i < m ; ++i)
				for (size_t j = 0 ; j < n ; ++j)
					if (!(GF.areEqual(D.getEntry(0,0),C.getEntry(0,0))))
						return 1;
		}
		else {
			int r =100 ;
			while (--r)  {
				size_t i = rand()%m;
				size_t j = rand()%n;
				if (!(GF.areEqual(D.getEntry(i,j),C.getEntry(i,j))))
					return 1;
			}
			// TODO check with apply !
		}
	}

	{
		std::cout << "ToomCook high mem" << std::endl;
		LinBox::BlasMatrix<GFpe> D(GF,m,n);
		Tim.clear(); Tim.start();
		LinBox::BLAS3::mul(D,A,B,LinBox::BLAS3::mulMethod::ToomCook<GFpe>(GF,true));
		Tim.stop();
		std::cout << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;
		if (m*n < 100) {
			for (size_t i = 0 ; i < m ; ++i)
				for (size_t j = 0 ; j < n ; ++j)
					if (!(GF.areEqual(D.getEntry(0,0),C.getEntry(0,0))))
						return 1;
		}
		else {
			int r =100 ;
			while (--r)  {
				size_t i = rand()%m;
				size_t j = rand()%n;
				if (!(GF.areEqual(D.getEntry(i,j),C.getEntry(i,j))))
					return 1;
			}
			// TODO check with apply !
		}
	}
	return 0;
}
