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
#include "linbox-config.h"
#include "linbox/field/modular.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/random-matrix.h"
// #include <fflas-ffpack/fflas/fflas.h>
#include "linbox/field/givaro.h"
#include "linbox/util/timer.h"

#include "linbox/algorithms/matrix-blas3/mul.h"

#include "test-common.h"



namespace LinBox { namespace Protected {
	struct IntegerSparseCraMatMul {


		typedef SparseMatrixFormat::CSR spfmt ;
		typedef Modular<double>         Field;
		typedef Field::Element          Element;
		typedef SparseMatrix<Field,spfmt>       ModularMatrix ;
		typedef BlasVector<Field>               ModularVector;
		typedef BlasVector<PID_integer>         IntegerVector;
		typedef SparseMatrix<PID_integer,spfmt> IntegerMatrix ;

#ifdef _LB_MM_TIMING
#ifdef _OPENMP
		typedef LinBox::OMPTimer Mytime;
#else
		typedef LinBox::Timer    Mytime;
#endif
#endif

		const IntegerMatrix &_A_ ;
		const IntegerVector &_B_ ;

#ifdef _LB_MM_TIMING
		mutable Mytime chrono;
#endif

		IntegerSparseCraMatMul(const IntegerMatrix& A, const IntegerVector& B) :
			_A_(A), _B_(B)
		{
#ifdef _LB_MM_TIMING
			chrono.clear();
#endif
			linbox_check(A.getPointer() == _A_.getPointer());
		}

		IntegerSparseCraMatMul(IntegerMatrix& A, IntegerVector& B) :
			_A_(A), _B_(B)
		{
#ifdef _LB_MM_TIMING
			chrono.clear();
#endif
			linbox_check(A.getPointer() == _A_.getPointer());
		}

		ModularVector& operator()(ModularVector& Cp, const Field& F) const
		{
			// BlasMatrixDomain<Field>   BMD(F);

			/*  intialisation */
			// ModularMatrix Cpp(_A_.rowdim(),_B_.coldim());
			// Cp = Cpp ;
			ModularMatrix Ap(_A_, F);
			ModularVector Bp(F, _B_);
			Cp.resize(Ap.rowdim());

			/*  multiplication mod p */

#ifdef _LB_MM_TIMING
			Mytime matmul; matmul.clear(); matmul.start();
#endif
			// BMD.mul(Cp,Ap,Bp);
			Ap.apply(Cp,Bp);
			// BMD.axpyin(Cp,Ap,Bp);
#if 0
			// BMD.mul( static_cast<BlasMatrix<double>&>(Cp),Ap,Bp);
			if (FAM_TYPE == _axpy)
				BMD.axpyin(Cp,Ap,Bp);
			else if (FAM_TYPE == _axmy)
				BMD.axmyin(Cp,Ap,Bp);
			else if (FAM_TYPE == _maxpy)
				BMD.maxpyin(Cp,Ap,Bp);
#endif
#ifdef _LB_MM_TIMING
			matmul.stop();
			this->chrono+=matmul;
#endif
#if 0
			if (Ap.rowdim() <= 20 && Ap.coldim() <= 20) {
				Integer chara;
				F.characteristic(chara);
				F.write(cout) << endl;
				cout << "p:=" << chara << ';' << std::endl;
				A.write(cout<< "A:=",true) << ';' << std::endl;
				Ap.write(cout << "Ap:=", F, true) << ';' << endl;
				Bp.write(cout << "Bp:=", F, true) << ';' << endl;
				Cp.write(cout<< "Cp:=", F, true) << ';' << endl;
			}
#endif
			return Cp;
		}



	};
} // Protected
} // LinBox

namespace LinBox { namespace BLAS2 {

	template<class _anyVector>
	_anyVector & mul (_anyVector& C,
			  const SparseMatrix<typename _anyVector::Field, SparseMatrixFormat::CSR> & A,
			  const _anyVector& B,
			  const BLAS3::mulMethod::CRA &)
	{

		size_t PrimeSize = 22; //! @todo pourqoi ?

		integer mA, mB ;
		// MatrixDomain<typename _anyMatrix::Field> MD(A.field());
		// VectorDomain<typename _anyMatrix::Field> VD(A.field());
		// MD.Magnitude(mA,A);
		mA = A.magnitude();
		// VD.Magnitude(mB,B);
		mB = B.magnitude();
		integer cA = (integer) A.maxrow();
		double logC = Givaro::naturallog(mA*mB*cA);

		typedef Modular<double> ModularField ;
		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);


		{

			RandomPrimeIterator genprime( (unsigned int)PrimeSize );
			ChineseRemainder< FullMultipBlasMatCRA< ModularField > > cra( std::pair<size_t,double>(C.size(), logC) );
			Protected::IntegerSparseCraMatMul iteration(A,B);

			cra(C, iteration, genprime);

#ifdef _LB_DEBUG
#ifdef _LB_MM_TIMING
			// report << "Sole modular matrix multiplications: " << iteration.chrono << std::endl;
#endif

			Integer mC;
			// VD.Magnitude(mC, C);
			mC = C.magnitude();
			report << "C max: " << logtwo(mC) <<  " (" << LinBox::naturallog(mC) << ')' << std::endl;
#endif

		}

		report << mA << ',' << mB << ',' <<  C.magnitude() << std::endl;

		return C;

	}
} // BLAS2
} // LinBox

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
	LinBox::commentator().start("toom-cook suite", "toom");
	ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);


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

		GF.write(report << "This is the field with " << (LinBox::Integer)pow((LinBox::Integer)p,e) << " elements: ") << ", using: "   << GF.irreducible() << " as irreducible polynomial" << std::endl;
		report << "matrices are " << m << 'x' << k << " and " << k << 'x' << n <<  std::endl;

		LinBox::BlasMatrix<GFpe> A(GF,m,k);
		LinBox::BlasMatrix<GFpe> B(GF,k,n);
		LinBox::BlasMatrix<GFpe> C(GF,m,n);

		typedef GFpe::RandIter Randiter;
		Randiter R(GF);
		LinBox::RandomDenseMatrix<Randiter,GFpe> randomizer(GF,R) ;
		randomizer.random(A);
		// report << "A[0,0] = " << A.getEntry(0,0) << std::endl;
		randomizer.random(B);

		report << "naive over GFq" << std::endl;
		Tim.clear(); Tim.start();
		LinBox::BLAS3::mul(C,A,B,LinBox::BLAS3::mulMethod::naive());
		Tim.stop();
		report << Tim << '(' << C.getEntry(0,0) << ')' << std::endl;

		// for (size_t i = 0 ; i< m ; ++i)
		// for (size_t j = 0 ; j< n ; ++j)
		// C.setEntry(i,j,GF.zero);

		// report << "A[0,0] = " << A.getEntry(0,0) << std::endl;
		{
			report << "ToomCook low mem" << std::endl;
			LinBox::BlasMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			LinBox::BLAS3::mul(D,A,B,LinBox::BLAS3::mulMethod::ToomCook<GFpe>(GF,false));
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "error" << std::endl;
				return 1;
			}
		}

		{
			report << "ToomCook high mem" << std::endl;
			LinBox::BlasMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			LinBox::BLAS3::mul(D,A,B,LinBox::BLAS3::mulMethod::ToomCook<GFpe>(GF,true));
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "error" << std::endl;
				return 1;
			}
		}

		{
			report << "Matrix Domain" << std::endl;
			LinBox::BlasMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			MD.mul(D,A,B);
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "error" << std::endl;
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

		report << "Naïve " << std::endl ;
		Tim.clear() ; Tim.start() ;
		LinBox::BLAS3::mul(C,A,B,LinBox::BLAS3::mulMethod::naive());
		Tim.stop();
		report << Tim << '(' << C.getEntry(0,0) << ')' << std::endl;

#ifdef __LINBOX_HAVE_FLINT
		{
			report << "FLINT " << std::endl;
			LinBox::BlasMatrix<LinBox::PID_integer> D(ZZ,m,n);
			Tim.clear(); Tim.start();
			LinBox::BLAS3::mul(D,A,B,LinBox::BLAS3::mulMethod::FLINT());
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "error" << std::endl;
				return 1;
			}
		}
#endif // __LINBOX_HAVE_FLINT

		{
			report << "Matrix Domain" << std::endl;
			LinBox::BlasMatrix<LinBox::PID_integer> D(ZZ,m,n);
			Tim.clear(); Tim.start();
			MD.mul(D,A,B);
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "error" << std::endl;
				return 1;
			}
		}

		{
			report << "CRA " << std::endl;
			LinBox::BlasMatrix<LinBox::PID_integer> D(ZZ,m,n);
			Tim.clear(); Tim.start();
			LinBox::BLAS3::mul(D,A,B,LinBox::BLAS3::mulMethod::CRA());
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				// report << D << std::endl;
				// report << C << std::endl;
				report << "error" << std::endl;
				return 1;
			}
		}
	}


	{ /* ZZ spmat mul */

		typedef LinBox::PID_integer Field;
		typedef LinBox::SparseMatrix<Field,LinBox::SparseMatrixFormat::CSR> BlackBox;
		Field ZZ ;

		LinBox::VectorDomain<LinBox::PID_integer> MD(ZZ);

		// typename Field::RandIter ri (ZZ,b);
		LinBox::RandomIntegerIter<false> ri((unsigned int)b);
		double sparsity = 0.05;
		LinBox::RandomSparseStream<Field, typename BlackBox::Row, LinBox::RandomIntegerIter<false> > stream (ZZ, ri, sparsity, k, m);

		BlackBox A (ZZ, stream);

		typedef LinBox::BlasVector<Field> Vector ;
		Vector x(ZZ,k);
		Vector y(ZZ,m);

		// size_t iter = 1 ;
		// LinBox::RandomDenseStream<Field, Vector> vs (ZZ, k, iter);
		x.random(ri);

		report << "Naïve " << std::endl ;
		Tim.clear() ; Tim.start() ;
		A.apply(y,x);
		Tim.stop();
		report << Tim << '(' << y[0] << ')' << std::endl;


		report << "CRA " << std::endl;
		{


			LinBox::BlasVector<LinBox::PID_integer> z(ZZ,m);
			Tim.clear(); Tim.start();
			LinBox::BLAS2::mul(z,A,x,LinBox::BLAS3::mulMethod::CRA());
			Tim.stop();
			report << Tim << '(' << z[0] << ')' << std::endl;

			if (!MD.areEqual(y,z)) {
				// report << D << std::endl;
				// report << C << std::endl;
				report << "error" << std::endl;
				return 1;
			}
		}
	}
	LinBox::commentator().stop("toom-cook suite");

	return 0;
}
