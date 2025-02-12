/*
 * Copyright (C) 2012 the LinBox group
 *
 * written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

#include <linbox/linbox-config.h>

#include <iostream>
#include "linbox/integer.h"
#include "linbox/ring/modular.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/solutions/hadamard-bound.h"
// #include <fflas-ffpack/fflas/fflas.h>
#include "givaro/modular.h"
#include "linbox/util/timer.h"

#include "linbox/algorithms/matrix-blas3/mul.h"

#include "test-common.h"



namespace LinBox { namespace Protected {
	struct IntegerSparseCraMatMul {


		typedef SparseMatrixFormat::CSR spfmt ;
		typedef Givaro::Modular<double>         Field;
		typedef Field::Element          Element;
		typedef SparseMatrix<Field,spfmt>       ModularMatrix ;
		typedef DenseVector<Field>               ModularVector;
		typedef DenseVector<Givaro::ZRing<Integer> >         IntegerVector;
		typedef SparseMatrix<Givaro::ZRing<Integer>,spfmt> IntegerMatrix ;

#ifdef _LB_MM_TIMING
#ifdef __LINBOX_USE_OPENMP
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
			// linbox_check(A.getPointer() == _A_.getPointer());
		}

		IntegerSparseCraMatMul(IntegerMatrix& A, IntegerVector& B) :
			_A_(A), _B_(B)
		{
#ifdef _LB_MM_TIMING
			chrono.clear();
#endif
			// linbox_check(A.getPointer() == _A_.getPointer());
		}

		IterationResult operator()(ModularVector& Cp, const Field& F) const
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
                        return IterationResult::CONTINUE;
		}



	};
} // Protected
} // LinBox


namespace LinBox {

template<typename Container>
Integer& magnitude(Integer& max_elt, const Container& v) {
    max_elt = integer(0);
    for(const auto& iter: v)
			if (max_elt < Givaro::abs(iter))
				max_elt = Givaro::abs(iter) ;
    return max_elt;
}



namespace BLAS2 {

	template<class _anyVector>
	_anyVector & mul (_anyVector& C,
			  const SparseMatrix<typename _anyVector::Field, SparseMatrixFormat::CSR> & A,
			  const _anyVector& B,
			  const BLAS3::mulMethod::CRA &)
	{


		integer mA, mB, mC ;
        mA = A.magnitude();
		magnitude(mB, B);
		integer cA = uint64_t(A.maxrow());
		double logC = Givaro::naturallog(mA*mB*cA);

		typedef Givaro::Modular<double> ModularField ;
		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);


		{

            PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<ModularField>::bestBitSize(A.coldim()));
            ChineseRemainder< CRABuilderFullMultipMatrix< ModularField > > cra( std::pair<size_t,double>(C.size(), logC) );
            Protected::IntegerSparseCraMatMul iteration(A,B);

            cra(C, iteration, genprime);


		}

        magnitude(mC, C);
#ifdef _LB_DEBUG
            report << "C max: " << logtwo(mC) <<  " (" << LinBox::naturallog(mC) << ')' << std::endl;
#endif

		report << mA << ',' << mB << ',' << mC << std::endl;

		return C;

	}
} // BLAS2
} // LinBox

using namespace LinBox;
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


	parseArguments (ac, av, as);
	commentator().start("toom-cook suite", "toom");
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);


	Timer Tim ;
	{ /* Toom Cook over GivarorExtension */
		//typedef Givaro::Modular<int64_t> Zpz;
		typedef Givaro::Modular<double> Zpz;
		typedef Givaro::Extension< Zpz > GFpe ;

		// Z/pZ
		Zpz F(p);
		// GF(p^e) ;
		GFpe GF (F, 2);
		MatrixDomain<GFpe> MD(GF);

		GF.write(report << "This is the field with " << (Integer)pow((Integer)p,e) << " elements: ") << ", using: "   << GF.irreducible() << " as irreducible polynomial" << std::endl;
		report << "matrices are " << m << 'x' << k << " and " << k << 'x' << n <<  std::endl;

		DenseMatrix<GFpe> A(GF,m,k);
		DenseMatrix<GFpe> B(GF,k,n);
		DenseMatrix<GFpe> C(GF,m,n);

		typedef GFpe::RandIter Randiter;
		Randiter R(GF);
		RandomDenseMatrix<Randiter,GFpe> randomizer(GF,R) ;
		randomizer.random(A);
		// report << "A[0,0] = " << A.getEntry(0,0) << std::endl;
		randomizer.random(B);

		report << "naive over GFq" << std::endl;
		Tim.clear(); Tim.start();
		BLAS3::mul(C,A,B,BLAS3::mulMethod::naive());
		Tim.stop();
		report << Tim << '(' << C.getEntry(0,0) << ')' << std::endl;

		{
			report << "ToomCook low mem" << std::endl;
			DenseMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			BLAS3::mul(D,A,B,BLAS3::mulMethod::ToomCook<GFpe>(GF,false));
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "low mem error" << std::endl;
//				return 1;
			}
		}

		{
			report << "ToomCook high mem" << std::endl;
			DenseMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			BLAS3::mul(D,A,B,BLAS3::mulMethod::ToomCook<GFpe>(GF,true));
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "high mem error" << std::endl;
//				return 1;
			}
		}

		{
			report << "Matrix Domain" << std::endl;
			DenseMatrix<GFpe> D(GF,m,n);
			Tim.clear(); Tim.start();
			MD.mul(D,A,B);
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "Matrix Domain error" << std::endl;
//				return 1;
			}
		}
	}

	{ /* ZZ mat mul */

            Givaro::ZRing<Integer> ZZ ;
		MatrixDomain<Givaro::ZRing<Integer> > MD(ZZ);
		DenseMatrix<Givaro::ZRing<Integer> > A(ZZ,m,k) ;
		DenseMatrix<Givaro::ZRing<Integer> > B(ZZ,k,n) ;
		DenseMatrix<Givaro::ZRing<Integer> > C(ZZ,m,n) ;

		//A.random((unsigned)b);
		//B.random((unsigned)b);
		A.random(); // BUG: b is ignored but was already the case before
		B.random();

		report << "Naïve " << std::endl ;
		Tim.clear() ; Tim.start() ;
		BLAS3::mul(C,A,B,BLAS3::mulMethod::naive());
		Tim.stop();
		report << Tim << '(' << C.getEntry(0,0) << ')' << std::endl;

#ifdef __LINBOX_HAVE_FLINT
		{
			report << "FLINT " << std::endl;
			DenseMatrix<Givaro::ZRing<Integer> > D(ZZ,m,n);
			Tim.clear(); Tim.start();
			BLAS3::mul(D,A,B,BLAS3::mulMethod::FLINT());
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "FLINT error" << std::endl;
				return 1;
			}
		}
#endif // __LINBOX_HAVE_FLINT

		{
			report << "Matrix Domain" << std::endl;
			DenseMatrix<Givaro::ZRing<Integer> > D(ZZ,m,n);
			Tim.clear(); Tim.start();
			MD.mul(D,A,B);
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				report << "Matrix Domain error" << std::endl;
//				return 1;
			}
		}

		{
			report << "CRA " << std::endl;
			DenseMatrix<Givaro::ZRing<Integer> > D(ZZ,m,n);
			Tim.clear(); Tim.start();
			BLAS3::mul(D,A,B,BLAS3::mulMethod::CRA());
			Tim.stop();
			report << Tim << '(' << D.getEntry(0,0) << ')' << std::endl;

			if (!MD.areEqual(D,C)) {
				// report << D << std::endl;
				// report << C << std::endl;
				report << "CRA error" << std::endl;
//				return 1;
			}
		}
	}

	{ /* ZZ spmat mul */

		typedef Givaro::ZRing<Integer> Field;
		typedef SparseMatrix<Field,SparseMatrixFormat::CSR> BlackBox;
		Field ZZ ;

		VectorDomain<Givaro::ZRing<Integer> > MD(ZZ);

		// typename Field::RandIter ri (ZZ,b);
                Givaro::RandomIntegerIterator<false> ri(ZZ,(size_t)b);
		double sparsity = 0.05;
		RandomSparseStream<Field, typename BlackBox::Row, Givaro::RandomIntegerIterator<false> > stream (ZZ, ri, sparsity, k, m);

		BlackBox A (ZZ, stream);

		typedef DenseVector<Field> Vector ;
		Vector x(ZZ,k);
		Vector y(ZZ,m);

		// size_t iter = 1 ;
		// RandomDenseStream<Field, Vector> vs (ZZ, ri, k, iter);
		x.random(ri);

		report << "Naïve " << std::endl ;
		Tim.clear() ; Tim.start() ;
		A.apply(y,x);
		Tim.stop();
		report << Tim << '(' << y[0] << ')' << std::endl;


		report << "CRA " << std::endl;
		{


			DenseVector<Givaro::ZRing<Integer> > z(ZZ,m);
			Tim.clear(); Tim.start();
			BLAS2::mul(z,A,x,BLAS3::mulMethod::CRA());
			Tim.stop();
			report << Tim << '(' << z[0] << ')' << std::endl;

			if (!MD.areEqual(y,z)) {
				// report << D << std::endl;
				// report << C << std::endl;
				report << "CRA error" << std::endl;
//				return 1;
			}
		}
	}
	commentator().stop("toom-cook suite");

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
