/*  Copyright (C) 2014 the members of the LinBox group
 *
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 * This file is part of the LinBox library.
 *
 * ========LICENCE========
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * LinBox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

#ifndef __LINBOX_matrix_blas3_mul_cra_INL
#define __LINBOX_matrix_blas3_mul_cra_INL

#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-builder-full-multip-fixed.h"

#include "givaro/random-integer.h"
#include "linbox/randiter/random-prime.h"


namespace LinBox { namespace BLAS3 { namespace Protected {

	struct IntegerCraMatMul {


		typedef Givaro::Modular<double>     Field;
		typedef Field::Element      Element;
		typedef BlasMatrix<Field>   ModularMatrix ;
		typedef BlasMatrix<Givaro::ZRing<Integer> > IntegerMatrix ;

#ifdef _LB_MM_TIMING
#ifdef _OPENMP
		typedef LinBox::OMPTimer Mytime;
#else
		typedef LinBox::Timer    Mytime;
#endif
#endif

		const IntegerMatrix &_A_, &_B_;

#ifdef _LB_MM_TIMING
		mutable Mytime chrono;
#endif

		IntegerCraMatMul(const IntegerMatrix& A, const IntegerMatrix& B) :
			_A_(A), _B_(B)
		{
#ifdef _LB_MM_TIMING
			chrono.clear();
#endif
			linbox_check(A.getPointer() == _A_.getPointer());
		}

		IntegerCraMatMul(IntegerMatrix& A, IntegerMatrix& B) :
			_A_(A), _B_(B)
		{
#ifdef _LB_MM_TIMING
			chrono.clear();
#endif
			linbox_check(A.getPointer() == _A_.getPointer());
		}

		IterationResult operator()(ModularMatrix& Cp, const Field& F) const
		{
			BlasMatrixDomain<Field>   BMD(F);

			/*  intialisation */
			// ModularMatrix Cpp(_A_.rowdim(),_B_.coldim());
			// Cp = Cpp ;
			ModularMatrix Ap(_A_, F);
			ModularMatrix Bp(_B_, F);
			Cp.resize(Ap.rowdim(),Bp.coldim());

			/*  multiplication mod p */

#ifdef _LB_MM_TIMING
			Mytime matmul; matmul.clear(); matmul.start();
#endif
			BMD.mul(Cp,Ap,Bp);
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
} // BLAS3
} // LinBox

namespace LinBox { namespace BLAS3 {
	template<class _anyMatrix>
	_anyMatrix & mul (_anyMatrix& C,
			  const _anyMatrix& A,
			  const _anyMatrix& B,
			  const mulMethod::CRA &)
	{

		integer mA, mB ;
		BlasMatrixDomain<typename _anyMatrix::Field> BMD(A.field());
		BMD.Magnitude(mA,A);
		BMD.Magnitude(mB,B);
		double logC = Givaro::naturallog(mA*mB*uint64_t(A.coldim()));

		typedef Givaro::Modular<double> ModularField ;

		{

                        PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<ModularField>::bestBitSize(A.coldim()));
			ChineseRemainder< CRABuilderFullMultipMatrix< ModularField > > cra( std::pair<size_t,double>(C.rowdim()*C.coldim(), logC) );
			Protected::IntegerCraMatMul iteration(A,B);

			cra(C, iteration, genprime);

#ifdef _LB_DEBUG
#ifdef _LB_MM_TIMING
			std::cout << "Sole modular matrix multiplications: " << iteration.chrono << std::endl;
#endif

			Integer mC; BMD.Magnitude(mC, C);
			std::cout << "C max: " << logtwo(mC) <<  " (" << LinBox::naturallog(mC) << ')' << std::endl;
#endif

		}

		return C;

	}
} // BLAS3
} // LinBox

#endif // __LINBOX_matrix_blas3_mul_cra_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
