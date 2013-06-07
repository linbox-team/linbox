/* Copyright (C) 2010 LinBox
 * Written by JG Dumas
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file algorithms/matrix-hom.h
 * @ingroup algorithms
 * @brief Matrix Homomorphism
 * A map function  converts a matrix on a field/ring
 * to its natural image in another field/ring.
 * @sa
 * \c rebind operator.
 */

#ifndef __LINBOX_matrix_hom_H
#define __LINBOX_matrix_hom_H

//! @bug it is dangerous to include matrices defs that include hom for their rebind...
#include "linbox/integer.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/polynomial.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/field/hom.h"
#include "linbox/matrix/matrix-category.h"

namespace LinBox {
	template<class A, class B, class C>
	class SparseMatrix ;

	template<class A>
	class BlasMatrix;
}

namespace LinBox
{

	// try to map a blackbox over a homorphic ring
	// The most suitable type
	template <class Blackbox, class Field>
	struct MatrixHomTrait {
		//typedef ... FBlackbox
		// donot know
		typedef Blackbox value_type;
	};

#if 0
	// Vector<Ring>::Sparse = Vector<Ring>::SparsePar
	template <class Ring, class Field>
	struct Convert<SparseMatrix<Ring, typename Vector<Ring>::Sparse>, Field> {
		typedef SparseMatrix<Field, typename Vector<Field>::Sparse> value_type;
	};
#endif

	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix<Ring, typename Vector<Ring>::SparseSeq>, Field> {
		typedef SparseMatrix<Field, typename Vector<Field>::SparseSeq> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix<Ring, typename Vector<Ring>::SparsePar>, Field> {
		typedef SparseMatrix<Field, typename Vector<Field>::SparsePar> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix<Ring, typename Vector<Ring>::SparseMap>, Field> {
		typedef SparseMatrix<Field, typename Vector<Field>::SparseMap> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<BlasMatrix<Ring>, Field> {
		typedef BlasMatrix<Field> value_type;
	};

	/// \brief Limited doc so far. Used in RationalSolver.
	namespace MatrixHom
	{

		//public:

		template<class FMatrix, class IMatrix>
		void map (FMatrix & Ap, const IMatrix& A)
		{
			typename IMatrix::template rebind<typename FMatrix::Field>()( Ap, A);
		}

		// construct a sparse matrix over finite field, such that Ap = A mod p, where F = Ring / <p>
		template<class Field, class Vect, class IMatrix>
		void map (SparseMatrix<Field, Vect> &Ap, const IMatrix& A, const Field &F);

		// construct a sparse matrix over finite field, such that Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Vect1, class Field, class Vect2>
		void map (SparseMatrix<Field, Vect2>& Ap, const SparseMatrix<Ring, Vect1>& A)
		{
			typename SparseMatrix<Ring,Vect1>::template rebind<Field,Vect2>()( Ap, A);
		}



		// function class to hanle map to BlasMatrix (needed to allow partial specialization)
		template< class Field, class IMatrix, class Type>
		class BlasMatrixMAP {
		public:
			void operator() (BlasMatrix<Field> &Ap, const IMatrix& A, const Field& F, Type type);
		};

		// construct a BlasMatrix over finite fiel, such that Ap - A mod p, where F = Ring / <p>

		template<class Ring, class Field>
		void map (BlasMatrix<Field> &Ap, const BlasMatrix<Ring>& A )
		{
			typename BlasMatrix<Ring>::template rebind<Field>()( Ap, A);
		}


		template <class Field, class IMatrix>
		void map (BlasMatrix<Field> &Ap, const IMatrix &A, const Field &F)
		{
			BlasMatrixMAP<Field, IMatrix, typename MatrixContainerTrait<IMatrix>::Type> ()(Ap, A, F, typename MatrixContainerTrait<IMatrix>::Type());
		}

		template <class Field, class IPoly, class IMatrix>
		void map (PolynomialBB< typename IMatrix::template rebind<Field>::other,
			 typename IPoly::template rebind<Field>::other> &Ap,
			  const PolynomialBB<IMatrix, IPoly> &A, const Field & F)
		{
			typename PolynomialBB<IMatrix,IPoly>::template rebind<Field>() (Ap, A, F);
		}

		template <class Field, class Ring>
		void map (ScalarMatrix<Field> &Ap,
			  const ScalarMatrix<Ring> &A)
		{
			typename ScalarMatrix<Ring>::template rebind<Field>() (Ap, A);
		}

	}

	template <class Field, class Vect, class IMatrix>
	void MatrixHom::map (SparseMatrix<Field, Vect> &Ap, const IMatrix& A, const Field &F)
	{

		// 		Ap = new SparseMatrix<Field, Vect>(F, A.rowdim(), A.coldim());

		typedef typename IMatrix::Field Ring;

		Ring r = A.field();

		typename Ring::Element one, zero;

		r. init(one, 1);

		r. init(zero, 0);

		std::vector<typename Ring::Element> e(A.coldim(), zero), tmp(A.rowdim());

		typename std::vector<typename Ring::Element>::iterator iter, e_p;

		typename Field::Element val;

		int i = 0;

		Hom<Ring, Field> hom(A. field(), F);

		for (e_p=e.begin();e_p != e.end(); ++e_p,i++){
			r.assign(*e_p, one);
			A.apply(tmp,e);
			int j;
			for (iter=tmp.begin(),j=0; iter != tmp.end(); ++iter,j++) {
				hom. image (val, *iter);
				if (!F.isZero(val))
					Ap.setEntry ((size_t)j,(size_t)i, val);

			}
			r.assign(*e_p, zero);
		}

	}

	namespace MatrixHom
	{

		template<class Field, class IMatrix>
		class BlasMatrixMAP<Field, IMatrix, MatrixContainerCategory::Blackbox> {
		public:
			void operator() (BlasMatrix<Field> &Ap, const IMatrix &A,  const Field &F, MatrixContainerCategory::Blackbox type)
			{

				// 				Ap = new BlasMatrix<Field>(F, A.rowdim(), A.coldim());

				typedef typename IMatrix::Field Ring;
				Ring r = A.field();

				typename Ring::Element one, zero;
				r. init(one, 1);
				r. init(zero, 0);

				std::vector<typename Ring::Element> e(A.coldim(), zero), tmp(A.rowdim());

				typename BlasMatrix<Field>::ColIterator col_p;

				typename BlasMatrix<Field>::Col::iterator elt_p;

				typename std::vector<typename Ring::Element>::iterator e_p, tmp_p;

				Hom<Ring, Field> hom(A. field(), F);

				for (col_p = Ap.colBegin(), e_p = e.begin();
				     e_p != e.end(); ++ col_p, ++ e_p) {

					r.assign(*e_p, one);
					A.apply (tmp, e);

					for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
					     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)
						hom.image (*elt_p, *tmp_p);
					r.assign(*e_p, zero);
				}
			}
		};


		template<class Field, class IMatrix>
		class BlasMatrixMAP<Field, IMatrix, MatrixContainerCategory::Container> {
		public:
			void operator() (BlasMatrix<Field> &Ap, const IMatrix &A, const Field &F, MatrixContainerCategory::Container type)
			{
				// 				Ap = new BlasMatrix<Field>(F, A.rowdim(), A.coldim());
				Hom<typename IMatrix::Field , Field> hom(A.field(), F);
				typename Field::Element e, zero;
				F.init(zero,0UL);
				for( typename IMatrix::ConstIndexedIterator indices = A.IndexedBegin();
				     (indices != A.IndexedEnd()) ;
				     ++indices ) {

					hom. image (e, A.getEntry(indices.rowIndex(),indices.colIndex()) );

					if (!F.isZero(e))
						Ap.setEntry (indices.rowIndex(),
							     indices.colIndex(), e);
					else
						Ap.setEntry (indices.rowIndex(),
							     indices.colIndex(), zero);
				}

			}
		};

		template<class Field, class IMatrix>
		class BlasMatrixMAP<Field, IMatrix, MatrixContainerCategory::BlasContainer> {
		public:
			void operator() (BlasMatrix<Field> &Ap, const IMatrix &A, const Field &F, MatrixContainerCategory::BlasContainer type)
			{
				//Ap = new BlasMatrix<Field>(F, A.rowdim(), A.coldim());
				Hom<typename IMatrix::Field , Field> hom(A.field(), F);

				typename IMatrix::ConstIterator        iterA  = A.Begin();
				typename BlasMatrix<Field>::Iterator iterAp = Ap.Begin();

				for(; iterA != A.End(); iterA++, iterAp++)
					hom. image (*iterAp, *iterA);
			}
		};

		template<class Field>
		class BlasMatrixMAP<Field, BlasMatrix<PID_integer>, MatrixContainerCategory::BlasContainer> {
		public:
			void operator() (BlasMatrix<Field> &Ap, const BlasMatrix<PID_integer> &A, const Field &F, MatrixContainerCategory::BlasContainer type)
			{
				//Ap = new BlasMatrix<Field>(F, A.rowdim(), A.coldim());
				PID_integer ZZ ;
				Hom<PID_integer , Field> hom(ZZ, F);

				typename BlasMatrix<PID_integer>::ConstIterator        iterA  = A.Begin();
				typename BlasMatrix<Field>::Iterator iterAp = Ap.Begin();

				for(; iterA != A.End(); iterA++, iterAp++)
					hom. image (*iterAp, *iterA);
			}
		};

#ifdef __LINBOX_blas_matrix_multimod_H
		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::BlasContainer > {
		public:
			void operator() (BlasMatrix<MultiModDouble> &Ap, const IMatrix &A, const MultiModDouble &F,  MatrixContainerCategory::BlasContainer type)
			{
				// Ap = new BlasMatrix<MultiModDouble>(F, A.rowdim(), A.coldim());
				for (size_t i=0; i<F.size();++i)
					MatrixHom::map(Ap.getMatrix(i), A, F.getBase(i));
			}
		};

		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::Container > {
		public:
			void operator() (BlasMatrix<MultiModDouble> &Ap, const IMatrix &A, const MultiModDouble &F,  MatrixContainerCategory::Container type)
			{
				// 				Ap = new BlasMatrix<MultiModDouble>(F, A.rowdim(), A.coldim());
				for (size_t i=0; i<F.size();++i)
					MatrixHom::map(Ap.getMatrix(i), A, F.getBase(i));
			}
		};

		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::Blackbox > {
		public:
			void operator() (BlasMatrix<MultiModDouble> &Ap, const IMatrix &A, const MultiModDouble &F,  MatrixContainerCategory::Blackbox type)
			{
				// 				Ap = new BlasMatrix<MultiModDouble>(F, A.rowdim(), A.coldim());
				for (size_t i=0; i<F.size();++i)
					MatrixHom::map(Ap.getMatrix(i), A, F.getBase(i));
			}
		};
#endif
	}
}

#endif //__LINBOX_matrix_hom_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

