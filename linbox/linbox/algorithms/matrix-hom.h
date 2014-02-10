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
#include "linbox/field/hom.h"
#include "linbox/matrix/matrix-category.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/polynomial.h"
#include "linbox/blackbox/scalar-matrix.h"

// namespace LinBox {
	// template<class A, class B>
	// class SparseMatrix ;

	// template<class A, class R>
	// class BlasMatrix;
// } // LinBox

namespace LinBox
{

	/// \brief Limited doc so far. Used in RationalSolver.
	namespace MatrixHom
	{

		template<class Field, class Vect>
		struct SparseVectorTranslate {
			typedef Vect other_t;
		};

		template<class Field>
		struct SparseVectorTranslate<Field,SparseMatrixFormat::SparseSeq> {
			typedef typename Vector<Field>::SparseSeq other_t;
		};

		template<class Field>
		struct SparseVectorTranslate<Field,SparseMatrixFormat::SparsePar> {
			typedef typename Vector<Field>::SparsePar other_t;
		};
		template<class Field>
		struct SparseVectorTranslate<Field,SparseMatrixFormat::SparseMap> {
			typedef typename Vector<Field>::SparseMap other_t;
		};



		//public:

		template<class FMatrix, class IMatrix>
		void map (FMatrix & Ap, const IMatrix& A)
		{
			typename IMatrix::template rebind<typename FMatrix::Field>()( Ap, A);
		}

		// construct a sparse matrix over finite field, such that Ap = A mod p, where F = Ring / <p>
		template<class Field, class Vect, class IMatrix>
		void map (SparseMatrix<Field, Vect> &Ap, const IMatrix& A);

		// construct a sparse matrix over finite field, such that Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Vect1, class Field, class Vect2>
		void map (SparseMatrix<Field, Vect2>& Ap, const SparseMatrix<Ring, Vect1>& A)
		{
			// typedef typename SparseVectorTranslate<Field,Vect1>::other_t Vect_1 ;
			typedef typename SparseVectorTranslate<Field,Vect2>::other_t Vect_2 ;
			typename SparseMatrix<Ring,Vect1>::template rebind<Field,Vect_2>()( Ap, A);
		}



		// function class to hanle map to BlasMatrix (needed to allow partial specialization)
		template< class Field, class IMatrix, class Type>
		class BlasMatrixMAP {
		public:
			template<class _Rep>
			void operator() (BlasMatrix<Field,_Rep> &Ap, const IMatrix& A, Type type);
		};

		// construct a BlasMatrix over finite fiel, such that Ap - A mod p, where F = Ring / <p>

		template<class Ring, class Field, class _Rep>
		void map (BlasMatrix<Field,_Rep> &Ap, const BlasMatrix<Ring,_Rep>& A )
		{
			typename BlasMatrix<Ring,_Rep>::template rebind<Field>()( Ap, A);
		}


		template <class Field, class IMatrix, class _Rep>
		void map (BlasMatrix<Field,_Rep> &Ap, const IMatrix &A)
		{
			BlasMatrixMAP<Field, IMatrix, typename MatrixContainerTrait<IMatrix>::Type> ()(Ap, A, typename MatrixContainerTrait<IMatrix>::Type());
		}

		template <class Field, class IPoly, class IMatrix>
		void map (PolynomialBB< typename IMatrix::template rebind<Field>::other,
			 typename IPoly::template rebind<Field>::other> &Ap,
			  const PolynomialBB<IMatrix, IPoly> &A)
		{
			typename PolynomialBB<IMatrix,IPoly>::template rebind<Field>() (Ap, A);
		}

		template <class Field, class Ring>
		void map (ScalarMatrix<Field> &Ap,
			  const ScalarMatrix<Ring> &A)
		{
			typename ScalarMatrix<Ring>::template rebind<Field>() (Ap, A);
		}

	}

	template <class Field, class Vect, class IMatrix>
	void MatrixHom::map (SparseMatrix<Field, Vect> &Ap, const IMatrix& A)
	{

		typedef typename IMatrix::Field Ring;

		Ring r = A.field();


		std::vector<typename Ring::Element> e(A.coldim(), r.zero), tmp(A.rowdim());

		typename std::vector<typename Ring::Element>::iterator iter, e_p;

		typename Field::Element val;

		int i = 0;

		Hom<Ring, Field> hom(A. field(), Ap.field());

		for (e_p=e.begin();e_p != e.end(); ++e_p,++i){
			r.assign(*e_p, r.one);
			A.apply(tmp,e);
			int j;
			for (iter=tmp.begin(),j=0; iter != tmp.end(); ++iter,++j) {
				hom. image (val, *iter);
				if (!Ap.field().isZero(val))
					Ap.setEntry ((size_t)j,(size_t)i, val);

			}
			r.assign(*e_p, r.zero);
		}

	}

	namespace MatrixHom
	{

		template<class Field, class IMatrix>
		class BlasMatrixMAP<Field, IMatrix, MatrixContainerCategory::Blackbox> {
		public:
			template<class _Rep>
			void operator() (BlasMatrix<Field,_Rep> &Ap, const IMatrix &A, MatrixContainerCategory::Blackbox type)
			{


				typedef typename IMatrix::Field Ring;
				Ring r = A.field();


				std::vector<typename Ring::Element> e(A.coldim(), r.zero), tmp(A.rowdim());

				typename BlasMatrix<Field,_Rep>::ColIterator col_p;

				typename BlasMatrix<Field,_Rep>::Col::iterator elt_p;

				typename std::vector<typename Ring::Element>::iterator e_p, tmp_p;

				Hom<Ring, Field> hom(A. field(), Ap.field());

				for (col_p = Ap.colBegin(), e_p = e.begin();
				     e_p != e.end(); ++ col_p, ++ e_p) {

					r.assign(*e_p, r.one);
					A.apply (tmp, e);

					for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
					     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)
						hom.image (*elt_p, *tmp_p);
					r.assign(*e_p, r.zero);
				}
			}
		};


		template<class Field, class IMatrix>
		class BlasMatrixMAP<Field, IMatrix, MatrixContainerCategory::Container> {
		public:
			template<class _Rep>
			void operator() (BlasMatrix<Field,_Rep> &Ap, const IMatrix &A,  MatrixContainerCategory::Container type)
			{
				Hom<typename IMatrix::Field , Field> hom(A.field(), Ap.field());
				typename Field::Element e;
				for( typename IMatrix::ConstIndexedIterator indices = A.IndexedBegin();
				     (indices != A.IndexedEnd()) ;
				     ++indices ) {

					hom. image (e, A.getEntry(indices.rowIndex(),indices.colIndex()) );

					if (!Ap.field().isZero(e))
						Ap.setEntry (indices.rowIndex(),
							     indices.colIndex(), e);
					else
						Ap.setEntry (indices.rowIndex(),
							     indices.colIndex(), Ap.field().zero);
				}

			}
		};

		template<class Field, class IMatrix>
		class BlasMatrixMAP<Field, IMatrix, MatrixContainerCategory::BlasContainer> {
		public:
			template<class _Rep>
			void operator() (BlasMatrix<Field,_Rep> &Ap, const IMatrix &A, MatrixContainerCategory::BlasContainer type)
			{
				Hom<typename IMatrix::Field , Field> hom(A.field(), Ap.field());

				typename IMatrix::ConstIterator        iterA  = A.Begin();
				typename BlasMatrix<Field,_Rep>::Iterator iterAp = Ap.Begin();

				for(; iterA != A.End(); ++iterA, ++iterAp)
					hom. image (*iterAp, *iterA);
			}
		};

		template<class Field, class _Rep>
		class BlasMatrixMAP<Field, BlasMatrix<PID_integer,_Rep>, MatrixContainerCategory::BlasContainer> {
		public:
			template<class _Rep2>
			void operator() (BlasMatrix<Field,_Rep2> &Ap, const BlasMatrix<PID_integer,_Rep> &A, MatrixContainerCategory::BlasContainer type)
			{
				PID_integer ZZ ;
				Hom<PID_integer , Field> hom(ZZ, Ap.field());

				typename BlasMatrix<PID_integer,_Rep>::ConstIterator        iterA  = A.Begin();
				typename BlasMatrix<Field,_Rep2>::Iterator iterAp = Ap.Begin();

				for(; iterA != A.End(); ++iterA, ++iterAp)
					hom. image (*iterAp, *iterA);
			}
		};

#ifdef __LINBOX_blas_matrix_multimod_H
		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::BlasContainer > {
		public:
			template<class _Rep>
			void operator() (BlasMatrix<MultiModDouble,_Rep> &Ap, const IMatrix &A,  MatrixContainerCategory::BlasContainer type)
			{
				for (size_t i=0; i<Ap.field().size();++i)
					MatrixHom::map(Ap.getMatrix(i), A);
			}
		};

		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::Container > {
		public:
			template<class _Rep>
			void operator() (BlasMatrix<MultiModDouble,_Rep> &Ap, const IMatrix &A,  MatrixContainerCategory::Container type)
			{
				for (size_t i=0; i<Ap.field().size();++i)
					MatrixHom::map(Ap.getMatrix(i), A);
			}
		};

		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::Blackbox > {
		public:
			template<class _Rep>
			void operator() (BlasMatrix<MultiModDouble,_Rep> &Ap, const IMatrix &A,  MatrixContainerCategory::Blackbox type)
			{
				for (size_t i=0; i<Ap.field().size();++i)
					MatrixHom::map(Ap.getMatrix(i), A);
			}
		};
#endif
	}
} // LinBox

#endif //__LINBOX_matrix_hom_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

