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
#include "linbox/matrix/sparse.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/polynomial.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/field/hom.h"
#include "linbox/matrix/matrix-category.h"

namespace LinBox {
	template<class A, class B, class C>
	class SparseMatrix ;

	template<class A, class R>
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


	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix2<Ring, SparseMatrixFormat::SparseSeq>, Field> {
		typedef SparseMatrix2<Field, SparseMatrixFormat::SparseSeq> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix2<Ring, SparseMatrixFormat::SparsePar>, Field> {
		typedef SparseMatrix2<Field, SparseMatrixFormat::SparsePar> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<SparseMatrix2<Ring, SparseMatrixFormat::SparseMap>, Field> {
		typedef SparseMatrix2<Field, SparseMatrixFormat::SparseMap> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<BlasMatrix<Ring>, Field> {
		typedef BlasMatrix<Field> value_type;
	};

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
		void map (SparseMatrix2<Field, Vect> &Ap, const IMatrix& A);

		// construct a sparse matrix over finite field, such that Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Vect1, class Field, class Vect2>
		void map (SparseMatrix2<Field, Vect2>& Ap, const SparseMatrix2<Ring, Vect1>& A)
		{
			// typedef typename SparseVectorTranslate<Field,Vect1>::other_t Vect_1 ;
			typedef typename SparseVectorTranslate<Field,Vect2>::other_t Vect_2 ;
			typename SparseMatrix2<Ring,Vect1>::template rebind<Field,Vect_2>()( Ap, A);
		}



		// function class to hanle map to BlasMatrix (needed to allow partial specialization)
		template< class Field, class IMatrix, class Type>
		class BlasMatrixMAP {
		public:
			void operator() (BlasMatrix<Field> &Ap, const IMatrix& A, Type type);
		};

		// construct a BlasMatrix over finite fiel, such that Ap - A mod p, where F = Ring / <p>

		template<class Ring, class Field>
		void map (BlasMatrix<Field> &Ap, const BlasMatrix<Ring>& A )
		{
			typename BlasMatrix<Ring>::template rebind<Field>()( Ap, A);
		}


		template <class Field, class IMatrix>
		void map (BlasMatrix<Field> &Ap, const IMatrix &A)
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
	void MatrixHom::map (SparseMatrix2<Field, Vect> &Ap, const IMatrix& A)
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
			void operator() (BlasMatrix<Field> &Ap, const IMatrix &A, MatrixContainerCategory::Blackbox type)
			{


				typedef typename IMatrix::Field Ring;
				Ring r = A.field();


				std::vector<typename Ring::Element> e(A.coldim(), r.zero), tmp(A.rowdim());

				typename BlasMatrix<Field>::ColIterator col_p;

				typename BlasMatrix<Field>::Col::iterator elt_p;

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
			void operator() (BlasMatrix<Field> &Ap, const IMatrix &A,  MatrixContainerCategory::Container type)
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
			void operator() (BlasMatrix<Field> &Ap, const IMatrix &A, MatrixContainerCategory::BlasContainer type)
			{
				Hom<typename IMatrix::Field , Field> hom(A.field(), Ap.field());

				typename IMatrix::ConstIterator        iterA  = A.Begin();
				typename BlasMatrix<Field>::Iterator iterAp = Ap.Begin();

				for(; iterA != A.End(); ++iterA, ++iterAp)
					hom. image (*iterAp, *iterA);
			}
		};

		template<class Field>
		class BlasMatrixMAP<Field, BlasMatrix<PID_integer>, MatrixContainerCategory::BlasContainer> {
		public:
			void operator() (BlasMatrix<Field> &Ap, const BlasMatrix<PID_integer> &A, MatrixContainerCategory::BlasContainer type)
			{
				PID_integer ZZ ;
				Hom<PID_integer , Field> hom(ZZ, Ap.field());

				typename BlasMatrix<PID_integer>::ConstIterator        iterA  = A.Begin();
				typename BlasMatrix<Field>::Iterator iterAp = Ap.Begin();

				for(; iterA != A.End(); ++iterA, ++iterAp)
					hom. image (*iterAp, *iterA);
			}
		};

#ifdef __LINBOX_blas_matrix_multimod_H
		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::BlasContainer > {
		public:
			void operator() (BlasMatrix<MultiModDouble> &Ap, const IMatrix &A,  MatrixContainerCategory::BlasContainer type)
			{
				for (size_t i=0; i<Ap.field().size();++i)
					MatrixHom::map(Ap.getMatrix(i), A);
			}
		};

		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::Container > {
		public:
			void operator() (BlasMatrix<MultiModDouble> &Ap, const IMatrix &A,  MatrixContainerCategory::Container type)
			{
				for (size_t i=0; i<Ap.field().size();++i)
					MatrixHom::map(Ap.getMatrix(i), A);
			}
		};

		template< class IMatrix>
		class BlasMatrixMAP<MultiModDouble, IMatrix, MatrixContainerCategory::Blackbox > {
		public:
			void operator() (BlasMatrix<MultiModDouble> &Ap, const IMatrix &A,  MatrixContainerCategory::Blackbox type)
			{
				for (size_t i=0; i<Ap.field().size();++i)
					MatrixHom::map(Ap.getMatrix(i), A);
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

