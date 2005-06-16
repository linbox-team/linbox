/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// Time-stamp: <14 Jun 05 18:59:39 Jean-Guillaume.Dumas@imag.fr> 
#ifndef __LINBOX_MATRIX_HOM_H__
#define __LINBOX_MATRIX_HOM_H__

#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/compose.h>
#include <linbox/integer.h>
#include <linbox/field/hom.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/matrix-category.h>


namespace LinBox {

	// try to map a blackbox over a homorphic ring
	// The most suitable type
	template <class Blackbox, class Field>
	struct MatrixHomTrait{
		//typedef ... FBlackbox
		// donot know
		typedef Blackbox value_type;
	};

	/*
	// Vector<Ring>::Sparse = Vector<Ring>::SparsePar
	template <class Ring, class Field>
	struct Convert<SparseMatrix<Ring, typename Vector<Ring>::Sparse>, Field> {
		typedef SparseMatrix<Field, typename Vector<Field>::Sparse> value_type;
	};
	*/

	template <class Ring, class Field>
	struct MatrixHomTrait<DenseMatrixBase<typename Ring::Element>, Field> {
		typedef DenseMatrixBase<typename Field::Element> value_type;
	};
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
	struct MatrixHomTrait<DenseMatrix<Ring>, Field> {
		typedef DenseMatrix<Field> value_type;
	};

	template <class Ring, class Field>
	struct MatrixHomTrait<BlasBlackbox<Ring>, Field> {
		typedef BlasBlackbox<Field> value_type;
	};

	/// \brief Limited doc so far. Used in RationalSolver.
	namespace MatrixHom {
		
	//public:
	
		//template<class FMatrix, class IMatrix, class Field>
		//void map (FMatrix* & Ap, const IMatrix& A, const Field& F) {
                //    typename IMatrix::template rebind<Field>()( Ap, A, F);
                // }
		
		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Field, class IMatrix>
		void map (DenseMatrix<Field>* &Ap, const IMatrix& A, const Field& F);
		
		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		void map (DenseMatrixBase<typename Field::Element>* &Ap, const DenseMatrixBase<typename Ring::Element>& A, const Field& F){
                    typename DenseMatrixBase<typename Ring::Element>::template rebind<Field>()( Ap, A, F);
                }
		
		// construct a sparse matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Field, class IMatrix>
		void map (SparseMatrix<Field>* &Ap, const IMatrix& A, const Field &F);

		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		void map (DenseMatrix<Field>* &Ap, const DenseMatrix<Ring>& A, const Field &F){
                    typename DenseMatrix<Ring>::template rebind<Field>()( Ap, A, F);
                }

		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		void map (DenseMatrix<Field>* &Ap, const SparseMatrix<Ring>& A, const Field &F);
		
		// construct a sparse matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Vect1, class Field, class Vect2>
		void map (SparseMatrix<Field, Vect2>*& Ap, const SparseMatrix<Ring, Vect1>& A, const Field& F) {
                    typename SparseMatrix<Ring,Vect1>::template rebind<Field,Vect2>()( Ap, A, F);
                }

		template <class Field, class IMatrix>
		void MatrixHom::map (BlasBlackbox<Field> *&Ap, const IMatrix &A, const Field &F) {
			MatrixHom::map(Ap, A, F, typename MatrixContainerTrait<IMatrix>::Type());
		}

		template <class Field, class IMatrix>
		void MatrixHom::map (BlasBlackbox<Field> *&Ap, const IMatrix &A, 
				     const Field &F, MatrixContainerCategory::Container type);

		template <class Field, class IMatrix>
		void MatrixHom::map (BlasBlackbox<Field> *&Ap, const IMatrix &A, 
				     const Field &F, MatrixContainerCategory::Blackbox type);

	}		

	template <class Field, class IMatrix>
	void MatrixHom::map (DenseMatrix<Field>* &Ap, const IMatrix& A, const Field &F) {

		Ap = new DenseMatrix<Field>(F, A.rowdim(), A.coldim());

		typedef typename IMatrix::Field Ring;

		Ring r = A.field();

		typename Ring::Element one, zero;

		r. init(one, 1);

		r. init(zero, 0);

		std::vector<typename Ring::Element> e(A.coldim(), zero), tmp(A.rowdim());

		typename DenseMatrix<Field>::ColIterator col_p;

		typename DenseMatrix<Field>::Col::iterator elt_p;

		typename std::vector<typename Ring::Element>::iterator e_p, tmp_p;

		Hom<Ring, Field> hom(A. field(), F);

		for (col_p = Ap -> colBegin(), e_p = e.begin();
		     e_p != e.end(); ++ col_p, ++ e_p) {
			
			r.assign(*e_p, one);
			
			A.apply (tmp, e);

			for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
			     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)

				hom.image (*elt_p, *tmp_p);

			r.assign(*e_p, zero);
		}
	}


	template <class Field, class IMatrix>
	void MatrixHom::map (SparseMatrix<Field>* &Ap, const IMatrix& A, const Field &F) {

		Ap = new SparseMatrix<Field>(F, A.rowdim(), A.coldim());

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
					Ap -> setEntry (j,i, val);		
			
			}
			r.assign(*e_p, zero);
		}
		
	}




	template <class Ring, class Field>
	void MatrixHom::map (DenseMatrix<Field>*& Ap, const SparseMatrix<Ring>& A, const Field &F) {
	
		Ap = new DenseMatrix<Field>(F, A.rowdim(), A.coldim());
		
		typename DenseMatrix<Field>::Element zero; F. init (zero, 0);
		typename DenseMatrix<Field>::RawIterator raw_p;
		for (raw_p = Ap -> rawBegin(); raw_p != Ap -> rawEnd(); ++ raw_p)
			F. assign (*raw_p, zero);
		
		typename SparseMatrix<Ring>::ConstRowIterator row_p;
	
		std::vector<size_t>::const_iterator j_p;
		
		typename std::vector<typename Ring::Element>::const_iterator e_p;
		
		typename Field::Element e;
		
		int i = 0;
		Hom<Ring, Field> hom(A. field(), F);
		
		for (row_p = A.rowBegin(); row_p != A.rowEnd(); ++ row_p, ++ i)
			for (j_p = row_p -> first. begin(), e_p = row_p -> second. begin(); 
			     j_p != row_p -> first. end(); ++ e_p, ++ j_p) {
				
				//F.init (e, *e_p);
				hom. image (e, *e_p);
				
				if (!F.isZero(e)) 
					Ap -> setEntry (i, *j_p, e);		
				
			}
		
	}


	
	template <class Field, class IMatrix>
	void MatrixHom::map (BlasBlackbox<Field> *&Ap, const IMatrix &A, 
			     const Field &F, MatrixContainerCategory::Container type) {

		Ap = new BlasBlackbox<Field>(F, A.rowdim(), A.coldim());
		Hom<typename IMatrix::Field , Field> hom(A.field(), F);

		typename IMatrix::ConstRawIterator          it = A.rawBegin();
		typename BlasBlackbox<Field>::RawIterator it_p = Ap->rawBegin();
		for (; it != A.rawEnd(); ++it, ++it_p)
			hom.image(*it_p, *it);
	}

	template <class Field, class IMatrix>
	void MatrixHom::map (BlasBlackbox<Field> *&Ap, const IMatrix &A, 
			     const Field &F, MatrixContainerCategory::Blackbox type) {

		Ap = new BlasBlackbox<Field>(F, A.rowdim(), A.coldim());

		typedef typename IMatrix::Field Ring;
		Ring r = A.field();

		typename Ring::Element one, zero;
		r. init(one, 1);
		r. init(zero, 0);

		std::vector<typename Ring::Element> e(A.coldim(), zero), tmp(A.rowdim());

		typename BlasBlackbox<Field>::ColIterator col_p;

		typename BlasBlackbox<Field>::Col::iterator elt_p;

		typename std::vector<typename Ring::Element>::iterator e_p, tmp_p;

		Hom<Ring, Field> hom(A. field(), F);

		for (col_p = Ap -> colBegin(), e_p = e.begin();
		     e_p != e.end(); ++ col_p, ++ e_p) {
			
			r.assign(*e_p, one);			
			A.apply (tmp, e);

			for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
			     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)
				hom.image (*elt_p, *tmp_p);
			r.assign(*e_p, zero);
		}
	}


}

#endif
