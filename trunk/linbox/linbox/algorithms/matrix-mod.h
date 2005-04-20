/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
#ifndef __LINBOX_MATRIX_MOD_H__
#define __LINBOX_MATRIX_MOD_H__

#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/compose.h>
#include <linbox/integer.h>
#include <linbox/field/hom.h>
#include <linbox/matrix/dense.h>
namespace LinBox {

	// try to map a blackbox over a homorphic ring
	// The most suitable type
	template <class Blackbox, class Field>
	struct MatrixModTrait{
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
	struct MatrixModTrait<DenseMatrixBase<typename Ring::Element>, Field> {
		typedef DenseMatrixBase<typename Field::Element> value_type;
	};
	template <class Ring, class Field>
	struct MatrixModTrait<SparseMatrix<Ring, typename Vector<Ring>::SparseSeq>, Field> {
		typedef SparseMatrix<Field, typename Vector<Field>::SparseSeq> value_type;
	};

	template <class Ring, class Field>
	struct MatrixModTrait<SparseMatrix<Ring, typename Vector<Ring>::SparsePar>, Field> {
		typedef SparseMatrix<Field, typename Vector<Field>::SparsePar> value_type;
	};

	template <class Ring, class Field>
	struct MatrixModTrait<SparseMatrix<Ring, typename Vector<Ring>::SparseMap>, Field> {
		typedef SparseMatrix<Field, typename Vector<Field>::SparseMap> value_type;
	};

	template <class Ring, class Field>
	struct MatrixModTrait<DenseMatrix<Ring>, Field> {
		typedef DenseMatrix<Field> value_type;
	};

	template <class Ring, class Field>
	struct MatrixModTrait<BlasBlackbox<Ring>, Field> {
		typedef BlasBlackbox<Field> value_type;
	};

	/// @memo Limited doc so far. Used in RationalSolver.
	namespace MatrixMod {
		
	//public:
	
		// general case, I do not how to do it.
		template<class FMatrix, class IMatrix, class Field>
		void mod (FMatrix* & Ap, const IMatrix& A, Field F);
		
		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Field, class IMatrix>
		void mod (DenseMatrix<Field>* &Ap, const IMatrix& A, Field F);
		
		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		void mod (DenseMatrixBase<typename Field::Element>* &Ap, const DenseMatrixBase<typename Ring::Element>& A, const Field& F){
			Ap = new DenseMatrixBase<typename Field::Element>(A.rowdim(), A.coldim());
			
			typedef DenseMatrixBase<typename Ring::Element> IMatrix;
			typename IMatrix::ConstRawIterator         iter_value = A.rawBegin();
			typename IMatrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
			typename Field::Element tmp;
			for (;iter_value != A.rawEnd(); ++iter_value,++iter_index){
				F.init(  tmp, *iter_value ); 
				Ap->setEntry(iter_index.rowIndex(), iter_index.colIndex(),tmp);
			}
		}
		
		// construct a sparse matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Field, class IMatrix>
		void mod (SparseMatrix<Field>* &Ap, const IMatrix& A, Field F);

		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		void mod (DenseMatrix<Field>* &Ap, const DenseMatrix<Ring>& A, Field F);

		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		void mod (DenseMatrix<Field>* &Ap, const SparseMatrix<Ring>& A, Field F);
		
		// construct a sparse matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		void mod (SparseMatrix<Field>*& Ap, const SparseMatrix<Ring>& A, Field F);
	}		

	template <class Field, class IMatrix>
	void MatrixMod::mod (DenseMatrix<Field>* &Ap, const IMatrix& A, Field F) {

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
	void MatrixMod::mod (SparseMatrix<Field>* &Ap, const IMatrix& A, Field F) {

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
	void MatrixMod::mod (DenseMatrix<Field>*& Ap, const DenseMatrix<Ring>& A, Field F) {
		
		Ap = new DenseMatrix<Field>(F, A.rowdim(), A.coldim());
		typename DenseMatrix<Ring>::ConstRawIterator A_p;
		typename DenseMatrix<Field>::RawIterator Ap_p;
		Hom<Ring, Field> hom(A. field(), F);
		for (A_p = A. rawBegin(), Ap_p = Ap -> rawBegin();
		     A_p != A. rawEnd(); ++ A_p, ++ Ap_p) 
			hom.image (*Ap_p, *A_p);
	}
	

	template <class Ring, class Field>
	void MatrixMod::mod (BlasBlackbox<Field>*& Ap, const BlasBlackbox<Ring>& A, Field F) {
		
		Ap = new BlasBlackbox<Field>(F, A.rowdim(), A.coldim());
		typename BlasBlackbox<Ring>::ConstRawIterator A_p;
		typename BlasBlackbox<Field>::RawIterator Ap_p;
		Hom<Ring, Field> hom(A. field(), F);
		for (A_p = A. rawBegin(), Ap_p = Ap -> rawBegin();
		     A_p != A. rawEnd(); ++ A_p, ++ Ap_p) 
			hom.image (*Ap_p, *A_p);
	}

	template <class Ring, class Field>
	void MatrixMod::mod (DenseMatrix<Field>*& Ap, const SparseMatrix<Ring>& A, Field F) {
	
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

	template <class Ring, class Field>
	void MatrixMod::mod (SparseMatrix<Field>*& Ap, const SparseMatrix<Ring>& A, Field F) {
	
		Ap = new SparseMatrix<Field>(F, A.rowdim(), A.coldim());
		
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

}

#endif
