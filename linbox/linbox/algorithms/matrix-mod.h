/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
#ifndef __LINBOX_MATRIX_MOD_H__
#define __LINBOX_MATRIX_MOD_H__

#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/compose.h>
#include <linbox/integer.h>

namespace LinBox {
	
	/// @memo Limited doc so far. Used in RationalSolver.
	class MatrixMod {
		
	public:
	
		// general case, I donot how to do it.
		template<class FMatrix, class IMatrix, class Field>
		static void mod (FMatrix* & Ap, const IMatrix& A, Field F);
		
		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Field, class IMatrix>
		static void mod (DenseMatrix<Field>* &Ap, const IMatrix& A, Field F);
		
		// construct a sparse matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Field, class IMatrix>
		static void mod (SparseMatrix<Field>* &Ap, const IMatrix& A, Field F);

		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		static void mod (DenseMatrix<Field>* &Ap, const DenseMatrix<Ring>& A, Field F);

		// construct a dense matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		static void mod (DenseMatrix<Field>* &Ap, const SparseMatrix<Ring>& A, Field F);
		
		// construct a sparse matrix over finite field, such that *Ap = A mod p, where F = Ring / <p>
		template<class Ring, class Field>
		static void mod (SparseMatrix<Field>*& Ap, const SparseMatrix<Ring>& A, Field F);
	};		

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

		for (col_p = Ap -> colBegin(), e_p = e.begin();
		     e_p != e.end(); ++ col_p, ++ e_p) {
			
			r.assign(*e_p, one);
			
			A.apply (tmp, e);

			for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
			     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)

				F. init (*elt_p, *tmp_p);

			r.assign(*e_p, zero);
		}
	}


	template <class Field, class IMatrix>
	void MatrixMod::mod (SparseMatrix<Field>* &Ap, const IMatrix& A, Field F) {

		Ap = new SparseMatrix<Field>(F, A.rowdim(), A.coldim());

		typedef typename IMatrix::Field Ring;

		Ring r = A.field();
		integer buff;

		typename Ring::Element one, zero;

		r. init(one, 1);

		r. init(zero, 0);

		std::vector<typename Ring::Element> e(A.coldim(), zero), tmp(A.rowdim());
		
		typename std::vector<typename Ring::Element>::iterator iter, e_p;
		
		typename Field::Element val;
		
		int i = 0;
		
		for (e_p=e.begin();e_p != e.end(); ++e_p,i++){
			r.assign(*e_p, one);
			A.apply(tmp,e);
			int j;
			for (iter=tmp.begin(),j=0; iter != tmp.end(); ++iter,j++) {
				F.init (val, r.convert(buff, *iter));	       	
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
		
		integer tmp;
		for (A_p = A. rawBegin(), Ap_p = Ap -> rawBegin();
		     A_p != A. rawEnd(); ++ A_p, ++ Ap_p) 
			//F.init (*Ap_p, *A_p);
			{A. field(). convert (tmp, *A_p); F. init (*Ap_p, tmp);}
	}

	template <class Ring, class Field>
	void MatrixMod::mod (DenseMatrix<Field>*& Ap, const SparseMatrix<Ring>& A, Field F) {
	
		Ap = new DenseMatrix<Field>(F, A.rowdim(), A.coldim());
		
		typename SparseMatrix<Ring>::ConstRowIterator row_p;
	
		std::vector<size_t>::const_iterator j_p;
		
		typename std::vector<typename Ring::Element>::const_iterator e_p;
		
		typename Field::Element e;
		
		int i = 0;
		
		for (row_p = A.rowBegin(); row_p != A.rowEnd(); ++ row_p, ++ i)
			for (j_p = row_p -> first. begin(), e_p = row_p -> second. begin(); 
			     j_p != row_p -> first. end(); ++ e_p, ++ j_p) {
				
				F.init (e, *e_p);
				
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
		
		for (row_p = A.rowBegin(); row_p != A.rowEnd(); ++ row_p, ++ i)
			for (j_p = row_p -> first. begin(), e_p = row_p -> second. begin(); 
			     j_p != row_p -> first. end(); ++ e_p, ++ j_p) {
				
				F.init (e, *e_p);
				
				if (!F.isZero(e)) 
					Ap -> setEntry (i, *j_p, e);		
				
			}
		
	}

}

#endif
