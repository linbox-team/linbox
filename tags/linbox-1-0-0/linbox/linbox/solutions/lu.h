/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/LU.h
 *
 * writtend by zhendong wan <wan@mail.eecis.udel.edu>
 */

#ifndef LU_H
#define LU_H

#include <iostream>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/submatrix.h>

namespace LinBox
{
	/**
	 * M <-- the LU decomposition of M.
	 *
	 *@param M is a dense matrix on input.  
	 * require all principle minor are non-singular.
	 * Every leading principal minor must be nonzero.
	 * M is modified to represent
	 * the LU decomposition of the original M.
	 * L is unit lower triangular and occupies the strictly lower triangular
	 * part of M.  The diagonal one's are implicit.
	 * U is upper triangular and occupies the rest of M.
	 */
	template<class Field>
	void LU(DenseMatrix<Field>& M);
  
	// M <-- LU decomp of M (just as for DenseMatrix argument).
	template<class Field>
	void LU(Submatrix<DenseMatrix<Field> >& M);
  
	// M <-- M L^{-1}, where L is unit lower triangular with implicit diagonal.
	template<class Matrix>
	void LL_MULIN(Matrix& M, const Matrix& L);
  
	// R <-- RU, where U is upper triangular.
	template<class Matrix>
	void RU_MULIN(Matrix& R, const Matrix& U);
  
	template<class Matrix>
	void  AXMYIN(Matrix&, const Matrix&, const Matrix&);
  
	template<class Field>
	void LU(DenseMatrix<Field>& M) 
	{
		linbox_check (M.rowdim() == M.coldim());

		if (M. rowdim() <= 1)

			return;

		Submatrix<DenseMatrix<Field> > sub (M, 0, 0, M. rowdim(), M. coldim());      
	
		LU (sub);

	}
  
	template<class Field>
	void LU(Submatrix<DenseMatrix<Field> >& M) {

		int dim = M. rowdim ();

		if ( dim <= 1)

			return;

		int HALF = dim / 2;

		int rest = dim - HALF;

		Submatrix<DenseMatrix<Field> > M00 (M, 0, 0, HALF, HALF);

		Submatrix<DenseMatrix<Field> > M01 (M, 0, HALF, HALF, rest);

		Submatrix<DenseMatrix<Field> > M10 (M, HALF, 0, rest, HALF);

		Submatrix<DenseMatrix<Field> > M11 (M, HALF, HALF, rest, rest);

		LU (M00);

		LL_MULIN (M01,M00);

		RU_MULIN (M10,M00);

		AXMYIN (M11,M10,M01);

		LU (M11);     

	}

	/*
	 *@param L is a unit lower triangle matrix.
	 * All diagonal entries in L are one.
	 * $ M <- L^{-1} M$
	 */
	template<class Matrix>
	void LL_MULIN(Matrix& M, const Matrix& L) {

		typename Matrix::ColIterator colp;

		typename Matrix::ConstRowIterator crowp;

		typename Matrix::Col::iterator ep, epo;

		typename Matrix::ConstRow::const_iterator cep;

		typename Matrix::Element e;

		for(colp = M. colBegin();colp != M. colEnd(); ++colp) {

			crowp = L. rowBegin();

			for(ep = colp -> begin(); ep != colp -> end(); ++ ep, ++ crowp)

				for(cep = crowp -> begin(), epo = colp -> begin(); epo !=ep; ++ cep, ++ epo) {

					M.field().mul(e,*epo,*cep);

					M.field().subin(*ep,e);
				}
		}
	}
  
	/*
	 *param U is an upper triangle matrix.
	 *$M <- M U^{-1}$
	 */
	template<class Matrix>
	void RU_MULIN (Matrix& M, const Matrix& U) {

		typename Matrix::RowIterator rowp;

		typename Matrix::ConstColIterator ccolp;

		typename Matrix::Element e;

		typename Matrix::Row::iterator ep,epo;

		typename Matrix::ConstCol::const_iterator cep;
      
		for (rowp = M. rowBegin(); rowp != M. rowEnd(); ++ rowp) {

			ccolp = U. colBegin();

			for (ep = rowp -> begin (); ep != rowp -> end(); ++ ep, ++ ccolp) {

				for ( cep= ccolp -> begin (), epo = rowp -> begin (); epo != ep; ++ cep,++ epo) {

					M. field(). mul (e, *epo, *cep);

					M. field(). subin (*ep, e);
				}

				M. field(). divin (*ep, *cep);

			}
		}
	}

	/*
	 *@M1, M2, dense matrix.
	 *M<-M-M1*M2
	 */
	template<class Matrix>
	void  AXMYIN(Matrix& M, const Matrix& M1, const Matrix& M2) {

		typename Matrix::ColIterator colp;

		typename Matrix::ConstRowIterator crowp1;

		typename Matrix::ConstColIterator ccolp2; 

		typename Matrix::Element e;

		typename Matrix::Col::iterator ep;

		typename Matrix::ConstRow::const_iterator cep1;

		typename Matrix::ConstCol::const_iterator cep2;

		for (colp = M. colBegin(), ccolp2 = M2. colBegin(); colp != M. colEnd(); ++colp, ++ccolp2)

			for (ep = colp -> begin(), crowp1 = M1. rowBegin(); ep != colp->end(); ++ep, ++crowp1) {
				M. field(). init(e,0);

				for(cep1 = crowp1 -> begin(), cep2 = ccolp2 -> begin(); cep1 != crowp1 -> end(); ++cep1, ++cep2)
					M. field(). axpyin (e, *cep1, *cep2);
	
				M. field(). subin (*ep, e);
			}
	}
}

#endif
