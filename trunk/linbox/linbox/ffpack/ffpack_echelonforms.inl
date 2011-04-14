/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* ffpack_echelon.h
 * Copyright (C) 2009, 2010 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ffpack_echelon_forms_INL
#define __LINBOX_ffpack_echelon_forms_INL

template <class Field>
size_t FFPACK::ColumnEchelonForm (const Field& F, const size_t M, const size_t N,
				  typename Field::Element * A, const size_t lda,
				  size_t* P, size_t* Qt, const bool transform)
{

	typename Field::Element one, mone;
	F.init (one, 1.0);
	F.neg (mone, one);
	size_t r;
	r = LUdivine (F, FflasNonUnit, FflasNoTrans, M, N, A, lda, P, Qt);

	if (transform){
		ftrtri (F, FflasUpper, FflasNonUnit, r, A, lda);
		ftrmm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit, r, N-r,
		       mone, A, lda, A+r, lda);
	}

	return r;
}

template <class Field>
size_t FFPACK::RowEchelonForm (const Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       size_t* P, size_t* Qt, const bool transform)
{

	typename Field::Element one, mone;
	F.init (one, 1.0);
	F.neg (mone, one);
	size_t r;
	r = LUdivine (F, FflasNonUnit, FflasTrans,  M, N, A, lda, P, Qt);

	if (transform){

		ftrtri (F, FflasLower, FflasNonUnit, r, A, lda);
		ftrmm (F, FflasRight, FflasLower, FflasNoTrans, FflasNonUnit, M-r, r,
		       mone, A, lda, A+r*lda, lda);
	}

	return r;
}

template <class Field>
size_t
FFPACK::ReducedColumnEchelonForm (const Field& F, const size_t M, const size_t N,
				  typename Field::Element * A, const size_t lda,
				  size_t* P, size_t* Qt, const bool transform)
{

	typename Field::Element one, mone, zero;
	F.init (one, 1.0);
	F.init (zero, 0.0);
	F.neg (mone, one);
	size_t r;
	r = ColumnEchelonForm (F, M, N, A, lda, P, Qt, transform);
	// M = Q^T M
	for (size_t i=0; i<r; ++i){
		if ( Qt[i]> (size_t) i ){
			fswap( F, i,
			       A + Qt[i]*lda, 1,
			       A + i*lda, 1 );
		}
	}
	if (transform){
		ftrtri (F, FflasLower, FflasUnit, r, A, lda);
		ftrmm (F, FflasRight, FflasLower, FflasNoTrans, FflasUnit, M-r, r,
		       one, A, lda, A+r*lda, lda);
		ftrtrm (F, FflasNonUnit, r, A, lda);
	} else {
		ftrsm (F, FflasRight, FflasLower, FflasNoTrans, FflasUnit, M-r, r,
		       one, A, lda, A+r*lda, lda);
		for (size_t i=0; i<r; i++){
			for (size_t j=0; j<N; j++)
				F.assign (*(A+i*lda+j),zero);
			F.assign (*(A + i*(lda+1)), one);
		}
		applyP(F, FflasLeft, FflasTrans, r, 0, r, A, lda, Qt);
	}

	return r;
}

template <class Field>
size_t
FFPACK::ReducedRowEchelonForm (const Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       size_t* P, size_t* Qt, const bool transform)
{

	typename Field::Element one, mone, zero;
	F.init (one, 1.0);
	F.init (zero, 0.0);
	F.neg (mone, one);
	size_t r;
	r = RowEchelonForm (F, M, N, A, lda, P, Qt, transform);
	// M = M Q
	for (size_t i=0; i<r; ++i){
		if ( Qt[i]>  i ){
			fswap( F, i,
			       A + Qt[i], lda,
			       A + i, lda );
		}
	}
	if (transform){
		ftrtri (F, FflasUpper, FflasUnit, r, A, lda);
		ftrmm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasUnit, r, N-r,
		       one, A, lda, A+r, lda);
		ftrtrm (F, FflasUnit, r, A, lda);
	} else {
		ftrsm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasUnit, r, N-r,
		       one, A, lda, A+r, lda);
		for (size_t i=0; i<r; i++){
			for (size_t j=0; j<M; j++)
				F.assign (*(A+j*lda+i),zero);
			F.assign (*(A + i*(lda+1)), one);
		}
		applyP(F, FflasRight, FflasNoTrans, r, 0, r, A, lda, Qt);
	}
	return r;
}

/*
 * Warning, this implementation is currently broken:
 * the LAPACK permutation mechanism can not be used here as is
 * More work required on the construction of the permutation P...
 */
template <class Field>
size_t
FFPACK::REF (const Field& F, const size_t M, const size_t N,
	     typename Field::Element * A, const size_t lda,
	     const size_t colbeg, const size_t rowbeg, const size_t colsize,
	     size_t* Qt, size_t* P)
{

	typedef typename Field::Element Element;
	Element one, mone, zero;
	F.init(zero, 0.0);
	F.init(one, 1.0);
	F.neg(mone, one);

	if (colsize == 1){
		for (size_t i=rowbeg; i<M; ++i){
			if (!F.isZero(*(A+i*lda+colbeg))){
				Qt[rowbeg] = i;
				if (i!= rowbeg){
					F.assign(*(A+rowbeg*lda+colbeg),*(A+i*lda+colbeg));
					F.assign(*(A+i*lda+colbeg), zero);
				}
				Element invpiv;
				F.inv(invpiv, *(A+rowbeg*lda + colbeg));
				F.assign(*(A+rowbeg*lda+colbeg), invpiv);
				F.negin(invpiv);
				for (size_t j=0; j<rowbeg; ++j)
					F.mulin (*(A+j*lda+colbeg), invpiv);
				for (size_t j=rowbeg+1; j<M; ++j)
					F.mulin (*(A+j*lda+colbeg), invpiv);
				return 1;
			}
		}
		Qt[rowbeg]=colbeg;
		return 0;
	}
	size_t recsize = colsize / 2;

	// Recurive call on slice A*1
	size_t r1 = REF(F, M, N, A, lda, colbeg, rowbeg, recsize, Qt, P);

	Element* A11 = A+colbeg;
	Element* A12 = A11+recsize;
	Element* A22 = A12+rowbeg*lda;
	Element* A21 = A11+rowbeg*lda;
	Element* A31 = A21+r1*lda;
	Element* A32 = A22+r1*lda;

	/**
	 *  ---------------------
	 * | I  | A11 | A12 |    |
	 * |----|-----|-----|----|
	 * |    |I | *| A22 |    |
	 * |    |0 | 0| A22 |    |
	 * |----|-----|-----|----|
	 * |    | 0   | A32 |    |
	 * |----|-----|-----|----|
	 *
	 * where the transformation matrix is stored at the pivot column position
	 */
	// Apply row permutation on A*2
	applyP (F, FflasLeft, FflasNoTrans, colsize - recsize, rowbeg, rowbeg+r1, A12, lda, Qt);

	// A12 <- A12 - A11 * A22
	fgemm (F, FflasNoTrans, FflasNoTrans, rowbeg, colsize - recsize, r1,
	       one, A11, lda, A22, lda, one, A12, lda);

	// A32 <- A32 - A31 * A22
	fgemm (F, FflasNoTrans, FflasNoTrans, M-rowbeg-r1, colsize - recsize, r1,
	       one, A31, lda, A22, lda, one, A32, lda);

	// A22 <- A21*A22
	Element* tmp = new Element [r1*(colsize-recsize)];
	for (size_t i = 0; i < r1; ++i)
		fcopy (F, colsize-recsize, tmp+i*(colsize-recsize), 1, A22+i*lda, 1);
	fgemm (F, FflasNoTrans, FflasNoTrans, r1, colsize-recsize, r1,
	       one, A21, lda, tmp, colsize-recsize, zero, A22, lda);
	delete[] tmp;

	// Recurive call on slice A*2
	size_t r2 = REF(F, M, N, A, lda, colbeg + recsize, rowbeg + r1,
			colsize - recsize, Qt, P);

	// Apply permutation on A*1
	applyP (F, FflasLeft, FflasNoTrans, r1, rowbeg+r1, rowbeg+r1+r2, A11, lda, Qt);

	Element * U11 = A11;
	Element * U12 = A12;
	Element * U21 = A31;
	Element * U22 = A32;
	Element * U31 = U21+r2*lda;
	Element * U32 = U31+recsize;

	// U11 <- U11 + U12 * U21
	fgemm (F, FflasNoTrans, FflasNoTrans, rowbeg+r1, r1, r2,
	       one, U12, lda, U21, lda, one, U11, lda);

	// U31 <- U31 + U32 * U21
	fgemm (F, FflasNoTrans, FflasNoTrans, M-rowbeg-r1-r2, r1, r2,
	       one, U32, lda, U21, lda, one, U31, lda);

	// U21 <- U22*U21
	tmp = new Element [r2*r1];
	for (size_t i = 0; i < r2; ++i)
		fcopy (F, r1, tmp+i*r1, 1, U21+i*lda, 1);

	fgemm (F, FflasNoTrans, FflasNoTrans, r2, r1, r2,
	       one, U22, lda, tmp, r1, zero, U21, lda);
	delete[] tmp;

	//Permute the non pivot columns to the end
	if (r1 < recsize){
		size_t ncol = recsize -r1;
		size_t nrow = rowbeg + r1;
		Element * NZ1 = A11+r1;

		tmp = new Element [nrow*ncol];
		for (size_t i=0; i < nrow; ++i)
			fcopy (F, ncol, tmp+i*ncol, 1, NZ1 + i*lda, 1);
		for (size_t i=0; i < M; ++i)
			// Risky copy with overlap, but safe with the naive
			// implementation of fcopy
			fcopy (F, r2, NZ1+i*lda, 1, A12 + i*lda, 1);
		NZ1 +=  r2;
		for (size_t i=0; i<nrow; ++i)
			fcopy (F, ncol, NZ1 + i*lda, 1, tmp+i*ncol,1);
		delete[] tmp;

		for (size_t i=rowbeg+r1; i<M; ++i)
			for (size_t j=0; j<recsize-r1; ++j)
				F.assign(*(NZ1+i*lda+j), zero);
		// size_t * temp = new size_t[recsize-r1];
		// for (size_t i=0,j = colbeg+r1; j<colbeg+recsize; ++i,++j)
		//  	temp[i] = P[j];
		// for (size_t  i = colbeg+recsize, j = colbeg+r1; i<colbeg+recsize+r2; ++i,++j)
		// 	P[j] = P[i];
		// for (size_t i=0,j = colbeg+r1+r2; i<recsize-r1; ++i,++j)
		// 	P[j] = temp[i];
		// delete temp;
		for (size_t  i = colbeg+recsize, j = colbeg+r1; i<colbeg+recsize+r2; ++i,++j){
			size_t t = P[i];
			P[i] = P[j];
			P[j] = t;
			//P[j]=P[i];
		}

	}

	return r1+r2;
}

#endif  // __LINBOX_ffpack_echelon_forms_INL
