/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflapack/fflapack.h
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __FFLAPACK_H
#define __FFLAPACK_H

#include "linbox/fflas/fflas.h"

class FFLAPACK : public FFLAS{
	
	
public:
	enum FFLAPACK_LUDIVINE_TAG { FflapackLUP=1, FflapackRank=2, 
				  FflapackLSP=3, FflapackSingular=4 };

	//---------------------------------------------------------------------
	// Rank: Rank for dense matrices based on LUP factorisation of A 
	//---------------------------------------------------------------------
	
	template <class Field>
	static size_t 
	Rank( const Field& F, const size_t M, const size_t N,
	      typename Field::element * A, const size_t lda){
		size_t *P = new size_t[N];
		size_t R = LUdivine( F, FflasUnit, M, N, A, lda, P, FflapackRank);
		delete[] P;
		return R;
 	}

	//---------------------------------------------------------------------
	// IsSingular: return true if A is singular
	// ( using rank computation with early termination ) 
	//---------------------------------------------------------------------
	
	template <class Field>
	static bool 
	IsSingular( const Field& F, const size_t M, const size_t N,
		    typename Field::element * A, const size_t lda){
		size_t *P = new size_t[N];
		return ( (bool) !LUdivine( F, FflasUnit, M, N, A, lda,
					  P, FflapackSingular));
 	}
	//---------------------------------------------------------------------
	// LUdivine: LUP factorisation of A 
	// P is the permutation matrix stored in an array of indexes
	//---------------------------------------------------------------------
	
	template <class Field>
	static size_t 
	LUdivine( const Field& F, const enum FFLAS_DIAG Diag,
		  const size_t M, const size_t N,
		  typename Field::element * A, const size_t lda,
		  size_t* P, const enum FFLAPACK_LUDIVINE_TAG LuTag=FflapackLUP );
	
	//---------------------------------------------------------------------
	// flaswp: apply the permutation P to the matrix A
	// P is stored in an array of indexes
	//---------------------------------------------------------------------
	template <class Field>
	static void 
	flaswp(const Field& F, const size_t N ,
	       typename Field::element * A, const size_t lda, 
	       const size_t K1, const size_t K2, size_t *piv, int inci);

	//---------------------------------------------------------------------
	// CharPoly: Compute the characteristic polynomial of A using Krylov
	// Method, and LUP factorization of the Krylov Base
	//---------------------------------------------------------------------
	template <class Field, class Polynomial>
	static list<Polynomial>&
	CharPoly( const Field& F, list<Polynomial>& charp, const size_t N,
		  const typename Field::element * A, const size_t lda,
		  typename Field::element * U, const size_t ldu);

	//---------------------------------------------------------------------
	// MinPoly: Compute the minimal polynomial of (A,v) using an LUP 
	// factorization of the Krylov Base (v, Av, .., A^kv)
	// U,X must be (n+1)*n
	// U contains the Krylov matrix and X, its LSP factorization
	//---------------------------------------------------------------------
	template <class Field, class Polynomial>
	static Polynomial&
	MinPoly( const Field& F, Polynomial& minP, const size_t N,
		 const typename Field::element *A, const size_t lda,
		 typename Field::element* U, size_t ldu,typename Field::element* X, size_t ldx,
		 size_t* P);


protected:
	
	// Compute the new d after a LSP ( d[i] can be zero )
	template<class Field>
	static size_t 
	newD( const Field& F, size_t * d, bool& KeepOn,
	      const size_t l, const size_t N, 
	      const typename Field::element * X,
	      typename Field::element ** minpt);

	template<class Field>
	static void 
	updateD(const Field& F, size_t * d, size_t& k,
		typename Field::element** minpt);
	
	//---------------------------------------------------------------------
	// TriangleCopy: copy a semi-upper-triangular matrix A to its triangular
	//               form in T, by removing the zero rows of A.
	//               cf Ibara, Moran, Hui 1982
	//               T is R*R
	//---------------------------------------------------------------------
	template <class Field>
	static void
	TriangleCopy( const Field& F, const enum FFLAS_UPLO Side,
		      const enum FFLAS_DIAG Diag, const size_t R, 
		      typename Field::element * T, const size_t ldt, 
		      const typename Field::element * A, const size_t lda ){

		const typename Field::element * Ai = A;
		typename Field::element * Ti = T;
		size_t j ;
		if ( Side == FflasUpper ){
			j=R-1;
			for (; Ti<T+R*ldt; Ti+=ldt, Ai+=lda+1, --j){
				while (F.iszero(*Ai))
					Ai+=lda;
				if ( Diag == FflasUnit ){
					*(Ti++) = F.one;
					fcopy( F, j, Ti, 1, Ai+1, 1);
				}
				else{
					fcopy( F, j+1, Ti++, 1, Ai, 1);
				}
			}
		}
		else{
			j=0;
			for (; Ti<T+R*ldt; Ti+=ldt+1, Ai+=lda+1, ++j){
				while (F.iszero(*Ai)){
					Ai+=lda;
				}
				if ( Diag == FflasUnit ){
					*Ti = F.one;
					fcopy( F, j, Ti-j, 1, Ai-j, 1);
				}
				else{
					fcopy( F, j+1, Ti-j, 1, Ai-j, 1);
				}
			}
		}
	}
	
	//---------------------------------------------------------------------
	// RectangleCopy: copy a rectangular matrix A to its reduced
	//               form in T, by removing the zero rows of A.
	//               T is M*N.
	//---------------------------------------------------------------------
	template <class Field>
	static void
	RectangleCopy( const Field& F, const size_t M, const size_t N, 
		       typename Field::element * T, const size_t ldt, 
		       const typename Field::element * A, const size_t lda ){

		const typename Field::element * Ai = A;
		typename Field::element * Ti = T;
		size_t x = M;
		for (; Ti<T+M*ldt; Ti+=ldt, Ai+=lda){
			while (F.iszero(*(Ai-x))){ // test if the pivot is 0
				Ai+=lda;
			}
			fcopy( F, N, Ti, 1, Ai, 1);
			x--;
		}
	}

	//---------------------------------------------------------------------
	// LUdivine_construct: (Specialisation of LUdivine)
	// LUP factorisation of X, the Krylov base matrix of A^t and v, in A.
	// X contains the nRowX first vectors v, vA, .., vA^{nRowX-1}
	// A contains the LUP factorisation of the nUsedRowX first row of X.
	// When all rows of X have been factorized in A, and rank is full,
	// then X is updated by the following scheme: X <= ( X; X.B ), where
	// B = A^2^i.
	// This enables to make use of Matrix multiplication, and stop computing
	// Krylov vector, when the rank is not longer full.
	// P is the permutation matrix stored in an array of indexes
	//---------------------------------------------------------------------
	
	template <class Field>
	static size_t
	LUdivine_construct( const Field& F, const enum FFLAS_DIAG Diag,
				   const size_t M, const size_t N,
				   typename Field::element * B, const size_t ldb,
				   typename Field::element * X, const size_t ldx,
				   typename Field::element * A, const size_t lda,
				   size_t* P, size_t* nRowX, const size_t nRowXMax,
				   size_t* nUsedRowX);
};


#include "linbox/fflapack/fflapack_flaswp.inl"
#include "linbox/fflapack/fflapack_LUdivine.inl"
//#ifndef CONSTRUCT
//#include "linbox/fflapack/fflapack_Minpoly.inl"
//#else
#include "linbox/fflapack/fflapack_Minpoly_construct.inl"
//#endif
#ifndef KGLU
#include "linbox/fflapack/fflapack_Charpoly.inl"
#else
#include "linbox/fflapack/fflapack_Charpoly_KGLU.inl"
#endif

#endif // __FFLAPACK_H
