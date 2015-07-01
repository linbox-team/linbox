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

class FFLAPACK : public FFLAS {
	
	
public:
	enum FFLAPACK_LUDIVINE_TAG { FflapackLQUP=1, FflapackRank=2, 
                                     FflapackLSP=3, FflapackSingular=4,
                                     FflapackTURBO=5};

	//---------------------------------------------------------------------
	// Rank: Rank for dense matrices based on LUP factorisation of A 
	//---------------------------------------------------------------------
	
	template <class Field>
	static size_t 
	Rank( const Field& F, const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda){
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
		    typename Field::Element * A, const size_t lda){
		size_t *P = new size_t[N];
		return ( (bool) !LUdivine( F, FflasUnit, M, N, A, lda,
					  P, FflapackSingular));
 	}
	
	//---------------------------------------------------------------------
	// TURBO: rank computation algorithm 
	//---------------------------------------------------------------------
	template <class Field>
	static size_t 
	FFLAPACK::TURBO( const Field& F, const size_t M, const size_t N,		
			 typename Field::Element * A, const size_t lda );

	//---------------------------------------------------------------------
	// LUdivine: LUP factorisation of A 
	// P is the permutation matrix stored in an array of indexes
	//---------------------------------------------------------------------
	
	template <class Field>
	static size_t 
	LUdivine( const Field& F, const enum FFLAS_DIAG Diag,
		  const size_t M, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  size_t* P, const enum FFLAPACK_LUDIVINE_TAG LuTag=FflapackLQUP,
		  size_t* rowP = NULL);
        
	
	//---------------------------------------------------------------------
	// flaswp: apply the permutation P to the matrix A
	// P is stored in an array of indexes
	//---------------------------------------------------------------------
	template <class Field>
	static void 
	flaswp(const Field& F, const size_t N ,
	       typename Field::Element * A, const size_t lda, 
	       const size_t K1, const size_t K2, size_t *piv, int inci);

	//---------------------------------------------------------------------
	// CharPoly: Compute the characteristic polynomial of A using Krylov
	// Method, and LUP factorization of the Krylov Base
	//---------------------------------------------------------------------
	template <class Field, class Polynomial>
	static list<Polynomial>&
	CharPoly( const Field& F, list<Polynomial>& charp, const size_t N,
		  const typename Field::Element * A, const size_t lda,
		  typename Field::Element * U, const size_t ldu);
	
	//---------------------------------------------------------------------
	// MinPoly: Compute the minimal polynomial of (A,v) using an LUP 
	// factorization of the Krylov Base (v, Av, .., A^kv)
	// U,X must be (n+1)*n
	// U contains the Krylov matrix and X, its LSP factorization
	//---------------------------------------------------------------------
	template <class Field, class Polynomial>
	static Polynomial&
	MinPoly( const Field& F, Polynomial& minP, const size_t N,
		 const typename Field::Element *A, const size_t lda,
		 typename Field::Element* U, size_t ldu,typename Field::Element* X, size_t ldx,
		 size_t* P);


protected:
	
	// Solve L X = B in place
	// L is M*M, B is M*N.
	// Only the R non trivial column of L are stored in the M*R matrix L
	template<class Field>
	static void
	solveLB( const Field& F, const size_t M, const size_t N, const size_t R, 
		 typename Field::Element * L, const size_t ldl, 
		 const size_t * Q,
		 typename Field::Element * B, const size_t ldb ){
		

		for (int i=R-1; i>=0; --i){
			if (Q[i] > i){
				//for (size_t j=0; j<=Q[i]; ++j)
				//F.init( *(L+Q[i]+j*ldl), 0 );
				fcopy( F, M-Q[i]-1, L+Q[i]*(ldl+1)+ldl,ldl, L+(Q[i]+1)*ldl+i, ldl );
				for ( size_t j=Q[i]; j<M; ++j)
					F.init( *(L+i+j*ldl), 0 );
			}
		} 
		ftrsm( F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit, M, N, F.one, L, ldl , B, ldb);
// 		static typename Field::Element Mone;
// 		F.init( Mone, -1 );
// 		typename Field::Element * Lcurr,* Rcurr,* Bcurr;
// 		size_t ib, k, Ldim,j=0;
// 		//cerr<<"In solveLB"<<endl;
// 		while ( j<R ) {
// 			k = ib = Q[j];
// 			//cerr<<"j avant="<<j<<endl;
// 			while ( (Q[j] == k) && (j<R) ) {k++;j++;}
		
// 			Ldim = k-ib;
// 			//cerr<<"k, ib, j, R "<<k<<" "<<ib<<" "<<j<<" "<<R<<endl;
// 			//cerr<<"M,k="<<M<<" "<<k<<endl;
// 			//cerr<<" ftrsm with M, N="<<Ldim<<" "<<N<<endl;
			
// 			Lcurr = L + j-Ldim + ib*ldl;
// 			Bcurr = B + ib*ldb;
// 			Rcurr = Lcurr + Ldim*ldl;

// 			ftrsm( F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit, Ldim, N, F.one, Lcurr, ldl , Bcurr, ldb );
			
// 			//cerr<<"M,k="<<M<<" "<<k<<endl;
// 			//cerr<<" fgemm with M, N, K="<<M-k<<" "<<N<<" "<<Ldim<<endl;

// 			fgemm( F, FflasNoTrans, FflasNoTrans, M-k, N, Ldim, Mone, Rcurr , ldl, Bcurr, ldb, F.one, Bcurr+Ldim*ldb, ldb);
// 		}
	}
	
	// Solve L X = B in place
	// L is M*M, B is M*N.
	// Only the R non trivial column of L are stored in the M*R matrix L
	template<class Field>
	static void
	solveLB2( const Field& F, const size_t M, const size_t N, const size_t R, 
		  typename Field::Element * L, const size_t ldl, 
		  const size_t * Q,
		  typename Field::Element * B, const size_t ldb ){
		

		
		static typename Field::Element Mone;
		F.init( Mone, -1 );
		typename Field::Element * Lcurr,* Rcurr,* Bcurr;
		size_t ib, k, Ldim,j=0;
		//cerr<<"In solveLB"<<endl;
		while ( j<R ) {
			k = ib = Q[j];
			//cerr<<"j avant="<<j<<endl;
			while ( (Q[j] == k) && (j<R) ) {k++;j++;}
		
			Ldim = k-ib;
			//cerr<<"k, ib, j, R "<<k<<" "<<ib<<" "<<j<<" "<<R<<endl;
			//cerr<<"M,k="<<M<<" "<<k<<endl;
			//cerr<<" ftrsm with M, N="<<Ldim<<" "<<N<<endl;
			
			Lcurr = L + j-Ldim + ib*ldl;
			Bcurr = B + ib*ldb;
			Rcurr = Lcurr + Ldim*ldl;

			ftrsm( F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit, Ldim, N, F.one, Lcurr, ldl , Bcurr, ldb );
			
			//cerr<<"M,k="<<M<<" "<<k<<endl;
			//cerr<<" fgemm with M, N, K="<<M-k<<" "<<N<<" "<<Ldim<<endl;

			fgemm( F, FflasNoTrans, FflasNoTrans, M-k, N, Ldim, Mone, Rcurr , ldl, Bcurr, ldb, F.one, Bcurr+Ldim*ldb, ldb);
		}
	}
	
	// Apply a row permutation matrix Q to an M*N matrix A
	template<class Field>
	static void 
	applyrowp( const Field& F, const size_t M, const size_t N, 
		   typename Field::Element * A, const size_t lda, size_t * Q ){
		
		typename Field::Element * temp = new typename Field::Element[N];
		//cerr<<"Q[i]=";
		for (size_t i=0; i<M; ++i){
			//cerr<<Q[i];
			if ( Q[i]>i )
				fswap( F, N, A + Q[i]*lda, 1, A +i*lda, 1 );
		}
		delete[] temp;
	}
	
	// Compute the new d after a LSP ( d[i] can be zero )
	template<class Field>
	static size_t 
	newD( const Field& F, size_t * d, bool& KeepOn,
	      const size_t l, const size_t N, 
	      const typename Field::Element * X,
	      typename Field::Element ** minpt);

	template<class Field>
	static void 
	updateD(const Field& F, size_t * d, size_t& k,
		typename Field::Element** minpt);
	
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
		      typename Field::Element * T, const size_t ldt, 
		      const typename Field::Element * A, const size_t lda ){

		const typename Field::Element * Ai = A;
		typename Field::Element * Ti = T;
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
	//---------------------------------------------------------------------
	// RectangleCopy: copy a rectangular matrix A to its reduced
	//               form in T, by removing the zero rows of A.
	//               T is M*N.
	//---------------------------------------------------------------------
	template <class Field>
	static void
	RectangleCopy( const Field& F, const size_t M, const size_t N, 
		       typename Field::Element * T, const size_t ldt, 
		       const typename Field::Element * A, const size_t lda ){

		const typename Field::Element * Ai = A;
		typename Field::Element * Ti = T;
		size_t x = M;
		for (; Ti<T+M*ldt; Ti+=ldt, Ai+=lda){
			while (F.iszero(*(Ai-x))){ // test if the pivot is 0
				Ai+=lda;
			}
			fcopy( F, N, Ti, 1, Ai, 1);
			x--;
		}
	}
	/*	template <class Field>
	static void
	RectangleCopy( const Field& F, const size_t M, const size_t N, 
		       const long dist2pivot,
		       typename Field::Element * T, const size_t ldt, 
		       const typename Field::Element * A, const size_t lda ){

		const typename Field::Element * Ai = A;
		typename Field::Element * Ti = T;
	        long x = dist2pivot;
	        for (; Ti<T+M*ldt; Ti+=ldt, Ai+=lda){
			while (F.iszero(*(Ai-x))){ // test if the pivot is 0
	        		Ai+=lda;
	        	}
			fcopy( F, N, Ti, 1, Ai, 1);
			x--;
		}
	}
	*/
	//---------------------------------------------------------------------
	// RectangleCopy2: copy a rectangular matrix A to its reduced
	//               form in T, by removing the rows of A corresponding to 
	//               triangular part of the lsp located in A-dist2pivot.
	//               T is M*N.
	//---------------------------------------------------------------------
	template <class Field>
	static void
	RectangleCopy2( const Field& F, const size_t M, const size_t N, 
		       const long dist2pivot,
		       typename Field::Element * T, const size_t ldt, 
		       const typename Field::Element * A, const size_t lda ){

		const typename Field::Element * Ai = A;
		typename Field::Element * Ti = T;
		long x = dist2pivot;
		for (; Ti<T+M*ldt; Ai+=lda){
			while (F.iszero(*(Ai-x))){ // test if the pivot is 0
				fcopy( F, N, Ti, 1, Ai, 1);
				Ai += lda;
				Ti += ldt;
			}
			x--;
		}
	}
	//---------------------------------------------------------------------
	// RectangleCopyTURBO: Copy A to T, with respect to the row permutation 
	//                     defined by the lsp factorization of located in 
	//                     A-dist2pivot
	//---------------------------------------------------------------------
	template <class Field>
	static void
	RectangleCopyTURBO( const Field& F, const size_t M, const size_t N, 
		       const size_t dist2pivot, const size_t rank,
		       typename Field::Element * T, const size_t ldt, 
		       const typename Field::Element * A, const size_t lda ){

		const typename Field::Element * Ai = A;
		typename Field::Element * T1i = T, T2i = T + rank*ldt;
		size_t x = dist2pivot;
		for (; Ai<A+M*lda; Ai+=lda){
			while ( F.iszero(*(Ai-x)) ) { // test if the pivot is 0
				fcopy( F, N, T2i, 1, Ai, 1);
				Ai += lda;
				T2i += ldt;
			}
			fcopy( F, N, T1i, 1, Ai, 1);
			T1i += ldt;
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
				   typename Field::Element * B, const size_t ldb,
				   typename Field::Element * X, const size_t ldx,
				   typename Field::Element * A, const size_t lda,
				   size_t* P, size_t* nRowX, const size_t nRowXMax,
				   size_t* nUsedRowX);
};


#include "linbox/fflapack/fflapack_flaswp.inl"
#include "fflapack_ludivine.inl"
//#ifndef CONSTRUCT
//#include "linbox/fflapack/fflapack_Minpoly.inl"
//#else
#include "linbox/fflapack/fflapack_minpoly_construct.inl"
//#endif
#ifndef KGLU
#include "linbox/fflapack/fflapack_charpoly.inl"
#else
#include "linbox/fflapack/fflapack_charpoly_kglu.inl"
#endif

#endif // __FFLAPACK_H
