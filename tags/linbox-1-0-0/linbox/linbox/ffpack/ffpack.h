/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/ffpack/ffpack.h
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __FFPACK_H
#define __FFPACK_H
//#define DEBUG 2
template<class Field>
std::ostream& write_field (const Field& F, std::ostream& os, const typename Field::Element* A, size_t M, size_t N, size_t lda){
	os<<"----------------------------\n";
	for (size_t i=0;i<M;i++){
		for (size_t j=0;j<N;j++)
			F.write(os,*(A+j+i*lda))<<",";
		os<<std::endl;
	}
	return os<<"----------------------------\n";
}


#include "linbox/fflas/fflas.h"
//#include <linbox/../examples/Matio.h>
#include <list>
#include <vector>

namespace LinBox{

	/**
	 * \brief Set of elimination based routines for dense linear algebra with matrices over finite prime field of characteristic less than 2^26.
	 *
	 *  This class only provides a set of static member functions. No instantiation is allowed.
	 *
	 * It enlarges the set of BLAS routines of the class FFLAS, with higher 
	 * level routines based on elimination.
	 \ingroup ffpack
	 */
class FFPACK : public FFLAS {
	
	
public:
	enum FFPACK_LUDIVINE_TAG { FfpackLQUP=1,
				   FfpackSingular=2};
	
	enum FFPACK_CHARPOLY_TAG { FfpackLUK=1,
				   FfpackKG=2,
				   FfpackHybrid=3,
				   FfpackKGFast=4,
				   FfpackHybrid2=5,
				   FfpackDanilevski=6 };
	
	enum FFPACK_MINPOLY_TAG { FfpackDense=1,
				    FfpackKGF=2 };

	/** 
	 * Computes the rank of the given matrix using a LQUP factorization.
	 * The input matrix is modified.
	 * @param M row dimension of the matrix
	 * @param N column dimension of the matrix
	 * @param A input matrix
	 * @param lda leading dimension of A
	 */
	
	template <class Field>
	static size_t 
	Rank( const Field& F, const size_t M, const size_t N,
	      typename Field::Element * A, const size_t lda){
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		size_t R = LUdivine( F, FflasUnit, M, N, A, lda, P, Q, FfpackLQUP);
		delete[] Q;
		delete[] P;
		return R;
 	}

	//---------------------------------------------------------------------
	// IsSingular: return true if A is singular ( A is modified)
	// ( using rank computation with early termination ) 
	//---------------------------------------------------------------------
	
	/// \brief  using rank computation with early termination. 
	template <class Field>
	static bool 
	IsSingular( const Field& F, const size_t M, const size_t N,
		    typename Field::Element * A, const size_t lda){
		
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		bool singular  = !LUdivine( F, FflasNonUnit, M, N, A, lda, P, Q, FfpackSingular);
		
		delete[] P;
		delete[] Q;
		return singular;
 	}
	
	//---------------------------------------------------------------------
	// Det: return det(A)
	// ( using LQUP factorization  with early termination ) 
	//---------------------------------------------------------------------
	
	///  using LQUP factorization  with early termination. 
	template <class Field>
	static typename Field::Element
	Det( const Field& F, const size_t M, const size_t N,
	     typename Field::Element * A, const size_t lda){
		
		typename Field::Element det;
		bool singular;
		size_t *P = new size_t[N];
		size_t *Q = new size_t[M];
		singular  = !LUdivine( F, FflasNonUnit, M, N, A, lda, P, Q, FfpackSingular);
		if (singular){
			F.init(det,0);
			delete[] P;
			delete[] Q; 
			return det;
		}
		else{
			F.init(det,1);
			typename Field::Element *Ai=A;
			for (; Ai < A+ M*lda+N; Ai+=lda+1 )
				F.mulin( det, *Ai );

			int count=0;
			for (size_t i=0;i<N;++i)
				if (P[i] != i) ++count;
				
			if ((count&1) == 1)
				F.negin(det);
		}
		delete[] P; 
		delete[] Q; 
		return det;
 	}
	//---------------------------------------------------------------------
	// Solve ( using LQUP factorization  ) 
	//---------------------------------------------------------------------
	
	/// Solve linear system using LQUP factorization. 
	template <class Field>
	static typename Field::Element*
	Solve( const Field& F, const size_t M,
		typename Field::Element * A, const size_t lda,
		typename Field::Element * x, const int incx,
		const typename Field::Element * b, const int incb ){
		
		typename Field::Element one, zero;
		F.init(one,1);
		F.init(zero,0);

		size_t *P = new size_t[M];
		size_t *rowP = new size_t[M];
		
		if (LUdivine( F, FflasNonUnit, M, M, A, lda, P, rowP, FfpackLQUP) < M){
			std::cerr<<"SINGULAR MATRIX"<<std::endl;
			delete[] P; 
			delete[] rowP; 
			return x;
		}
		else{
			fcopy( F, M, x, incx, b, incb );
			
			ftrsv(F,  FflasLower, FflasNoTrans, FflasUnit, M, 
			      A, lda , x, incx);
			ftrsv(F,  FflasUpper, FflasNoTrans, FflasNonUnit, M, 
			      A, lda , x, incx);
			applyP( F, FflasRight, FflasTrans, M, 0, M, x, incx, P );
			delete[] rowP; 
			delete[] P; 

			return x;
		
		}
 	}
	
	//---------------------------------------------------------------------
	// Invert ( using LQUP factorization  ) 
	//---------------------------------------------------------------------
	
	template <class Field>
	static typename Field::Element*
	Invert( const Field& F, const size_t M,
		typename Field::Element * A, const size_t lda,
		typename Field::Element * X, const size_t ldx,
		int& nullity) 
	{
		typename Field::Element one, zero;
		F.init (one,1);
		F.init (zero,0);

		size_t *P = new size_t[M];
		size_t *rowP = new size_t[M];
		
		nullity = M - LUdivine( F, FflasNonUnit, M, M, A, lda, P, rowP, FfpackLQUP);
		if (nullity > 0){
			delete[] P; 
			delete[] rowP; 
			return NULL;
		}
		// Improvement: construct X=P^1 directly
		for (size_t i=0;i<M;++i)
			for (size_t j=0; j<M;++j)
				if (i==j)
					F.assign( *(X+i*ldx+j), one);
				else
					F.assign( *(X+i*ldx+j), zero);
		
		applyP( F, FflasRight, FflasTrans, M, 0, M, X, ldx, P );
		ftrsm(F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, M, M, one, 
		      A, lda , X, ldx);
		ftrsm(F, FflasRight, FflasLower, FflasNoTrans, FflasUnit, M, M, one, 
		      A, lda , X, ldx);
		delete[] P;
		delete[] rowP;
		return X;
	}
	
	template <class Field>
	static typename Field::Element*
	Invert2( const Field& F, const size_t M,
		 typename Field::Element * A, const size_t lda,
		 typename Field::Element * X, const size_t ldx,
		 int& nullity){
		
		typename Field::Element one, zero;
		F.init(one,1);
		F.init(zero,0);

		size_t *P = new size_t[M];
		size_t *rowP = new size_t[M];
		
		nullity = M - LUdivine( F, FflasNonUnit, M, M, A, lda, P, rowP, FfpackLQUP);
		if (nullity > 0){
			delete[] P;
			delete[] rowP; 
			return NULL;
		} else {
			// Initializing X to 0
			for (size_t i=0; i<M; ++i)
				for (size_t j=0; j<M;++j)
					F.assign(*(X+i*ldx+j), zero);
			// X = L^-1 in n^3/3
			invL( F, M, A, lda, X, ldx );

			// X = Q^-1.X is not necessary since Q = Id
			
			// X = U^-1.X
			ftrsm( F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit, 
			       M, M, one, A, lda , X, ldx);

			// X = P^-1.X
			applyP( F, FflasLeft, FflasTrans, M, 0, M, X, ldx, P ); 
			
			delete[] P;
			delete[] rowP;
			return X;
		}	
 	}

	/** 
	 *  LQUPtoInverseOfFullRankMinor
	 * Suppose A has been factorized as L.Q.U.P, with rank r.
	 * Then Qt.A.Pt has an invertible leading principal r x r submatrix
	 * This procedure efficiently computes the inverse of this minor and puts it into X.
	 * NOTE: It changes the lower entries of A_factors in the process (NB: unless A was nonsingular and square)
	 *
	 * @param rank:       rank of the matrix.
	 * @param A_factors:  matrix containing the L and U entries of the factorization
	 * @param QtPointer:  theLQUP->getQ()->getPointer() (note: getQ returns Qt!)
	 * @param X:          desired location for output
	 */
	template <class Field>
	static typename Field::Element*
	LQUPtoInverseOfFullRankMinor( const Field& F, const size_t rank,
				      typename Field::Element * A_factors, const size_t lda,
				      const size_t* QtPointer,
				      typename Field::Element * X, const size_t ldx){
		
		typename Field::Element one, zero;
		F.init(one,1);
		F.init(zero,0);

		// upper entries are okay, just need to move up bottom ones
		const size_t* srcRow = QtPointer;
		for (size_t row=0; row<rank; row++, srcRow++) 
			if (*srcRow != row) {
				typename Field::Element* oldRow = A_factors + (*srcRow) * lda;
				typename Field::Element* newRow = A_factors + row * lda;
				for (size_t col=0; col<row; col++, oldRow++, newRow++) 
					F.assign(*newRow, *oldRow); 
			}
		
		// X <- (Qt.L.Q)^(-1)
		invL( F, rank, A_factors, lda, X, ldx); 

		// X = U^-1.X
		ftrsm( F, FflasLeft, FflasUpper, FflasNoTrans, 
		       FflasNonUnit, rank, rank, one, A_factors, lda, X, ldx); 

		return X;
		
 	}
	

	//---------------------------------------------------------------------
	// TURBO: rank computation algorithm 
	//---------------------------------------------------------------------
	template <class Field>
	static size_t 
	TURBO( const Field& F, const size_t M, const size_t N,		
	       typename Field::Element * NW, const size_t ld1,
	       typename Field::Element * NE, const size_t ld2,
	       typename Field::Element * SW, const size_t ld3,
	       typename Field::Element * SE, const size_t ld4	 );

	//---------------------------------------------------------------------
	// LUdivine: LUP factorisation of A 
	// P is the permutation matrix stored in an array of indexes
	//---------------------------------------------------------------------
	
	template <class Field>
	static size_t 
	LUdivine( const Field& F, const enum FFLAS_DIAG Diag,
		  const size_t M, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  size_t* P, size_t* Q, const enum FFPACK_LUDIVINE_TAG LuTag=FfpackLQUP);
        
	
	// Apply a permutation submatrix of P (between ibeg and iend) to a matrix
	// to (iend-ibeg) vectors of size M stored in A (as column for NoTrans 
	// and rows for Trans)
	// Side==FflasLeft for row permutation Side==FflasRight for a column 
	// permutation
	// Trans==FflasTrans for the inverse permutation of P
	template<class Field>
	static void 
	applyP( const Field& F, 
		const enum FFLAS_SIDE Side,
		const enum FFLAS_TRANSPOSE Trans,
		const size_t M, const int ibeg, const int iend,
		typename Field::Element * A, const size_t lda, const size_t * P ){
		
		if ( Side == FflasRight )
			if ( Trans == FflasTrans ){
				for ( size_t i=ibeg; i<(size_t) iend; ++i){
					if ( P[i]> i )
						fswap( F, M, 
						       A + P[i]*1, lda, 
						       A + i*1, lda );
				}
			}
			else{ // Trans == FflasNoTrans
				for (int i=iend-1; i>=ibeg; --i){
					if ( P[i]>(size_t)i ){
						fswap( F, M, 
						       A + P[i]*1, lda, 
						       A + i*1, lda );
					}
				}
			}
		else // Side == FflasLeft
			if ( Trans == FflasNoTrans ){
				for (size_t i=ibeg; i<(size_t)iend; ++i){
					if ( P[i]> (size_t) i )
						fswap( F, M, 
						       A + P[i]*lda, 1, 
						       A + i*lda, 1 );
				}
			}
			else{ // Trans == FflasTrans
				for (int i=iend-1; i>=ibeg; --i){
					if ( P[i]> (size_t) i ){
						fswap( F, M, 
						       A + P[i]*lda, 1, 
						       A + i*lda, 1 );
					}
				}
			}
			
	}
	//---------------------------------------------------------------------
	// CharPoly: Compute the characteristic polynomial of A using Krylov
	// Method, and LUP factorization of the Krylov Base
	//---------------------------------------------------------------------
	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	CharPoly( const Field& F, std::list<Polynomial>& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  const enum FFPACK_CHARPOLY_TAG CharpTag= FfpackHybrid);
	
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
		 typename Field::Element* X, const size_t ldx, size_t* P,
		 const enum FFPACK_MINPOLY_TAG MinTag,
		 const size_t kg_mc, const size_t kg_mb, const size_t kg_j );


	// Solve L X = B or X L = B in place
	// L is M*M if Side == FflasLeft and N*N if Side == FflasRigt, B is M*N.
	// Only the R non trivial column of L are stored in the M*R matrix L
	// Requirement :  so that L could  be expanded in-place
	template<class Field>
	static void
	solveLB( const Field& F, const enum FFLAS_SIDE Side,
		 const size_t M, const size_t N, const size_t R, 
		 typename Field::Element * L, const size_t ldl, 
		 const size_t * Q,
		 typename Field::Element * B, const size_t ldb ){
		
		 typename Field::Element one, zero;
		F.init(one, 1);
		F.init(zero, 0);
		for (int i=R-1; i>=0; --i){
			if (  Q[i] > (size_t) i){
				//for (size_t j=0; j<=Q[i]; ++j)
				//F.init( *(L+Q[i]+j*ldl), 0 );
				fcopy( F, M-Q[i]-1, L+Q[i]*(ldl+1)+ldl,ldl, L+(Q[i]+1)*ldl+i, ldl );
				for ( size_t j=Q[i]*ldl; j<M*ldl; j+=ldl)
					F.assign( *(L+i+j), zero );
			}
		}
		ftrsm( F, Side, FflasLower, FflasNoTrans, FflasUnit, M, N, one, L, ldl , B, ldb);
		
		// Undo the permutation of L
		for (size_t i=0; i<R; ++i){
			if ( Q[i] > (size_t) i){
				//for (size_t j=0; j<=Q[i]; ++j)
				//F.init( *(L+Q[i]+j*ldl), 0 );
				fcopy( F, M-Q[i]-1, L+(Q[i]+1)*ldl+i, ldl, L+Q[i]*(ldl+1)+ldl,ldl );
				for ( size_t j=Q[i]*ldl; j<M*ldl; j+=ldl)
					F.assign( *(L+Q[i]+j), zero );
			}
		} 
	}
	
	// Solve L X = B in place
	// L is M*M or N*N, B is M*N.
	// Only the R non trivial column of L are stored in the M*R matrix L
	template<class Field>
	static void
	solveLB2( const Field& F, const enum FFLAS_SIDE Side,
		  const size_t M, const size_t N, const size_t R, 
		  typename Field::Element * L, const size_t ldl, 
		  const size_t * Q,
		  typename Field::Element * B, const size_t ldb ){
		

		
		typename Field::Element Mone, one;
		F.init( Mone, -1 );
		F.init( one, 1 );
		typename Field::Element * Lcurr,* Rcurr,* Bcurr;
		size_t ib, k, Ldim;
		//cerr<<"In solveLB"<<endl;
		if ( Side == FflasLeft ){
			size_t j = 0;
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
				ftrsm( F, Side, FflasLower, FflasNoTrans, FflasUnit, Ldim, N, one, Lcurr, ldl , Bcurr, ldb );
				//cerr<<"M,k="<<M<<" "<<k<<endl;
				//cerr<<" fgemm with M, N, K="<<M-k<<" "<<N<<" "<<Ldim<<endl;
				fgemm( F, FflasNoTrans, FflasNoTrans, M-k, N, Ldim, Mone, Rcurr , ldl, Bcurr, ldb, one, Bcurr+Ldim*ldb, ldb);
			}
		}
		else{ // Side == FflasRight
			int j=R-1;
			while ( j >=0 ) {
				k = ib = Q[j];
				while ( (j>=0) &&  (Q[j] == k)  ) {--k;--j;}
				Ldim = ib-k;
				Lcurr = L + j+1 + (k+1)*ldl;
				Bcurr = B + ib;
				Rcurr = Lcurr + Ldim*ldl;
				fgemm( F, FflasNoTrans, FflasNoTrans, M, N-ib-1, Ldim, Mone, Bcurr, ldb, Rcurr, ldl,  one, Bcurr-Ldim, ldb);
				//cerr<<"j avant="<<j<<endl;
				//cerr<<"k, ib, j, R "<<k<<" "<<ib<<" "<<j<<" "<<R<<endl;
				//cerr<<"M,k="<<M<<" "<<k<<endl;
				//cerr<<" ftrsm with M, N="<<Ldim<<" "<<N<<endl;
				ftrsm( F, Side, FflasLower, FflasNoTrans, FflasUnit, M, Ldim, one, Lcurr, ldl , Bcurr-Ldim, ldb );
				//cerr<<"M,k="<<M<<" "<<k<<endl;
				//cerr<<" fgemm with M, N, K="<<M-k<<" "<<N<<" "<<Ldim<<endl;
			}
		}
	}

	template<class Field>
	static void trinv_left( const Field& F, const size_t N, const typename Field::Element * L, const size_t ldl,
				typename Field::Element * X, const size_t ldx ){
		invL(F,N,L,ldl,X,ldx);
	}
	

protected:
	

	// Inversion of a lower triangular matrix with a unit diagonal
	template<class Field>
	static void 
	invL( const Field& F, const size_t N, const typename Field::Element * L, const size_t ldl,
	      typename Field::Element * X, const size_t ldx ){
		//assumes X2 is initialized to 0
		typename Field::Element mone, one;
		F.init(one,1UL);
		F.init(mone,-1);
		
		if (N == 1){
			F.assign(*X, one);
		}
		else{
			size_t N1 = N >> 1;
			size_t N2 = N-N1;
			typename Field::Element * X11 = X;
			const typename Field::Element * L11 = L;
			typename Field::Element * X21 = X+N1*ldx;
			const typename Field::Element * L21 = L+N1*ldl;
			typename Field::Element * X22 = X21+N1;
			const typename Field::Element * L22 = L21+N1;
			// recursive call for X11
			// X11 = L11^-1
			invL( F, N1, L11, ldl, X11, ldx );
#if DEBUG==2
			cerr<<"Apres invl1 L11^-1="<<endl;
			write_field(F,cerr,X11,N1,N1,N);
#endif

			// recursive call for X11
			// X22 = L22^-1
			invL( F, N2, L22, ldl, X22, ldx );
#if DEBUG==2
			cerr<<"Apres invl2 L22^-1="<<endl;
			write_field(F,cerr,X22,N1,N1,N);
#endif
			
			// Copy L21 into X21
			for ( size_t i=0; i<N2; ++i)
				fcopy( F, N1, X21+i*ldx, 1, L21+i*ldl, 1 );

			// X21 = X21 . -X11^-1 (pascal 2004-10-12, make the negation
			// after the multiplication, problem in ftrmm)
			ftrmm (F, FflasRight, FflasLower, FflasNoTrans, FflasUnit, 
			       N2, N1, one, X11, ldx, X21, ldx );
			for (size_t i=0; i<N2; ++i)
				for (size_t j=0; j<N1; ++j)
					F.negin(*(X21+i*ldx+j));

#if DEBUG==2
			cerr<<"Apres trmm1 X21^-1="<<endl;
			write_field(F,cerr,X21,N2,N1,N);
#endif
			// X21 = X22^-1 . X21
			ftrmm( F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit, N2, N1, one, X22, ldx, X21, ldx );
#if DEBUG==2
			cerr<<"Apres trmm2 X21^-1="<<endl;
			write_field(F,cerr,X21,N2,N1,N);
#endif
		}
	}
		
		// Compute the new d after a LSP ( d[i] can be zero )
	template<class Field>
	static size_t 
	newD( const Field& F, size_t * d, bool& KeepOn,
	      const size_t l, const size_t N, 
	      typename Field::Element * X,
	      const size_t* Q,
	      std::vector<std::vector<typename Field::Element> >& minpt);

	template<class Field>
	static size_t
	updateD(const Field& F, size_t * d, size_t k,
		std::vector<std::vector<typename Field::Element> >& minpt );
	
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
			while ( F.isZero(*(Ai-x)) ) { // test if the pivot is 0
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
			    const typename Field::Element * A, const size_t lda,
			    typename Field::Element * X, const size_t ldx,
			    typename Field::Element * u, size_t* P,
			    bool computeX, const enum FFPACK_MINPOLY_TAG MinTag,
			    const size_t kg_mc, const size_t kg_mb, const size_t kg_j );
		
	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	KellerGehrig( const Field& F, std::list<Polynomial>& charp, const size_t N,
		      const typename Field::Element * A, const size_t lda );

	template <class Field, class Polynomial>
	static int
	KGFast ( const Field& F, std::list<Polynomial>& charp, const size_t N,
		 typename Field::Element * A, const size_t lda, 
		 size_t * kg_mc, size_t* kg_mc, size_t* kg_j );

	template<class Field>
	static void 
	fgemv_kgf( const Field& F,  const size_t N, 
		   const typename Field::Element * A, const size_t lda,
		   const typename Field::Element * X, const size_t incX,
		   typename Field::Element * Y, const size_t incY, 
		   const size_t kg_mc, const size_t kg_mb, const size_t kg_j );

	template <class Field, class Polynomial>
	static std::list<Polynomial>& 
	LUKrylov( const Field& F, std::list<Polynomial>& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  typename Field::Element * U, const size_t ldu,
		  const enum FFPACK_CHARPOLY_TAG CharpTag);

	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	Danilevski (const Field& F, std::list<Polynomial>& charp, 
		    const size_t N, typename Field::Element * A, const size_t lda);
		
	template <class Field, class Polynomial>
	static std::list<Polynomial>&
	LUKrylov_KGFast( const Field& F, std::list<Polynomial>& charp, const size_t N,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * X, const size_t ldx);
		
};

}
#include "linbox/ffpack/ffpack_ludivine.inl"
#include "linbox/ffpack/ffpack_minpoly.inl"
#include "linbox/ffpack/ffpack_charpoly_kglu.inl"
#include "linbox/ffpack/ffpack_charpoly_kgfast.inl"
#include "linbox/ffpack/ffpack_charpoly_danilevski.inl"
#include "linbox/ffpack/ffpack_charpoly.inl"

#endif // __FFPACK_H