/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflapack/fflapack_charpoly.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

//---------------------------------------------------------------------
// CharPoly: Compute the characteristic polynomial of A using Krylov
// Method, and LUP factorization of the Krylov Base
// Create a (N-k)*(N-k) matrix for a recursive call, were k is the degree
// of the minpoly(A,v)
//---------------------------------------------------------------------
template <class Field, class Polynomial>
std::list<Polynomial>&
LinBox::FFLAPACK::CharPoly( const Field& F, std::list<Polynomial>& charp, const size_t N,
			    const typename Field::Element * A, const size_t lda,
			    typename Field::Element * U, const size_t ldu,
			    const enum FFLAPACK_CHARPOLY_TAG CharpTag ){
	switch ( CharpTag ) {
	case FflapackLUK: 
		return LUKrylov( F, charp, N, A, lda, U, ldu, FflapackLUK );
	case FflapackHybrid: 
		return LUKrylov( F, charp, N, A, lda, U, ldu, FflapackHybrid );
	case FflapackKG:
		return KellerGehrig( F, charp, N, A, lda, U, ldu );
	default:
		return LUKrylov( F, charp, N, A, lda, U, ldu, FflapackHybrid );
	}
}

template <class Field, class Polynomial>
std::list<Polynomial>&
LinBox::FFLAPACK::LUKrylov( const Field& F, std::list<Polynomial>& charp, const size_t N,
		    const typename Field::Element * A, const size_t lda,
		    typename Field::Element * U, const size_t ldu,
		    const enum FFLAPACK_CHARPOLY_TAG CharpTag){
	
	typedef typename Field::Element elt;
	Polynomial *minP = new Polynomial();
	const elt* Ai;
	elt* A2i, *Xi;
	static elt Mone, one, zero;
	F.init(zero,0UL);
	F.init(one, 1UL);
	F.neg(Mone,one);

	size_t P[N];
#ifdef __MINP_CONSTRUCT
	// X contains the LSP factorisation of U, the Krylov Matrix
	elt* X = new elt[N*(N+1)]; 
	size_t ldx = N;
	MinPoly( F, *minP, N, A, lda, U, ldu, X, N, P );
#else
	elt* X = U;
	size_t ldx = ldu;
	MinPoly( F, *minP, N, A, lda, X, ldx, P );
#endif
	
	size_t k = minP->size()-1; // degre of minpoly
	if ( k==1 && F.isZero( (*minP)[0] ) ){ // minpoly is X
		Ai = A;
		int j = N*N;
		while ( j-- && F.isZero(*(Ai++)) );
		if ( !j ){ // A is 0, CharPoly=X^n
			minP->resize(N+1);
			(*minP)[1] = zero;
			(*minP)[N] = one;
			k=N;
		}
	}
	
	if ( k==N ){
		charp.clear();
		charp.push_back(*minP); // CharPoly = MinPoly
		return charp;
	}
	
	size_t Nrest = N-k;
	elt * X21 = X + k*ldx;
        elt * X22 = X21 + k;
	
	//  Compute the n-k last rows of A' = PA^tP^t in X2_
	
	// A = A . P^t
	applyP( F, FflasRight, FflasTrans, N, 0, k, 
		const_cast<typename Field::Element* &>(A), lda, P);
	//flaswp( F, N, const_cast<typename Field::Element* &>(A), N, 0, k, P, 1);
	
	// Copy X2_ = (A'_2)^t
	for ( Xi = X21, Ai = A+k; Xi != X21 + Nrest*ldx; Ai++, Xi+=ldx-N ){
		for ( size_t jj=0; jj<N*lda; jj+=lda ){
			*(Xi++) = *(Ai+jj);
		}
	}

	// A = A . P : Undo the permutation on A
	applyP( F, FflasRight, FflasNoTrans, N, 0, k, 
		const_cast<typename Field::Element* &>(A), lda, P);
	//flaswp( F, N,const_cast<typename Field::Element* &>( A), N, 0, k, P, -1);
	
	// X2_ = X2_ . P^t (=  ( P A^t P^t )2_ ) 
	applyP( F, FflasRight, FflasTrans, Nrest, 0, k, X21, ldx, P);
	//flaswp( F, Nrest, X21, N, 0, k, P, 1);  
	

	// X21 = X21 . S1^-1
	ftrsm(F, FflasRight, FflasUpper, FflasNoTrans, FflasUnit, Nrest, k,
	      one, X, ldx, X21, ldx);  
	
	// Creation of the matrix A2 for recurise call 
	elt * A2 = new elt[Nrest*Nrest];
	
	for ( Xi = X22,  A2i = A2;
	      Xi != X22 + Nrest*ldx;
	      Xi += (ldx-Nrest) ){
		for ( size_t jj=0; jj<Nrest; ++jj ){
			*(A2i++) = *(Xi++);
		}
	}
	fgemm( F, FflasNoTrans, FflasNoTrans, Nrest, Nrest, k, Mone,
	       X21, ldx, X+k, ldx, one, A2, Nrest, 0);
	
#ifdef __MINP_CONSTRUCT
	delete[] X;
#endif
	// Recursive call on X22
	charp.push_back( *minP );
	
	if ( (CharpTag == FflapackHybrid) && (k < (N>>3) ) )
		KellerGehrig( F, charp, Nrest, A2, Nrest, U+k*(ldu+1), ldu );
	else
		LUKrylov( F, charp, Nrest, A2, Nrest, U+k*(ldu+1), ldu, CharpTag );

	delete[] A2;
	return charp;
}

