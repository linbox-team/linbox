/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflapack/fflapack_minpoly.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

//---------------------------------------------------------------------
// MinPoly: Compute the minimal polynomial of (A,v) using an LUP 
// factorization of the Krylov Base (v, Av, .., A^kv)
// U must be (n+1)*n
//---------------------------------------------------------------------
template <class Field, class Polynomial>
Polynomial&
LinBox::FFLAPACK::MinPoly( const Field& F, Polynomial& minP, const size_t N,
			   const typename Field::Element *A, const size_t lda,
			   typename Field::Element* U, size_t ldu, size_t* P){

	typedef typename Field::Element elt;
	// nRow is the number of row in the krylov base already computed
	size_t  nNewRow, nRow = 2;
	elt* B = new elt[ N*N ];
	elt* Ui, *Bi=B;
	static elt one, zero;
	F.init(one, 1UL);
	F.init(zero, 0UL);
	typename Field::RandIter g (F);

	// Picking a vector in F\{0}
	for (Ui=U; Ui<U+N; ++Ui)
		while ( F.isZero (g.random (*Ui)) );
	
	Ui = U+ldu;
	
	// Try memcopy here
	const elt* Ai=A;
	for (; Ai<A+lda*N; Ai+=lda-N)
		for ( size_t j=0; j<N; ++j){
			*(Bi++) = *(Ai++);
		}
	// Computing vA^t
	// Try gemv here
	fgemv( F, FflasNoTrans, N, N, one, B, N, U, 1, zero, Ui, 1 );
	// 	fgemm(F, FflasNoTrans, FflasTrans, 1, N, N, one,
	// 	      U, ldu, B, N, zero, Ui, ldu, 0);
	Ui += ldu;

	while (nRow < N+1){ 
		// All available rows have been used, compute new ones 
		// number of new rows to be computed:
		nNewRow = MIN( nRow,  N+1 - nRow);

	        // B <= B*B
		fsquare(F, FflasNoTrans, N, one, B, N, zero, B, N);
			
		// ( U <= ( U ; U.A^2^i )
		fgemm(F, FflasNoTrans, FflasTrans, nNewRow, N, N, one,
		      U, ldu, B, N, zero, Ui, ldu, 0);
		
		Ui += nNewRow*ldu;
		nRow += nNewRow;
	}
	
	delete[] B;
	// LUP factorization of the Krylov Basis Matrix
	size_t Q[N]; 
	size_t k = LUdivine(F, FflasUnit, N+1, N, U, ldu, P, FflapackLQUP, Q );
	minP.resize(k+1);
	minP[k] = one;
	if (k==1 && F.isZero(*(U+ldu))){ // minpoly is X
		return minP;
	}

	elt* m= new elt[k];
	fcopy( F, k, m, 1, U+k*ldu, 1);
	ftrsv( F, FflasLower, FflasTrans, FflasNonUnit, k, U, ldu, m, 1);
	typename Polynomial::iterator it = minP.begin();
	for (size_t j=0; j<k; ++j, it++){
		F.neg(*it, m[j]);
	}
	delete[] m;
	return minP;
}
