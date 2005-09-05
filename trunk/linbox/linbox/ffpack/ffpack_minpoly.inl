/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/ffpack/ffpack_minpoly_construct.inl
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
LinBox::FFPACK::MinPoly( const Field& F, Polynomial& minP, const size_t N,
			 const typename Field::Element *A, const size_t lda,
			 typename Field::Element* X, const size_t ldx,
			 size_t* P, const enum FFPACK_MINPOLY_TAG MinTag = FfpackDense,
			 const size_t kg_mc =0, const size_t kg_mb=0, const size_t kg_j=0 ){

	typedef typename Field::Element elt;
	static elt one,zero;
	F.init( one, 1UL );
	F.init( zero, 0UL );
	// nRow is the number of row in the krylov base already computed
	size_t j, k, nRow = 2;
	typename Polynomial::iterator it;
	elt* Xi, *Ui;
	typename Field::RandIter g (F);
	bool KeepOn=true;

#if DEBUG==2
	for (j=0;j<(N+1)*N;j++)
		X[j] = zero;
#endif
	
	// Vector to store the Krylov iterate A^i.v
	elt* U = new elt[N];
	// Picking a non zero vector
	do{
		for (Ui=U, Xi = X; Ui<U+N; ++Ui, ++Xi){
			g.random (*Ui);
		 	*Xi = *Ui;
			if (!F.isZero(*Ui))
				KeepOn = false;
		}
	}while(KeepOn);
	
	nRow = 1;
	
	// LUP factorization of the Krylov Base Matrix
	k = LUdivine_construct (F, FflasUnit, N+1, N, A, lda, X, ldx, U, P, true,
				MinTag, kg_mc, kg_mb, kg_j);

	//delete[] U;
	minP.resize(k+1);
	minP[k] = one;
	if ( (k==1) && F.isZero(*(X+ldx))){ // minpoly is X
		delete[] U;
		for (size_t i=0; i<k; ++i)
			minP[i] = zero;
		return minP;
	}
	// U contains the k first coefs of the minpoly
	//elt* m= new elt[k];
	fcopy( F, k, U, 1, X+k*ldx, 1);
	ftrsv( F, FflasLower, FflasTrans, FflasNonUnit, k, X, ldx, U, 1);
	it = minP.begin();
	for (j=0; j<k; ++j, it++){
		F.neg(*it, U[j]);
	}
	delete[] U;
	return minP;
}
