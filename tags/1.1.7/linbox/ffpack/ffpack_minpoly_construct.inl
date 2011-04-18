/* linbox/ffpack/ffpack_minpoly_construct.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ffpack_minpoly_construct_INL
#define __LINBOX_ffpack_minpoly_construct_INL

//#define LB_DEBUG
//---------------------------------------------------------------------
// MinPoly: Compute the minimal polynomial of (A,v) using an LUP 
// factorization of the Krylov Base (v, Av, .., A^kv)
// U must be (n+1)*n
//---------------------------------------------------------------------
template <class Field, class Polynomial>
Polynomial&
LinBox::FFPACK::MinPoly( const Field& F, Polynomial& minP, const size_t N,
		const typename Field::Element *A, const size_t lda,
		typename Field::Element* U, size_t ldu,typename Field::Element* X, size_t ldx,
		size_t* P)
{

	typedef typename Field::Element elt;
	static elt one,zero;
	F.init( one, 1UL );
	F.init( zero, 0UL );
	// nRow is the number of row in the krylov base already computed
	size_t j, k, nRow = 2;
	elt* B = new elt[ N*N ];
	typename Polynomial::iterator it;
	elt* Xi, *Ui;
	typename Field::RandIter g (F);
	bool KeepOn=true;
	// Creating the Krylov Base copy matrix X where to factorize 
	//elt * X = new elt[(N+1)*N];
#ifdef LB_DEBUG
	for (j=0;j<(N+1)*N;j++)
		X[j] = zero;
#endif
	// Creating the copy of A, where to compute A^2^i
	// Try memcopy here
	for (size_t i=0; i<N; ++i)
		for ( size_t j=0; j<N; ++j)
			*(B+i*N+j) = *(A+i*N+j);
	
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
	size_t nUsedRow = 0;
	// LUP factorization of the Krylov Base Matrix
	
	k = LUdivine_construct(F, FflasUnit, N+1, N, B, N, U, ldu, X, N, P,
			       &nRow, N+1, &nUsedRow );
	delete[] B;
	minP.resize(k+1);
	minP[k] = one;
	if (k==1 && F.isZero(*(X+N))){ // minpoly is X
		return minP;
	}
	// m contains the k first coefs of the minpoly
	elt* m= new elt[k];
	fcopy( F, k, m, 1, X+k*N, 1);
	ftrsv( F, FflasLower, FflasTrans, FflasNonUnit, k, X, N, m, 1);
	//delete[] X;
	it = minP.begin();
	for (j=0; j<k; ++j, it++){
		F.neg(*it, m[j]);
	}
	delete[] m;
	return minP;
}

#undef LB_DEBUG

#endif //__LINBOX_ffpack_minpoly_construct_INL
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
