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
FFLAPACK::MinPoly( const Field& F, Polynomial& minP, const size_t N,
		const typename Field::Element *A, const size_t lda,
		typename Field::Element* U, size_t ldu, size_t* P){

	typedef typename Field::Element elt;
	// nRow is the number of row in the krylov base already computed
	size_t j, k, nNewRow, nRow = 2;
	elt* B = new elt[ N*N ];
	const elt* Ai=A;
	typename Polynomial::iterator it;
	elt* Ui, *Bi=B, *L, *Li, *y;
	typename Field::randIter g (F);
	bool KeepOn=true;
	// Picking a non zero vector
	do{
		for (Ui=U; Ui<U+N; ++Ui){
			g.random (*Ui);
		 	if (!F.iszero(*Ui))
				KeepOn = false;
				}
	}while(KeepOn);
	Ui = U+ldu;

#if DEBUG==2	
	cerr<<"U="<<endl;
	write_field(F,cerr,U,N,N,ldu);
	cerr<<"A="<<endl;
	write_field(F,cerr,A,N,N,lda);
	
#endif
	// Try memcopy here
	for (; Ai<A+lda*N; Ai+=lda-N)
		for ( j=0; j<N; ++j){
			*(Bi++) = *(Ai++);
		}
	// Computing vA
	// Try gemv here
	fgemm(F, FflasNoTrans, FflasTrans, 1, N, N, F.one,
	      U, ldu, B, N, F.zero, Ui, ldu, 0);
	Ui += ldu;

#if DEBUG==2	
	cerr<<"after first BX U="<<endl;
	write_field(F,cerr,U,N+1,N,ldu);
#endif
	while (nRow < N+1){ 
		// All available rows have been used, compute new ones 
		// number of new rows to be computed:
		nNewRow = MIN( nRow,  N+1 - nRow);
	        // B <= B*B
#if DEBUG==2
		cerr<<"nNewRow ="<<nNewRow<<" nRow = "
		    <<(nRow)<<endl;
#endif

		fsquare(F, FflasNoTrans, N, F.one, B, N, F.zero, B, N);
			
		// ( U <= ( U ; U.A^2^i )
		fgemm(F, FflasNoTrans, FflasTrans, nNewRow, N, N, F.one,
		      U, ldu, B, N, F.zero, Ui, ldu, 0);
#if DEBUG==2
		cerr<<"apres (U,BU), U="<<endl;
		write_field(F,cerr,U,N+1,N,ldu);
#endif				
		Ui += nNewRow*ldu;
		nRow += nNewRow;
	}
	
	delete[] B;
	// LUP factorization of the Krilov Base Matrix
	k = LUdivine(F, FflasUnit, N+1, N, U, ldu, P );
	minP.resize(k+1);
	minP[k] = F.one;
	if (k==1 && F.isZero(*(U+ldu))){ // minpoly is X
		return minP;
	}
	// m contains the k first coefs of the minpoly
	elt* mi,* m= new elt[k];
	fcopy( F, k, m, 1, U+k*ldu, 1);
	ftrsv( F, FflasLower, FflasTrans, FflasNonUnit, k, U, ldu, m, 1);
	it = minP.begin();
	for (j=0; j<k; ++j, it++){
		F.neg(*it, m[j]);
	}
	delete[] m;
	return minP;
}
