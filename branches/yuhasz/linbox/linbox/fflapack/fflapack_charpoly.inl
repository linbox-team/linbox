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
list<Polynomial>&
FFLAPACK::CharPoly( const Field& F, list<Polynomial>& charp, const size_t N,
		 const typename Field::Element * A, const size_t lda,
		 typename Field::Element * U, const size_t ldu){
	
	typedef typename Field::Element elt;
	Polynomial minP;
	const elt* Ai;
	elt* A2i, *Xi;
	int  j;
	size_t k;
	size_t P[N];
	static elt Mone;
	F.neg(Mone,F.one);

#if DEBUG==2	
	cerr<<"U="<<endl;
	write_field(F,cerr,U,N,N,ldu);
	cerr<<"A="<<endl;
	write_field(F,cerr,A,N,N,lda);
#endif

#if DEBUG==2
	cerr<<"Computing the Minpoly...";
#endif
#if DEBUG
	cerr<<".";
#endif
	// X contains the LSP factorisation of U, the Krylov Matrix
	elt* X = new elt[N*(N+1)]; 
	MinPoly( F, minP, N, A, lda, U, ldu, X, N, P );	
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif
	k = minP.size()-1; // degre of minpoly
	if ( k==1 && F.iszero( minP[0] ) ){ // minpoly is X
		Ai = A;
		j = N*N;
		while ( j-- && F.isZero(*(Ai++)) ){}
		if ( j<0 ){ // A is 0, CharPoly=X^n
#if DEBUG==2
			cerr<<"Matrix is 0"<<endl;
#endif
			minP.resize(N+1);
			minP[1] = F.zero;
			minP[N] = F.one;
			
			k=N;
		}
	}

	if ( k==N ){
		charp.clear();
		charp.push_front(minP); // CharPoly = MinPoly
#if DEBUG==2	
		cerr<<"Charpoly==Minpoly"<<endl;
		cerr<<"k="<<k<<endl;
#endif			
		return charp;
	}
	
	size_t Nrest = N-k;
	elt * X21 = X + k*N;
        elt * X22 = X21 + k;
	
	// Apply P on rows and on columns of A12 and A22: X2_=(((PA_2)^tP^-1))
	//  X2_ <- ((A_2^t.P))
#if DEBUG==2
	cerr<<"Applying first permutation and copy...";
#endif			
	//	cerr<<"Before Flaswp copy X="<<endl;
	//	write_field(F,cerr,X,N+1,N,N);

	// Copy X2_ <- (A_2)^t
	flaswp( F, N, const_cast<typename Field::Element* &>(A), N, 0, k, P, 1);
	for ( Xi = X21, Ai = A+k;
	      Xi != X21 + N*Nrest;
	      Ai++ ){
		for ( j=0; j<N*lda; j+=lda ){
			*(Xi++) = *(Ai+j);
		}
	}
	// Undo the permutation
	flaswp( F, N,const_cast<typename Field::Element* &>( A), N, 0, k, P, -1);
	
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif			

	//	cerr<<"After flaswp and copy X="<<endl;
	//	write_field(F,cerr,X,N,N,N);
	
	// X2_ <- X2 . P^t 
#if DEBUG==2
	cerr<<"Applying second permutation...";
#endif			
	flaswp( F, Nrest, X21, N, 0, k, P, 1);  
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif			
	
	//	cerr<<"After flaswp X="<<endl;
	//	write_field(F,cerr,X,N,N,ldu);
	
	// A12^t <= A12^t*L1^-1^t ( in X21 )

	//	cerr<<"Before Ftrsm X="<<endl;
	//	write_field(F,cerr,X,N+1,N,N);


#if DEBUG==2
	cerr<<"Applying Ftrsm...";
#endif			
	ftrsm(F, FflasRight, FflasUpper, FflasNoTrans, FflasUnit, Nrest, k,
	      F.one, X, N, X21, N);  
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif			
	
	//	cerr<<"After Ftrsm X="<<endl;
	//	write_field(F,cerr,X,N+1,N,N);
	
	// Creation of the matrix A2 for recurise call 
#if DEBUG==2
	cerr<<"Computing A2 for recursive call...";
#endif			
	//X22 = X+(N+1)*(k);
	elt * A2 = new elt[Nrest*Nrest];
	for ( Xi = X22,  A2i = A2;
	      Xi != X22 + ldu*Nrest;
	      Xi += N-Nrest )
		for ( j=Nrest; j; --j ){
			*(A2i++) = *(Xi++);
		}
	
	fgemm( F, FflasNoTrans, FflasNoTrans, Nrest, Nrest, k, Mone,
		     X21, N, X+k, N, F.one, A2, Nrest, 0);
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif			
#if DEBUG==2	
	cerr<<"A2="<<endl;
	write_field(F,cerr,A2,Nrest,Nrest,Nrest);
#endif
	 // Recursive call on X22
	CharPoly(F, charp, Nrest, A2, Nrest, U+k*(ldu+1), ldu );
	delete[] A2;

	charp.push_front( minP ); 
	return charp;			
}

