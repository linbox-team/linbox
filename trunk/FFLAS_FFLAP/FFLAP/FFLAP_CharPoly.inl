/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//----------------------------------------------------------------------
//                FFLAP: Finite Field Fast Linear Algebra Package
//                FFLAP_CharPoly: Characteristic Polynomial of A
//----------------------------------------------------------------------
// by Clement PERNET ( clement.pernet@imag.fr )
// 30/04/2003
//----------------------------------------------------------------------

//---------------------------------------------------------------------
// CharPoly: Compute the characteristic polynomial of A using Krylov
// Method, and LUP factorization of the Krylov Base
// Create a (N-k)*(N-k) matrix for a recursive call, were k is the degree
// of the minpoly(A,v)
//---------------------------------------------------------------------
template <class Field, class Polynomial>
list<Polynomial>&
FFLAP::CharPoly( const Field& F, list<Polynomial>& charp, const size_t N,
	  const typename Field::element * A, const size_t lda,
	  typename Field::element * U, const size_t ldu){
	
	typedef typename Field::element elt;
	Polynomial minP;
	const elt* Ai;
	elt* A2i, *Ui;
	int  j;
	size_t k;
	size_t P[N];
	static elt Mone;
	F.neg(Mone,F.one);
	Timer ti;

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
	ti.clear();
	ti.start();
	MinPoly( F, minP, N, A, lda, U, ldu, P );	
	ti.stop();
	cerr<<"Time for  Minpoly():"<<ti.usertime()<<endl;
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
	elt * U21 = U + k*ldu;
	typename Field::element * U22 = U21 + k;
	
	// Apply P on rows and on columns of A12 and A22: U2_=(((PA_2)^tP^-1))
	//  U2_ <- ((A_2^t.P))
#if DEBUG==2
	cerr<<"Applying first permutation and copy...";
#endif			
	//	cerr<<"Before Flaswp copy U="<<endl;
	//	write_field(F,cerr,U,N+1,N,ldu);

	// Copy U2_ <- (A_2)^t
	flaswp( F, N, const_cast<typename Field::element* &>(A), N, 0, k, P, 1);
	for ( Ui = U21, Ai = A+k;
	      Ui != U21 + ldu*Nrest;
	      Ui += ldu-N, Ai++ ){
		for ( j=0; j<N*lda; j+=lda ){
			*(Ui++) = *(Ai+j);
		}
	}
	// Undo the permutation
	flaswp( F, N,const_cast<typename Field::element* &>( A), N, 0, k, P, -1);
	
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif			

	//	cerr<<"After flaswp and copy U="<<endl;
	//	write_field(F,cerr,U,N,N,ldu);
	
	// U2_ <- U2 . P^t 
#if DEBUG==2
	cerr<<"Applying second permutation...";
#endif			
	flaswp( F, Nrest, U21, ldu, 0, k, P, 1);  
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif			
	
	//	cerr<<"After flaswp U="<<endl;
	//	write_field(F,cerr,U,N,N,ldu);
	
	// A12^t <= A12^t*L1^-1^t ( in U21 )

	//	cerr<<"Before Ftrsm U="<<endl;
	//	write_field(F,cerr,U,N+1,N,ldu);


#if DEBUG==2
	cerr<<"Applying Ftrsm...";
#endif			
	ftrsm(F, FflasRight, FflasUpper, FflasNoTrans, FflasUnit, Nrest, k,
	      F.one, U, ldu, U21, ldu);  
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif			
	
	//	cerr<<"After Ftrsm U="<<endl;
	//	write_field(F,cerr,U,N+1,N,ldu);
	
	// Creation of the matrix A2 for recurise call 
#if DEBUG==2
	cerr<<"Computing A2 for recursive call...";
#endif			
	U22 = U+(ldu+1)*(k);
	elt * A2 = new elt[Nrest*Nrest];
	for ( Ui = U22,  A2i = A2;
	      Ui != U22 + ldu*Nrest;
	      Ui += ldu-Nrest )
		for ( j=Nrest; j; --j ){
			*(A2i++) = *(Ui++);
		}
	
	fgemm( F, FflasNoTrans, FflasNoTrans, Nrest, Nrest, k, Mone,
		     U21, ldu, U+k, ldu, F.one, A2, Nrest, 0);
#if DEBUG==2
	cerr<<"Ok"<<endl;
#endif			
#if DEBUG==2	
	cerr<<"A2="<<endl;
	write_field(F,cerr,A2,Nrest,Nrest,Nrest);
#endif
	 // Recursive call on U22
	CharPoly(F, charp, Nrest, A2, Nrest, U22, ldu );
	delete[] A2;

	charp.push_front( minP ); 
	return charp;			
}

