/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//----------------------------------------------------------------------
//                FFLAP: Finite Field Fast Linear Algebra Package
//                FFLAP_CharPoly_KGLU: Characteristic Polynomial of A
//		  using Kelleg-Gehrig M(n)logn branching algorithm
//----------------------------------------------------------------------
// by Clement PERNET ( clement.pernet@imag.fr )
// 09/05/2003
//----------------------------------------------------------------------

//---------------------------------------------------------------------
// CharPoly: Compute the characteristic polynomial of A 
//---------------------------------------------------------------------
template <class Field, class Polynomial>
list<Polynomial>&
FFLAP::CharPoly( const Field& F, list<Polynomial>& charp, const size_t N,
	  const typename Field::element * A, const size_t lda,
	  typename Field::element * U, const size_t ldu){
	
	typedef typename Field::element elt;
	const elt * Ai;
	elt *Ui, *Uj, *Uk, *Ukp1, *Ukp1new, *Bi, *Vi, *Vk, *Xi, *Xj;
	elt * B = new elt[N*N];     // to store A^2^i
	elt * V = new elt[N*N];     // to store A^2^i.U
	elt * X = new elt[2*N*N];   // to compute the LSP factorization
	size_t * P = new size_t[N]; // Permutation matrix for LSP
	size_t * d = new size_t[N]; // dimensions of Vect(ei, Aei...)
	size_t * dold = new size_t[N]; // copy of d
	
	size_t i=0, l=1, j, k,s, newRowNb, nrowX;
	bool  KeepOn;

	cerr<<"A="<<endl;
	write_field( F,cerr,A,N,N,lda);
	Ai = A; 
	Xi = X;
	for ( i=0; i<N; ++i)
		dold[i] = d[i] = 1;
	// Computing the first X: (e1; e1A; e2; e2A;...;en;enA)
	for ( i=0, Ui=U, Vi=V, Bi=B; i<N; ++i, Ai -= N*lda-1, Ui+=ldu  ){
		for ( Xj=Xi, Uj=Ui; Xj<Xi+N; ++Xj, ++Uj){
			*Xj = F.zero;
			*Uj = F.zero;
		}
		*(Ui+i) = *(Xi+i) = F.one;
		while ( Xj<Xi+2*N) {
			*(Bi++) = *(Xj++) = *(Vi++) = *Ai;
			Ai+=lda;
		}
		Xi = Xj;
	}
	cerr<<"U="<<endl;
	write_field(F,cerr,U,N,N,N);
	cerr<<"V="<<endl;
	write_field(F,cerr,V,N,N,N);
	cerr<<"B="<<endl;
	write_field(F,cerr,B,N,N,N);
	cerr<<"X="<<endl;
	write_field(F,cerr,X,2*N,N,N);
	// step form elimination using LSP
	for ( i=0;i<N;++i) 
		P[i]=0;
	LUdivine( F, FflasUnit, 2*N, N, X, N, P, FflapLSP );

	cerr<<"After LSP X="<<endl;
	write_field(F,cerr,X,2*N,N,N);

	// updating d
	Xi = X;
	KeepOn = false;
	for ( i=0; i<N; ++i){
		j = 0;
		s = ( d[i]==l )? 2*l : d[i];  
		while ( !F.iszero(*Xi) && j<s){
			Xi += N+1; 
			j++;
		}
		Xi += (s-j)*N;
		d[i] = j;
		if ( j==2*l)
			KeepOn = true;
	}

	cerr<<"d=";
	for (i=0;i<N;++i) cerr<<d[i]<<" ";
	cerr<<endl;
	// Updating U:
	Uk = U;
	// Firstly, removing extra rows
	for ( i = 0; i<N; ++i){
		if (d[i] < dold[i]){
			Ukp1new = Uk +d[i];  // new position of first row of Uk+1
			Ukp1 = Uk + dold[i]; // first row of Uk+1
			Ui = Ukp1new;
			Uj = Ukp1;
			while ( Uj < U + N*ldu ){
				for ( j=N; j; --j)
					*(Ui++) = *(Uj++);
				Uj+=ldu-N;
				Ui+=ldu-N;
			}
			Uk = Ukp1new;
			dold[i] = d[i];
		}
		else 
			Uk += dold[i];
	}
	cerr<<"After removing rows U="<<endl;
	write_field(F,cerr,U,N,N,N);
	
	// Then inserting the duplicated blocks
	Uk = U;
	Vk = V;
	for ( i = 0; Uk<U+N*ldu; ++i){
		Ukp1 = Uk + dold[i]; // first row of Uk+1
		newRowNb = d[i] - dold[i];
		if ( newRowNb > 0){
			Ui = U+N*ldu-1; // last row of U
			Uj = U+(N-newRowNb)*ldu-1; // last row future position
			while ( Uj >= Ukp1+N-1){
				for ( j=N;j;--j)
					*(Ui--) = *(Uj--);   // moving block
				cerr<<"After transfert of 1 row:"<<endl;
				write_field(F,cerr,U,N,N,N);
				Ui-=ldu-N;
				Uj-=ldu-N;
			}
			// copying the d[i]-dold[i] new rows from Vi
			Vi = Vk;
			Uj++;
			cerr<<"Then copy of V block"<<endl;
			while ( Uj < Ukp1+ldu*newRowNb ){
				for ( j=N;j;--j)
					*(Uj++) = *(Vi++);
				Uj+=ldu-N;
			}
			dold[i] = d[i];
		}
		Uk = Uk + d[i]*ldu;
		Vk += l*N;
	}
	cerr<<"After inserting rows U="<<endl;
	write_field(F,cerr,U,N,N,N);
	while(KeepOn){ // Main loop, until each subspace dimension has been found

		// block size of X, U, V is l=2^i
		l*=2;
		// B = A^2^i
		fsquare( F, FflasNoTrans, N, F.one, B, N, F.zero, B, N );

		// V = U.B^t
		fgemm( F, FflasNoTrans, FflasNoTrans, N, N, N, F.one, 
		       U, ldu, B, N, F.zero, V, N);
		
		
		// X = ( U1, V1, U2, V2, ... )
		Xi = X;
		Ui = U;
		Vi = V;
		k=0;
		nrowX = 0;
		cerr<<"Loop for l="<<l<<" k="<<k<<endl;
		cerr<<"V="<<endl;
		write_field(F,cerr,V,N,N,N);
		cerr<<"B="<<endl;
		write_field(F,cerr,B,N,N,N);
		while ( Vi < V + N*N ){
			// Copying Uk
			cerr<<"d["<<k<<"]="<<d[k]<<endl;
			for ( i = d[k]; i; --i, nrowX++){ 
				for ( j = N; j; --j )
					*(Xi++) = *(Ui++);
				Ui += ldu - N;
			}
			cerr<<"After copying Uk to X="<<endl;
			write_field(F,cerr,X,2*N,N,N);

			// Copying Vk
			if ( d[k] == l )
				for ( i=d[k]; i; i--, nrowX++)
					for ( j=N; j; j--)
						*(Xi++) = *(Vi++);
			
			else
				Vi = Vi + N*d[k];
			cerr<<"After copying Vk to X="<<endl;
			write_field(F,cerr,X,2*N,N,N);
			k++;
		} // at this point, k is the number of blocks ( Ui, Vi ) in X
		
		cerr<<"U="<<endl;
		write_field(F,cerr,U,N,N,N);
	
		cerr<<"X="<<endl;
		write_field(F,cerr,X,2*N,N,N);


		// step form elimination of X using LSP
		for ( i=0;i<N;++i) 
			P[i]=0;
		LUdivine( F, FflasUnit, nrowX, N, X, N, P, FflapLSP );
		
		// Updating d
		Xi = X;
		KeepOn = false;
		size_t dtot = 0;
		for ( i=0; i<k; ++i){
			j = 0;
			s = ( d[i]==l )? 2*l : d[i];
			if (s > N-dtot) 
				s= N-dtot; 
			cerr<<"computing d["<<i<<"], s="<<s<<endl;
			
			while ( !F.iszero(*Xi) && j<s){
				Xi += N+1; 
				j++;
			}
			Xi += (s-j)*N;
			d[i] = j;
			dtot += j;
			if ( j==2*l)
				KeepOn = true;
		} // Invariant: sum(d[i],i=0..k-1) = n
		
		cerr<<"d=";
		for (i=0;i<N;++i) cerr<<d[i]<<" ";
		cerr<<endl;
		// Updating U:
		Uk = U;
		// Firstly, removing extra rows
		for ( i = 0; i<N; ++i){
			if (d[i] < dold[i]){
				Ukp1new = Uk +d[i];  // new position of Uk+1
				Ukp1 = Uk + dold[i]; // first row of Uk+1
				Ui = Ukp1new;
				Uj = Ukp1;
				while ( Uj < U + N*ldu ){
					for ( j=N; j; --j)
						*(Ui++) = *(Uj++);
					Uj+=ldu-N;
					Ui+=ldu-N;
				}
				Uk = Ukp1new;
				dold[i] = d[i];
			}
			else 
				Uk += dold[i];
		}
		// Then inserting the duplicated blocks
		Uk = U;
		Vk = V;
		for ( i = 0; i<N; ++i){
			Ukp1 = Uk + dold[i]; // first row of Uk+1
			newRowNb = d[i] - dold[i];
			if ( newRowNb > 0){
				Ui = U+N*ldu-1; // last row of U
				Uj = U+(N-newRowNb)*ldu-1; // last row future pos
				while ( Uj >= Ukp1){
					for ( j=N;j;--j)
						*(Ui--) = *(Uj--);// moving block
					Ui-=ldu+N;
					Uj-=ldu+N;
				}
				// copying the d[i]-dold[i] new rows from Vi
				Vi = Vk;
				while ( Uj < Ukp1+ldu*newRowNb ){
					for ( j=N;j;--j)
						*(Uj++) = *(Vi++);
					Uj+=ldu-N;
				}
				dold[i] = d[i];
			}
			Uk = Uk + d[i]*ldu;
			Vk += l*N;
		}
		l<<1;
	}
	delete[] X;
	delete[] B;
	delete[] V;
	delete[] P;
	delete[] dold;
	// Constructing the CharPoly
	charp.clear();
	elt * m = new elt[N];
	for ( i=0; i<k; ++i){
		Polynomial minP(d[i]+1);
		minP[d[i]] = F.one;
		fcopy( F, d[i], m, 1, U+d[i]*ldu, 1);
		ftrsv( F, FflasLower, FflasTrans, FflasNonUnit, d[i], U, ldu, m, 1);
		typename Polynomial::iterator it = minP.begin();
		for (j=0; j<k; ++j, it++){
			F.neg(*it, m[j]);
		}
		charp.push_front( minP );
	}
	delete[] d;
	delete[] m;
}
