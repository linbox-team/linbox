/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflapack/fflapack_charpoly_kglu.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

// removes zero dimension blocks
template<class Field>
void FFLAPACK::updateD(const Field& F, size_t * d, size_t& k,typename Field::Element** minpt){
	size_t ind=0, i=0;
	while(i<k){
		if (d[i]){
			d[ind] = d[i];
			minpt[ind++]= minpt[i];
		}
		i++;
	}
	k = ind;
}

// Compute the new d after a LSP ( d[i] can be zero )
template<class Field>
size_t FFLAPACK::newD( const Field& F, size_t * d, bool& KeepOn, 
		    const size_t l, const size_t N, 
		    const typename Field::Element * X,
		    typename Field::Element ** minpt){
	typedef typename Field::Element elt;
	const elt * Xi = X,*Li, *Xminp;
	KeepOn = false;
	size_t s,m,i,j, ind=0, nr, dtot = 0;
	for ( i=0; dtot<N; ++i){ // for each block
		j = 0;
		nr = s = ( d[i]==l )? 2*l : d[i];
		if (s > N-dtot) 
			s= N-dtot; 
#if DEBUG==2
		cerr<<"computing d["<<i<<"], s="<<s<<endl;
#endif	       
		Li = Xi;
		for ( m=0; m<nr; ++m, Xi+=N ){
#if DEBUG==2
			cerr<<"*Xi="<<*Xi<<" *(Xi+1)="<<*(Xi+1)<<endl;
#endif		
			if (!F.iszero(*Xi) &&  j<s ){
				Xminp = ++Xi;
				++j;
			}
		}
		// Xi points to the next block
		d[i] = j;
		if (j<nr){ // this block has just been completed, store its minpoly
#if DEBUG==2
			cerr<<"block complete, storing minpoly"<<endl;
#endif
			// Xminp-j is the beginning of the minpoly vector
			if (minpt[i]==NULL){ // No polynomial already computed for this block
#if DEBUG==2
				cerr<<"First time: creation of minpt["<<i<<"] of size j="<<j<<endl;
#endif
				// try resize
				minpt[i] = new elt[j];
			}
#if DEBUG==2
			cerr<<"Filling the vector N,i,j,nr="<<N<<" "<<i<<" "<<j<<" "<<nr<<endl;
			cerr<<"Xminp-j+N="<<*(Xminp-j+N)<<" "<<*(Xminp-j+N+1)<<endl;
#endif
			fcopy(F,j,minpt[i],1,Xminp-j+N,1);
			elt * T=new elt[ j*j ];
			TriangleCopy( F, FflasLower, FflasNonUnit, j, T, j, Li,N);
			ftrsv( F, FflasLower, FflasTrans, FflasNonUnit, j, T, j, minpt[i],1);
			delete[] T;

		}
		if (j){
			ind++;
			dtot += j;
#if DEBUG==2
			cerr<<"d["<<i<<"]="<<j<<endl;
#endif
			if ( j==2*l)
				KeepOn = true;
		}
	} // Invariant: sum(d[i],i=0..k-1) = n
	
#if DEBUG==2
	cerr<<"k="<<ind<<endl;
	cerr<<"d=";
	for (j=0;j<i;++j) cerr<<d[j]<<" ";
	cerr<<endl;
#endif
	return i;
}

//---------------------------------------------------------------------
// CharPoly: Compute the characteristic polynomial of A 
//---------------------------------------------------------------------
template <class Field, class Polynomial>
list<Polynomial>&
FFLAPACK::CharPoly( const Field& F, list<Polynomial>& charp, const size_t N,
		 const typename Field::Element * A, const size_t lda,
		 typename Field::Element * U, const size_t ldu){
	
	typedef typename Field::Element elt;
	const elt * Ai;
	elt *Ui, *Uj, *Uk, *Ukp1, *Ukp1new, *Bi, *Vi, *Vk, *Xi, *Xj;
	elt * B = new elt[N*N];     // to store A^2^i
	elt * V = new elt[N*N];     // to store A^2^i.U
	elt * X = new elt[2*N*N];   // to compute the LSP factorization
	size_t * P = new size_t[N]; // Permutation matrix for LSP
	size_t * d= new size_t[N];
	size_t * dv = new size_t[N]; // dimensions of Vect(ei, Aei...)
	size_t * dold = new size_t[N]; // copy of d
	elt ** m = new elt*[N]; // table of each vector storing block minpolys
	typename Polynomial::iterator it;
	size_t i=0, l=1, j, k=N, s, cpt, newRowNb, nrowX, ind, v;
	bool  KeepOn;
#if DEBUG==2
	cerr<<"A="<<endl;
	write_field( F,cerr,A,N,N,lda);
#endif
	for (i=0;i<N;++i)
		m[i]=0;
	Ai = A; 
	Xi = X;
	for ( i=0, v=0; i<N; ++i){
		dv[i] = dold[i] = d[i] = 1;
		v+=2;
	}
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
#if DEBUG==2
	cerr<<"U="<<endl;
	write_field(F,cerr,U,N,N,N);
	cerr<<"V="<<endl;
	write_field(F,cerr,V,N,N,N);
	cerr<<"B="<<endl;
	write_field(F,cerr,B,N,N,N);
	cerr<<"X="<<endl;
	write_field(F,cerr,X,2*N,N,N);
#endif
	// step form elimination using LSP
	for ( i=0;i<N;++i) 
		P[i]=0;
	LUdivine( F, FflasUnit, 2*N, N, X, N, P, FflapackLSP );
#if DEBUG==2
	cerr<<"After LSP X="<<endl;
	write_field(F,cerr,X,2*N,N,N);
#endif

	k = newD( F,d, KeepOn, l, N, X, m);
	while(KeepOn){ // Main loop, until each subspace dimension has been found

#if DEBUG==2
		cerr<<"Before removing extra rows U="<<endl;
		write_field(F,cerr,U,N,N,N);
#endif		
		// Updating U:
		Uk = U;
		// Firstly, removing extra rows
		for ( i = 0, cpt = 0; cpt<N; ++i){
			if (d[i] < dold[i]){
				Ukp1new = Uk +d[i]*ldu;  // new position of Uk+1
				Ukp1 = Uk + dold[i]*ldu; // first row of Uk+1
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
			else {
				Uk += dold[i]*ldu;
			}
			cpt += d[i];
		}

#if DEBUG==2
		cerr<<"After removing extra rows U="<<endl;
		write_field(F,cerr,U,N,N,N);
#endif		
		// Then inserting the duplicated blocks
		Uk = U;
		Vk = V;
		for ( i = 0; i<k; ++i){
			Ukp1 = Uk + dold[i]*ldu; // first row of Uk+1
			newRowNb = d[i] - dold[i];
			if ( newRowNb > 0){
				Ui = U+N*ldu-1; // last row of U
				Uj = U+(N-newRowNb)*ldu-1; // last row future pos
				while ( Uj >= Ukp1){
					for ( j=N;j;--j)
						*(Ui--) = *(Uj--);// moving block
					Ui-=ldu-N;
					Uj-=ldu-N;
				}
				Uj++;
				// copying the d[i]-dold[i] new rows from Vi
				Vi = Vk;
				while ( Uj < Ukp1+ldu*newRowNb ){
					for ( j=N;j;--j)
						*(Uj++) = *(Vi++);
					Uj+=ldu-N;
				}
			}
			Uk = Uk + d[i]*ldu;
			Vk += dv[i]*N;
		}

#if DEBUG==2
		cerr<<"After inserting duplicated rows U="<<endl;
		write_field(F,cerr,U,N,N,N);
#endif
		// block size of X, U, V is l=2^i
		l*=2;
		// B = A^2^i
		fsquare( F, FflasNoTrans, N, F.one, B, N, F.zero, B, N );
		// V = U.B^t
		fgemm( F, FflasNoTrans, FflasNoTrans, N, N, N, F.one, 
		       U, ldu, B, N, F.zero, V, N);		
#if DEBUG==2
		cerr<<"V=U.B^t="<<endl;
		write_field(F,cerr,V,N,N,N);
#endif
		// X = ( U1, V1, U2, V2, ... )
		Xi = X; Ui = U; Vi = V;
		ind=0; cpt=0; nrowX = 0;
		while ( Vi < V + N*N ){
			// Copying Uk
			for ( i = d[ind]; i; --i, nrowX++){ 
				for ( j = N; j; --j )
					*(Xi++) = *(Ui++);
				Ui += ldu - N;
			}
			// Copying Vk
			if ( d[ind] == l || ind==k-1 ){
				cpt+=2*d[ind];
				for ( i=d[ind]; i; i--, nrowX++)
					for ( j=N; j; j--)
						*(Xi++) = *(Vi++);
			}
			else{
				cpt += d[ind];
				Vi = Vi + N*d[ind];
			}
			ind++;
		} 
		updateD( F, d, k, m);
		for (i=0;i<k;++i)
			dv[i] = dold[i] = d[i];
#if DEBUG==2
		cerr<<"U="<<endl;
		write_field(F,cerr,U,N,N,N);
		cerr<<"X="<<endl;
		write_field(F,cerr,X,2*N,N,N);
#endif
		// step form elimination of X using LSP
		for ( i=0;i<N;++i) 
			P[i]=0;
		LUdivine( F, FflasUnit, nrowX, N, X, N, P, FflapackLSP );
		
#if DEBUG==2
		cerr<<"After LSP"<<endl;
		cerr<<"X="<<endl;
		write_field(F,cerr,X,2*N,N,N);
#endif
		k = newD(F, d, KeepOn, l, N, X, m);
	}
	delete[] V;
	delete[] B;
	delete[] P;
	delete[] dv;
	delete[] dold;
	updateD( F, d, k, m);
	
	// Constructing the CharPoly
	charp.clear();
	for ( i=0; i<k; ++i){
		Polynomial * minP = new Polynomial(d[i]+1);
		minP->operator[](d[i]) = F.one;
		it = minP->begin();
		for ( j=0; j<d[i]; ++j, it++)
			F.neg(*it, m[i][j]);
		charp.push_back( *minP );
	}
	delete[] X;
	delete[] d;
	for (i=0; i<k;++i)
		delete[] m[i];
	delete[] m;
}
