
/* ffpack/ffpack_charpoly_kglu.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

template<class Field>
size_t FFPACK::updateD(const Field& F, size_t * d, size_t k,
			       std::vector<std::vector<typename Field::Element> >& minpt){
	size_t ind=0, i=0;
	while(i<k){
		if (d[i]){
			if (ind<i){
				d[ind] = d[i];
				minpt[ind++] = minpt[i];
			}
			else ind++;
		}
		i++;
	}
	for (i=ind; i<k; ++i)
		minpt[i].resize(0);
	minpt.resize(ind);
	return ind;
}

// Compute the new d after a LSP ( d[i] can be zero )
template<class Field>
size_t FFPACK::newD( const Field& F, size_t * d, bool& KeepOn, 
			       const size_t l, const size_t N, 
			       typename Field::Element * X,
			       const size_t * Q,
			       std::vector<std::vector<typename Field::Element> >& minpt){
	typedef typename Field::Element elt;
	//const elt * Xi = X; // Xi points to the begining of each block
	elt *Li=X, *Xminp=X;
	KeepOn = false;
	size_t nr, s, i, j, jtot=0, dtot = 0, nrtot=0;
	
	for ( i=0; dtot<N; ++i){ // for each block
		j = 0;
		nr = s = ( d[i]==l )? 2*l : d[i];
		if (s > N-dtot) 
			s= N-dtot; 
		nrtot += nr;
		
		while ( (Q[j+jtot] <nrtot) && (j+jtot<N) )
			j++;
		
		Xminp = X+Q[j+jtot-1]*N+jtot+j ;
		d[i] = j;
		jtot+=j;
		dtot += j;

		if (j<nr){ 
			minpt[i].resize(j);
			ftrsv( F, FflasLower, FflasTrans, FflasUnit, 
			       j, Li, N, Xminp-j+N,1);
			elt* Xi = Xminp-j+N;
			for (size_t ii = 0; ii<j; ++ii, ++Xi)
				minpt[i][ii] = *Xi;
		}
		Li += nr*N+j;
		if ( j==2*l )
			KeepOn = true;
	} // Invariant: sum(d[i],i=0..k-1) = n
	
	return i;
}

//---------------------------------------------------------------------
// CharPoly: Compute the characteristic polynomial of A using 
// Keller-Gehrig's algorithm
//---------------------------------------------------------------------
template <class Field, class Polynomial>
std::list<Polynomial>&
FFPACK::KellerGehrig( const Field& F, std::list<Polynomial>& charp, const size_t N,
				const typename Field::Element * A, const size_t lda ){
	
	
	typedef typename Field::Element elt;
	const elt * Ai=A;
	static elt one, zero;
	elt * U = new elt[N*N];     // to store A^2^i
	elt * B = new elt[N*N];     // to store A^2^i
	elt * V = new elt[N*N];     // to store A^2^i.U
	elt * X = new elt[2*N*N];   // to compute the LSP factorization
	elt *Ui, *Uj, *Uk, *Ukp1, *Ukp1new, *Bi, *Vi, *Vk, *Xi=X, *Xj;
	size_t * P = new size_t[N]; // Column Permutation for LQUP
	size_t * Q = new size_t[2*N]; // Row Permutation for LQUP

	size_t * d= new size_t[N];   // dimensions of Vect(ei, Aei...)
	size_t * dv = new size_t[N];
	size_t * dold = new size_t[N]; // copy of d
	// vector of the opposite of the coefficient of computed minpolys
	std::vector< std::vector< elt > > m(N); 
	typename Polynomial::iterator it;
	size_t i=0, l=1, j, k=N,  cpt, newRowNb, nrowX, ind;
	bool  KeepOn;
	F.init( one, 1.0);
	F.init( zero, 0.0);

	for ( i=0; i<N; ++i)
		dv[i] = dold[i] = d[i] = 1;
	
	// Computing the first X: (e1; e1A^t; e2; e2A^t;...;en;enA^t)
	for ( i=0, Ui=U, Vi=V, Bi=B; i<N; ++i, Ai -= N*lda-1  ){
		for ( Xj=Xi, Uj=Ui; Xj<Xi+N; ++Xj, ++Uj){
			F.assign(*Xj, zero);
			F.assign(*Ui, zero);
		}
		F.assign(*(Ui+i), one);
		F.assign(*(Xi+i), one);
		while ( Xj<Xi+2*N) {
			*(Bi++) = *(Xj++) = *(Vi++) = *Ai;
			Ai+=lda;
		}
		Xi = Xj;
		Ui = Uj;
	}
	
	// step form elimination using LQUP factorization
	for ( i=0;i<N;++i) 
		P[i]=0;
	for ( i=0;i<2*N;++i) 
		Q[i]=0;
	LUdivine( F, FflasNonUnit, FflasNoTrans, 2*N, N, X, N, P, Q, FfpackLQUP);
	
	k = newD( F,d, KeepOn, l, N, X, Q, m);
	
	while(KeepOn){ // Main loop, until each subspace dimension has been found
		// Updating U:
		Uk = U;
		// Firstly, removing extra rows
		for ( i = 0, cpt = 0; i<N; ++i){
			if (d[i] < dold[i]){
				Ukp1new = Uk + d[i]*N;  // new position of Uk+1
				Ukp1 = Uk + dold[i]*N; // first row of Uk+1
				Ui = Ukp1new;
				Uj = Ukp1;
				while ( Uj < U + N*N ){
					for ( j=N; j; --j)
						*(Ui++) = *(Uj++);
				}
				Uk = Ukp1new;
				dold[i] = d[i];
			}
			else {
				Uk += dold[i]*N;
			}
			cpt += d[i];
		}

		// Then inserting the duplicated blocks
		Uk = U;
		Vk = V;
		for ( i = 0; i<k; ++i){
			Ukp1 = Uk + dold[i]*N; // first row of Uk+1
			newRowNb = d[i] - dold[i];
			if ( newRowNb > 0){
				Ui = U+N*N-1; // last row of U
				Uj = U+(N-newRowNb)*N-1; // last row future pos
				while ( Uj > Ukp1-1){
					for ( j=N;j;--j)
						*(Ui--) = *(Uj--);// moving block
				}
				Uj++;
				Vi = Vk;
				while ( Uj < Ukp1+N*newRowNb ){
					for ( j=N;j;--j)
						*(Uj++) = *(Vi++);
				}
			}
			Uk = Uk + d[i]*N;
			Vk += dv[i]*N;
		}

		// max block size of X, U, V is l=2^i
		l*=2;
		// B = A^2^i
		fsquare( F, FflasNoTrans, N, one, B, N, zero, B, N );
		// V = U.B^t
		fgemm( F, FflasNoTrans, FflasNoTrans, N, N, N, one, 
		       U, N, B, N, zero, V, N);		
		// X = ( U1, V1, U2, V2, ... )
		Xi = X; Ui = U; Vi = V;
		ind=0; cpt=0; nrowX = 0;
		while ( Vi < V + N*N ){
			// Copying Uk
			for ( i = d[ind]; i; --i, nrowX++){ 
				for ( j = N; j; --j )
					*(Xi++) = *(Ui++);
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
		// removes factors of degree 0 in m
		k = updateD( F, d, k, m);
		
		for (i=0;i<k;++i)
			dv[i] = dold[i] = d[i];
		
		// step form elimination of X using LSP
		for ( i=0;i<N;++i) 
			P[i]=0;
		for ( i=0;i<2*N;++i) 
			Q[i]=0;
		LUdivine( F, FflasNonUnit, FflasNoTrans, nrowX, N, X, N, P, Q, FfpackLQUP);
		
		// Recompute the degrees of the list factors
		k = newD(F, d, KeepOn, l, N, X,Q, m);
	}
	delete[] U;
	delete[] V;
	delete[] B;
	delete[] P;
	delete[] Q;
	delete[] dv;
	delete[] dold;
	
	k = updateD( F, d, k, m);
	// Constructing the CharPoly
	for ( i=0; i<k; ++i){
		Polynomial * minP = new Polynomial(d[i]+1);
		minP->operator[](d[i]) = one;
		it = minP->begin();
		for ( j=0; j<d[i]; ++j, it++)
			F.neg(*it, m[i][j]);
		charp.push_back( *minP );
	}
	delete[] X;
	delete[] d;
	return charp;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
