
/* ffpack/ffpack_charpoly_danilevski.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif

//---------------------------------------------------------------------
// CharPoly: Compute the characteristic polynomial of A using 
// Danilevski's algorithm.
//---------------------------------------------------------------------

template <class Field, class Polynomial>
std::list<Polynomial>&
FFPACK::Danilevski (const Field& F, std::list<Polynomial>& charp, 
			    const size_t N, typename Field::Element * A, 
			    const size_t lda){	
	static typename Field::Element one, mone;
	F.init(one, 1.0);
	F.neg(mone, one);
	charp.clear();
	size_t dtot=0;
	typename Field::Element *pivot,*e,*u1,invp;
	for (size_t k=0; k<N; ++k){
		size_t i = k+1; 
		size_t d;
		e = pivot = A + (k+1) * lda + k; // coef
		while ((i<N) && F.isZero(*e)) { e += lda; i++; }
		if (i < N){
			if (i > k + 1) {
				fswap (F, N-k, e, 1, pivot, 1);
				fswap (F, N, A+i, lda, A+k+1, lda);
			}
			F.inv (invp, *pivot);
			fscal (F, N-k-1, invp, pivot+1, 1);
			fscal (F, N-dtot, *pivot, A+dtot*lda+k+1, lda);
			// X <- X - uw
			fger (F, k + 1-dtot, N - k -1, mone, 
			      A + dtot*lda + k, lda, pivot+1, 1, 
			      A+k+1+dtot*lda, lda);
			if (k<N-2){
				
				// Y <- Y - vw
				fger (F, N-k-2, N-k-1, mone, pivot+lda, lda, pivot+1, 1, 
				      pivot+lda+1,lda);
				//6
				fgemv (F, FflasNoTrans, N-dtot, N-k-2, 
				       one, A+dtot*lda+k+2, lda, pivot+lda, lda, one, 
				       A+dtot*lda+k+1,lda);
			}
			//5
			u1 = A+dtot*lda+k;
			for (i = dtot; i <= k; ++i){
				F.addin( *(u1+lda+1), *u1);
				u1+=lda;
			}
		}
		if (i==N){// completed one companion block
			d = k+1-dtot;
			typename Field::Element *Ai = A+k+dtot*lda;
			Polynomial * P = new Polynomial(d+1);
			for (i = 0; i < d; ++i){
				F.neg (P->operator[](i), *(Ai+i*lda));
			}
			F.assign( (*P)[d], one);
			charp.push_front(*P);
			dtot+=d;
		} 
	}
	return charp;
}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
