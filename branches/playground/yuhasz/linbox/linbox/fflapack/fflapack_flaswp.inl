/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflapack/fflapack_flaswp.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

//---------------------------------------------------------------------
// flaswp: apply the permutation P to the matrix A
// P is stored in an array of indexes
//---------------------------------------------------------------------
template <class Field>
void FFLAPACK::flaswp( const Field& F, const size_t N ,
		       typename Field::Element * A, const size_t lda, 
		       const size_t K1, const size_t K2, size_t *piv, int inci){
	
	//	const size_t n=K2-K1;
	size_t nb = N>>5;
	const size_t mr = N - (nb<<5);
	const size_t incA = lda <<5;
	size_t *ipiv;
	size_t i, ip, i1, i2;
	int KeepOn;
	typename Field::Element *a0;
	typename Field::Element *a1;
	typename Field::Element r;	
	register size_t h;
	
	if (K2 <= K1) return;
	if (inci < 0){
		piv -= (K2 - 1) * inci;
		i1 = K2-1;
		i2 = K1;
	}
	else{
		piv += K1*inci;
		i1 = K1;
		i2 = K2-1;
	}
	
	if (nb){
		do{
			ipiv = piv;
			i = i1;
			do{
				ip = *ipiv; ipiv += inci;
				if ( ip!= i ){
					a0 = A + i;
					a1 = A + ip;
					for (h=32;h;h--){
						r = *a0;
						*a0 = *a1;
						*a1 = r;
						a0 += lda;
						a1 += lda;
					}
				}
				if (inci > 0) KeepOn = (i++ < i2);
				else KeepOn = (i-- > i2);
			}
			while (KeepOn);
			A += incA;
		}
		while(--nb);
	}
	if (mr){
		ipiv = piv;
		i = i1;
		do{
			ip = *ipiv; ipiv += inci;
			if ( ip!= i ){
				a0 = A + i;
				a1 = A + ip;
				for (h=mr;h;h--){
					r = *a0;
					*a0 = *a1;
					*a1 = r;
					a0 += lda;
					a1 += lda;
				}
			}
			if (inci > 0) KeepOn = (i++<i2);
			else{
				KeepOn = (i-- > i2);
			}
		}
		while (KeepOn);
	}
}
