/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflapack/fflapack_charpoly_kgfast.inl
 * Copyright (C) 2004 Clement Pernet
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
// Keller-Gehrig's fast algorithm. A must be generic.
//---------------------------------------------------------------------
template <class Field, class Polynomial>
int
LinBox::FFLAPACK::KGFast ( const Field& F, std::list<Polynomial>& charp, const size_t N,
			   const typename Field::Element * A, const size_t lda ){
	
	size_t mc=N>>1; // Matrix A is transformed into a mc_Frobenius form
	size_t mb=N-mc;
	typename Field::Element * C, *B;

	while ( m > 1 ) do{
		
		size_t j=0;
		C = A + (N-mc);
		
		while ( (j+1)*mc < N ) do{
			
			mb = MIN ( mb, N-(j+1)*mc );
			B = A + (N-mc-mb);
		
			// B1 <- C1^-1.B1
			typename Field::Element * LUP = new typename Field::Element[mc*mc];
			for (size_t i=0; i<mc; ++i)
				fcopy( F, mc, LUP+i*mc, 1, C+i*lda, 1);
			if (LUdivine( F, FflasNonUnit, mc, mc, LUP, mc, P, FflapackLQUP, Q ) < mc)
				return -1;
			ftrsm(F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit, mc, mb, one, 
			      A, lda , X, ldx);
			ftrsm(F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit, mc, mb, one, 
			      A, lda , X, ldx);
			delete[] LUP;
			applyP( F, FflasLeft, FflasTrans, mb, 0, mc, B, lda, P );

			// B2 <- B2 - C2.B1
			fgemm(F, FflasNoTrans, FflasNoTrans, N-mc, mb, mc, 
			      Mone, C+mc*lda, lda, B, lda, 
			      one, B+mc*lda, lda);

			// Shifting B: B1;B2 -> B2;B1
			typename Field::Element * tmp = new typename Field::Element[mc*mb];
			for (size_t i=0; i<mc; ++i)
				fcopy( F, mb, tmp+i*mb, 1, B+i*lda, 1);
			for (size_t i=mc; i<N; ++i)		
				fcopy( F, mb, B+(i-mc)*lda, 1, B+i*lda, 1);
			for (size_t i=0; i<mc; ++i)
				fcopy( F, mb, B+(i+N-mc)*lda, 1, tmp+i*mb, 1);
			delete[] tmp;

			int lambda = N - mb - (j+1)*mc;
			if ( mb < lambda ){

				typename Field::Element * tmp = new typename Field::Element[lambda*mc];

				// tmp <- C1
				for (size_t i=0; i<lambda; ++i)
					fcopy( F, mc, tmp+i*mc, 1, C+i*lda, 1);	
				
				// C1' <- B1.C2
				fgemm(F, FflasNoTrans, FflasNoTrans, mb, mc, mb, 
				      one, B, lda, C+lambda*lda, lda,
				      zero, C, lda);

				// tmp <- B2.C2 + tmp
				fgemm(F, FflasNoTrans, FflasNoTrans, lambda, mc, mb, 
				      one, B+mb*lda, lda, C+lambda*lda, lda,
				      one, tmp, mc);
				
				// C2' <- tmp
				for (size_t i=0; i<lambda; ++i)
					fcopy( F, mc, C+mb*lda+i*lda, 1, tmp+i*mc, 1);	
				delete[] tmp;
			}
			else if ( lambda > 0 ){
				typename Field::Element * tmp = new typename Field::Element[mb*mc];
				// C1 <- B2.C2 + C1
				fgemm(F, FflasNoTrans, FflasNoTrans, lambda, mc, mb, 
				      one, B+mb*lda, lda, C+lambda*lda, lda,
				      one, C, lda);

				// tmp <-B1.C2
				fgemm(F, FflasNoTrans, FflasNoTrans, mb, mc, mb, 
				      one, B, lda, C+lambda*lda, lda,
				      zero, tmp, mc);
				
				// C2' <- C1
				for (size_t i=0; i<lambda; ++i)
					fcopy( F, mc, C+mb*lda+i*lda, 1, C+i*lda, 1);

				// C1' <- tmp
				for (size_t i=0; i<mb; ++i)
					fcopy( F, mc, C+i*lda, 1, tmp+i*mc, 1);
				delete[] tmp;
			}
			else{
				mb = N - (j+1)*mc;
				typename Field::Element * tmp = new typename Field::Element[mb*mc];

				// tmp <-B1.C1
				fgemm(F, FflasNoTrans, FflasNoTrans, mb, mc, mb, 
				      one, B, lda, C, lda,
				      zero, tmp, mc);

				// C1' <- tmp
				for (size_t i=0; i<mb; ++i)
					fcopy( F, mc, C+i*lda, 1, tmp+i*mc, 1);
				delete[] tmp;
			}
			
			// C3 <- B3.C2 + C3
			fgemm(F, FflasNoTrans, FflasNoTrans, (j+1)*mc, mc, mb, 
			      Mone, B+(N-(j+1)*mc)*lda, lda, C+(N-(j+1)*mc-mb)*lda, lda,
			      one, C+(N-(j+1)*mc)*lda, lda);
			j++;
		}
		mc/=2;
	}
	Polynomial *minP = new Polynomial();
	minP.resize(N+1);
	minP[N] = one;
	typename Polynomial::iterator it = minP.begin();
	for (int j=0; j<N; ++j, it++){
		F.neg(*it, A+N-1+j*lda);
	}
	charp.clear();
	charp.push_front(minP);
	return 0;
}
