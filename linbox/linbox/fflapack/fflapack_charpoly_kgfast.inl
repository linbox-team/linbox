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
using namespace std;
template <class Field, class Polynomial>
int
LinBox::FFLAPACK::KGFast ( const Field& F, std::list<Polynomial>& charp, 
			   const size_t N,
			   typename Field::Element * A, const size_t lda ){
	
	static typename Field::Element one, zero, mone;
	F.init(one, 1UL);
	F.neg(mone, one);
	F.init(zero, 0UL);
	size_t mc=N>>1; // Matrix A is transformed into a mc_Frobenius form
        size_t mb=N-mc;
	size_t r;
	typename Field::Element * C, *B;

	
	while ( mc > 0 ) {
// 		cerr<<"Boucle1: mc,mb,N="<<mc<<" "<<mb<<" "<<N<<endl;
// 		write_field( F, cerr, A, N, N, lda );
		size_t j=0;
		C = A + (N-mc);
		cerr<<endl<<"mc="<<mc<<":";
		while ( (j+1)*mc < N ) {
			mb = MIN ( mb, N-(j+1)*mc );
// 			cerr<<"Boucle2: j,mb="<<j<<" "<<mb<<endl;
// 			write_field( F, cerr, A, N, N, lda );
			B = A + (N-mc-mb);
		
			// B1 <- C1^-1.B1
			typename Field::Element * LUP = new typename Field::Element[mc*mc];
			for (size_t i=0; i<mc; ++i)
				fcopy( F, mc, LUP+i*mc, 1, C+i*lda, 1);
			size_t * P = new size_t[mc];
			size_t * Q = new size_t[mc];

			if ( r = LUdivine( F, FflasNonUnit, mc, mc, 
				      LUP, mc, P, FflapackLQUP, Q ) < mc )
				return -1;
// 			cerr<<"LUP="<<endl;
// 			write_field( F, cerr, LUP, mc, mc, mc );
                        cerr<<" "<<r;
			ftrsm(F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit,
			      mc, mb, one, LUP, mc , B, lda);
			ftrsm(F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit, 
			      mc, mb, one, LUP, mc , B, lda);
			delete[] LUP;
			applyP( F, FflasLeft, FflasTrans, mb, 0, mc, B, lda, P );

// 			cerr<<"Apres B1<-C1^-1"<<endl;
// 			write_field( F, cerr, A, N, N, lda );
                        
			// B2 <- B2 - C2.B1
			fgemm(F, FflasNoTrans, FflasNoTrans, N-mc, mb, mc, 
			      mone, C+mc*lda, lda, B, lda, 
			      one, B+mc*lda, lda);

// 			cerr<<"Apres B2<-B2-C2.B1"<<endl;
//                         write_field( F, cerr, A, N, N, lda );

			// Shifting B: B1;B2 -> B2;B1
			typename Field::Element * tmp = new typename Field::Element[mc*mb];
			for (size_t i=0; i<mc; ++i)
				fcopy( F, mb, tmp+i*mb, 1, B+i*lda, 1);
			for (size_t i=mc; i<N; ++i)		
				fcopy( F, mb, B+(i-mc)*lda, 1, B+i*lda, 1);
			for (size_t i=0; i<mc; ++i)
				fcopy( F, mb, B+(i+N-mc)*lda, 1, tmp+i*mb, 1);
			delete[] tmp;
// 			cerr<<"Apres shift de B"<<endl;
//                         write_field( F, cerr, A, N, N, lda );

			
			// C3 <- B3.C1 + C3
			fgemm(F, FflasNoTrans, FflasNoTrans, (j+1)*mc, mc, mb, 
			      one, B+(N-(j+1)*mc)*lda, lda, C+(N-(j+1)*mc-mb)*lda, lda,
			      one, C+(N-(j+1)*mc)*lda, lda);
			// cerr<<"C3 <- B3.C1 + C3: B3="<<endl;
//                         write_field( F, cerr, B+(N-(j+1)*mc)*lda, (j+1)*mc, mb, lda );
// 			cerr<<"C3 <- B3.C1 + C3: C1"<<endl;
//                         write_field( F, cerr,  C+(N-(j+1)*mc-mb)*lda, mb, mc, lda );
// 			cerr<<"C3 <- B3.C1 + C3: C3="<<endl;
//                         write_field( F, cerr, C+(N-(j+1)*mc)*lda, (j+1)*mc, mc, lda );

			int lambda = N - mb - (j+1)*mc;
			if ( int(mb) < lambda ){

// 				cerr<<"mb<lambda"<<endl;
				typename Field::Element * tmp = new typename Field::Element[lambda*mc];

				// tmp <- C1
				for (int i=0; i<lambda; ++i)
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
				for (int i=0; i<lambda; ++i)
					fcopy( F, mc, C+mb*lda+i*lda, 1, tmp+i*mc, 1);	
				delete[] tmp;
			}
			else if ( lambda > 0 ){
// 				cerr<<"lambda>0"<<endl;
                                
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
				for (int i=0; i<lambda; ++i)
					fcopy( F, mc, C+mb*lda+i*lda, 1, C+i*lda, 1);

				// C1' <- tmp
				for (size_t i=0; i<mb; ++i)
					fcopy( F, mc, C+i*lda, 1, tmp+i*mc, 1);
				delete[] tmp;
			}
			else{
// 				cerr<<"lambda<0"<<endl;
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
			
			j++;
		}
		mb = mc;
	        mc>>=1;
		mb -= mc;
	}

	Polynomial *minP = new Polynomial();
	minP->resize(N+1);
	minP->operator[](N) = one;
	typename Polynomial::iterator it = minP->begin();
	for (size_t j=0; j<N; ++j, it++){
		F.neg(*it, *(A+N-1+j*lda));
	}
	charp.clear();
	charp.push_front(*minP);
	return 0;
}
