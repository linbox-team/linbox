/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflas/fflas_fgemv.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */


#include "linbox/field/givaro-zpz.h"
#include "linbox/field/modular.h"

//---------------------------------------------------------------------
// fgemv: GEneral Matrix Vector Multiplication
// Computes  Y <- alpha.op(A).X + beta.Y
// A is M*N, 
//---------------------------------------------------------------------
template<class Field>
inline void 
LinBox::FFLAS::fgemv( const Field& F, const enum FFLAS_TRANSPOSE TransA, 
	      const size_t M, const size_t N,
	      const typename Field::Element alpha, 
	      const typename Field::Element * A, const size_t lda,
	      const typename Field::Element * X, const size_t incX,
	      const typename Field::Element beta,
	      typename Field::Element * Y, const size_t incY){

// 	static typename Field::Element Mone, one;
// 	 F.init(Mone,-1);
// 	F.init(one,1);
// 	typename Field::Element tmp;
	

	double* Ad = new double[M*N];
	double* Xd = new double[N];
	double* Yd = new double[M];	
	double alphad, betad;
	F.convert( alphad, alpha );
	F.convert( betad, beta );

	MatF2MatD( F, Ad, A, lda, M, N );
	
	double *Xdi=Xd;	
	for (const typename Field::Element* Xi=X; Xi != X+N*incX; Xi+=incX, Xdi++)
		F.convert( *(Xdi), *Xi);
	

	double  *Ydi=Yd;
	if (!F.isZero(beta))
		for (typename Field::Element* Yi = Y; Yi != Y+M*incY; Yi+=incY, Ydi++)
			F.convert( *(Ydi), *Yi);

	cblas_dgemv( CblasRowMajor, (enum CBLAS_TRANSPOSE) TransA, M, N, alphad, Ad, N, Xd, 1, 
		     betad, Yd, 1);

	Ydi=Yd;
	for ( typename Field::Element* Yi = Y; Yi != Y+M*incY; Yi+=incY, Ydi++)
		F.init( *Yi, *(Ydi));
	

	delete[] Ad;
	delete[] Xd;
	if (!F.isZero(beta))
		delete[] Yd;

// 	for (size_t i=0; i<M; ++i){
// 		//tmp = fdot( F, i, A+i*lda, 1, X, incX);
// 		const size_t NB_blocks; i/DIM_block;
		
// 		for(unsigned int k = 0; k < NB_blocks; ++k) {
// 			for(unsigned int j = 0; j<DIM_block; ++j)
// 				    res += v1[j] * v2[j];
// 			    Element dbl = res;
// 			    dbl *= _F.inv_modulus;

// 			    dbl = (double)( (long long)dbl );
// 			    dbl *= _F.modulus;
// 			    res -= dbl;
// 		    }
// 		    for(unsigned int i = NB_blocks*DIM_block; i < DIM; ++i)
// 			    res += v1[i] * v2[i];
// 		    Element dbl(res);
// 		    dbl *= _F.inv_modulus;
// 		    dbl = (double)( (long long)dbl );
// 		    dbl *= _F.modulus;
// 		    return res -= dbl;
// 	    }
// 		if (!F.isOne(beta))
// 			F.mulin(*(Y+i*lda), beta );
// 		if ( F.areEqual(alpha, Mone))
// 			F.subin( *(Y+i*lda), tmp);
// 		else if (F.areEqual(alpha, one))
// 			F.addin( *(Y+i*lda), tmp);
// 		else
// 			F.axpyin(alpha, tmp, *(Y+i*lda));
// 		//A finir... optimiser pour beta=0
// 	}
}

template<>
inline void
LinBox::FFLAS::fgemv( const Modular<double>& F, const enum FFLAS_TRANSPOSE TransA, 
	      const size_t M, const size_t N,
	      const double alpha, 
	      const double * A, const size_t lda,
	      const double * X, const size_t incX,
	      const double beta,
	      double * Y, const size_t incY){

	cblas_dgemv( CblasRowMajor, (enum CBLAS_TRANSPOSE) TransA, M, N, alpha, A, lda, X, incX, 
		     beta, Y, incY);

	for ( double * Yi = Y; Yi != Y+M*incY; Yi+=incY)
		F.init( *Yi, *Yi);
}

#if __LINBOX_HAVE_GIVARO
template<>
inline void
LinBox::FFLAS::fgemv( const GivaroZpz<Std32>& F, const enum FFLAS_TRANSPOSE TransA,
	      const size_t M, const size_t N,
	      const GivaroZpz<Std32>::Element alpha,
	      const GivaroZpz<Std32>::Element * A, const size_t lda,
	      const GivaroZpz<Std32>::Element * X, const size_t incX,
	      const  GivaroZpz<Std32>::Element beta,
	      GivaroZpz<Std32>::Element * Y, const size_t incY) {
           // Suppose alpha == -1
           // beta == 1



       static const long MOD = F.characteristic();
       static const long MODmun = MOD-1;
       static unsigned long   temp = ( (unsigned long)(-MODmun) / MODmun) / MODmun;
       static const unsigned long Kbest = (temp<1?1:temp);
       static const unsigned long CORR = (( (2147483648UL) % (unsigned long)MOD)<<1) % (unsigned long)MOD;
       const unsigned long DIMoKbest = (N / Kbest);
       const unsigned long DIMoKfKbest = DIMoKbest*Kbest;

       GivaroZpz<Std32>::Element * Yit = Y;
       const GivaroZpz<Std32>::Element * Ait = A;//, * Xit = X;
       GivaroZpz<Std32>::Element best, inter;

       if (DIMoKfKbest == 0)
           for(unsigned long i=0; i<M; ++i) {
               best = 0;
               for(unsigned long j=0; j<N; ++j)
		       //    best += Ait[j] * X[j];
		       //	       if incX != 1, Change to : 
		       best += Ait[j] * ( *(X+j*incX));

               best %= MOD;
       
                   // Suppose alpha == -1, beta == 1 !!!
	       if (*Yit<best)  *Yit+=MOD;
	       *Yit-= best;

               Yit += incY;
               Ait += lda;
           }
       else
           for(unsigned long i=0; i<M; ++i) {
               best = inter = 0;
               unsigned long k = 0;
               for(unsigned long  j=0 ; j < Kbest; ++j, ++k)
		       //                   best += Ait[k] * X[k];
// 		      if incX != 1, Change to : 
		       best += Ait[k] * ( *(X+k*incX));
               for(inter = best ; k < DIMoKfKbest; inter=best) {
                   for(unsigned long  j=0 ; j < Kbest; ++j, ++k)
			   //                       best += Ait[k] * X[k];
// 		      if incX != 1, Change to : 
			   best += Ait[k] * ( *(X+k*incX));

                   if (inter > best) best += CORR;
               
               }
               inter = best;
               for(; k<N; ++k)
		       //                   best += Ait[k] * X[k];
// 		      if incX != 1, Change to : 
		       best += Ait[k] * ( *(X+k*incX));
               if (inter > best) best += CORR;
               best %= MOD;
               
                    // Suppose alpha == -1, beta == 1 !!!
	       if (*Yit<best)  *Yit+=MOD;
               *Yit -= best;

               Yit += incY;
               Ait += lda;
           }
   }


template<>
inline void
LinBox::FFLAS::fgemv( const GivaroZpz<Std64>& F, const enum FFLAS_TRANSPOSE TransA,
	      const size_t M, const size_t N,
	      const GivaroZpz<Std64>::Element alpha,
	      const GivaroZpz<Std64>::Element * A, const size_t lda,
	      const GivaroZpz<Std64>::Element * X, const size_t incX,
	      const  GivaroZpz<Std64>::Element beta,
	      GivaroZpz<Std64>::Element * Y, const size_t incY) {

	
	static const long MOD=F.size();
	GivaroZpz<Std64>::Element * Yit = Y;
	const GivaroZpz<Std64>::Element * Ait = A;//, * Xit = X;
	GivaroZpz<Std64>::Element acc;
	for(unsigned long i=0; i<M; ++i) {
               acc = 0;
               for(unsigned long j=0; j<N; ++j)
		       acc += Ait[j] * ( *(X+j*incX));
               acc %= MOD;
               
                   // Suppose alpha == -1, beta == 1 !!!
               *Yit -= acc;

               Yit += incY;
               Ait += lda;

	}
}
#endif
