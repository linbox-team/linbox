/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflas/fflas_fgemv.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */


#if __LINBOX_HAVE_GIVARO
// Note: all use of givaro is commented out below, so I also comment out here. -bds
//#include "linbox/field/givaro-zpz.h"
#endif

#include "linbox/field/modular-double.h"

namespace LinBox {

//---------------------------------------------------------------------
// fgemv: GEneral Matrix Vector Multiplication
// Computes  Y <- alpha.op(A).X + beta.Y
// A is M*N, 
//---------------------------------------------------------------------
template<class Field>
inline void 
FFLAS::fgemv (const Field& F, const enum FFLAS_TRANSPOSE TransA, 
	      const size_t M, const size_t N,
	      const typename Field::Element alpha, 
	      const typename Field::Element * A, const size_t lda,
	      const typename Field::Element * X, const size_t incX,
	      const typename Field::Element beta,
	      typename Field::Element * Y, const size_t incY)
{
 	static typename Field::Element  one, mone, zero;
	F.init(one,1);
 	F.init(mone,-1);
 	F.init(zero,0);

	static Field G = F;
	static integer pig;
	G.characteristic(pig);
	integer pif;
	F.characteristic(pif);
	static typename Field::Element b = beta;
	static size_t kmax = FflasKmax (F, 0, beta)-1;
     
 
	if ((b != beta) || (pif != pig)){
		G = F;
		b = beta;
		kmax = FflasKmax (F, 0, beta)-1;
	}
	//	FflasKmax (kmax, F, 0, beta);
	//kmax--;
	if (kmax > 1) {
		if  (TransA == FflasNoTrans) {
			size_t nblock = N / kmax;
			size_t remblock = N % kmax;
			// To ensure the initial computation with beta
			if (!remblock){
				remblock = kmax;
				--nblock;
			}
			
			MatVectProd (F, FflasNoTrans, M, remblock, alpha, 
				     A+kmax*nblock, lda, X+kmax*nblock*incX, incX, beta, 
				     Y, incY);
			for  (size_t i = 0; i < nblock; ++i){
				MatVectProd (F, FflasNoTrans, M, kmax, alpha, 
					     A+i*kmax, lda, X+i*kmax*incX, incX, one, 
					     Y, incY);
			}
		}
		else{ // FflasTrans
			size_t nblock = M / kmax;
			size_t remblock = M % kmax;
			// To ensure the initial computation with beta
			if (!remblock){
				remblock = kmax;
				--nblock;
			}
			
			MatVectProd (F, FflasTrans, remblock, N, alpha, 
				     A+kmax*nblock*lda, lda, X+kmax*nblock*incX, incX, beta, 
				     Y, incY);
			for  (size_t i = 0; i < nblock; ++i){
				MatVectProd (F, FflasTrans, kmax, N, alpha, 
					     A+i*kmax*lda, lda, X+i*kmax*incX, incX, one, 
					     Y, incY);
			}
		
		}
	} else {
		if  (TransA == FflasNoTrans) {
			if (F.isZero (beta))
				for (size_t i = 0; i < M; ++i)
					F.assign( *(Y+i*incY), zero);
			else {
				typename Field::Element betadivalpha;
				F.div (betadivalpha, beta, alpha); 
				for (size_t i = 0; i < M; ++i)
					F.mulin( *(Y+i*incY), betadivalpha);
			}
			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j)
					F.axpyin (*(Y+i*incY), *(A+i*lda+j), *(X+j*incX));
			if (! F.isOne(alpha))
				for (size_t i = 0; i < M; ++i)
					F.mulin (*(Y+i*incY), alpha);
		} else {
			if (F.isZero (beta))
				for (size_t i = 0; i < N; ++i)
					F.assign( *(Y+i*incY), zero);
			else {
				typename Field::Element betadivalpha;
				F.div (betadivalpha, beta, alpha); 
				for (size_t i = 0; i < N; ++i)
					F.mulin( *(Y+i*incY), betadivalpha);
			}
			

			for (size_t i = 0; i < M; ++i)
				for (size_t j = 0; j < N; ++j){
					F.axpyin (*(Y+j*incY), *(A+i*lda+j), *(X+i*incX));
				}
			if (! F.isOne(alpha))
				for (size_t i = 0; i < N; ++i)
					F.mulin (*(Y+i*incY), alpha);
		}
	}
}
	
// MatVectProd: computes y <- alpha.op(A)*x +  beta.y. 
// Assumes that the condition k(p-1)^2 <2^53 is satisfied
template<class Field>
inline void 
FFLAS::MatVectProd (const Field& F, const enum FFLAS_TRANSPOSE TransA, 
			    const size_t M, const size_t N,
			    const typename Field::Element alpha, 
			    const typename Field::Element * A, const size_t lda,
			    const typename Field::Element * X, const size_t incX,
			    const typename Field::Element beta,
			    typename Field::Element * Y, const size_t incY){

 	static typename Field::Element  one, mone;
 	typename Field::Element tmp;
	F.init(one,1);
	F.init(mone,-1);
	
	size_t Xl, Yl;
	if  (TransA == FflasNoTrans){
		Xl = N;
		Yl = M;
	}
	else{
		Xl = M;
		Yl = N;
	}
	double* Ad = new double[M*N];
	double* Xd = new double[Xl];
	double* Yd = new double[Yl];
	double alphad, betad;

	if (F.areEqual (mone, alpha)){
		alphad = -1.0;
		F.convert (betad, beta);
	}
	else{
		if (! F.areEqual (one, alpha)){
			// Compute C = A*B + beta/alpha.C
			// and after C *= alpha
			F.div (tmp, beta, alpha);
			F.convert (betad, tmp);
		}
		else
			F.convert (betad, beta);
		alphad = 1.0;
	}
	
	MatF2MatD (F, Ad, N, A, lda, M, N);
	
	double *Xdi=Xd;	
	for (const typename Field::Element* Xi=X; Xi != X+Xl*incX; Xi+=incX, Xdi++)
		F.convert (*(Xdi), *Xi);
	double  *Ydi=Yd;
	if (!F.isZero(beta))
		for (typename Field::Element* Yi = Y; Yi != Y+Yl*incY; Yi+=incY, Ydi++)
			F.convert (*(Ydi), *Yi);

	cblas_dgemv (CblasRowMajor, (enum CBLAS_TRANSPOSE) TransA, M, N, alphad, Ad, N, Xd, 1, 
		     betad, Yd, 1);

	Ydi=Yd;
	for  (typename Field::Element* Yi = Y; Yi != Y+Yl*incY; Yi+=incY, Ydi++)
		F.init (*Yi, *(Ydi));
	
	if  (!F.areEqual (one, alpha) && !F.areEqual (mone, alpha)){
		// Fix-up: compute Y *= alpha
		for (typename Field::Element* Yi = Y; Yi != Y+Yl*incY; Yi += incY)
			F.mulin (*Yi , alpha);
	}
	
	delete[] Ad;
	delete[] Xd;
	delete[] Yd;
}

template<>
inline void
FFLAS::MatVectProd (const Modular<double>& F, const enum FFLAS_TRANSPOSE TransA, 
			    const size_t M, const size_t N,
			    const double alpha, 
			    const double * A, const size_t lda,
			    const double * X, const size_t incX,
			    const double beta,
			    double * Y, const size_t incY){

	static  double Mone, one;
	F.init(Mone, -1);
	F.init(one, 1);	
	double _alpha, _beta;
	
	if  (F.areEqual (Mone, beta))
		_beta = -1.0;
	else
		_beta = beta;

	if (F.areEqual (Mone, alpha)){
		_alpha = -1.0;
	}
	else{
		_alpha = 1.0;
		if (! F.areEqual (one, alpha)){
			// Compute y = A*x + beta/alpha.y
			// and after y *= alpha
			F.divin (_beta, alpha);
		}
	}
	
	cblas_dgemv (CblasRowMajor, (enum CBLAS_TRANSPOSE) TransA, M, N, 
		     _alpha, A, lda, X, incX, _beta, Y, incY);
	
	for  (double * Yi = Y; Yi != Y+((TransA == FflasNoTrans)?M:N)*incY; Yi+=incY)
		F.init (*Yi, *Yi);
	
	if ( (!F.areEqual (one, alpha)) && (!F.areEqual (Mone, alpha))){
		// Fix-up: compute y *= alpha
		for (double* Yi = Y; Yi != Y+((TransA == FflasNoTrans)?M:N)*incY; Yi += incY)
			F.mulin (*Yi , alpha);
	}
}


template<>
inline void
FFLAS::fgemv (const DoubleDomain& D, const enum FFLAS_TRANSPOSE TransA, 
		      const size_t M, const size_t N,
		      const DoubleDomain::Element  alpha, 
		      const DoubleDomain::Element * A, const size_t lda,
		      const DoubleDomain::Element * X, const size_t incX,
		      const DoubleDomain::Element beta,
		      DoubleDomain::Element * Y, const size_t incY)
{
	cblas_dgemv (CblasRowMajor, (enum CBLAS_TRANSPOSE) TransA, M, N, 
		     alpha, A, lda, X, incX, beta, Y, incY);
}


#if __LINBOX_HAVE_GIVARO
// template<>
// inline void
// FFLAS::fgemv (const GivaroZpz<Std32>& F, const enum FFLAS_TRANSPOSE TransA,
// 	      const size_t M, const size_t N,
// 	      const GivaroZpz<Std32>::Element alpha,
// 	      const GivaroZpz<Std32>::Element * A, const size_t lda,
// 	      const GivaroZpz<Std32>::Element * X, const size_t incX,
// 	      const  GivaroZpz<Std32>::Element beta,
// 	      GivaroZpz<Std32>::Element * Y, const size_t incY) {
// 	/* y := beta*y + alpha*A*x */
// 	// dpritcha -- for some reason this used to assume alpha=-1, beta=1. fixed now.
// 	// todo: optimize for alpha=0

//        static const long MOD = F.characteristic();
//        static const long MODmun = MOD-1;
//        static unsigned long   temp =  ((unsigned long)(-MODmun) / MODmun) / MODmun;
//        static const unsigned long Kbest = (temp<1?1:temp);
//        static const unsigned long CORR = ( ((2147483648UL) % (unsigned long)MOD)<<1) % (unsigned long)MOD;
//        const unsigned long DIMoKbest = (N / Kbest);
//        const unsigned long DIMoKfKbest = DIMoKbest*Kbest;

//        GivaroZpz<Std32>::Element * Yit = Y;
//        const GivaroZpz<Std32>::Element * Ait = A;//, * Xit = X;
//        GivaroZpz<Std32>::Element best, inter;

//        if (TransA == FflasTrans) 
// 	       std::cerr << "Warning, FflasTrans fgemv not implemented over Givaro<std32>\n";

//        if (DIMoKfKbest == 0)
//            for(unsigned long i=0; i<M; ++i) {
//                best = 0;
//                for(unsigned long j=0; j<N; ++j)
// 		       //    best += Ait[j] * X[j];
// 		       //	       if incX != 1, Change to : 
// 		       best += Ait[j] *  (*(X+j*incX));
	       
//                best %= MOD;
// 	       if (alpha != 1) { best *= alpha; best %= MOD;}

// 	       if (beta == 0) *Yit = best;
// 	       else {
// 		       if (beta != 1) { 
// 			       if (beta == -1) {*Yit -= beta; if (*Yit < 0) *Yit += MOD;}
// 			       else {*Yit *= beta; *Yit %= MOD;}
// 		       }
// 		       *Yit += best;
// 		       if (*Yit >= MOD) *Yit -= MOD;
// 	       }

//                Yit += incY;
//                Ait += lda;
//            }
//        else
//            for(unsigned long i=0; i<M; ++i) {
//                best = inter = 0;
//                unsigned long k = 0;
//                for(unsigned long  j=0 ; j < Kbest; ++j, ++k)
// 		       //                   best += Ait[k] * X[k];
// // 		      if incX != 1, Change to : 
// 		       best += Ait[k] *  (*(X+k*incX));
//                for(inter = best ; k < DIMoKfKbest; inter=best) {
//                    for(unsigned long  j=0 ; j < Kbest; ++j, ++k)
// 			   //                       best += Ait[k] * X[k];
// // 		      if incX != 1, Change to : 
// 			   best += Ait[k] *  (*(X+k*incX));

//                    if (inter > best) best += CORR;
               
//                }
//                inter = best;
//                for(; k<N; ++k)
// 		       //                   best += Ait[k] * X[k];
// // 		      if incX != 1, Change to : 
// 		       best += Ait[k] *  (*(X+k*incX));
//                if (inter > best) best += CORR;

//                best %= MOD;
// 	       if (alpha != 1) { best *= alpha; best %= MOD;}

// 	       if (beta == 0) *Yit = best;
// 	       else {
// 		       if (beta != 1) { 
// 			       if (beta == -1) {*Yit -= beta; if (*Yit < 0) *Yit += MOD;}
// 			       else {*Yit *= beta; *Yit %= MOD;}
// 		       }
// 		       *Yit += best;
// 		       if (*Yit >= MOD) *Yit -= MOD;
// 	       }


//                Yit += incY;
//                Ait += lda;
//            }
//    }


// template<>
// inline void
// FFLAS::fgemv (const GivaroZpz<Std64>& F, const enum FFLAS_TRANSPOSE TransA,
// 	      const size_t M, const size_t N,
// 	      const GivaroZpz<Std64>::Element alpha,
// 	      const GivaroZpz<Std64>::Element * A, const size_t lda,
// 	      const GivaroZpz<Std64>::Element * X, const size_t incX,
// 	      const  GivaroZpz<Std64>::Element beta,
// 	      GivaroZpz<Std64>::Element * Y, const size_t incY) {
// 	/* y := beta*y + alpha*A*x */
// 	// dpritcha -- for some reason this used to assume alpha=-1, beta=1. fixed now.
// 	// todo: optimize for alpha=0

// 	if (TransA == FflasTrans) 
// 		std::cerr << "Warning, FflasTrans fgemv not implemented over Givaro<std64>\n";

// 	static const long MOD=F.size();
// 	GivaroZpz<Std64>::Element * Yit = Y;
// 	const GivaroZpz<Std64>::Element * Ait = A;//, * Xit = X;
// 	GivaroZpz<Std64>::Element acc;
// 	for(unsigned long i=0; i<M; ++i) {
//                acc = 0;
//                for(unsigned long j=0; j<N; ++j)
// 		       acc += Ait[j] *  (*(X+j*incX));

//                acc %= MOD;
// 	       if (alpha != 1) { acc *= alpha; acc %= MOD;}

// 	       if (beta == 0) *Yit = acc;
// 	       else {
// 		       if (beta != 1) { 
// 			       if (beta == -1) {*Yit -= beta; if (*Yit < 0) *Yit += MOD;}
// 			       else {*Yit *= beta; *Yit %= MOD;}
// 		       }
// 		       *Yit += acc;
// 		       if (*Yit >= MOD) *Yit -= MOD;
// 	       }

//                Yit += incY;
//                Ait += lda;

// 	}
// }
#endif
}//namespace LinBox
