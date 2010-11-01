
/* fflas/fflas_fgemv.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */


//---------------------------------------------------------------------
// fgemv: GEneral Matrix Vector Multiplication
// Computes  Y <- alpha.op(A).X + beta.Y
// A is M*N, 
//---------------------------------------------------------------------
template<class Field>
inline void 
FFLAS::fgemv (const Field& F, const FFLAS_TRANSPOSE TransA, 
	      const size_t M, const size_t N,
	      const typename Field::Element alpha, 
	      const typename Field::Element * A, const size_t lda,
	      const typename Field::Element * X, const size_t incX,
	      const typename Field::Element beta,
	      typename Field::Element * Y, const size_t incY)
{
 	static typename Field::Element  one, mone, zero;
	F.init(one,1UL);
 	F.neg(mone,one);
 	F.init(zero,0UL);

	if (F.isZero (alpha)){
		for  (typename Field::Element * Yi = Y; Yi != Y+((TransA == FflasNoTrans)?M:N)*incY; Yi+=incY)
			F.mulin(*Yi, beta);
		return;
	}
	
	FFLAS_BASE base = BaseCompute (F, 0);
	size_t kmax = DotProdBound (F, 0, beta, base);

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
FFLAS::MatVectProd (const Field& F, const FFLAS_TRANSPOSE TransA, 
		    const size_t M, const size_t N,
		    const typename Field::Element alpha, 
		    const typename Field::Element * A, const size_t lda,
		    const typename Field::Element * X, const size_t incX,
		    const typename Field::Element beta,
		    typename Field::Element * Y, const size_t incY){
	typename Field::Element  one, mone;
	typename Field::Element tmp;
	F.init(one,1.0);
	F.init(mone,-1.0);
	
	size_t Xl, Yl;
	if  (TransA == FflasNoTrans){Xl = N;Yl = M;}
	else {Xl = M; Yl = N;}
	double* Ad = new double[M*N];
	double* Xd = new double[Xl];
	double* Yd = new double[Yl];
	double alphad, betad;
	
	if (F.areEqual (mone, alpha)){
		alphad = -1.0;
		F.convert (betad, beta);
	} else {
		if (! F.areEqual (one, alpha)){
			// Compute C = A*B + beta/alpha.C
			// and after C *= alpha
			F.div (tmp, beta, alpha);
			F.convert (betad, tmp);
		} else
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
	
	cblas_dgemv (CblasRowMajor, (CBLAS_TRANSPOSE) TransA, M, N, alphad,
		     Ad, N, Xd, 1, betad, Yd, 1);
	
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
inline void FFLAS::MatVectProd (const ModularBalanced<double>& F,
				const FFLAS_TRANSPOSE TransA, 
				const size_t M, const size_t N,
				const double alpha, 
				const double * A, const size_t lda,
				const double * X, const size_t incX,
				const double beta,
				double * Y, const size_t incY){
	
	double Mone, one, _alpha, _beta;
	F.init(one, 1UL);	
	F.neg(Mone, one);
	if (F.areEqual (Mone, beta)) _beta = -1.0;
	else _beta = beta;
	
	if (F.areEqual (Mone, alpha)) _alpha = -1.0;
	else{
		_alpha = 1.0;
		if (! F.areEqual (one, alpha))
			// Compute y = A*x + beta/alpha.y
			// and after y *= alpha
			F.divin (_beta, alpha);
	}
	
	cblas_dgemv (CblasRowMajor, (CBLAS_TRANSPOSE) TransA, M, N, 
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
inline void FFLAS::MatVectProd (const ModularBalanced<float>& F,
				const FFLAS_TRANSPOSE TransA, 
				const size_t M, const size_t N,
				const float alpha, 
				const float * A, const size_t lda,
				const float * X, const size_t incX,
				const float beta,
				float * Y, const size_t incY){
	
	float Mone, one, _alpha, _beta;
	F.init(one, 1UL);	
	F.neg(Mone, one);
	if  (F.areEqual (Mone, beta)) _beta = -1.0;
	else _beta = beta;

	if (F.areEqual (Mone, alpha)) _alpha = -1.0;
	else{
		_alpha = 1.0;
		if (! F.areEqual (one, alpha)){
			// Compute y = A*x + beta/alpha.y
			// and after y *= alpha
			F.divin (_beta, alpha);
		}
	}
	cblas_sgemv (CblasRowMajor, (CBLAS_TRANSPOSE) TransA, M, N, 
		     _alpha, A, lda, X, incX, _beta, Y, incY);
	for  (float * Yi = Y; Yi != Y+((TransA == FflasNoTrans)?M:N)*incY; Yi+=incY)
		F.init (*Yi, *Yi);
	if ( (!F.areEqual (one, alpha)) && (!F.areEqual (Mone, alpha))){
		// Fix-up: compute y *= alpha
		for (float* Yi = Y; Yi != Y+((TransA == FflasNoTrans)?M:N)*incY; Yi += incY)
			F.mulin (*Yi , alpha);
	}
}

template<>
inline void FFLAS::MatVectProd (const Modular<double>& F,
				const FFLAS_TRANSPOSE TransA, 
				const size_t M, const size_t N,
				const double alpha, 
				const double * A, const size_t lda,
				const double * X, const size_t incX,
				const double beta,
				double * Y, const size_t incY){
	
	double Mone, one, _alpha, _beta;
	F.init(one, 1UL);	
	F.neg(Mone, one);
	if (F.areEqual (Mone, beta)) _beta = -1.0;
	else _beta = beta;
	
	if (F.areEqual (Mone, alpha)) _alpha = -1.0;
	else{
		_alpha = 1.0;
		if (! F.areEqual (one, alpha))
			// Compute y = A*x + beta/alpha.y
			// and after y *= alpha
			F.divin (_beta, alpha);
	}
	
	cblas_dgemv (CblasRowMajor, (CBLAS_TRANSPOSE) TransA, M, N, 
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
inline void FFLAS::MatVectProd (const Modular<float>& F,
				const FFLAS_TRANSPOSE TransA, 
				const size_t M, const size_t N,
				const float alpha, 
				const float * A, const size_t lda,
				const float * X, const size_t incX,
				const float beta,
				float * Y, const size_t incY){
	
	float Mone, one, _alpha, _beta;
	F.init(one, 1UL);	
	F.neg(Mone, one);
	if  (F.areEqual (Mone, beta)) _beta = -1.0;
	else _beta = beta;

	if (F.areEqual (Mone, alpha)) _alpha = -1.0;
	else{
		_alpha = 1.0;
		if (! F.areEqual (one, alpha)){
			// Compute y = A*x + beta/alpha.y
			// and after y *= alpha
			F.divin (_beta, alpha);
		}
	}
	cblas_sgemv (CblasRowMajor, (CBLAS_TRANSPOSE) TransA, M, N, 
		     _alpha, A, lda, X, incX, _beta, Y, incY);
	for  (float * Yi = Y; Yi != Y+((TransA == FflasNoTrans)?M:N)*incY; Yi+=incY)
		F.init (*Yi, *Yi);
	if ( (!F.areEqual (one, alpha)) && (!F.areEqual (Mone, alpha))){
		// Fix-up: compute y *= alpha
		for (float* Yi = Y; Yi != Y+((TransA == FflasNoTrans)?M:N)*incY; Yi += incY)
			F.mulin (*Yi , alpha);
	}
}

template<>
inline void
FFLAS::fgemv (const DoubleDomain& , const FFLAS_TRANSPOSE TransA, 
	      const size_t M, const size_t N,
	      const DoubleDomain::Element  alpha, 
	      const DoubleDomain::Element * A, const size_t lda,
	      const DoubleDomain::Element * X, const size_t incX,
	      const DoubleDomain::Element beta,
	      DoubleDomain::Element * Y, const size_t incY) {
	cblas_dgemv (CblasRowMajor, (CBLAS_TRANSPOSE) TransA, M, N, 
		     alpha, A, lda, X, incX, beta, Y, incY);
}

template<>
inline void
FFLAS::fgemv (const FloatDomain& , const FFLAS_TRANSPOSE TransA, 
	      const size_t M, const size_t N,
	      const FloatDomain::Element  alpha, 
	      const FloatDomain::Element * A, const size_t lda,
	      const FloatDomain::Element * X, const size_t incX,
	      const FloatDomain::Element beta,
	      FloatDomain::Element * Y, const size_t incY)
{
	cblas_sgemv (CblasRowMajor, (CBLAS_TRANSPOSE) TransA, M, N, 
		     alpha, A, lda, X, incX, beta, Y, incY);
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
