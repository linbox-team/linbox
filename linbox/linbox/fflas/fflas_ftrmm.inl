/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflas/fflas_ftrmm.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

//---------------------------------------------------------------------
// ftrmm: TRiangular Matrix Multiply
// Computes  B <- alpha.op(A).B,  B <- alpha.B.op(A)
// B is M*N, A is M*M if Side==FflasLeft, N*N if Side==FflasRight
// Warning :Assumes alpha==1
//---------------------------------------------------------------------
#ifndef DOUBLE_MANTISSA
#define DOUBLE_MANTISSA 53
#endif

// dpritcha: added functions to compute bound
// had overflow issues and static problems before; needs to work over multiple moduli
size_t mul_bound_compute(const long long pi) {
	long long p=pi;//,p1=1,p2=1;
	long long nmax = ( (  1ULL<<(DOUBLE_MANTISSA) )/((p-1)*(p-1)))-1;
	//if overflow, return maxint32. if people ever use matrices bigger than 2 billion to a side this
	//will be trouble. 
	if (nmax >= 2147483647) return 2147483647; 
	return nmax;

}
size_t mul_bound(const long long pi) {
	static long long p=pi;
	static size_t nmax=mul_bound_compute(pi);
	if (p == pi) 
		return nmax;
	else 
		return nmax=mul_bound_compute(p=pi);
}


template<class Field>
inline void
LinBox::FFLAS::ftrmm(const Field& F, const enum FFLAS_SIDE Side,
	     const enum FFLAS_UPLO Uplo, 
	     const enum FFLAS_TRANSPOSE TransA,
	     const enum FFLAS_DIAG Diag, 
	     const size_t M, const size_t N,
	     const typename Field::Element alpha,
	     typename Field::Element * A, const size_t lda,
	     typename Field::Element * B, const size_t ldb){
	
	if (!M || !N ) return; 
	
	integer pi;
	F.characteristic(pi);
	long long p = pi;
	size_t nmax = mul_bound(p);
	//cout << "nmax="<<nmax<<"\n";
	
	if ( Side==FflasLeft ){
		if ( Uplo==FflasUpper){
			if (TransA == FflasNoTrans){
				ftrmmLeftUpNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
			else{
				ftrmmLeftUpTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
		else{
			if (TransA == FflasNoTrans){
				ftrmmLeftLowNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb, nmax);
			}
			else{
				ftrmmLeftLowTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
	}
	else{
	if ( Uplo==FflasUpper){
			if (TransA == FflasNoTrans){
				ftrmmRightUpNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb,nmax);
			}
			else{
				ftrmmRightUpTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
		else{
			if (TransA == FflasNoTrans){
				ftrmmRightLowNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb, nmax);
			}
			else{
				ftrmmRightLowTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
	}
}

template<class Field>
inline void 
LinBox::FFLAS::ftrmmLeftUpNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element alpha,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb){
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( M==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::Element inv;
			F.inv(inv, *A);
			fscal(F, N, inv, B, 1);
		}
	}
	else{
		size_t Mup=M>>1;
		size_t Mdown = M-Mup;
		ftrmmLeftUpNoTrans( F, Diag, Mdown, N, alpha, 
				    A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
		fgemm( F, FflasNoTrans, FflasNoTrans, Mup, N, Mdown, Mone,
		       A+Mup, lda, B+Mup*ldb, ldb, alpha, B, ldb);
		ftrmmLeftUpNoTrans( F, Diag, Mup, N, one, A, lda, B, ldb);
	}
}

template<class Field>
inline void
LinBox::FFLAS::ftrmmLeftUpTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		      const size_t M, const size_t N,
		      const typename Field::Element alpha,
		      const typename Field::Element * A, const size_t lda,
		      typename Field::Element * B, const size_t ldb){

	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( M==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::Element inv;
			F.inv(inv, *A);
			fscal(F, N, inv, B, 1);
		}
	}
	else{
		size_t Mup=M>>1;
		size_t Mdown = M-Mup;
		ftrmmLeftUpTrans( F, Diag, Mdown, N, alpha, 
				  A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
		fgemm( F, FflasTrans, FflasNoTrans, Mup, N, Mdown, Mone, 
		       A+Mup*lda, lda, B+Mup*ldb, ldb, alpha, B, ldb);
		ftrmmLeftUpTrans( F, Diag, Mup, N, one, A, lda, B, ldb);
	}
}

template<class Field>
inline void
LinBox::FFLAS::ftrmmLeftLowNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element alpha,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, const size_t nmax){

	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( M <= nmax ){
		double alphad;
		if (F.areEqual(alpha, Mone))
			alphad = -1.0;
		else
			F.convert( alphad, alpha );
		DoubleDomain::Element * Ad = new DoubleDomain::Element[M*M];
		DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
		MatF2MatD( F, Ad, A, lda, M, M );
		MatF2MatD( F, Bd, B, ldb, M, N );
		cblas_dtrmm(  CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
			      (CBLAS_DIAG) Diag, M, N, alphad, Ad, M, Bd, N );
		delete[] Ad;
		MatD2MatF( F, B, ldb, Bd, M, N );
		delete[] Bd;
	}
	else{
		size_t Mup=M>>1;
		size_t Mdown = M-Mup;
		ftrmmLeftLowNoTrans( F, Diag, Mdown, N, one, 
				     A+Mup*(lda+1), lda, B+Mup*ldb, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, Mdown, N, Mup,
		       one, A+Mup*lda, lda, B, ldb, alpha, B+Mup*ldb, ldb);
		ftrmmLeftLowNoTrans( F, Diag, Mup, N, alpha, A, lda, B, ldb, nmax);
	}
}

template<>
inline void
LinBox::FFLAS::ftrmmLeftLowNoTrans(const Modular<double>& F, const enum FFLAS_DIAG Diag, 
			   const size_t M, const size_t N,
			   const double alpha,
			   double * A, const size_t lda,
			   double * B, const size_t ldb, const size_t nmax){
	
	static double Mone;
	static double one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( M <= nmax ){
		cblas_dtrmm(  CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
			      (CBLAS_DIAG) Diag, M, N, alpha, A, lda, B, ldb );
		for (size_t i=0; i< M; ++i)
			for (size_t j=0; j<N; ++j)
				F.init(*(B+i*ldb+j),*(B+i*ldb+j));
	}
	else{
		size_t Mup=M>>1;
		size_t Mdown = M-Mup;
		ftrmmLeftLowNoTrans( F, Diag, Mdown, N, one, 
				     A+Mup*(lda+1), lda, B+Mup*ldb, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, Mdown, N, Mup,
		       one, A+Mup*lda, lda, B, ldb, alpha, B+Mup*ldb, ldb);
		ftrmmLeftLowNoTrans( F, Diag, Mup, N, alpha, A, lda, B, ldb, nmax);
	}
}

template<class Field>
inline void 
LinBox::FFLAS::ftrmmLeftLowTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		       const size_t M, const size_t N,
		       const typename Field::Element alpha,
		       const typename Field::Element * A, const size_t lda,
		       typename Field::Element * B, const size_t ldb){

	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( M==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::Element inv;
			F.inv(inv, *A);
			fscal(F, N, inv, B, 1);
		}
	}
	else{
		size_t Mup=M>>1;
		size_t Mdown = M-Mup;
		ftrmmLeftLowTrans( F, Diag, Mup, N, alpha, A, lda, B, ldb);
		fgemm( F, FflasTrans, FflasNoTrans, Mdown, N, Mup,
		       Mone, A+Mup, lda, B, ldb, alpha, B+Mup*ldb, ldb);
		ftrmmLeftLowTrans( F, Diag, Mdown, N, one, 
				   A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
	}
}

template<class Field>
inline void 
LinBox::FFLAS::ftrmmRightUpNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::Element alpha,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb, const size_t nmax){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);

	if ( N <= nmax ){
		
		double alphad;
		if (F.areEqual(alpha, Mone))
			alphad = -1.0;
		else
			F.convert( alphad, alpha );
		DoubleDomain::Element * Ad = new DoubleDomain::Element[N*N];
		DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
		MatF2MatD( F, Ad, A, lda, N, N );
		MatF2MatD( F, Bd, B, ldb, M, N );
		cblas_dtrmm(  CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
			      (CBLAS_DIAG) Diag, M, N, alphad, Ad, N, Bd, N );
		delete[] Ad;
		MatD2MatF( F, B, ldb, Bd, M, N );
		delete[] Bd;
	}
	else{
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrmmRightUpNoTrans( F, Diag, M, Ndown, one, 
				     A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, M, Ndown, Nup,
		       one, B, ldb, A+Nup, lda, alpha, B+Nup, ldb);
		ftrmmRightUpNoTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb, nmax);
	}
}

template<>
inline void 
LinBox::FFLAS::ftrmmRightUpNoTrans(const Modular<double>& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const double alpha,
			 double * A, const size_t lda,
			 double * B, const size_t ldb, const size_t nmax){
	
	static double Mone;
	static double one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( N <= nmax ){
		cblas_dtrmm(  CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
			      (CBLAS_DIAG)Diag, M, N, alpha, A, lda, B, ldb );
		for (size_t i=0; i< M; ++i)
			for (size_t j=0; j<N; ++j)
				F.init(*(B+i*ldb+j));
	}
	else{
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrmmRightUpNoTrans( F, Diag, M, Ndown, one, 
				     A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, M, Ndown, Nup,
		       one, B, ldb, A+Nup, lda, alpha, B+Nup, ldb);
		ftrmmRightUpNoTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb, nmax);
		}
}


template<class Field>
inline void
LinBox::FFLAS::ftrmmRightUpTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		       const size_t M, const size_t N,
		       const typename Field::Element alpha,
		       const typename Field::Element * A, const size_t lda,
		       typename Field::Element * B, const size_t ldb){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( N==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::Element inv;
			F.inv(inv, *A);
			fscal(F, M, inv, B, ldb);
		}
	}
	else{	
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrmmRightUpTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb);
		fgemm( F, FflasNoTrans, FflasTrans, M, Ndown, Nup, Mone, 
		       B, ldb, A+Nup*lda, lda, alpha, B+Nup, ldb);
		ftrmmRightUpTrans( F, Diag, M, Ndown, one, 
				   A+Nup*(lda+1), lda, B+Nup, ldb);
	}
}


template<class Field>
inline void
LinBox::FFLAS::ftrmmRightLowNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			    const size_t M, const size_t N,
			    const typename Field::Element alpha,
			    const typename Field::Element * A, const size_t lda,
			    typename Field::Element * B, const size_t ldb, const size_t nmax){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( N <= nmax ){
		double alphad;
		if (F.areEqual(alpha, Mone))
			alphad = -1.0;
		else
			F.convert( alphad, alpha );
		DoubleDomain::Element * Ad = new DoubleDomain::Element[N*N];
		DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
		MatF2MatD( F, Ad, A, lda, N, N );
		MatF2MatD( F, Bd, B, ldb, M, N );
		cblas_dtrmm(  CblasRowMajor, CblasRight, CblasLower, CblasNoTrans,
			      (CBLAS_DIAG) Diag, M, N, alphad, Ad, N, Bd, N );
		delete[] Ad;
		MatD2MatF( F, B, ldb, Bd, M, N );
		delete[] Bd;
	}
	else{
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrmmRightLowNoTrans( F, Diag, M, Ndown, alpha, 
				      A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, M, Nup, Ndown,
		       one, B+Nup, ldb, A+Nup*lda, lda, alpha, B, ldb);
		ftrmmRightLowNoTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb, nmax);
	}
}
template<>
inline void
LinBox::FFLAS::ftrmmRightLowNoTrans(const Modular<double>& F, const enum FFLAS_DIAG Diag, 
			    const size_t M, const size_t N,
			    const double alpha,
			    const double * A, const size_t lda,
			    double * B, const size_t ldb, const size_t nmax){
	
	static double Mone;
 	static double one;
 	F.init(Mone, -1);
 	F.init(one, 1);
	if ( N <= nmax ){
		// pascal 2004-10-12, 
		// there is a problem here if alpha > 1, the bound is wrong, need to be fix
		cblas_dtrmm(  CblasRowMajor, CblasRight, CblasLower, CblasNoTrans,
			      (CBLAS_DIAG) Diag, M, N, alpha, A, lda, B, ldb );
	
		for (size_t i=0; i< M; ++i)
			for (size_t j=0; j<N; ++j)
				F.init(*(B+i*ldb+j),*(B+i*ldb+j));			
			
	}
	else{
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrmmRightLowNoTrans( F, Diag, M, Ndown, alpha, 
				      A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, M, Nup, Ndown,
		       one, B+Nup, ldb, A+Nup*lda, lda, alpha, B, ldb);
		ftrmmRightLowNoTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb, nmax);	
	}	
}

template<class Field>
inline void
LinBox::FFLAS::ftrmmRightLowTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element alpha,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	if ( N==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::Element inv;
			F.inv(inv, *A);
			fscal(F, M, inv, B, ldb);
		}
	}
	else{
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrmmRightLowTrans( F, Diag, M, Ndown, alpha, 
				    A+Nup*(lda+1), lda, B+Nup, ldb);
		fgemm( F, FflasNoTrans, FflasTrans, M, Nup, Ndown, Mone, 
		       B+Nup, ldb, A+Nup, lda, alpha, B, ldb);
		ftrmmRightLowTrans( F, Diag, M, Nup, one, A, lda, B, ldb);
	}
}


