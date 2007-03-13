/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* fflas/fflas_ftrsm.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */


//---------------------------------------------------------------------
// ftrsm: TRiangular System solve with matrix
// Computes  B <- alpha.op(A^-1).B,  B <- alpha.B.op(A^-1)
// B is M*N, A is M*M if Side==FflasLeft, N*N if Side==FflasRight
// Warning :Assumes alpha==1
//---------------------------------------------------------------------
template<class Field>
inline void
FFLAS::ftrsm (const Field& F, const enum FFLAS_SIDE Side,
	      const enum FFLAS_UPLO Uplo, 
	      const enum FFLAS_TRANSPOSE TransA,
	      const enum FFLAS_DIAG Diag, 
	      const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      typename Field::Element * A, const size_t lda,
	      typename Field::Element * B, const size_t ldb) 
{
	if (!M || !N ) return; 

	size_t nmax = TRSMBound<Field> (F);

	if ( Side==FflasLeft ){
		if ( Uplo==FflasUpper){
			
			if (TransA == FflasNoTrans){
				ftrsmLeftUpNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
			else{
				ftrsmLeftUpTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
		else{
			if (TransA == FflasNoTrans){
				ftrsmLeftLowNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb, nmax);
			}
			else{
				ftrsmLeftLowTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
	}
	else{
	if ( Uplo==FflasUpper){
			if (TransA == FflasNoTrans){
				ftrsmRightUpNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb,nmax);
			}
			else{
				ftrsmRightUpTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
		else{
			if (TransA == FflasNoTrans){
				ftrsmRightLowNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
			else{
				ftrsmRightLowTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
	}
	
}

template<class Field>
inline void 
FFLAS::ftrsmLeftUpNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element alpha,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb){
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
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
		ftrsmLeftUpNoTrans( F, Diag, Mdown, N, alpha, 
				    A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
		fgemm( F, FflasNoTrans, FflasNoTrans, Mup, N, Mdown, Mone,
		       A+Mup, lda, B+Mup*ldb, ldb, alpha, B, ldb);
		ftrsmLeftUpNoTrans( F, Diag, Mup, N, one, A, lda, B, ldb);
	}

}

template<class Field>
inline void
FFLAS::ftrsmLeftUpTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		      const size_t M, const size_t N,
		      const typename Field::Element alpha,
		      const typename Field::Element * A, const size_t lda,
		      typename Field::Element * B, const size_t ldb){

	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
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
		ftrsmLeftUpTrans( F, Diag, Mdown, N, alpha, 
				  A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
		fgemm( F, FflasTrans, FflasNoTrans, Mup, N, Mdown, Mone, 
		       A+Mup*lda, lda, B+Mup*ldb, ldb, alpha, B, ldb);
		ftrsmLeftUpTrans( F, Diag, Mup, N, one, A, lda, B, ldb);
	}
}
template<class Field>
inline void
FFLAS::ftrsmLeftLowNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			   const size_t M, const size_t N,
			   const typename Field::Element alpha,
			   typename Field::Element * A, const size_t lda,
			   typename Field::Element * B, const size_t ldb, const size_t nmax){

	if (AreEqual<typename Field::Element, double>::value)
		ftrsmLeftLowNoTrans_dbl(F,Diag,M,N,alpha,A,lda,B,ldb,nmax);
	else
		ftrsmLeftLowNoTrans_gen(F,Diag,M,N,alpha,A,lda,B,ldb,nmax);
}

template<class Field>
inline void
FFLAS::ftrsmLeftLowNoTrans_gen(const Field& F, const enum FFLAS_DIAG Diag, 
			       const size_t M, const size_t N,
			       const typename Field::Element alpha,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * B, const size_t ldb, const size_t nmax){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
	if ( M <= nmax ){
		typename Field::Element inv;
		if (Diag == FflasNonUnit ){
			//Normalization of A and correction of B
			// A<-DA, B<-DB
			typename Field::Element * Ai = A;
			typename Field::Element * Bi = B;
			for (size_t i=0; i<M; ++i){
				F.inv( inv, *(Ai+i) );
				fscal(F, i, inv, Ai, 1 );
				fscal(F, N, inv, Bi, 1 );
				Ai += lda; Bi+=ldb;
			}
		}
		double alphad;
		if (F.areEqual(alpha, Mone))
			alphad = -1.0;
		else
			F.convert( alphad, alpha );
		DoubleDomain::Element * Ad = new DoubleDomain::Element[M*M];
		DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
		MatF2MatD( F, Ad, N, A, lda, M, M );
		MatF2MatD( F, Bd, N, B, ldb, M, N );
		cblas_dtrsm(  CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
			      CblasUnit, M, N, alphad, Ad, M, Bd, N );
		delete[] Ad;
		MatD2MatF( F, B, ldb, Bd, N, M, N );
		delete[] Bd;
		if (Diag == FflasNonUnit ){
			//Denormalization of A
			// A-> D^(-1)A
			typename Field::Element *  Ai=A;
			for (size_t i=0; i<M; ++i){
				fscal( F, i, *(Ai+i), Ai, 1 );
				Ai += lda;
			}
		}
	}
	else{
		size_t Mup=M>>1;
		size_t Mdown = M-Mup;
		ftrsmLeftLowNoTrans( F, Diag, Mup, N, alpha, A, lda, B, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, Mdown, N, Mup,
		       Mone, A+Mup*lda, lda, B, ldb, alpha, B+Mup*ldb, ldb);
		ftrsmLeftLowNoTrans( F, Diag, Mdown, N, one, 
				     A+Mup*(lda+1), lda, B+Mup*ldb, ldb, nmax);
	}
}

template<class Field>
inline void
FFLAS::ftrsmLeftLowNoTrans_dbl(const Field& F, const enum FFLAS_DIAG Diag, 
			       const size_t M, const size_t N,
			       const double alpha,
			       double * A, const size_t lda,
			       double * B, const size_t ldb, const size_t nmax){
	static double Mone;
	static double one;
	F.init(one, 1.0);
	F.neg(Mone, one);

	if ( M <= nmax ){ 
		double inv;
		if (Diag == FflasNonUnit ){
			//Normalization of A and correction of B
			double * Ai = A;
			double * Bi = B;
			for (size_t i=0; i<M; ++i){
				F.inv( inv, *(Ai+i) );
				fscal(F, i, inv, Ai, 1 );
				fscal(F, N, inv, Bi, 1 );
				Ai += lda; Bi+=ldb;

			}
		}
		
		cblas_dtrsm(  CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
			      CblasUnit, M, N, alpha, A, lda, B, ldb );
		for (size_t i=0; i< M; ++i)
			for (size_t j=0; j<N; ++j)
				F.init(*(B+i*ldb+j),*(B+i*ldb+j));
	
	if (Diag == FflasNonUnit ){
			//Denormalization of A
			double *  Ai=A;
			for (size_t i=0; i<M; ++i){
				fscal( F, i, *(Ai+i), Ai, 1 );
				Ai += lda;
			}
		}
	}
	else{
		size_t Mup=M>>1;
		size_t Mdown = M-Mup;
		ftrsmLeftLowNoTrans( F, Diag, Mup, N, alpha, A, lda, B, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, Mdown, N, Mup,
		       Mone, A+Mup*lda, lda, B, ldb, alpha, B+Mup*ldb, ldb);
		ftrsmLeftLowNoTrans( F, Diag, Mdown, N, one, 
				     A+Mup*(lda+1), lda, B+Mup*ldb, ldb, nmax);
	}

}

template<class Field>
inline void 
FFLAS::ftrsmLeftLowTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		       const size_t M, const size_t N,
		       const typename Field::Element alpha,
		       const typename Field::Element * A, const size_t lda,
		       typename Field::Element * B, const size_t ldb){

	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
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
		ftrsmLeftLowTrans( F, Diag, Mup, N, alpha, A, lda, B, ldb);
		fgemm( F, FflasTrans, FflasNoTrans, Mdown, N, Mup,
		       Mone, A+Mup, lda, B, ldb, alpha, B+Mup*ldb, ldb);
		ftrsmLeftLowTrans( F, Diag, Mdown, N, one, 
				   A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
	}
}

template<class Field>
inline void 
FFLAS::ftrsmRightUpNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			   const size_t M, const size_t N,
			   const typename Field::Element alpha,
			   typename Field::Element * A, const size_t lda,
			   typename Field::Element * B, const size_t ldb, const size_t nmax){
	if (AreEqual<typename Field::Element,double>::value)
		ftrsmRightUpNoTrans_dbl(F,Diag,M,N,alpha,A,lda,B,ldb,nmax);
	else
		ftrsmRightUpNoTrans_gen(F,Diag,M,N,alpha,A,lda,B,ldb,nmax);

}
template<class Field>
inline void 
FFLAS::ftrsmRightUpNoTrans_gen(const Field& F, const enum FFLAS_DIAG Diag, 
			       const size_t M, const size_t N,
			       const typename Field::Element alpha,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * B, const size_t ldb, const size_t nmax){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
		
	if ( N <= nmax ){
		typename Field::Element inv;
		if (Diag == FflasNonUnit){
			//Normalization of A and B
			typename Field::Element *  Ai = A, * Bi = B;
			for (size_t i=0; i<N; ++i){
				F.inv( inv, *(Ai+i*lda) );
				fscal( F, i, inv, Ai, lda );
				fscal( F, M, inv, Bi, ldb );
				Ai++;
				Bi++;
			}
		}
		double alphad;
		if (F.areEqual(alpha, Mone))
			alphad = -1.0;
		else
			F.convert( alphad, alpha );
		DoubleDomain::Element * Ad = new DoubleDomain::Element[N*N];
		DoubleDomain::Element * Bd = new DoubleDomain::Element[M*N];
		MatF2MatD( F, Ad, N, A, lda, N, N );
		MatF2MatD( F, Bd, N, B, ldb, M, N );
		cblas_dtrsm(  CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
			      CblasUnit, M, N, alphad, Ad, N, Bd, N );
		delete[] Ad;
		MatD2MatF( F, B, ldb, Bd, N, M, N );
		delete[] Bd;
		if (Diag == FflasNonUnit ){
			//Denormalization of A
			typename Field::Element *  Ai=A;
			for (size_t i=0; i<N; ++i){
				fscal( F, i, *(Ai+i*lda), Ai, lda );
				Ai++;
			}	
		}
	}
	else{
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrsmRightUpNoTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, M, Ndown, Nup,
		       Mone, B, ldb, A+Nup, lda, alpha, B+Nup, ldb);
		ftrsmRightUpNoTrans( F, Diag, M, Ndown, one, 
				     A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
	}

}

template<class Field>
inline void 
FFLAS::ftrsmRightUpNoTrans_dbl (const Field& F, const enum FFLAS_DIAG Diag, 
			    const size_t M, const size_t N, const double alpha,
			    double * A, const size_t lda,
			    double * B, const size_t ldb, const size_t nmax){
	
	
	static double Mone;
	static double one;
	F.init(one, 1.0);
	F.neg(Mone,one);
	if ( N <= nmax ){
		double inv;
		if (Diag == FflasNonUnit ){
			//Normalization of A and B
			double *  Ai = A, * Bi = B;
			for (size_t i=0; i<N; ++i){
				F.inv( inv, *(Ai+i*lda) );
				fscal( F, i, inv, Ai, lda );
				fscal( F, M, inv, Bi, ldb );
				Ai++;
				Bi++;
			}
		}
		cblas_dtrsm(  CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
			      CblasUnit, M, N, alpha, A, lda, B, ldb );
		for (size_t i=0; i< M; ++i)
			for (size_t j=0; j<N; ++j){
				F.init(*(B+i*ldb+j),*(B+i*ldb+j));
				
			}
		if (Diag == FflasNonUnit ){
			//Denormalization of A
			double *  Ai=A;
			for (size_t i=0; i<N; ++i){
				fscal( F, i, *(Ai+i*lda), Ai, lda );
				Ai++;
			}
			//Correction on B
			// Ai =A;
// 			double *Bi=B;
// 			for (size_t i=0; i<N; ++i){
// 				F.inv( inv, *Ai);
// 				fscal( F, M, inv, Bi, ldb );
// 				Ai += lda+1; Bi++;
// 			}
		}
	}
	else{
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrsmRightUpNoTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb, nmax);
		fgemm( F, FflasNoTrans, FflasNoTrans, M, Ndown, Nup,
		       Mone, B, ldb, A+Nup, lda, alpha, B+Nup, ldb);
		ftrsmRightUpNoTrans( F, Diag, M, Ndown, one, 
				     A+Nup*(lda+1), lda, B+Nup, ldb, nmax);
	}
}


template<class Field>
inline void
FFLAS::ftrsmRightUpTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		       const size_t M, const size_t N,
		       const typename Field::Element alpha,
		       const typename Field::Element * A, const size_t lda,
		       typename Field::Element * B, const size_t ldb){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
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
		ftrsmRightUpTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb);
		fgemm( F, FflasNoTrans, FflasTrans, M, Ndown, Nup, Mone, 
		       B, ldb, A+Nup*lda, lda, alpha, B+Nup, ldb);
		ftrsmRightUpTrans( F, Diag, M, Ndown, one, 
				   A+Nup*(lda+1), lda, B+Nup, ldb);
	}
}


template<class Field>
inline void
FFLAS::ftrsmRightLowNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element alpha,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
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
		ftrsmRightLowNoTrans( F, Diag, M, Ndown, alpha, 
				      A+Nup*(lda+1), lda, B+Nup, ldb);
		fgemm( F, FflasNoTrans, FflasNoTrans, M, Nup, Ndown,
		       Mone, B+Nup, ldb, A+Nup*lda, lda, alpha, B, ldb);
		ftrsmRightLowNoTrans( F, Diag, M, Nup, one, A, lda, B, ldb);
	}
}

template<class Field>
inline void
FFLAS::ftrsmRightLowTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::Element alpha,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb){
	
	static typename Field::Element Mone;
	static typename Field::Element one;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
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
		ftrsmRightLowTrans( F, Diag, M, Ndown, alpha, 
				    A+Nup*(lda+1), lda, B+Nup, ldb);
		fgemm( F, FflasNoTrans, FflasTrans, M, Nup, Ndown, Mone, 
		       B+Nup, ldb, A+Nup, lda, alpha, B, ldb);
		ftrsmRightLowTrans( F, Diag, M, Nup, one, A, lda, B, ldb);
	}
}
