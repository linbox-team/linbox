/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflas/fflas_ftrsm.inl
 * Copyright (C) 2003 Clement Pernet
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
FFLAS::ftrsm(const Field& F, const enum FFLAS_SIDE Side,
	     const enum FFLAS_UPLO Uplo, 
	     const enum FFLAS_TRANSPOSE TransA,
	     const enum FFLAS_DIAG Diag, 
	     const size_t M, const size_t N,
	     const typename Field::element alpha,
	     const typename Field::element * A, size_t lda,
	     typename Field::element * B, size_t ldb){
	
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
				ftrsmLeftLowNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
			else{
				ftrsmLeftLowTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
			}
		}
	}
	else{
	if ( Uplo==FflasUpper){
			if (TransA == FflasNoTrans){
				ftrsmRightUpNoTrans(F,Diag,M,N,alpha,A,lda,B,ldb);
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
			  const typename Field::element alpha,
			  const typename Field::element * A, size_t lda,
			  typename Field::element * B, size_t ldb){
	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( M==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::element inv;
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
		ftrsmLeftUpNoTrans( F, Diag, Mup, N, F.one, A, lda, B, ldb);
	}
}

template<class Field>
inline void
FFLAS::ftrsmLeftUpTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		      const size_t M, const size_t N,
		      const typename Field::element alpha,
		      const typename Field::element * A, size_t lda,
		      typename Field::element * B, size_t ldb){

	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( M==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::element inv;
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
		ftrsmLeftUpTrans( F, Diag, Mup, N, F.one, A, lda, B, ldb);
	}
}

template<class Field>
inline void
FFLAS::ftrsmLeftLowNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::element alpha,
			 const typename Field::element * A, size_t lda,
			 typename Field::element * B, size_t ldb){

	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( M==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::element inv;
			F.inv(inv, *A);
			fscal(F, N, inv, B, 1);
		}
	}
	else{
		size_t Mup=M>>1;
		size_t Mdown = M-Mup;
		ftrsmLeftLowNoTrans( F, Diag, Mup, N, alpha, A, lda, B, ldb);
		fgemm( F, FflasNoTrans, FflasNoTrans, Mdown, N, Mup,
		       Mone, A+Mup*lda, lda, B, ldb, alpha, B+Mup*ldb, ldb);
		ftrsmLeftLowNoTrans( F, Diag, Mdown, N, F.one, 
				     A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
	}
}

template<class Field>
inline void 
FFLAS::ftrsmLeftLowTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		       const size_t M, const size_t N,
		       const typename Field::element alpha,
		       const typename Field::element * A, size_t lda,
		       typename Field::element * B, size_t ldb){

	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( M==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::element inv;
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
		ftrsmLeftLowTrans( F, Diag, Mdown, N, F.one, 
				   A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
	}
}

template<class Field>
inline void 
FFLAS::ftrsmRightUpNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			 const size_t M, const size_t N,
			 const typename Field::element alpha,
			 const typename Field::element * A, size_t lda,
			 typename Field::element * B, size_t ldb){
	
	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( N==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::element inv;
			F.inv(inv, *A);
			fscal(F, M, inv, B, ldb);
		}
	}
	else{
		size_t Nup=N>>1;
		size_t Ndown = N-Nup;
		ftrsmRightUpNoTrans( F, Diag, M, Nup, alpha, A, lda, B, ldb);
		fgemm( F, FflasNoTrans, FflasNoTrans, M, Ndown, Nup,
		       Mone, B, ldb, A+Nup, lda, alpha, B+Nup, ldb);
		ftrsmRightUpNoTrans( F, Diag, M, Ndown, F.one, 
				     A+Nup*(lda+1), lda, B+Nup, ldb);
	}
}



template<class Field>
inline void
FFLAS::ftrsmRightUpTrans(const Field& F, const enum FFLAS_DIAG Diag, 
		       const size_t M, const size_t N,
		       const typename Field::element alpha,
		       const typename Field::element * A, size_t lda,
		       typename Field::element * B, size_t ldb){
	
	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( N==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::element inv;
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
		ftrsmRightUpTrans( F, Diag, M, Ndown, F.one, 
				   A+Nup*(lda+1), lda, B+Nup, ldb);
	}
}


template<class Field>
inline void
FFLAS::ftrsmRightLowNoTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::element alpha,
			  const typename Field::element * A, size_t lda,
			  typename Field::element * B, size_t ldb){
	
	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( N==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::element inv;
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
		ftrsmRightLowNoTrans( F, Diag, M, Nup, F.one, A, lda, B, ldb);
	}
}

template<class Field>
inline void
FFLAS::ftrsmRightLowTrans(const Field& F, const enum FFLAS_DIAG Diag, 
			  const size_t M, const size_t N,
			  const typename Field::element alpha,
			  const typename Field::element * A, size_t lda,
			  typename Field::element * B, size_t ldb){
	
	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( N==1 ){
		if (Diag == FflasNonUnit ){
			typename Field::element inv;
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
		ftrsmRightLowTrans( F, Diag, M, Nup, F.one, A, lda, B, ldb);
	}
}


