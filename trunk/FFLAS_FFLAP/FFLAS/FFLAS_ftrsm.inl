/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//----------------------------------------------------------------------
//                FFLAS: Finite Field Fast Linear Algebra Subroutines 
//                FFLAS_trsm: Triangular Matrix Solve with Matrix
//----------------------------------------------------------------------
// by Clement PERNET ( clement.pernet@imag.fr )
// 25/04/2003
//----------------------------------------------------------------------


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
	
	static typename Field::element Mone;
	F.neg(Mone, F.one);
	if ( Side==FflasLeft ){
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
			if ( Uplo==FflasUpper){
				ftrsm( F, Side, Uplo, TransA, Diag, Mdown, N, alpha, 
				       A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
				if (TransA == FflasNoTrans)
					fgemm( F, TransA, FflasNoTrans, Mup, N, 
					       Mdown, Mone, A+Mup, lda, B+Mup*ldb, ldb,
					       alpha, B, ldb);
				else  // FflasTrans
					fgemm( F, TransA, FflasNoTrans, Mup, N, 
					       Mdown, Mone, A+Mup*lda, lda,
					       B+Mup*ldb, ldb, alpha, B, ldb);
				ftrsm( F, Side, Uplo, TransA, Diag, Mup, N, F.one, 
				       A, lda, B, ldb);
			}
			else{ // FflasLower
				ftrsm( F, Side, Uplo, TransA, Diag, Mup, N, alpha, 
				       A, lda, B, ldb);
				if (TransA == FflasNoTrans)
					fgemm( F, TransA, FflasNoTrans, Mdown, N, Mup,
					       Mone, A+Mup*lda, lda, B, ldb,
					       alpha, B+Mup*ldb, ldb);
				else  // FflasTrans
					fgemm( F, TransA, FflasNoTrans, Mdown, N, Mup,
					       Mone, A+Mup, lda, B, ldb,
					       alpha, B+Mup*ldb, ldb);
				ftrsm( F, Side, Uplo, TransA, Diag, Mdown, N, F.one, 
				       A+Mup*(lda+1), lda, B+Mup*ldb, ldb);
			}
		}
	}
	else{ // FflasRight
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
			if ( Uplo==FflasUpper){ 
				ftrsm( F, Side, Uplo, TransA, Diag, M, Nup, alpha, 
				       A, lda, B, ldb);
				if (TransA == FflasNoTrans)
					fgemm( F, FflasNoTrans, TransA, M, Ndown, Nup,
					       Mone, B, ldb, A+Nup, lda,
					       alpha, B+Nup, ldb);
				else  // FflasTrans
					fgemm( F, FflasNoTrans, TransA, M, Ndown, Nup,
					       Mone, B, ldb, A+Nup*lda, lda,
					       alpha, B+Nup, ldb);
				ftrsm( F, Side, Uplo, TransA, Diag, M, Ndown, F.one, 
				       A+Nup*(lda+1), lda, B+Nup, ldb);
			}
			else{ // FflasLower
				ftrsm( F, Side, Uplo, TransA, Diag, M, Ndown, alpha, 
				       A+Nup*(lda+1), lda, B+Nup, ldb);
				if (TransA == FflasNoTrans)
					fgemm( F, FflasNoTrans, TransA, M, Nup, Ndown,
					       Mone, B+Nup, ldb, A+Nup*lda, lda,
					       alpha, B, ldb);
				else  // FflasTrans
					fgemm( F, FflasNoTrans, TransA, M, Nup, Ndown,
					       Mone, B+Nup, ldb, A+Nup, lda,
					       alpha, B, ldb);
				ftrsm( F, Side, Uplo, TransA, Diag, M, Nup, F.one, 
				       A, lda, B, ldb);
			}
		}
	}
}

