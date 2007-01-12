/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/lifting-container-base.h
 * Copyright (C) 2005  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __LINBOX_CONFIG_BLAS
#define __LINBOX_CONFIG_BLAS

#ifndef __LINBOX_CONFIGURATION
#include <linbox-config.h>
#endif

 

//#define __LINBOX_HAVE_CBLAS


#ifndef __LINBOX_HAVE_CBLAS

// CBLAS are not available define our own wrapper


// define external link to BLAS function
extern "C" {
  
	// level 1 routines
	void   daxpy_   (const int*, const double*, const double*, const int*, double*, const int*);
	double ddot_    (const int*, const double*, const int*, const double*, const int*);
	double dasum_   (const int*, const double*, const int*);
	int    idamax_  (const int*, const double*, const int*);
	double dnrm2_   (const int*, const double*, const int*);

	// level 2 routines
	void dgemv_ (const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
	void dger_  (const int*, const int*, const double*, const double*, const int*, const double*, const int*, double*, const int*);

	// level 3 routines
	void dtrsm_ (const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
	void dtrmm_ (const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
	void dgemm_ (const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
}

// define external link to LAPACK routines
extern "C" {
  
#if HAVE_DGETRF
        void dgetrf_(const int *, const int *, double *, const int *, int *, int *);
#endif
#if HAVE_DGETRI
        void dgetri_(const int *, double *, const int *, const int *, double *, const int *, int *);
#endif
}

// define C wrappers
extern "C" {

#define CBLAS_ENUM_DEFINED_H
	enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
	enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, AtlasConj=114};
	enum CBLAS_UPLO  {CblasUpper=121, CblasLower=122};
	enum CBLAS_DIAG  {CblasNonUnit=131, CblasUnit=132};
	enum CBLAS_SIDE  {CblasLeft=141, CblasRight=142};
	
	static const char* EXT_BLAS_TRANSPOSE    (CBLAS_TRANSPOSE t) { if (t == CblasNoTrans) return "N"; else if (t == CblasTrans) return "T"; else return "";}
	static const char* EXT_BLAS_TRANSPOSE_tr (CBLAS_TRANSPOSE t) { if (t == CblasNoTrans) return "T"; else if (t == CblasTrans) return "N"; else return "";}
	
	static const char* EXT_BLAS_UPLO         (CBLAS_UPLO t)      { if (t == CblasUpper) return "U"; else return "L";}
	static const char* EXT_BLAS_UPLO_tr      (CBLAS_UPLO t)      { if (t == CblasUpper) return "L"; else return "U";}
	
	static const char* EXT_BLAS_DIAG         (CBLAS_DIAG t)      { if (t == CblasUnit)  return "U"; else return "N";}
	
	static const char* EXT_BLAS_SIDE         (CBLAS_SIDE t)      { if (t == CblasLeft)  return "L"; else return "R";}
	static const char* EXT_BLAS_SIDE_tr      (CBLAS_SIDE t)      { if (t == CblasLeft)  return "R"; else return "L";}
	
#define CBLAS_INDEX int

	// level 1 routines

	void cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY)
	{
		daxpy_ (&N,&alpha, X, &incX, Y, &incY);
	}
  
	double cblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY)
	{
		return ddot_ (&N, X, &incX, Y, &incY);
	}
  
	double cblas_dasum(const int N, const double *X, const int incX){
		return dasum_ (&N, X, &incX);
	}
  
	int cblas_idamax(const int N, const double *X, const int incX){
		return idamax_ (&N, X, &incX);
	}

	double cblas_dnrm2(const int N, const double *X, const int incX){
		return dnrm2_(&N, X, &incX);
	}


	// level 2 routines

	void cblas_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, 
			 const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY)
	{
		if (Order == CblasRowMajor)
			dgemv_ ( EXT_BLAS_TRANSPOSE_tr(TransA), &N, &M, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
		else
			dgemv_ ( EXT_BLAS_TRANSPOSE(TransA), &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
	}
  
	void cblas_dger(const enum CBLAS_ORDER Order, const int M, const int N, const double alpha, const double *X, const int incX,
			const double *Y, const int incY, double *A, const int lda)
	{  
		if (Order == CblasRowMajor)
			dger_ (&N, &M, &alpha, Y, &incY, X, &incX, A, &lda);
		else
			dger_ (&M, &N, &alpha, X, &incX, Y, &incY, A, &lda);
	}



	// level 3 routines

	void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
			 double *B, const int ldb)
	{  
		if (Order == CblasRowMajor) 
			dtrsm_ ( EXT_BLAS_SIDE_tr(Side), EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &N, &M, &alpha, A, &lda, B, &ldb); 
		else
			dtrsm_ ( EXT_BLAS_SIDE(Side), EXT_BLAS_UPLO(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &M, &N, &alpha, A, &lda, B, &ldb);
	}
  
	void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
			 double *B, const int ldb)
	{  
		if (Order == CblasRowMajor)
			dtrmm_ ( EXT_BLAS_SIDE_tr(Side), EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &N, &M, &alpha, A, &lda, B, &ldb);
		else
			dtrmm_ ( EXT_BLAS_SIDE(Side), EXT_BLAS_UPLO(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &M, &N, &alpha, A, &lda, B, &ldb);
	}
  
	void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			 const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb,
			 const double beta, double *C, const int ldc) 
	{   
		if (Order == CblasRowMajor)
			dgemm_ ( EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_TRANSPOSE(TransB), &N, &M, &K, &alpha, B, &ldb, A, &lda, &beta, C, &ldc);
		else
			dgemm_ ( EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_TRANSPOSE(TransB), &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}

	// LAPACK routines

#if HAVE_DGETRF
	int clapack_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
			   double *A, const int lda, int *ipiv) 
        {
            int info;
            dgetrf_ ( &M, &N, A, &lda, ipiv, &info);
            return info;
        }
#endif

#if HAVE_DGETRI
	int clapack_dgetri(const enum CBLAS_ORDER Order, const int N, double *A,
			   const int lda, const int *ipiv);
#endif

} 

#else
extern "C" {
	

#define CBLAS_ENUM_DEFINED_H
	enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
	enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
			      AtlasConj=114};
	enum CBLAS_UPLO  {CblasUpper=121, CblasLower=122};
	enum CBLAS_DIAG  {CblasNonUnit=131, CblasUnit=132};
	enum CBLAS_SIDE  {CblasLeft=141, CblasRight=142};


#define CBLAS_INDEX int	
	

	int cblas_errprn(int ierr, int info, char *form, ...);

	// level 1 routines

	void   cblas_daxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);

	double cblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY);
  
	double cblas_dasum(const int N, const double *X, const int incX);
  
	int    cblas_idamax(const int N, const double *X, const int incX);

	double cblas_dnrm2(const int N, const double *X, const int incX);


	// level 2 routines

	void cblas_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, 
			 const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY);
  
	void cblas_dger(const enum CBLAS_ORDER Order, const int M, const int N, const double alpha, const double *X, const int incX,
			const double *Y, const int incY, double *A, const int lda);


	// level 3 routines

	void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
			 double *B, const int ldb);
  
	void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
			 double *B, const int ldb);
  
	void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			 const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb,
			 const double beta, double *C, const int ldc) ;

	// LAPACK routines

	int clapack_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
			   double *A, const int lda, int *ipiv);
	int clapack_dgetri(const enum CBLAS_ORDER Order, const int N, double *A,
			   const int lda, const int *ipiv);


}
#endif

#endif
