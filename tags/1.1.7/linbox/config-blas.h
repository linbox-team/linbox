/* config-blas.h
 * Copyright (C) 2005  Pascal Giorgi
 *               2007  Clement Pernet
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


#ifndef __LINBOX_config_blas_H
#define __LINBOX_config_blas_H

// #ifndef __LINBOX_CONFIGURATION
// #include <linbox-config.h>
// #endif

 

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
	void sgemv_ (const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
	void dger_  (const int*, const int*, const double*, const double*, const int*, const double*, const int*, double*, const int*);

	// level 3 routines
	void dtrsm_ (const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
	void strsm_ (const char*, const char*, const char*, const char*, const int*, const int*, const float*, const float*, const int*, float*, const int*);
	void dtrmm_ (const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
	void strmm_ (const char*, const char*, const char*, const char*, const int*, const int*, const float*, const float*, const int*, float*, const int*);
	void sgemm_ (const char*, const char*, const int*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
	void dgemm_ (const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
}

// define external link to LAPACK routines
extern "C" {
  
#ifdef  __LINBOX_HAVE_DGETRF
        void dgetrf_ (const int *, const int *, double *, const int *, int *, int *);
#endif
#ifdef  __LINBOX_HAVE_DGETRI
        void dgetri_ (const int *, double *, const int *, const int *, double *, const int *, int *);
#endif
#ifdef  __LINBOX_HAVE_DTRTRI
        void dtrtri_ (const char *, const char *, const int *, double *, const int *, int *);
#endif
	void dswap_ (const int *, double *, const int *, double *, const int *);
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
	void cblas_sgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const float alpha, 
			 const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY)
	{
		if (Order == CblasRowMajor)
			sgemv_ ( EXT_BLAS_TRANSPOSE_tr(TransA), &N, &M, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
		else
			sgemv_ ( EXT_BLAS_TRANSPOSE(TransA), &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
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
	void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float *A, const int lda,
			 float *B, const int ldb)
	{  
		if (Order == CblasRowMajor) 
			strsm_ ( EXT_BLAS_SIDE_tr(Side), EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &N, &M, &alpha, A, &lda, B, &ldb); 
		else
			strsm_ ( EXT_BLAS_SIDE(Side), EXT_BLAS_UPLO(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &M, &N, &alpha, A, &lda, B, &ldb);
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
	void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float *A, const int lda,
			 float *B, const int ldb)
	{  
		if (Order == CblasRowMajor)
			strmm_ ( EXT_BLAS_SIDE_tr(Side), EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &N, &M, &alpha, A, &lda, B, &ldb);
		else
			strmm_ ( EXT_BLAS_SIDE(Side), EXT_BLAS_UPLO(Uplo), EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_DIAG(Diag), &M, &N, &alpha, A, &lda, B, &ldb);
	}
  
	void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			 const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb,
			 const double beta, double *C, const int ldc) 
	{   
		if (Order == CblasRowMajor)
			dgemm_ ( EXT_BLAS_TRANSPOSE(TransB), EXT_BLAS_TRANSPOSE(TransA), &N, &M, &K, &alpha, B, &ldb, A, &lda, &beta, C, &ldc);
		else
			dgemm_ ( EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_TRANSPOSE(TransB), &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			 const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb,
			 const float beta, float *C, const int ldc) 
	{   
		if (Order == CblasRowMajor)
			sgemm_ ( EXT_BLAS_TRANSPOSE(TransB), EXT_BLAS_TRANSPOSE(TransA), &N, &M, &K, &alpha, B, &ldb, A, &lda, &beta, C, &ldc);
		else
			sgemm_ ( EXT_BLAS_TRANSPOSE(TransA), EXT_BLAS_TRANSPOSE(TransB), &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}

	// LAPACK routines

#ifdef  __LINBOX_HAVE_DGETRF

	// return A=P.L.U (L unitary) with ColMajor
	// return A=L.U.P (U unitary) with RowMajor
	int clapack_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
			   double *A, const int lda, int *ipiv) 
        {
            int info;
	    dgetrf_ ( &M, &N, A, &lda, ipiv, &info);
            return info;
        }
#endif

#ifdef  __LINBOX_HAVE_DGETRI
	int clapack_dgetri(const enum CBLAS_ORDER Order, const int N, double *A,
			   const int lda, const int *ipiv)
	{
		int info;
		double *work;	

#ifndef __LINBOX_AUTOIMPLEMENT_DGETRI
		// the optimum size of work can be determinted via the 
		// Lapack function ilaenv. 
		work= new double[N];
		dgetri_ (&N, A, &lda, ipiv, work, &N,  &info);
		delete[] work;
#else
		work= new double[N*N];
		dtrtri_("U","N", &N, A, &lda, &info);
		if (info > 0) 
			return 0;
		
		for (int i=0;i<N;++i){
			for(int j=i;j<N;++j){
				work[i*N+j]=A[i*N+j];
				if (j>i) A[i*N+j]=0.0;
			}	      
			work[i*N+i]=1.;
		}
		
		double cst=1.;
		dtrsm_ ("R", "L", "N", "U", &N, &N, &cst, work, &N, A, &N);
		
		int ip;
		const int incr=1;
		for (int i=0; i<N; ++i){
			ip = ipiv[i]-1;
			if (ip != i) 
				dswap_ (&N, &A[i*lda],&incr , &A[ip*lda], &incr);
		}
		
		delete[] work;
#endif
		return info;	
	}
#endif

#ifdef  __LINBOX_HAVE_DTRTRI
	int clapack_dtrtri(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
			   const enum CBLAS_DIAG Diag,const int N, double *A, const int lda)
	{
		int info;
		if (Order == CblasRowMajor)
			dtrtri_ (EXT_BLAS_UPLO_tr(Uplo), EXT_BLAS_DIAG(Diag), &N, A, &lda, &info);
		else
			dtrtri_ (EXT_BLAS_UPLO(Uplo), EXT_BLAS_DIAG(Diag), &N, A, &lda, &info);
		
		return info;		
	}
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

	void cblas_sgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const float alpha, 
			 const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY);
  
	void cblas_dger(const enum CBLAS_ORDER Order, const int M, const int N, const double alpha, const double *X, const int incX,
			const double *Y, const int incY, double *A, const int lda);


	// level 3 routines

	void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
			 double *B, const int ldb);
  
	void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float *A, const int lda,
			 float *B, const int ldb);
  
	void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda,
			 double *B, const int ldb);
  
	void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
			 const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float *A, const int lda,
			 float *B, const int ldb);
  
	void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			 const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb,
			 const double beta, double *C, const int ldc) ;
	void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			 const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb,
			 const float beta, float *C, const int ldc) ;

	// LAPACK routines

	int clapack_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
			   double *A, const int lda, int *ipiv);
	int clapack_dgetri(const enum CBLAS_ORDER Order, const int N, double *A,
			   const int lda, const int *ipiv);
	int clapack_dtrtri(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
			   const enum CBLAS_DIAG Diag,const int N, double *A, const int lda);
}
#endif

#endif //__LINBOX_config_blas_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
