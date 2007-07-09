/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* fflas.h
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */
#include <math.h>

#ifndef __FFLAS_H
#define __FFLAS_H
#include <math.h>

#ifndef MAX
#define MAX(a,b) ((a < b)?b:a)
#endif
#ifndef MIN
#define MIN(a,b) ((a > b)?b:a)
#endif

#ifdef _LINBOX_LINBOX_CONFIG_H 
#include "linbox/config-blas.h"
#include "linbox/field/unparametric.h"

namespace LinBox {
#else
#include "config-blas.h"
#include "unparametric.h"
#endif
	
#ifndef __LINBOX_STRASSEN_OPTIMIZATION
#define WINOTHRESHOLD 400
#else
#define WINOTHRESHOLD __LINBOX_WINOTHRESHOLD
#endif

// Thresholds determining which floating point representation to use,
// depending on the cardinality of the finite field. This is only used when
// the element representation is not a floating point type.
#define FLOAT_DOUBLE_THRESHOLD_0 430
#define FLOAT_DOUBLE_THRESHOLD_1 350
#define FLOAT_DOUBLE_THRESHOLD_2 175
	
#define DOUBLE_MANTISSA 53
#define FLOAT_MANTISSA 24
	
class FFLAS {

public:	
	enum FFLAS_TRANSPOSE { FflasNoTrans=111, FflasTrans=112};
	enum FFLAS_UPLO      { FflasUpper=121, FflasLower=122 };
	enum FFLAS_DIAG      { FflasNonUnit=131, FflasUnit=132 };
	enum FFLAS_SIDE      { FflasLeft=141, FflasRight = 142 };

	/* Determine the type of the element representation for Matrix Mult kernel
	 * FflasDouble: to use the double precision BLAS
	 * FflasFloat: to use the single precison BLAS
	 * FflasFloat: for any other domain, that can not be converted to floating point integers
	 */
	enum FFLAS_BASE      { FflasDouble = 151, FflasFloat = 152, FflasGeneric = 153};

	/* Representations of Z with floating point elements*/
	typedef UnparametricField<float> FloatDomain;
	typedef UnparametricField<double> DoubleDomain;


	
	
//---------------------------------------------------------------------
// Level 1 routines
//---------------------------------------------------------------------
	//---------------------------------------------------------------------
	// fscal: X <- alpha.X
	// X is a vector of size n
	//---------------------------------------------------------------------
	template<class Field>
	static void
	fscal (const Field& F, const size_t n, const typename Field::Element alpha, 
	       typename Field::Element * X, const size_t incX){
		
		typename Field::Element * Xi = X;
		for (; Xi < X+n*incX; Xi+=incX )
			F.mulin( *Xi, alpha );
	}
	//---------------------------------------------------------------------
	// fcopy: x <- y
	// x,y are vectors of size N
	//---------------------------------------------------------------------
	template<class Field>
	static void
	fcopy (const Field& F, const size_t N, 
	       typename Field::Element * X, const size_t incX,
	       const typename Field::Element * Y, const size_t incY );

	//---------------------------------------------------------------------
	// faxpy: y <- a.x + y
	// x,y are vectors of size N
	//---------------------------------------------------------------------
	template<class Field>
	static void
	faxpy (const Field& F, const size_t N, 
	       const typename Field::Element a,
	       const typename Field::Element * X, const size_t incX,
	       typename Field::Element * Y, const size_t incY );

	//---------------------------------------------------------------------
	// fdot: returns x^T . y
	// x and y are vectors of size N
	//---------------------------------------------------------------------
	template<class Field>
	static typename Field::Element
	fdot (const Field& F, const size_t N, 
	      const typename Field::Element * X, const size_t incX,
	      const typename Field::Element * Y, const size_t incY );

	//---------------------------------------------------------------------
	// fswap: X <-> Y
	// X,Y are vectors of size N
	//---------------------------------------------------------------------
	template<class Field>
	static void
	fswap (const Field& F, const size_t N, typename Field::Element * X, const size_t incX,
	       typename Field::Element * Y, const size_t incY ){
		
		typename Field::Element tmp;
		typename Field::Element * Xi = X;
		typename Field::Element * Yi=Y;
		for (; Xi < X+N*incX; Xi+=incX, Yi+=incY ){
			F.assign( tmp, *Xi );
			F.assign( *Xi, *Yi );
			F.assign( *Yi, tmp );
		}
	}
		
//---------------------------------------------------------------------
// Level 2 routines
//---------------------------------------------------------------------
	/**
	 *  @brief finite prime Field GEneral Matrix Vector multiplication
	 *
	 *  Computes  Y <- alpha op(A).X + beta.Y \\
	 *  A is m*n
	 */
	template<class Field>
	static void
	fgemv (const Field& F, const FFLAS_TRANSPOSE TransA, 
	       const size_t M, const size_t N,
	       const typename Field::Element alpha, 
	       const typename Field::Element * A, const size_t lda,
	       const typename Field::Element * X, const size_t incX, 
	       const  typename Field::Element beta,
	       typename Field::Element * Y, const size_t incY);

	/**
	 *  @brief fger: GEneral ?
	 *
	 *  Computes  A <- alpha x . y^T + A \\
	 *  A is m*n, x and y are vectors of size m and n
	 */
	template<class Field>
	static void
	fger (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element alpha, 
	      const typename Field::Element * x, const size_t incx,
	      const typename Field::Element * y, const size_t incy, 
	      typename Field::Element * A, const size_t lda);

	/**
	   @brief ftrsv: TRiangular System solve with Vector
	   Computes  X <- op(A^-1).X\\
	   size of X is N
	*/
	template<class Field>
	static void
	ftrsv (const Field& F, const FFLAS_UPLO Uplo, 
	       const FFLAS_TRANSPOSE TransA, const FFLAS_DIAG Diag,
	       const size_t N,const typename Field::Element * A, const size_t lda,
	       typename Field::Element * X, int incX);
	
//---------------------------------------------------------------------
// Level 3 routines
//---------------------------------------------------------------------

	//---------------------------------------------------------------------
	// ftrsm: TRiangular System solve with matrix
	// Computes  B <- alpha.op(A^-1).B,  B <- alpha.B.op(A^-1)
	// B is m*n
	//---------------------------------------------------------------------
	template<class Field>
	static void
	ftrsm (const Field& F, const FFLAS_SIDE Side,
	       const FFLAS_UPLO Uplo, 
	       const FFLAS_TRANSPOSE TransA,
	       const FFLAS_DIAG Diag, 
	       const size_t M, const size_t N,
	       const typename Field::Element alpha,
	       typename Field::Element * A, const size_t lda,
	       typename Field::Element * B, const size_t ldb);
	
	//---------------------------------------------------------------------
	// ftrmm: TRiangular Matrix Multiply
	// Computes  B <- alpha.op(A).B,  B <- alpha.B.op(A)
	// B is m*n
	//---------------------------------------------------------------------
	template<class Field>
	static void
	ftrmm (const Field& F, const FFLAS_SIDE Side,
	       const FFLAS_UPLO Uplo, 
	       const FFLAS_TRANSPOSE TransA,
	       const FFLAS_DIAG Diag, 
	       const size_t M, const size_t N,
	       const typename Field::Element alpha,
	       typename Field::Element * A, const size_t lda,
	       typename Field::Element * B, const size_t ldb);
	
	/** @brief  Field GEneral Matrix Multiply 
	 * 
	 * Computes C = alpha.op(A)*op(B) + beta.C ,
	 * op(A) = A, A<sup>T</sup>
	 * wl recursive levels of Winograd's algorithm are used 
	 */
	template<class Field>
	static typename Field::Element* 
	fgemm( const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m,
	       const size_t n,
	       const size_t k,
	       const typename Field::Element alpha,
	       const typename Field::Element* A, const size_t lda,
	       const typename Field::Element* B, const size_t ldb, 
	       const typename Field::Element beta,
	       typename Field::Element* C, const size_t ldc,
	       const size_t w){

		if (!(m && n && k)) return C;
		
		FFLAS_BASE base = BaseCompute<typename Field::Element> ()(F, w);

		WinoMain (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta,
				 C, ldc, DotProdBound (F, w, beta, base), w, base);
		return C;
		};
	
	/** @brief  Field GEneral Matrix Multiply 
	 * 
	 * Computes C = alpha.op(A)*op(B) + beta.C ,
	 * op(A) = A, A<sup>T</sup>
	 * Automitically set Winograd recursion level
	 */
	template<class Field>
	static typename Field::Element*
	fgemm (const Field& F,
	       const FFLAS_TRANSPOSE ta,
	       const FFLAS_TRANSPOSE tb,
	       const size_t m,
	       const size_t n,
	       const size_t k,
	       const typename Field::Element alpha,
	       const typename Field::Element* A, 
	       const size_t lda,
	       const typename Field::Element* B,
	       const size_t ldb, 
	       const typename Field::Element beta,
	       typename Field::Element* C, 
	       const size_t ldc){

		if (!(m && n && k)) return C;

		size_t w, kmax=0;
 		FFLAS_BASE base;

		setMatMulParam<typename Field::Element> ()(F, MIN(MIN(m,n),k), beta,
							   w, base, kmax);

		WinoMain (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta,
			  C, ldc, kmax, w, base);
		return C;
	}
	
	//---------------------------------------------------------------------
	// fsquare: 
	// compute C = alpha. op(A)*op(A) + beta.C over a Field
	// op(A) =A, A^T
	// Avoid the conversion of B 
	//---------------------------------------------------------------------
	template<class Field>
	static typename Field::Element* fsquare (const Field& F,
						 const FFLAS_TRANSPOSE ta,
						 const size_t n,
						 const typename Field::Element alpha,
						 const typename Field::Element* A, 
						 const size_t lda,
						 const typename Field::Element beta,
						 typename Field::Element* C, 
						 const size_t ldc);
protected:

	// Prevents the instantiation of the class
	FFLAS(){}
	template <class X,class Y>
	class AreEqual
	{
	public:
		static const bool value = false;
	};
	
	template <class X>
	class AreEqual<X,X>
	{
	public:
		static const bool value = true;
	};

	//-----------------------------------------------------------------------------
	// Some conversion functions
	//-----------------------------------------------------------------------------
	
	//---------------------------------------------------------------------
	// Finite Field matrix => double matrix
	//---------------------------------------------------------------------
	template<class Field>
	static void MatF2MatD (const Field& F,
			       DoubleDomain::Element* S, const size_t lds,
			       const typename Field::Element* E,
			       const size_t lde,const size_t m, const size_t n){
		
		const typename Field::Element* Ei = E;
		DoubleDomain::Element *Si=S;
		size_t j; 
		for (; Ei < E+lde*m; Ei+=lde, Si += lds)
			for ( j=0; j<n; ++j){
				F.convert(*(Si+j),*(Ei+j));
			}
		}
	//---------------------------------------------------------------------
	// Finite Field matrix => float matrix
	//---------------------------------------------------------------------
	template<class Field>
	static void MatF2MatFl (const Field& F,
				FloatDomain::Element* S, const size_t lds,
				const typename Field::Element* E,
				const size_t lde,const size_t m, const size_t n){
		
		const typename Field::Element* Ei = E;
		FloatDomain::Element *Si=S;
		size_t j; 
		for (; Ei < E+lde*m; Ei+=lde, Si += lds)
			for ( j=0; j<n; ++j){
				F.convert(*(Si+j),*(Ei+j));
			}
		}
	
	//---------------------------------------------------------------------
	// Finite Field matrix => double matrix
	// Special design for upper-triangular matrices
	//---------------------------------------------------------------------
	template<class Field>
	static void MatF2MatD_Triangular (const Field& F,
					  typename DoubleDomain::Element* S, const size_t lds,
					  const typename Field::Element* const E,
					  const size_t lde,
					  const size_t m, const size_t n){
		
		const typename Field::Element* Ei = E;
		typename DoubleDomain::Element* Si = S;
		size_t i=0, j;
		for ( ; i<m;++i, Ei+=lde, Si+=lds)
			for ( j=i; j<n;++j)
				F.convert(*(Si+j),*(Ei+j));
	}

	//---------------------------------------------------------------------
	// Finite Field matrix => float matrix
	// Special design for upper-triangular matrices
	//---------------------------------------------------------------------
	template<class Field>
	static void MatF2MatFl_Triangular (const Field& F,
					   typename FloatDomain::Element* S, const size_t lds,
					   const typename Field::Element* const E,
					   const size_t lde,
					   const size_t m, const size_t n){
		
		const typename Field::Element* Ei = E;
		typename FloatDomain::Element* Si = S;
		size_t i=0, j;
		for ( ; i<m;++i, Ei+=lde, Si+=lds)
			for ( j=i; j<n;++j)
				F.convert(*(Si+j),*(Ei+j));
	}
	
	//---------------------------------------------------------------------
	// double matrix => Finite Field matrix
	//---------------------------------------------------------------------
	template<class Field>
	static void MatD2MatF (const Field& F,
			       typename Field::Element* S, const size_t lds,
			       const typename DoubleDomain::Element* E, const size_t lde,
			       const size_t m, const size_t n){
		
		typename Field::Element* Si = S;
		const DoubleDomain::Element* Ei =E;
		size_t j;
		for ( ; Si < S+m*lds; Si += lds, Ei+= lde){
			for ( j=0; j<n;++j)
				F.init( *(Si+j), *(Ei+j) );
		}
	}

	//---------------------------------------------------------------------
	// float matrix => Finite Field matrix
	//---------------------------------------------------------------------
	template<class Field>
	static void MatFl2MatF (const Field& F,
				typename Field::Element* S, const size_t lds,
				const typename FloatDomain::Element* E, const size_t lde,
				const size_t m, const size_t n){
		
		typename Field::Element* Si = S;
		const FloatDomain::Element* Ei =E;
		size_t j;
		for ( ; Si < S+m*lds; Si += lds, Ei+= lde){
			for ( j=0; j<n;++j)
				F.init( *(Si+j), *(Ei+j) );
		}
	}

	/** @brief Bound for the delayed modulus matrix multiplication
	 *
	 *  @param F  - the finite field
	 *  @param w  - the number of recursive levels of Winograds algorithm being used
	 *  @param beta - for the computation of C <- AB + beta C
	 *
	 *  Compute the maximal dimension k, such that a matrix multiplication (m,k)x(k,n)
	 *  can be performed over Z without overflow of the 53 bits of the double
	 *  mantissa.
	 *  See [Dumas, Gautier, Pernet ISSAC'2002]
	 */
	template <class Field>
	static size_t DotProdBound (const Field& F, const size_t w,
				    const typename Field::Element& beta,
				    const FFLAS_BASE base);

	template <class Field>
	static size_t DotProdBoundCompute (const Field& F, const size_t w,
					   const typename Field::Element& beta,
					   const FFLAS_BASE base);
	

	static size_t WinoSteps (const size_t m);
	
	// 	template <class Element>
// 	class callDotProdBoundCompute;

	/** @brief Bound for the delayed modulus triangular system solving
	 *
	 *  @param F  - the finite field
	 *
	 *  Compute the maximal dimension k, such that a unit diagonal triangular
	 *  system of dimension k can be solved over Z without overflow of the
	 *  53 bits of the double mantissa.
  	 *  See [Dumas, Giorgi, Pernet ISSAC'2004]
	 */
	template <class Field>
	static size_t TRSMBound (const Field& F);

	template <class Element>
	class callTRSMBound;

	/** @brief Set the optimal parameters for the Matrix Multiplication
	 */
	template <class Element>
	class setMatMulParam;

	template <class Element>
	class BaseCompute;

	template <class Field>
	static void DynamicPealing( const Field& F, 
				    const FFLAS_TRANSPOSE ta,
				    const FFLAS_TRANSPOSE tb,
				    const size_t m, const size_t n, const size_t k,
				    const typename Field::Element alpha, 
				    const typename Field::Element* A, const size_t lda,
				    const typename Field::Element* B, const size_t ldb, 
				    const typename Field::Element beta,
				    typename Field::Element* C, const size_t ldc, 
				    const size_t kmax );

	template<class Field>
	static void MatVectProd (const Field& F, 
				 const FFLAS_TRANSPOSE TransA, 
				 const size_t M, const size_t N,
				 const typename Field::Element alpha, 
				 const typename Field::Element * A, const size_t lda,
				 const typename Field::Element * X, const size_t incX, 
				 const typename Field::Element beta,
				 typename Field::Element * Y, const size_t incY);

	template <class Field>
	static void ClassicMatmul(const Field& F,  
				  const FFLAS_TRANSPOSE ta,
				  const FFLAS_TRANSPOSE tb,
				  const size_t m, const size_t n, const size_t k,
				  const typename Field::Element alpha,
				  const typename Field::Element * A, const size_t lda,
				  const typename Field::Element * B, const size_t ldb,
				  const typename Field::Element beta,
				  typename Field::Element * C, const size_t ldc, 
				  const size_t kmax, const FFLAS_BASE base );
    
	// Winograd Multiplication  alpha.A(n*k) * B(k*m) + beta . C(n*m)
	// WinoCalc performs the 22 Winograd operations
	template<class Field>
	static void WinoCalc (const Field& F, 
			      const FFLAS_TRANSPOSE ta,
			      const FFLAS_TRANSPOSE tb,
			      const size_t mr, const size_t nr,const size_t kr,
			      const typename Field::Element alpha,
			      const typename Field::Element* A,const size_t lda,
			      const typename Field::Element* B,const size_t ldb,
			      const typename Field::Element beta,
			      typename Field::Element * C, const size_t ldc,
			      const size_t kmax, const size_t w, const FFLAS_BASE base);
	
	template<class Field>
	static void WinoMain (const Field& F, 
			      const FFLAS_TRANSPOSE ta,
			      const FFLAS_TRANSPOSE tb,
			      const size_t m, const size_t n, const size_t k,
			      const typename Field::Element alpha,
			      const typename Field::Element* A,const size_t lda,
			      const typename Field::Element* B,const size_t ldb,
			      const typename Field::Element beta,
			      typename Field::Element * C, const size_t ldc,
			      const size_t kmax, const size_t w, const FFLAS_BASE base);


	template<class Element>
	class callWinoMain;

	template<class Element>
	class callClassicMatmul;

	template<class Element>
	class callFsquare;

	template<class Element>
	class callMatVectProd;
		
	// Specialized routines for ftrsm
	template<class Field>
	static void ftrsmLeftUpNoTrans (const Field& F, const FFLAS_DIAG Diag, 
					const size_t M, const size_t N,
					const typename Field::Element alpha,
					const typename Field::Element * A, const size_t lda,
					typename Field::Element * B, const size_t ldb);
	
	template<class Field>
	static void ftrsmLeftUpTrans (const Field& F, const FFLAS_DIAG Diag, 
				      const size_t M, const size_t N,
				      const typename Field::Element alpha,
				      const typename Field::Element * A, const size_t lda,
				      typename Field::Element * B, const size_t ldb);
	
	template<class Field>
	static void ftrsmLeftLowNoTrans (const Field& F, const FFLAS_DIAG Diag, 
					 const size_t M, const size_t N,
					 const typename Field::Element alpha,
					 typename Field::Element * A, const size_t lda,
					 typename Field::Element * B, const size_t ldb, 
					 const size_t nmax);
	template<class Element>
	class callFtrsmLeftLowNoTrans;
	
	template<class Element>
	class callFtrsmRightUpNoTrans;
	
	template<class Element>
	class callFtrmmLeftUpNoTrans;
	
	template<class Element>
	class callFtrmmLeftUpTrans;

	template<class Element>
	class callFtrmmLeftLowNoTrans;		

	template<class Element>
	class callFtrmmLeftLowTrans;		

	template<class Element>
	class callFtrmmRightUpNoTrans;		

	template<class Element>
	class callFtrmmRightUpTrans;		

	template<class Element>
	class callFtrmmRightLowNoTrans;		

	template<class Element>
	class callFtrmmRightLowTrans;		
	
	template<class Field>
	static void ftrsmLeftLowTrans (const Field& F, const FFLAS_DIAG Diag, 
				       const size_t M, const size_t N,
				       const typename Field::Element alpha,
				       const typename Field::Element * A, const size_t lda,
				       typename Field::Element * B, const size_t ldb);
	
	template<class Field>
	static void ftrsmRightUpNoTrans (const Field& F, const FFLAS_DIAG Diag, 
					 const size_t M, const size_t N,
					 const typename Field::Element alpha,
					 typename Field::Element * A, const size_t lda,
					 typename Field::Element * B, const size_t ldb,
					 const size_t nmax);
	
	template<class Field>
	static void ftrsmRightUpTrans (const Field& F, const FFLAS_DIAG Diag, 
				       const size_t M, const size_t N,
				       const typename Field::Element alpha,
				       const typename Field::Element * A, const size_t lda,
				       typename Field::Element * B, const size_t ldb);

	template<class Field>
	static void ftrsmRightLowNoTrans (const Field& F, const FFLAS_DIAG Diag, 
					  const size_t M, const size_t N,
					  const typename Field::Element alpha,
					  const typename Field::Element * A, const size_t lda,
					  typename Field::Element * B, const size_t ldb);

	template<class Field>
	static void ftrsmRightLowTrans (const Field& F, const FFLAS_DIAG Diag, 
					const size_t M, const size_t N,
					const typename Field::Element alpha,
					const typename Field::Element * A, const size_t lda,
					typename Field::Element * B, const size_t ldb);

	// Specialized routines for ftrmm
	template<class Field>
	static void ftrmmLeftUpNoTrans (const Field& F, const FFLAS_DIAG Diag, 
					const size_t M, const size_t N,
					const typename Field::Element * A, const size_t lda,
					typename Field::Element * B, const size_t ldb, 
					const size_t nmax);

	template<class Field>
	static void ftrmmLeftUpTrans (const Field& F, const FFLAS_DIAG Diag, 
				      const size_t M, const size_t N,
				      const typename Field::Element * A, const size_t lda,
				      typename Field::Element * B, const size_t ldb,
				      const size_t nmax);

	template<class Field>
	static void ftrmmLeftLowNoTrans (const Field& F, const FFLAS_DIAG Diag, 
					 const size_t M, const size_t N,
					 const typename Field::Element * A, const size_t lda,
					 typename Field::Element * B, const size_t ldb, 
					 const size_t nmax);

	template<class Field>
	static void ftrmmLeftLowTrans (const Field& F, const FFLAS_DIAG Diag, 
				       const size_t M, const size_t N,
				       const typename Field::Element * A, const size_t lda,
				       typename Field::Element * B, const size_t ldb,
				       const size_t nmax);
	
	template<class Field>
	static void ftrmmRightUpNoTrans (const Field& F, const FFLAS_DIAG Diag, 
					 const size_t M, const size_t N,
					 const typename Field::Element * A, const size_t lda,
					 typename Field::Element * B, const size_t ldb, 
					 const size_t nmax);

	template<class Field>
	static void ftrmmRightUpTrans (const Field& F, const FFLAS_DIAG Diag, 
				       const size_t M, const size_t N,
				       const typename Field::Element * A, const size_t lda,
				       typename Field::Element * B, const size_t ldb, 
				       const size_t nmax);

	template<class Field>
	static void ftrmmRightLowNoTrans (const Field& F, const FFLAS_DIAG Diag, 
					  const size_t M, const size_t N,
					  const typename Field::Element * A, const size_t lda,
					  typename Field::Element * B, const size_t ldb,
					  const size_t nmax);
	template<class Field>
	static void ftrmmRightLowTrans (const Field& F, const FFLAS_DIAG Diag, 
					const size_t M, const size_t N,
					const typename Field::Element * A, const size_t lda,
					typename Field::Element * B, const size_t ldb, 
					const size_t nmax);
}; // class FFLAS

#include "fflas_bounds.inl"
#include "fflas_fgemm.inl"
#include "fflas_fgemv.inl"
#include "fflas_fger.inl"
#include "fflas_ftrsm.inl"
#include "fflas_ftrmm.inl"
#include "fflas_ftrsv.inl"
#include "fflas_faxpy.inl"
#include "fflas_fdot.inl"
#include "fflas_fcopy.inl"

#ifdef _LINBOX_LINBOX_CONFIG_H
}
#endif

#endif // __FFLAS_H

