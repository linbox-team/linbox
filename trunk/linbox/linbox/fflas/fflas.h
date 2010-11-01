/* fflas.h
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_fflas_H
#define __LINBOX_fflas_H
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
#include "linbox/field/modular-double.h"
#include "linbox/field/modular-float.h"
#include "linbox/field/modular-balanced-double.h"
#include "linbox/field/modular-balanced-float.h"
namespace LinBox 
{
#else
#include "config-blas.h"
#include "fflas-ffpack/unparametric.h"
#include "fflas-ffpack/modular-positive.h"
#include "fflas-ffpack/modular-balanced.h"
#endif
	
#ifndef __LINBOX_STRASSEN_OPTIMIZATION
#define WINOTHRESHOLD 1000
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

//#define LB_TRTR
	
class FFLAS 
{

public:	
	enum FFLAS_TRANSPOSE { FflasNoTrans=111, FflasTrans=112};
	enum FFLAS_UPLO      { FflasUpper=121, FflasLower=122 };
	enum FFLAS_DIAG      { FflasNonUnit=131, FflasUnit=132 };
	enum FFLAS_SIDE      { FflasLeft=141, FflasRight = 142 };

	/* Determine the type of the element representation for Matrix Mult kernel
	 * FflasDouble: to use the double precision BLAS
	 * FflasFloat: to use the single precison BLAS
	 * FflasGeneric: for any other domain, that can not be converted to floating point integers
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

		if (F.isZero (alpha)){
			for (size_t i = 0; i<m; ++i)
				fscal(F, n, beta, C + i*ldc, 1);
			return C;
		}
		
		size_t kmax = 0;
		size_t winolevel = w;
		FFLAS_BASE base;
		MatMulParameters (F, MIN(MIN(m,n),k), beta, kmax, base,
				  winolevel, true);
		WinoMain (F, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta,
				 C, ldc, kmax, winolevel, base);
		return C;
		}
	
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
	       const typename Field::Element* A, const size_t lda,
	       const typename Field::Element* B, const size_t ldb, 
	       const typename Field::Element beta,
	       typename Field::Element* C, const size_t ldc){

		if (!(m && n && k)) return C;
		if (F.isZero (alpha)){
			for (size_t i = 0; i<m; ++i)
				fscal(F, n, beta, C + i*ldc, 1);
			return C;
		}
		
		size_t w, kmax;
 		FFLAS_BASE base;

		MatMulParameters (F, MIN(MIN(m,n),k), beta, kmax, base, w);

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
#ifdef LB_TRTR
	// BB
	/* B \gets alpha op(A)*B (for FFLAS_SIDE::FflasLeft)
	* A and B are triangular, with B UpLo
	* and op(A) = A, A^T according to TransA
	* A and B can be (non)unit
	* 
	*/ 
	template<class Field>
	static typename Field::Element* ftrtr (const Field& F, const FFLAS_SIDE Side,
					       const FFLAS_UPLO Uplo, 
					       const FFLAS_TRANSPOSE TransA,
					       const FFLAS_DIAG ADiag, 
					       const FFLAS_DIAG BDiag, 
					       const size_t M, 
					       const typename Field::Element alpha,
					       typename Field::Element * A, const size_t lda,
					       typename Field::Element * B, const size_t ldb);
#endif
	

	/**
	 * MatCopy
	 * Makes a copy of the matrix M into a new allocated space.
	 */
	template<class Field>
	static typename Field::Element* MatCopy (const Field& F,
						 const size_t M, const size_t N,
						 const typename Field::Element * A,
						 const size_t lda){

		typename Field::Element * C = new typename Field::Element[M*N];
		for (size_t i = 0; i < N; ++i)
			for (size_t j = 0; j < N; ++j)
				F.assign(*(C + i*N + j),*(A + i*lda + j));
		return C;
	}
	
	/**
	 * Winosteps
	 *
	 * \brief Computes the number of recursive levels to perform
	 *
	 * \param m the common dimension in the product AxB
	 */
	static size_t WinoSteps (const size_t m);
	
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

	/**
	 * MatMulParameters
	 *
	 * \brief Computes the threshold parameters for the cascade
	 *        Matmul algorithm
	 *
	 * 
	 * \param F Finite Field/Ring of the computation.
	 * \param k Common dimension of A and B, in the product A x B
	 * \param bet Computing AB + beta C
	 * \param delayedDim Returns the size of blocks that can be multiplied
	 *                   over Z with no overflow
	 * \param base Returns the type of BLAS representation to use
	 * \param winoRecLevel Returns the number of recursion levels of
	 *                     Strassen-Winograd's algorithm to perform
	 * \param winoLevelProvided tells whether the user forced the number of
	 *                          recursive level of Winograd's algorithm
	 *
	 * See [Dumas, Giorgi, Pernet, arXiv cs/0601133]
	 * http://arxiv.org/abs/cs.SC/0601133
	 */
	template <class Field>
	static void MatMulParameters (const Field& F,
				      const size_t k,
				      const typename Field::Element& beta,
				      size_t& delayedDim,
				      FFLAS_BASE& base,
				      size_t& winoRecLevel,
				      bool winoLevelProvided=false);

	
	/**
	 * DotprodBound
	 *
	 * \brief  computes the maximal size for delaying the modular reduction
	 *         in a dotproduct
	 *
	 * This is the default version assuming a conversion to a positive modular representation
	 * 
	 * \param F Finite Field/Ring of the computation
	 * \param winoRecLevel Number of recusrive Strassen-Winograd levels (if any, 0 otherwise)
	 * \param beta Computing AB + beta C
	 * \param base Type of floating point representation for delayed modular computations
	 * 
	 */
	template <class Field>
	static size_t DotProdBound (const Field& F,
			     const size_t w, 
			     const typename Field::Element& beta,
			     const FFLAS_BASE base);
	

	/**
	 * Internal function for the bound computation
	 * Generic implementation for positive representations
	 */
	template <class Field>
	static double computeFactor (const Field& F, const size_t w);
	

	/**
	 * BaseCompute
	 *
	 * \brief Determines the type of floating point representation to convert to,
	 *        for BLAS computations
	 * \param F Finite Field/Ring of the computation
	 * \param w Number of recursive levels in Winograd's algorithm
	 */
	template <class Field>
	static FFLAS_BASE BaseCompute (const Field& F, const size_t w);
		
	/**
	 * TRSMBound
	 *
	 * \brief  computes the maximal size for delaying the modular reduction
	 *         in a triangular system resolution
	 *
	 *  Compute the maximal dimension k, such that a unit diagonal triangular
	 *  system of dimension k can be solved over Z without overflow of the
	 *  underlying floating point representation.
  	 *  See [Dumas, Giorgi, Pernet 06, arXiv:cs/0601133 ]
	 * 
	 * \param F Finite Field/Ring of the computation
	 * 
	 */
	template <class Field>
	static size_t TRSMBound (const Field& F);

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
				    const size_t  ); //kmax

	template<class Field>
	static void MatVectProd (const Field& F, 
				 const FFLAS_TRANSPOSE TransA, 
				 const size_t M, const size_t N,
				 const typename Field::Element alpha, 
				 const typename Field::Element * A, const size_t lda,
				 const typename Field::Element * X, const size_t incX, 
				 const typename Field::Element beta,
				 typename Field::Element * Y, const size_t incY);

	template<class Field>
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

	// Specialized routines for ftrsm
	template <class Element>
	class ftrsmLeftUpperNoTransNonUnit;
	template <class Element>
	class ftrsmLeftUpperNoTransUnit;
	template <class Element>
	class ftrsmLeftUpperTransNonUnit;
	template <class Element>
	class ftrsmLeftUpperTransUnit;
	template <class Element>
	class ftrsmLeftLowerNoTransNonUnit;
	template <class Element>
	class ftrsmLeftLowerNoTransUnit;
	template <class Element>
	class ftrsmLeftLowerTransNonUnit;
	template <class Element>
	class ftrsmLeftLowerTransUnit;
	template <class Element>
	class ftrsmRightUpperNoTransNonUnit;
	template <class Element>
	class ftrsmRightUpperNoTransUnit;
	template <class Element>
	class ftrsmRightUpperTransNonUnit;
	template <class Element>
	class ftrsmRightUpperTransUnit;
	template <class Element>
	class ftrsmRightLowerNoTransNonUnit;
	template <class Element>
	class ftrsmRightLowerNoTransUnit;
	template <class Element>
	class ftrsmRightLowerTransNonUnit;
	template <class Element>
	class ftrsmRightLowerTransUnit;

	// Specialized routines for ftrmm
	template <class Element>
	class ftrmmLeftUpperNoTransNonUnit;
	template <class Element>
	class ftrmmLeftUpperNoTransUnit;
	template <class Element>
	class ftrmmLeftUpperTransNonUnit;
	template <class Element>
	class ftrmmLeftUpperTransUnit;
	template <class Element>
	class ftrmmLeftLowerNoTransNonUnit;
	template <class Element>
	class ftrmmLeftLowerNoTransUnit;
	template <class Element>
	class ftrmmLeftLowerTransNonUnit;
	template <class Element>
	class ftrmmLeftLowerTransUnit;
	template <class Element>
	class ftrmmRightUpperNoTransNonUnit;
	template <class Element>
	class ftrmmRightUpperNoTransUnit;
	template <class Element>
	class ftrmmRightUpperTransNonUnit;
	template <class Element>
	class ftrmmRightUpperTransUnit;
	template <class Element>
	class ftrmmRightLowerNoTransNonUnit;
	template <class Element>
	class ftrmmRightLowerNoTransUnit;
	template <class Element>
	class ftrmmRightLowerTransNonUnit;
	template <class Element>
	class ftrmmRightLowerTransUnit;

	// BB : ça peut servir...
#ifdef LB_TRTR
	template <class Element>
	class ftrtrLeftUpperNoTransNonUnitNonUnit;
	template <class Element>
	class ftrtrLeftUpperNoTransUnitNonUnit;
	template <class Element>
	class ftrtrLeftUpperTransNonUnitNonUnit;
	template <class Element>
	class ftrtrLeftUpperTransUnitNonUnit;
	template <class Element>
	class ftrtrLeftLowerNoTransNonUnitNonUnit;
	template <class Element>
	class ftrtrLeftLowerNoTransUnitNonUnit;
	template <class Element>
	class ftrtrLeftLowerTransNonUnitNonUnit;
	template <class Element>
	class ftrtrLeftLowerTransUnitNonUnit;
	template <class Element>
	class ftrtrLeftUpperNoTransNonUnitUnit;
	template <class Element>
	class ftrtrLeftUpperNoTransUnitUnit;
	template <class Element>
	class ftrtrLeftUpperTransNonUnitUnit;
	template <class Element>
	class ftrtrLeftUpperTransUnitUnit;
	template <class Element>
	class ftrtrLeftLowerNoTransNonUnitUnit;
	template <class Element>
	class ftrtrLeftLowerNoTransUnitUnit;
	template <class Element>
	class ftrtrLeftLowerTransNonUnitUnit;
	template <class Element>
	class ftrtrLeftLowerTransUnitUnit;
	template <class Element>
	class ftrtrRightUpperNoTransNonUnitNonUnit;
	template <class Element>
	class ftrtrRightUpperNoTransUnitNonUnit;
	template <class Element>
	class ftrtrRightUpperTransNonUnitNonUnit;
	template <class Element>
	class ftrtrRightUpperTransUnitNonUnit;
	template <class Element>
	class ftrtrRightLowerNoTransNonUnitNonUnit;
	template <class Element>
	class ftrtrRightLowerNoTransUnitNonUnit;
	template <class Element>
	class ftrtrRightLowerTransNonUnitNonUnit;
	template <class Element>
	class ftrtrRightLowerTransUnitNonUnit;
	template <class Element>
	class ftrtrRightUpperNoTransNonUnitUnit;
	template <class Element>
	class ftrtrRightUpperNoTransUnitUnit;
	template <class Element>
	class ftrtrRightUpperTransNonUnitUnit;
	template <class Element>
	class ftrtrRightUpperTransUnitUnit;
	template <class Element>
	class ftrtrRightLowerNoTransNonUnitUnit;
	template <class Element>
	class ftrtrRightLowerNoTransUnitUnit;
	template <class Element>
	class ftrtrRightLowerTransNonUnitUnit;
	template <class Element>
	class ftrtrRightLowerTransUnitUnit;

#endif

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

//BB
#ifdef LB_TRTR
#include "fflas_ftrtr.inl"
#endif

#ifdef _LINBOX_LINBOX_CONFIG_H
}
#endif

#undef LB_TRTR

#endif // __LINBOX_fflas_H


/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
