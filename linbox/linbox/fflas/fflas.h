/* fflas.h
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

/** @file fflas.h
 * @author Clément Pernet.
 * @brief <b>F</b>inite <b>F</b>ield <b>L</b>inear <b>A</b>lgebra <b>S</b>ubroutines
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

/// @brief FFLAS: <b>F</b>inite <b>F</b>ield <b>L</b>inear <b>A</b>lgebra <b>S</b>ubroutines.
class FFLAS
{

public:
	enum FFLAS_TRANSPOSE
	{
		FflasNoTrans=111, /**< Matrix is not transposed */
		FflasTrans  =112  /**< Matrix is transposed */
	};
	enum FFLAS_UPLO
	{
		FflasUpper=121,  /**< Triangular matrix is Upper triangular (if \f$i>j\f$ then \f$T_{i,j} = 0\f$)*/
		FflasLower=122   /**< Triangular matrix is Lower triangular (if \f$i<j\f$ then \f$T_{i,j} = 0\f$)*/
	};
	enum FFLAS_DIAG
	{
		FflasNonUnit=131 ,  /**< Triangular matrix has an explicit general diagonal */
	       	FflasUnit   =132    /**< Triangular matrix has an implicit unit diagonal (\f$T_{i,i} = 1\f$)*//**< */
	};
	enum FFLAS_SIDE
	{
		FflasLeft  =141,  /**< Operator applied on the left */
		FflasRight = 142  /**< Operator applied on the rigth*/
	};

	/** \p FFLAS_BASE  determines the type of the element representation for Matrix Mult kernel.  */
	enum FFLAS_BASE
	{
		FflasDouble  = 151,  /**<  to use the double precision BLAS */
		FflasFloat   = 152,  /**<  to use the single precison BLAS */
		FflasGeneric = 153   /**< for any other domain, that can not be converted to floating point integers */
	};

	/* Representations of Z with floating point elements*/
	typedef UnparametricField<float> FloatDomain;
	typedef UnparametricField<double> DoubleDomain;


//---------------------------------------------------------------------
// Level 1 routines
//---------------------------------------------------------------------
	/** fscal.
	 * \f$x \gets a.x\f$
	 * \p X, is a vector of size \p N and stride \p incX.
	 */
	template<class Field>
	static void
	fscal (const Field& F, const size_t n, const typename Field::Element alpha,
	       typename Field::Element * X, const size_t incX)
	{

		typename Field::Element * Xi = X;
		for (; Xi < X+n*incX; Xi+=incX )
			F.mulin( *Xi, alpha );
	}

	/** \brief fcopy : \f$x \gets y \f$.
	 * \param X vector of size \p N and stride \p incX.
	 * \param Y vector of size \p N and stride \p incY.
	 */
	template<class Field>
	static void
	fcopy (const Field& F, const size_t N,
	       typename Field::Element * X, const size_t incX,
	       const typename Field::Element * Y, const size_t incY );

	/** \brief faxpy : \f$y \gets a.x + y\f$.
	 * \param X vector of size \p N and stride \p incX.
	 * \param Y vector of size \p N and stride \p incY.
	 */
	template<class Field>
	static void
	faxpy (const Field& F, const size_t N,
	       const typename Field::Element a,
	       const typename Field::Element * X, const size_t incX,
	       typename Field::Element * Y, const size_t incY );

	/** \brief fdot: dot product \f$x^T  y\f$.
	 * \param X vector of size \p N and stride \p incX.
	 * \param Y vector of size \p N and stride \p incY.
	 */
	template<class Field>
	static typename Field::Element
	fdot (const Field& F, const size_t N,
	      const typename Field::Element * X, const size_t incX,
	      const typename Field::Element * Y, const size_t incY );

	/** \brief fswap: \f$ X \leftrightarrow Y\f$.
	 * \param X vector of size \p N and stride \p incX.
	 * \param Y vector of size \p N and stride \p incY.
	 */
	template<class Field>
	static void
	fswap (const Field& F, const size_t N, typename Field::Element * X, const size_t incX,
	       typename Field::Element * Y, const size_t incY )
	{

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

	/** fsub : matrix addition.
	 * Computes \p C = \p A + \p B.
	 */
        template <class Field>
        static void
        fadd (const Field& F, const size_t M, const size_t N,
              const typename Field::Element* A, const size_t lda,
              const typename Field::Element* B, const size_t ldb,
              typename Field::Element* C, const size_t ldc)
        {
            const typename Field::Element *Ai = A, *Bi = B;
            typename Field::Element *Ci = C;
            for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
                for (size_t i=0; i<N; i++)
                    F.add (Ci[i], Ai[i], B[i]);
        }

	/** fsub : matrix subtraction.
	 * Computes \p C = \p A - \p B.
	 */
        template <class Field>
        static void
        fsub (const Field& F, const size_t M, const size_t N,
              const typename Field::Element* A, const size_t lda,
              const typename Field::Element* B, const size_t ldb,
              typename Field::Element* C, const size_t ldc)
        {
            const typename Field::Element * Ai = A, *Bi = B;
            typename Field::Element *Ci = C;
            for (; Ai < A+M*lda; Ai+=lda, Bi+=ldb, Ci+=ldc)
                for (size_t i=0; i<N; i++)
                    F.sub (Ci[i], Ai[i], B[i]);
        }

	/**  @brief finite prime Field GEneral Matrix Vector multiplication.
	 *
	 *  Computes  \f$Y \gets \alpha \mathrm{op}(A) X + \beta Y \f$.
	 * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
	 * \param A is \p M x \p N
	 * \param X is a vector of size \p N or \p M according to \p TransA
	 * \param Y is a vector of size \p M or \p N.
	 * \param incx stride in \p X
	 * \param incy stride in \p Y
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

	/**  @brief fger: GEneral ?
	 *
	 *  Computes  \f$A \gets \alpha x . y^T + A\f$
	 *  \p A is \p M x \p N
	 *  \p x and \p y are vectors of size \p M and \p N
	 */
	template<class Field>
	static void
	fger (const Field& F, const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      const typename Field::Element * x, const size_t incx,
	      const typename Field::Element * y, const size_t incy,
	      typename Field::Element * A, const size_t lda);

	/** @brief ftrsv: TRiangular System solve with Vector
	 *  Computes  \f$ X \gets \mathrm{op}(A^{-1}) X\f$
	 * @param X vector of size \p N on a field \p F
	 * @param A a matrix of leading dimension \p lda
	 * @param N number of rows or columns of \p A according to \p TransA
	 * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
	 * \param Diag if \c Diag==FflasUnit then \p A is unit.
	 * \param Uplo if \c Uplo==FflasUpper then \p A is upper triangular
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

	/** @brief ftrsm: TRiangular System solve with matrix.
	 * Computes  \f$ B \gets \alpha \mathrm{op}(A^{-1}) B\f$ or  \f$B \gets \alpha B \mathrm{op}(A^{-1})\f$.
	 * \param F field
	 * \param M rows of \p B
	 * \param N cols of \p B
	 * \param Side if \c Side==FflasLeft then  \f$ B \gets \alpha \mathrm{op}(A^{-1}) B\f$ is computed.
	 * \param Diag if \c Diag==FflasUnit then \p A is unit.
	 * \param Uplo if \c Uplo==FflasUpper then \p A is upper triangular
	 * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
	 * @warning unsafe with \c Trans==FflasTrans (debugging in progress)
	 */
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

	/** @brief ftrmm: TRiangular Matrix Multiply.
	 * Computes  \f$ B \gets \alpha \mathrm{op}(A) B\f$ or  \f$B \gets \alpha B \mathrm{op}(A)\f$
	 * B is M*N, A is M*M if
	 * \param M rows of \p B
	 * \param N cols of \p B
	 * \param A if \c Side==FflasLeft then \p A is \f$N\times N\f$, otherwise \p A is \f$M\times M\f$
	 * \param F field
	 * \param Side if \c Side==FflasLeft then  \f$ B \gets \alpha \mathrm{op}(A) B\f$ is computed.
	 * \param Diag if \c Diag==FflasUnit then \p A is unit.
	 * \param Uplo if \c Uplo==FflasUpper then \p A is upper triangular
	 * \param TransA if \c TransA==FflasTrans then \f$\mathrm{op}(A)=A^t\f$.
	 * @warning unsafe with \c Trans==FflasTrans (debugging in progress)
	 */
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

	/** @brief  Field GEneral Matrix Multiply.
	 *
	 * Computes \f$C = \alpha \mathrm{op}(A) \times \mathrm{op}(B) + \beta C\f$
	 * \param ta if \c ta==FflasTrans then \f$\mathrm{op}(A)=A^t\f$, else \f$\mathrm{op}(A)=A\f$,
	 * \param w recursive levels of Winograd's algorithm are used
	 * \param A \f$\mathrm{op}(A)\f$ is \f$m \times k\f$
	 * \param B \f$\mathrm{op}(B)\f$ is \f$k \times n\f$
	 * \param C \f$C\f$ is \f$m \times n\f$
	 * \param lda leading dimension of \p A
	 * \param ldb leading dimension of \p B
	 * \param ldc leading dimension of \p C
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
	       const size_t w)
	{

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

	/** @brief  Field GEneral Matrix Multiply.
	 *
	 * Computes \f$C = \alpha \mathrm{op}(A) \mathrm{op}(B) + \beta C\f$.
	 * Automatically set Winograd recursion level
	 * \param ta if \c ta==FflasTrans then \f$\mathrm{op}(A)=A^t\f$, else \f$\mathrm{op}(A)=A\f$,
	 * \param w recursive levels of Winograd's algorithm are used
	 * \param A \f$\mathrm{op}(A)\f$ is \f$m \times k\f$
	 * \param B \f$\mathrm{op}(B)\f$ is \f$k \times n\f$
	 * \param C \f$C\f$ is \f$m \times n\f$
	 * \param lda leading dimension of \p A
	 * \param ldb leading dimension of \p B
	 * \param ldc leading dimension of \p C

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
	       typename Field::Element* C, const size_t ldc)
	{

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

	/** @brief fsquare: Squares a matrix.
	 * compute \f$ C \gets \alpha \mathrm{op}(A) \mathrm{op}(A) + \beta C\f$ over a Field \p F
	 * Avoid the conversion of B
	 * @param ta  if \c ta==FflasTrans, \f$\mathrm{op}(A)=A^T\f$.
	 */
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
	/** @brief ftrtr: Triangular-Triangular matrix multiplication.
	 * B \gets alpha op(A)*B (for FFLAS_SIDE::FflasLeft)
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

	template<class Field>
	static void faddm(const Field & F,
			  const FFLAS_TRANSPOSE transA,
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb);

	template<class Field>
	static void faddm(const Field & F,
			  const FFLAS_TRANSPOSE transA,
			  const FFLAS_TRANSPOSE transB,
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  const typename Field::Element * B, const size_t ldb,
			  typename Field::Element * C, const size_t ldc );

	template<class Field>
	static void fsubm(const Field & F,
			  const FFLAS_TRANSPOSE transA,
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  typename Field::Element * B, const size_t ldb) ;

	template<class Field>
	static void fsubm(const Field & F,
			  const FFLAS_TRANSPOSE transA,
			  const FFLAS_TRANSPOSE transB,
			  const size_t M, const size_t N,
			  const typename Field::Element * A, const size_t lda,
			  const typename Field::Element * B, const size_t ldb,
			  typename Field::Element * C, const size_t ldc );

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
	template<class Element>
	class faddmTrans;
	template<class Element>
	class faddmNoTrans;
	template<class Element>
	class fsubmTrans;
	template<class Element>
	class fsubmNoTrans;
	template<class Element>
	class faddmTransTrans;
	template<class Element>
	class faddmNoTransTrans;
	template<class Element>
	class faddmTransNoTrans;
	template<class Element>
	class faddmNoTransNoTrans;
	template<class Element>
	class fsubmTransTrans;
	template<class Element>
	class fsubmNoTransTrans;
	template<class Element>
	class fsubmTransNoTrans;
	template<class Element>
	class fsubmNoTransNoTrans;



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

#include "fflas_faddm.inl"

#ifdef _LINBOX_LINBOX_CONFIG_H
}
#endif

#undef LB_TRTR

#endif // __LINBOX_fflas_H


/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
