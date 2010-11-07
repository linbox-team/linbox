/* fflas/fflas_ftrmm_src.inl
 * Copyright (C) 2010 LinBox
 *
 * Written by Brice Boyer <Brice.Boyer@imag.fr>
 *
 * See COPYING for license information.
 */

#define Mjoin(pre, nam) my_join(pre, nam)
#define my_join(pre, nam) pre ## nam

#ifdef  __FFLAS__BNOTRANSPOSE
#ifdef __FFLAS__ANOTRANSPOSE
#define __FFLAS_incA    1
#define __FFLAS_ldA     lda
#define __FFLAS_A_Trans NoTrans
#define __FFLAS_incB    1
#define __FFLAS_ldB     ldb
#define __FFLAS_B_Trans NoTrans
#endif
#ifdef __FFLAS__ATRANSPOSE
#define __FFLAS_incA    lda
#define __FFLAS_ldA     1   
#define __FFLAS_A_Trans Trans
#define __FFLAS_incB    1
#define __FFLAS_ldB     ldb
#define __FFLAS_B_Trans NoTrans
#endif
#else
#ifdef __FFLAS__ANOTRANSPOSE
#define __FFLAS_incA    1
#define __FFLAS_ldA     lda
#define __FFLAS_A_Trans NoTrans
#define __FFLAS_incB    ldb
#define __FFLAS_ldB     1
#define __FFLAS_B_Trans Trans
#endif
#ifdef __FFLAS__ATRANSPOSE
#define __FFLAS_incA    lda
#define __FFLAS_ldA     1   
#define __FFLAS_A_Trans Trans
#define __FFLAS_incB    ldb
#define __FFLAS_ldB     1
#define __FFLAS_B_Trans Trans
#endif
#endif


#ifdef __FFLAS__FLOAT
#define __FFLAS_Element float
#endif
#ifdef __FFLAS__DOUBLE
#define __FFLAS_Element double
#endif
#ifdef __FFLAS__GENERIC
#define __FFLAS_Element Element
#endif

#ifndef __FFLAS__GENERIC
template<>
class FFLAS::Mjoin(faddm, Mjoin(__FFLAS_A_Trans, __FFLAS_B_Trans))<__FFLAS_Element> 
{
public :
	template<class Field>
	void operator() (const Field & F,
			 const size_t M, const size_t N,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 typename Field::Element * C, const size_t ldc)
	{
		if (!M || !N ) return; // ne doit jamais arriver, déjà testé !

		// adding (precomputing tB ?)
		for (size_t i = 0 ; i < M ; ++i)
			for (size_t j = 0 ; j < N ; ++j)
				F.add(*(C+i*ldc+j), *(A+i*__FFLAS_ldA+j*__FFLAS_incA), *(B+i*__FFLAS_ldB+j*__FFLAS_incB)) ;

		// reducing :
		if (M == lda )
			for (size_t i = 0 ; i < M*N ; ++i)
				F.init(*A+i,*A+i);
		else
			for (size_t i = 0 ; i < M ; ++i)
				for (size_t j = 0 ; j < N ; ++j)
					F.init(*(A+i*lda+j), *(A+i*lda+j));


	}

};
#else
template<class Element>
class FFLAS::Mjoin(faddm, Mjoin(__FFLAS_A_Trans, __FFLAS_B_Trans))
{
public :
	template<class Field>
	void operator() (const Field & F,
			 const size_t M, const size_t N,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 typename Field::Element * C, const size_t ldc)
	{
		if (!M || !N ) return; // ne doit jamais arriver, déjà testé !

		// adding (precomputing tB ?)
		for (size_t i = 0 ; i < M ; ++i)
			for (size_t j = 0 ; j < N ; ++j)
				F.add(*(C+i*ldc+j), *(A+i*__FFLAS_ldA+j*__FFLAS_incA), *(B+i*__FFLAS_ldB+j*__FFLAS_incB)) ;

		// reducing :
		if (M == lda )
			for (size_t i = 0 ; i < M*N ; ++i)
				F.init(*A+i,*A+i);
		else
			for (size_t i = 0 ; i < M ; ++i)
				for (size_t j = 0 ; j < N ; ++j)
					F.init(*(A+i*lda+j), *(A+i*lda+j));


	}

};
#endif

#ifndef __FFLAS__GENERIC
template<>
class FFLAS::Mjoin(fsubm, Mjoin(__FFLAS_A_Trans, __FFLAS_B_Trans))<__FFLAS_Element > 
{
public :
	template<class Field>
	void operator() (const Field & F,
			 const size_t M, const size_t N,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 typename Field::Element * C, const size_t ldc)
	{
		if (!M || !N ) return; // ne doit jamais arriver, déjà testé !

		// adding (precomputing tB ?)
		for (size_t i = 0 ; i < M ; ++i)
			for (size_t j = 0 ; j < N ; ++j)
				*(C+i*ldc+j) = *(A+i*__FFLAS_ldA+j*__FFLAS_incA) - *(B+i*__FFLAS_ldB+j*__FFLAS_incB) ;

		// reducing :
		if (M == lda )
			for (size_t i = 0 ; i < M*N ; ++i)
				F.init(*A+i,*A+i);
		else
			for (size_t i = 0 ; i < M ; ++i)
				for (size_t j = 0 ; j < N ; ++j)
					F.init(*(A+i*lda+j), *(A+i*lda+j));


	}

};
#else
template<class Element>
class FFLAS::Mjoin(fsubm,Mjoin(__FFLAS_A_Trans,__FFLAS_B_Trans))
{
public :
	template<class Field>
	void operator() (const Field & F,
			 const size_t M, const size_t N,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * B, const size_t ldb,
			 typename Field::Element * C, const size_t ldc)
	{
		if (!M || !N ) return; // ne doit jamais arriver, déjà testé !

		// adding (precomputing tB ?)
		for (size_t i = 0 ; i < M ; ++i)
			for (size_t j = 0 ; j < N ; ++j)
				F.sub(*(C+i*ldc+j), *(A+i*__FFLAS_ldA+j*__FFLAS_incA), *(B+i*__FFLAS_ldB+j*__FFLAS_incB)) ;

		// reducing :
		if (M == lda )
			for (size_t i = 0 ; i < M*N ; ++i)
				F.init(*A+i,*A+i);
		else
			for (size_t i = 0 ; i < M ; ++i)
				for (size_t j = 0 ; j < N ; ++j)
					F.init(*(A+i*lda+j), *(A+i*lda+j));


	}

};
#endif


#undef Mjoin
#undef my_join
#undef __FFLAS_incB
#undef __FFLAS_incA
#undef __FFLAS_Element
#undef __FFLAS_ldA
#undef __FFLAS_ldB
#undef __FFLAS_Trans
#undef __FFLAS_A_Trans
#undef __FFLAS_B_Trans
#undef __FFLAS_B_inc
#undef __FFLAS_B_ld


/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
