/* Copyright (C) 2013 LinBox
 * Written by bb <boyer.brice@gmail.com>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 * Function definitions for block Lanczos iteration
 */

/** \file tests/test-matrix-utils.h
 * \brief Utilities for tests on matrices
 * @ingroup tests
 * @warning this should better written and shared among the tests
 * @bug put in \c matrix/random-matrix.h
 */


#ifndef __LINBOX_tests_test_matrix_utils_H
#define __LINBOX_tests_test_matrix_utils_H


#include "fflas-ffpack/utils/Matio.h"

#define RAPPORT(a) \
	if (pass) {   \
		commentator().report() << "\t \033[1;36m<<<\033[0;m \t" << (a) << " passed :)" << std::endl;  \
	} \
        else { \
		commentator().report() << "\t \033[1;    31m!!!\033[0;m \t " << (a) << " failed :(" << std::endl ;  \
		exit(-1); \
	}

#define TESTE(a) \
	commentator().report() << "\t \033[1;35m>>>\033[0;m \t testing " << (a) << " :" << std::endl ;


#define element_t(Field) \
typename Field ::Element


namespace LinBox { // namespace tests
/** @brief gives a random number such that \f$0 \leq RIII < s\f$.
 * @details basic..
 * @param [in]  s sup
 * \param [in]  seed seed. If \p 0 (default) we create a new one.
 * @param [out] RIII random integer in the interval  \f$[[0, s-1]]\f$.
 * @return a reference to \p RIII
 */
size_t & RandIntInInt ( const size_t & s, size_t & RIII, const int & seed = 0 )
{
	/*
	 * if (seed == 0)
	 *	srandom( (unsigned) time(NULL) );
	 *else
	 *	srandom ( seed );
	 */
	double alea = rand();
	RIII          = (size_t) ((double)s * (alea/(RAND_MAX+1.0)));
	assert(RIII<s);
	return RIII ;
}

/*!
 * Creates a random Lapack style Permutation \p P of size \p len.
 */
void RandomPermutation ( size_t * P, const size_t & len)
{
	size_t alea = 0 ;
	for (size_t i = 0 ; i < len ; ++i) {
		RandIntInInt(len-i, alea);
		*(P+i) = i + alea ;
	}
	return;
}

int permutationDet(size_t *P, const size_t len) {
	int d = 1 ;
	for(size_t i = 0 ; i < len ; ++i)
		if (P[i] != i) d = -d;
	return d;
}

/**
 * @brief Checks we got the right rank.
 *
 * @param F field
 * @param A matrix
 * @param m rows
 * @param n cols
 * @param lda leadin dimmension
 * @param alledged_rank supposedly correct rank.
 *
 * @return \c alledged_rank==rank(A)
 */
template <class Field>
bool CheckRank( const Field & F,
		const element_t(Field) * A,
		const size_t & m,
		const size_t & n,
		const size_t & lda, //!@bug not used
		const size_t & alledged_rank)
{
	//  std::cout << " is rank truely " << alledged_rank << " ?" << std::endl;
	element_t(Field) * Acopy = new element_t(Field)[m*lda] ;
	FFLAS::fcopy(F,m*lda,Acopy,1,A,1);
	size_t true_rank = FFPACK::Rank(F,m,n,Acopy,lda);
	delete[] Acopy ;
	//                std::cout << "It is " << true_rank << "." << std::endl;
	return (alledged_rank == true_rank);
}

template <class Field>
bool CheckRank( const Field & F,
		const BlasMatrix<Field> & A,
		const size_t & alledged_rank)
{
	return CheckRank(F,A.getPointer(),A.rowdim(),A.coldim(),A.stride(),alledged_rank);
}
template <class Field>
bool CheckRank( const Field & F,
		const BlasSubmatrix<Field> & A,
		const size_t & alledged_rank)
{
	return CheckRank(F,A.getPointer(),A.rowdim(),A.coldim(),A.stride(),alledged_rank);
}


template <class Field>
bool CheckDet( const Field & F,
	       const element_t(Field) * A,
	       const size_t & m,
	       const size_t & lda,
	       const element_t(Field) & alledged_det)
{
	 // std::cout << " is det truely " << alledged_det << " ?" << std::endl;
	element_t(Field) * Acopy = new element_t(Field)[m*m] ;
	FFLAS::fcopy(F,m*m,Acopy,1,A,1);
	element_t(Field) true_det = FFPACK::Det(F,m,m,Acopy,lda);
	delete[] Acopy ;
	// std::cout << "It is " << true_det << "." << std::endl;
	return (alledged_det == true_det);
}

/*!
 * Builds a \p m x \p n random matrix of rank \p rank over field \p F.
 */
template <class Field >
void RandomMatrixWithRank(const Field & F,
			  element_t(Field) * A,
			  const size_t & m,
			  const size_t & n,
			  const size_t & lda,
			  const size_t & rank)
{
	assert(rank <= m);
	assert(rank <= n);

	//                srandom( (unsigned) time(NULL)  ) ; // on met une nouvelle graine.
	typename Field::RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	typedef element_t(Field) Element;
	// element_t(Field) one,zero;
	// F.init(one,1UL);
	// F.init(zero,0UL);

	Element * B = new Element[m*m];
	Element * C = new Element[m*n];
	// Create B a random invertible matrix (m x m format)
	for (size_t j=0 ; j<m ; ++j){
		size_t i= 0;
		for ( ; i<j ; ++i)
			F.assign (*(B+i*m+j),F.zero); // triangulaire
		assert(i==j);
		Gn.random(*(B+j*m+j));
		for (++i; i<m;++i)
			G.random (*(B+i*m+j));
	}
	// Create C a random matrix of rank \p ( m x n format)
	for (size_t i = 0; i < rank; ++i){
		size_t j = 0;
		for ( ; j < std::min(i,n) ; ++j)
			F.assign (*(C+i*n+j),F.zero);
		for ( ; j < n ; ++j)
			Gn.random (*(C+i*n+j));
	}
	for (size_t i = n*rank; i < n*m; ++i){
		F.assign (*(C+i),F.zero);
	}
	linbox_check(CheckRank(F,C,m,n,n,rank));
	// create P a random permutation of size \p n
	size_t *P = new size_t[n];
	RandomPermutation(P,n);
	// create Q a random permutation of size \p m
	size_t *Q = new size_t[m];
	RandomPermutation(Q,m);
	FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
		       n, 0, (int)m, C, n, Q );
	FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
		       m, 0, (int)n, C, n, P );
	FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, m,
		      F.one, B, m, C, n, F.zero, A, lda );
	delete[] B;
	delete[] C;
	delete[] P;
	delete[] Q;
	linbox_check(CheckRank(F,A,m,n,lda,rank));
	return;

}

/*!
 * Builds a \p m x \p m random matrix of determinant \p det over field \p F.
 */
template <class Field >
void RandomMatrixWithDet(const Field & F,
			 element_t(Field) * A,
			 const size_t & m,
			 const size_t & lda, //!@bug not used
			 const element_t(Field) & det)
{
	typename Field::RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	typedef element_t(Field) Element;

	Element * B = new Element[m*m];
	Element * C = new Element[m*m];

	size_t *P = new size_t[m];
	//srandom( (unsigned) time(NULL) ) ; // on met une nouvelle graine.
	RandomPermutation(P,m);
	// create Q a random permutation of size \p m
	size_t *Q = new size_t[m];
	RandomPermutation(Q,m);

	int sgn = permutationDet(P,m);
	sgn *= permutationDet(Q,m);


	// Create B a random invertible matrix (m x m format) of determinant \p det
	for (size_t j=0 ; j<m ; ++j){
		for (size_t i = 0 ; i<j ; ++i)
			F.assign (*(B+i*m+j),F.zero); // triangulaire
		Gn.random(*(B+j*m+j) );
		for (size_t i=j+1; i<m;++i)
			G.random (*(B+i*m+j));
	}
	{
		Element dx = F.one;
		size_t j=0;
		for ( ; j<m-1 ; ++j){
			F.mulin(dx,B[j*m+j]);
		}
		F.invin(dx);
		F.mulin(dx,det);
		if (sgn < 0)
			F.negin(dx);
		B[j*m+j]=dx;
	}
	// write_field(F,std::cout << "L:=",B,m,m,m,true) << ';'<<std::endl;
	Element mdet ; F.neg(mdet,det);
	assert(CheckDet(F,B,m,m,(sgn<0)?mdet:det));

	// Create C a random matrix of rank \p ( m x m format)
	// for (size_t i = 0; i < std::min(rank,m); ++i)
	for (size_t i = 0; i < m; ++i){
		size_t j = 0;
		for ( ; j < i ; ++j)
			F.assign (*(C+i*m+j),F.zero);
		assert(i==j);
		F.assign (*(C+i*m+j),F.one);
		for ( ++j ; j < m ; ++j)
			G.random (*(C+i*m+j));
	}
	assert(CheckDet(F,C,m,m,F.one));
	// write_field(F,std::cout << "U:=",C,m,m,m,true) << ';'<<std::endl;
	// create P a random permutation of size \p m
	FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
		       m, 0, (int)m, C, m, Q );
	//PrintLapackPermutation(P,m,std::cout << "P == ");
	// write_field (F, std::cout<<"C_perm1="<<std::endl, C, m, m, m);
	FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans,
		       m, 0, (int)m, C, m, P );
	//PrintLapackPermutation(Q,m,std::cout << "Q == ");
	//write_field (F, std::cout<<"C_perm2="<<std::endl, C, m, m, m);
	// A = B*C (m x m format), of rank \p rank
	FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, m, m,
		      F.one, B, m, C, m, F.zero, A, m );
	delete[] B;
	delete[] C;
	delete[] Q;
	delete[] P;

	// write_field(F,std::cout << "A:=",A,m,m,m,true) << ';'<<std::endl;
	assert(CheckDet(F,A,m,m,det));
	return;

}

} // LinBox


#endif // __LINBOX_tests_test_matrix_utils_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
