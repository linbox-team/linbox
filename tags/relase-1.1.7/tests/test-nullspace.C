/* Copyright (C) 2009 LinBox
 * Written by <boyer.brice@gmail.com>
 * Inspired and adapted from test-ffpack.C
 * see COPYING for license details
 */

/** \file tests/nullspace.C
  \brief Tests the nullspace functions
  */

#include "../linbox/linbox-config.h"
#include <iostream>
#include "../linbox/integer.h"
#include "../linbox/matrix/matrix-domain.h"
//#include "linbox/field/givaro-zpz.h"
#include "../linbox/field/modular.h"
//#include "linbox/ffpack/ffpack.h"
#include "../linbox/solutions/nullspace.h"
#include <vector>
#include "./test-common.h"
//#include "Matio.h" // write_field ;

using namespace LinBox;


/** @brief gives a random number such that \f$0 \leq RIII < s\f$.
 * @details basic..
 * @param [in]  s sup
//         * /param [in]  seed is a seed. If \p 0 (default) we create a new one.
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
	RIII          = (size_t) s * (alea/(RAND_MAX+1.0));
	assert(RIII<s);
	return RIII ;
}

/*! Creates a random Lapack style Permutation \p P of size \p len.
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


/*! Builds a \p m x \p n random matrix of rank \p rank over field \p F.
*/
template <class Field >
void RandomMatrixWithRank(const Field & F,
			  typename Field::Element * A,
			  const size_t & m,
			  const size_t & n,
			  const size_t & rank)
{
	assert(rank <= m);
	assert(rank <= n);

	//                srandom( (unsigned) time(NULL)  ) ; // on met une nouvelle graine.
	typename Field::RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	typename Field::Element one,zero;
	F.init(one,1UL);
	F.init(zero,0UL);

	typename Field::Element * B = new typename Field::Element[m*m];
	typename Field::Element * C = new typename Field::Element[m*n];
	// Create B a random invertible matrix (m x m format)
	for (size_t j=0 ; j<m ; ++j){
		for (size_t i = 0 ; i<j ; ++i)
			F.assign (*(B+i*m+j),zero); // triangulaire
		F.assign(*(B+j*m+j),one  );
		for (size_t i=j+1; i<m;++i)
			Gn.random (*(B+i*m+j)); // random mais pas nul.. euh... et sur Z/2 ?? :/
	}
	// Create C a random matrix of rank \p ( m x n format)
	//                for (size_t i = 0; i < std::min(rank,m); ++i){
	for (size_t i = 0; i < rank; ++i){
		size_t j = 0;
		for ( ; j < std::min(i,n) ; ++j)
			F.assign (*(C+i*n+j),zero);
		for ( ; j < n ; ++j)
			Gn.random (*(C+i*n+j));
	}
	for (size_t i = n*rank; i < n*m; ++i){
		F.assign (*(C+i),zero);
	}
	assert(CheckRank(F,C,m,n,n,rank));
	// create P a random permutation of size \p n
	size_t *P = new size_t[n];
	//srandom( (unsigned) time(NULL) ) ; // on met une nouvelle graine.
	RandomPermutation(P,n);
	// create Q a random permutation of size \p m
	size_t *Q = new size_t[m];
	RandomPermutation(Q,m);
	FFPACK::applyP(F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, m, C, n, Q );
	//PrintLapackPermutation(P,n,std::cout << "P == ");
	//write_field (F, std::cout<<"C_perm1="<<std::endl, C, m, n, n);
	FFPACK::applyP(F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m, 0, n, C, n, P );
	//PrintLapackPermutation(Q,m,std::cout << "Q == ");
	//write_field (F, std::cout<<"C_perm2="<<std::endl, C, m, n, n);
	// A = B*C (m x n format), of rank \p rank
	FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, m,
		      one, B, m, C, n, zero, A, n );
	delete[] B;
	delete[] C;
	delete[] P;
	delete[] Q;
	assert(CheckRank(F,A,m,n,n,rank));
	return;

}


/** 
 * @brief Checks we got the right rank.
 * 
 * @param F 
 * @param A 
 * @param m 
 * @param n 
 * @param lda 
 * @param alledged_rank 
 * 
 * @return 
 */
template <class Field>
bool CheckRank( const Field & F, 
		const typename Field::Element * A, 
		const size_t & m, 
		const size_t & n, 
		const size_t & lda, 
		const size_t & alledged_rank)
{
	//                std::cout << " is rank truely " << alledged_rank << " ?" << std::endl;
	typename Field::Element * Acopy = new typename Field::Element[m*n] ;
	FFLAS::fcopy(F,m*n,Acopy,1,A,1);
	size_t true_rank = FFPACK::Rank(F,m,n,Acopy,lda);
	delete[] Acopy ;
	//                std::cout << "It is " << true_rank << "." << std::endl;
	return (alledged_rank == true_rank);
}

/**
 * @brief Tests the NullSpace routines
 * @param F field
 * @param m row
 * @param n col
 * @param rank \p n-rank is the size of the NullSpace
 * @param iterations number of its
 * @param a_droite \p true if.. \p false if on the left
 * @return \p true hopefully if test's passed!
 */
template <class Field >
static bool testNullSpaceBasis (const Field& F, size_t m, size_t n, size_t rank, int iterations, bool a_droite) 
{
	typedef typename Field::Element			Element;

	//Commentator commentator;
	//commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	//commentator.start (pretty("Testing NullSpace Decomposition"),"testNullSpace",iterations);
	typename Field::Element one,zero;
	F.init(one,1UL);
	F.init(zero,0UL);

	bool ret = true;
	{
		size_t min = std::min(m,n);
		if (rank > min ) {
			rank = min; // rank <= min(m,n)...
		}
	}
	for (int k=0; k<iterations; ++k) {

		//commentator.progress(k);
		Element * A = new Element[m*n];
		size_t ld_a =  n ;
		size_t wd_a =  m ;
		RandomMatrixWithRank(F,A,m,n,rank);

		Element * Abis = new Element[m*n]; // copie de A
		for (size_t i=0; i<m*n; ++i)
			*(Abis+i) = *(A+i);
		size_t ker_dim = 0 ; // or coker_dim
		Element * Kern  = NULL;
		size_t ld_k = 0 ;
		if (a_droite) { 
			NullSpaceBasis (F, FFLAS::FflasRight,m,n,A,ld_a,Kern,ld_k,ker_dim);
			if (ker_dim != (ld_a - rank)) {
				ret = false;
				cout << "faux : (1) mauvaises dim : " << ker_dim << " != " << (ld_a - rank) << endl;
				delete[] Kern;
				delete[] A;
				delete[] Abis;
				break ;
			}
		} else {
			NullSpaceBasis (F, FFLAS::FflasLeft,m,n,A,ld_a,Kern,ld_k,ker_dim);
			if (ker_dim != (wd_a - rank) ) {
				ret = false;
				cout << "faux : (1) mauvaises dim " << ker_dim << " != " << (wd_a - rank)  << endl;
				delete[] Kern;
				delete[] A;
				delete[] Abis;
				break ;
			}
		}
		size_t ld_ker = (a_droite)?ker_dim:m ;
		size_t wd_ker = (a_droite)?n:ker_dim ;
		assert(ld_ker == ld_k) ;
		size_t ld_n = (a_droite)?ker_dim:ld_a;
		size_t wd_n = (a_droite)?wd_a:ker_dim;
		assert(CheckRank(F,Kern,wd_ker,ld_ker,ld_ker,ker_dim)); // ...il est bien de rang plein...
		Element * NullMat = new Element[ld_n*wd_n] ;// ...et on s'attend à ce que ça soit nul !
		if ( a_droite){
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, wd_a, ld_ker, ld_a,
				     one, Abis, ld_a, Kern, ld_ker , zero, NullMat, ld_n);
		}else{
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, wd_ker,ld_a, ld_ker,
				     one,  Kern, ld_ker , Abis, ld_a, zero, NullMat, ld_n);
		}

		//write_field (F, std::cout<<"A="<<endl, A, m, n, n,true);
		//write_field (F, std::cout<<"Abis="<<endl, Abis, m, n, n, true);
		delete[] Abis ;
		delete[] A ;
		delete[] Kern ;
		for (size_t i = 0 ; i < wd_n ; ++i ){
			for (size_t j = 0 ; j < ld_n ; ++j ) {
				if (!F.isZero(*(NullMat + j+i*ld_n)) ){
					    //	write_field (F, std::cout<<"faux : (3) NullMat pas nulle. "<<std::endl, NullMat, wd_n, ld_n, ld_n, true);
					delete[] NullMat ;
					ret = false;
					break;
				}
			}
			if (!ret)
				break;
		}
		if (ret) delete[] NullMat ;
		else break;
		
		//delete[] Kern ;
	}

	//commentator.stop(MSG_STATUS (ret), (const char *) 0, "testNullSpace");

	return ret;
}


int main(int argc, char** argv)
{
	//-----------------------------------------------------------------------
	// Choice of the finite field representation
	//typedef GivaroZpz<Std32> Field;
	typedef Modular<double> Field;
	//typedef Modular<float> Field;
	//typedef Modular<LinBox::uint32> Field;
	//------------------------------------------------------------------------

	bool pass = true;

	static size_t n = 5;
	static size_t m = 4;
	static size_t r = 2;
	static integer q = 11;
	static int iterations =2;

	static Argument args[] = {
		{ 'n', "-n N", "Set width of test matrices.",			TYPE_INT,     &n },
		{ 'm', "-m M", "Set hight of test matrices.",			TYPE_INT,     &m },
		{ 'r', "-r R", "Set rank of test matrices.",			TYPE_INT,     &r },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].",		TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",           TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	Field F (q);

	//srand(time (NULL));

	//commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start("NullSpace test suite", "nullspace");

	std::ostream& report = commentator.report();
    
	report << "\t \033[1;35m>>>\033[0;m \t testing left kernel" << endl ;
	if (!testNullSpaceBasis (F, m,n,r, iterations, false))
		pass=false;
	if (pass) report << "\t \033[1;36m<<<\033[0;m \t left kernel passed :)" << endl; else {report << "\t \033[1;31m!!!\033[0;m \t left kernel failed :(" << endl ; exit(-1);}
	report << "\t \033[1;35m>>>\033[0;m \t testing right kernel" << endl ;
	if (!testNullSpaceBasis (F, m,n,r, iterations, true))
		pass=false;
	if (pass) report << "\t \033[1;36m<<<\033[0;m \t right kernel passed :)" << endl; else {report << "\t \033[1;31m!!!\033[0;m \t right kernel failed :(" << endl ; exit(-1);}

	report << "\033[1;32m +++ ALL MY TESTS PASSED +++\033[0;m" << endl;


	commentator.stop("NullSpace test suite");
	return (pass ? 0 : -1);
}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
