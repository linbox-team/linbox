/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//          Tests for the fflapack set of routines
//--------------------------------------------------------------------------
// usage: test-fflapack p A n, for n lsp factorization  
// of A over Z/pZ
//-------------------------------------------------------------------------

#include "linbox-config.h"
#include <iostream>
#include <linbox/integer.h>
#include <linbox/matrix/matrix-domain.h>
#include "linbox/field/givaro-zpz.h"
#include "linbox/field/modular-double.h"
#include "linbox/fflapack/fflapack.h"

#include <vector>
#include "test-common.h"

using namespace LinBox;

const int maxpretty = 35;

const char* pretty(const char* a) {

	string msg(a);
	string blank("     ");
	blank+=a;
	int msgsize= maxpretty - blank.size();
	string dot(".");
	for (int i=0;i<msgsize ;++i)
		blank+=dot;
	return blank.data();
}

/*
 *  Testing the rank of dense matrices 
 *  construct a n*n matrices of rank r and compute the rank
 */
template <class Field>
static bool testRank (const Field& F,size_t n, int iterations) {

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	Commentator mycommentator;
	mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing rank"),"testRank",iterations);

	Element one, zero;
	F.init( one, 1UL);
	F.init( zero, 0UL);
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 

	unsigned int r;
	bool ret = true;
	
	for (int k=0;k<iterations; ++k) {
    
		mycommentator.progress(k);
		Element * A = new Element[n*n];
		Element * S = new Element[n*n];
		Element * L = new Element[n*n];
     
		r = rand() % n;
		// create S as an upper triangular matrix with r nonzero rows
		for (size_t i=0;i<r;++i){
			Gn.random(*(S+i*n+i));
			for (size_t j=i+1;j<n;++j)     
				G.random(*(S+i*n+j));
		}
     
		// create L as a lower triangular matrix with nonzero elements on the diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				G.random(*(L+i*n+j));
			Gn.random(*(L+i*n+i));
		}
		
		//  compute A=LS
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, L, n, S, n, zero, A, n );
		
		// compute the rank of A
		unsigned int rank= FFLAPACK::Rank( F, n, n, A, n);
                
		if (rank!=r)
			ret=false;
	}
  	
	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testRank");
    
	return ret;
}

/*
 *  Testing the rank of dense matrices using TURBO algorithm
 *  construct a n*n matrices of rank r and compute the rank
 */
template <class Field>
static bool testTURBO (const Field& F,size_t n, int iterations) {

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	Commentator mycommentator;
	mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing TURBO"),"testTURBO",iterations);

	Element one, zero;
	F.init( one, 1UL);
	F.init( zero, 0UL);
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 

	unsigned int r;
	bool ret = true;
	size_t no2 = n>>1;
	
	for (int k=0;k<iterations; ++k) {
    
		mycommentator.progress(k);
		Element * A = new Element[n*n];
		Element * S = new Element[n*n];
		Element * L = new Element[n*n];
     
		r = rand() % n;
		// create S as an upper triangular matrix with r nonzero rows
		for (size_t i=0;i<r;++i){
			Gn.random(*(S+i*n+i));
			for (size_t j=i+1;j<n;++j)     
				G.random(*(S+i*n+j));
		}
     
		// create L as a lower triangular matrix with nonzero elements on the diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				G.random(*(L+i*n+j));
			Gn.random(*(L+i*n+i));
		}
		
		//  compute A=LS
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, L, n, S, n, zero, A, n );
		
		// compute the rank of A
		unsigned int rank= FFLAPACK::TURBO( F, n, n, 
						    A, n, A+no2,n, 
						    A+no2*n, n, A+no2*(n+1), n );
                
		if (rank!=r)
			ret=false;
	}
  	
	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testTURBO");
    
	return ret;
}


/*
 *  Testing the determinant of dense matrices using BlasDomain
 *  construct a n*n matrices of determinant d and compute the determinant
 */
template <class Field>
static bool testDet (const Field& F,size_t n, int iterations) {

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
  
	Commentator mycommentator;
	mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing determinant"),"testDet",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element one,zero,d;
	F.init(one,1UL);
	F.init(zero,0UL);

	bool ret = true;
	
	for (int k=0;k<iterations;++k) {
		
		mycommentator.progress(k);

		G.random(d);
		
		Element * A = new Element[n*n];
		Element * S = new Element[n*n];
		Element * L = new Element[n*n];

		

		// create S as an upper triangular matrix of full rank 
		// with diagonal's element equal to 1 except the first entry wich equals to d
		for (size_t i=0;i<n;++i){
			F.assign(*(S+i*n+i), one);
			for (size_t j=i+1;j<n;++j)      
				G.random(*(S+i*n+j));
		}
		F.assign(*S,d);
		
		// create L as a lower triangular matrix with only 1's on diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				G.random(*(L+i*n+j));
			F.assign(*(L+i*n+i),one);
		}

    
		//  compute A=LS
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, L, n, S, n, zero, A, n );

		// compute the determinant of A
		Element det= FFLAPACK::Det( F, n, n, A, n);
    
		if (!F.areEqual(det,d))
			ret=false;
	}
  
	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testDet");
    
	return ret;
}

/*
 * Test of the LQUP factorization routine
 */
template <class Field>
static bool testLUdivine (const Field& F, size_t m, size_t n, int iterations) {

	typedef typename Field::Element                  Element;
	typedef typename Field::RandIter                RandIter;

	Commentator mycommentator;
	mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing LQUP factorization"),"testLQUP",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element one,zero;
	F.init(one,1UL);
	F.init(zero,0UL);
  
	bool ret = true;
	MatrixDomain<Field> MD(F);
	for (int k=0;k<iterations;++k) {
    
		mycommentator.progress(k);    
		
		
		Element * A = new Element[m*n];
		Element * Abis = new Element[m*n];
		Element * B = new Element[m*m];
		Element * C = new Element[m*n];
		

		// Create B a random matrix of rank n/2
		for (size_t j=0;j<m;++j)
			if ( j % 2 )
				for (size_t i=0;i<m;++i)
					G.random(*(B+i*m+j));
			else
				for (size_t i=0;i<m;++i)
					F.assign(*(B+i*m+j),zero);
		// Create C a random matrix of rank n/2
		for (size_t i=0;i<m;++i)
			if ( i % 2 )
				for (size_t j=0;j<n;++j)
					G.random(*(C+i*n+j));
			else
				for (size_t j=0;j<n;++j)
					F.assign(*(C+i*n+j),zero);
		
		// A = B*C
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, m,
			      one, B, m, C, n, zero, A, n );
	
		
		for (size_t i=0; i<m*n; ++i)
			*(Abis+i) = *(A+i);

		size_t * P = new size_t[n];
		size_t * Q = new size_t[m];
		
		size_t r=FFLAPACK::LUdivine( F, FFLAS::FflasNonUnit, m, n, 
				    A, n, P, FFLAPACK::FflapackLQUP, Q);

		Element * L = new Element[m*m];
		Element * U = new Element[m*n];
	
		for (size_t i=0; i<m; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign(*(L+i*n+j), *(A+i*n+j) );
			for (size_t j=i; j<n; ++j)
				F.assign(*(L+i*n+j), zero );
		}
		FFLAPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m, 
				  0, m, L, m, Q);
		for (size_t i=0; i<m; ++i)
			F.assign( *(L+i*n+i), one);
		for (size_t i=0; i<m; ++i){
			for (size_t j=0; j<i; ++j)
				F.assign( *(U+i*n+j), zero );
			for (size_t j=i; j<n; ++j)
				F.assign( *(U+i*n+j), *(A+i*n+j) );
		}

		// C = U*P
		FFLAPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m, 
				  0, r, U, n, P);
		// C = Q*C
		FFLAPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 
				  0, m, U, n, Q);
		// A = L*C
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, m,
			      one, L, m, U, n, zero, A, n );
		
		
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<n;++j)
				if (!F.areEqual( *(A+i*n+j), *(Abis+i*n+j)))
					ret=false;
		

	}

	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testLQUP");
    
	return ret;
}


int main(int argc, char** argv){
	//-----------------------------------------------------------------------
	// Choice of the finite field representation
	//typedef GivaroZpz<Std32> Field;
	typedef Modular<double> Field;
	//------------------------------------------------------------------------

	bool pass = true;

	static size_t n = 200;
	static integer q = 101U;
	static int iterations =10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 256)",       TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q }, 
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
    
	Field F (q);

	srand(time (NULL));
	
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start("fflapack Test suite");
	std::cerr<<endl<<endl;

 	if (!testLUdivine (F, n,n, iterations)) pass=false;
 	if (!testRank (F, n, iterations))   pass = false;   
 	if (!testDet (F, n, iterations))   pass = false;   
  	if (!testTURBO (F, n, iterations))   pass = false;   
//  	if (!testapplyP  (F, n, iterations)) pass = false;
//  	if (!testInv  (F, n, iterations)) pass = false;
//  	if (!testMinPoly (F,n,iterations)) pass=false;
// 	if (!testCharPoly (F,n,iterations)) pass=false;
	
	std::cerr<<"\nfflapack Test suite...";
	commentator.stop(MSG_STATUS(pass),"fflapack Test suite");
    
	return pass ? 0 : -1;
}














