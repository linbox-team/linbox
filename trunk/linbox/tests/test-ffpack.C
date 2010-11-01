/* Copyright (C) LinBox
 *
 * written by C. Pernet
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


/*! @file tests/test-ffpack.C
 * @brief Tests for the ffpack set of routines
 * usage: test-ffpack p A n, for n lsp factorization  of A over Z/pZ
 */

#include "linbox/linbox-config.h"
#include <iostream>
#include <linbox/integer.h>
#include <linbox/matrix/matrix-domain.h>
//#include "linbox/field/givaro-zpz.h"
#include "linbox/field/modular.h"
#include "linbox/ffpack/ffpack.h"
#include <vector>
#include "test-common.h"

using namespace LinBox;

string blank;

const int maxpretty = 35;
const char* pretty(string a) 
{

	blank = "     " + a;
	int msgsize= maxpretty - blank.size();
	string dot(".");
	for (int i=0;i<msgsize ;++i)
		blank+=dot;
	return blank.c_str();
}

/*
 *  Testing the rank of dense matrices 
 *  construct a n*n matrices of rank r and compute the rank
 */
template <class Field>
static bool testRank (const Field& F,size_t n, int iterations) 
{

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	//Commentator commentator;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start (pretty("Testing rank"),"testRank",iterations);

	Element one, zero;
	F.init( one, 1UL);
	F.init( zero, 0UL);
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 

	unsigned int r;
	bool ret = true;
	
	for (int k=0;k<iterations; ++k) {
    
		commentator.progress(k);
		Element * A = new Element[n*n];
		Element * S = new Element[n*n];
		Element * L = new Element[n*n];
     
		r = rand() % n;
		// create S as an upper triangular matrix with r nonzero rows
		for (size_t i=0;i<r;++i){
			for (size_t j=0;j<i;++j)
				F.assign(*(S+j+i*n),zero);
			Gn.random(*(S+i*n+i));
			for (size_t j=i+1;j<n;++j)     
				G.random(*(S+i*n+j));
		}
		for (size_t i=r;i<n;++i){
			for (size_t j=0;j<n;++j)
				F.assign(*(S+j+i*n),zero);
		}			
		// create L as a lower triangular matrix with nonzero elements on the diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				G.random(*(L+i*n+j));
			Gn.random(*(L+i*n+i));
			for (size_t j=i+1;j<n;++j)
				F.assign(*(L+j+i*n),zero);

		}
		
		//  compute A=LS
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, L, n, S, n, zero, A, n );
		delete[] L;
		delete[] S;

		// compute the rank of A
		unsigned int rank= FFPACK::Rank( F, n, n, A, n);
                delete[] A;
		if (rank!=r)
			ret=false;
	}
  	
	commentator.stop(MSG_STATUS (ret), (const char *) 0, "testRank");
    
	return ret;
}

/*
 *  Testing the rank of dense matrices using TURBO algorithm
 *  construct a n*n matrices of rank r and compute the rank
 */
template <class Field>
static bool testTURBO (const Field& F,size_t n, int iterations) 
{

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	//Commentator commentator;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start (pretty("Testing TURBO"),"testTURBO",iterations);

	Element one, zero;
	F.init( one, 1UL);
	F.init( zero, 0UL);
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 

	unsigned int r;
	bool ret = true;
		
	for (int k=0;k<iterations; ++k) {
    
		commentator.progress(k);
		Element * A = new Element[n*n];
		Element * S = new Element[n*n];
		Element * L = new Element[n*n];
     
		r = rand() % n;
		// create S as an upper triangular matrix with r nonzero rows
		for (size_t i=0;i<r;++i){
			for (size_t j=0;j<i;++j)
				F.assign(*(S+j+i*n),zero);
			Gn.random(*(S+i*n+i));
			for (size_t j=i+1;j<n;++j)     
				G.random(*(S+i*n+j));
		}
		for (size_t i=r; i<n; ++i)
			for (size_t j=0;j<n;++j)
				F.assign( *(S+i*n+j),zero);
		// create L as a lower triangular matrix with nonzero elements on the diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				G.random(*(L+i*n+j));
			Gn.random(*(L+i*n+i));
			for (size_t j=i+1;j<n;++j)
				F.assign(*(L+j+i*n),zero);
		}
	
		
		//  compute A=LS
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, L, n, S, n, zero, A, n );
		
		delete[] L;
		delete[] S;
		// compute the rank of A
		size_t * P = new size_t[n];
		size_t * Q = new size_t[n];
		unsigned int rank= FFPACK::TURBO( F, n, n, 
						  A,n, P, Q, 100);
// 						  A, n, A+no2,n, 
// 						    A+no2*n, n, A+no2*(n+1), n );

		delete[] P;
		delete[] Q;
		delete[] A;

		if (rank!=r)
			ret=false;
	}
  	
	commentator.stop(MSG_STATUS (ret), (const char *) 0, "testTURBO");
    
	return ret;
}


/*
 *  Testing the determinant of dense matrices using BlasDomain
 *  construct a n*n matrices of determinant d and compute the determinant
 */
template <class Field>
static bool testDet (const Field& F,size_t n, int iterations) 
{

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
  
	//Commentator commentator;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start (pretty("Testing determinant"),"testDet",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element one,zero,d;
	F.init(one,1UL);
	F.init(zero,0UL);

	bool ret = true;
	
	for (int k=0;k<iterations;++k) {
		
		commentator.progress(k);

		G.random(d);
		
		Element * A = new Element[n*n];
		Element * S = new Element[n*n];
		Element * L = new Element[n*n];

		

		// create S as an upper triangular matrix of full rank 
		// with diagonal's element equal to 1 except the first entry wich equals to d
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				F.assign(*(S+j+i*n),zero);
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
			for (size_t j=i+1;j<n;++j)
				F.assign(*(L+j+i*n),zero);
		}

    
		//  compute A=LS
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, L, n, S, n, zero, A, n );
		delete[] L;
		delete[] S;


		// compute the determinant of A
		Element det= FFPACK::Det( F, n, n, A, n);
    		delete[] A;


		if (!F.areEqual(det,d))
			ret=false;
	}
  
	commentator.stop(MSG_STATUS (ret), (const char *) 0, "testDet");
    
	return ret;
}

/*
 * Test of the LQUP factorization routine
 */
template <class Field>
static bool testLUdivine (const Field& F, size_t m, size_t n, int iterations) 
{

	typedef typename Field::Element                  Element;
	typedef typename Field::RandIter                RandIter;

	//Commentator commentator;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start (pretty("Testing LQUP factorization"),"testLQUP",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element one,zero;
	F.init(one,1UL);
	F.init(zero,0UL);
  
	bool ret = true;
	
	for (int k=0;k<iterations;++k) {
    
		commentator.progress(k);    
		
		
		Element * A = new Element[m*n];
		Element * B = new Element[m*m];
		Element * C = new Element[m*n];
		

		// Create B a random matrix of rank n/2
		for (size_t j=0;j<m;++j)
			if ( j % 2 ){
				for (size_t i=0;i<j;++i)
					F.assign (*(B+i*m+j),zero);
				for (size_t i=j;i<m;++i)
					Gn.random (*(B+i*m+j));
			} else
				for (size_t i=0;i<m;++i)
					F.assign (*(B+i*m+j), zero);
		// Create C a random matrix of rank n/2
		for (size_t i = 0; i < m; ++i)
			if ( i % 2 ){
				for (size_t j = 0; j < i; ++j)
					F.assign (*(C+i*n+j),zero);
				for (size_t j = i; j < n; ++j)
					Gn.random (*(C+i*n+j));
			} else
				for (size_t j = 0; j < n; ++j)
					F.assign (*(C+i*n+j),zero);
		
		// A = B*C
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, m,
			      one, B, m, C, n, zero, A, n );
		delete[] B;
		delete[] C;

		Element * Abis = new Element[m*n];
		for (size_t i=0; i<m*n; ++i)
			*(Abis+i) = *(A+i);

		size_t * P = new size_t[n];
		size_t * Q = new size_t[m];
		
// 		write_field (F, cerr<<"A="<<endl, A, m, n, n);
		size_t r = FFPACK::LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans,
					     m, n, A, n, P, Q, FFPACK::FfpackLQUP);
// 		write_field (F, cerr<<"LQUP(A)="<<endl, A, m, n, n);

		Element * L = new Element[m*m];
		Element * U = new Element[m*n];
	
		for (size_t i=0; i<m; ++i){
			for (size_t j = 0; j < i; ++j)
				if (j<n)
					F.assign (*(L+i*m+j), *(A+i*n+j) );
				else
					F.assign (*(L+i*m+j), zero );
					    
			for (size_t j = i; j < m; ++j)
				F.assign (*(L+i*m+j), zero );
		}
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m, 
				  0, r, L, m, Q);
		for (size_t i=0; i<m; ++i)
			F.assign( *(L+i*m+i), one);
		for (size_t i=0; i<m; ++i){
			for (size_t j=0; j<i; ++j)
				if (j<n)
					F.assign( *(U+i*n+j), zero );
			for (size_t j=i; j<n; ++j)
				F.assign( *(U+i*n+j), *(A+i*n+j) );
		}
// 		write_field (F, cerr<<"L="<<endl, L, m, m, m);
// 		write_field (F, cerr<<"U"<<endl, U, m, n, n);
		// C = U*P
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, m, 
				  0, r, U, n, P);
		//		write_field (F, cerr<<"UP"<<endl, U, m, n, n);
		// C = Q*C
		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, n, 
				  0, r, U, n, Q);
		//		write_field (F, cerr<<"QUP"<<endl, U, m, n, n);

		delete[] P;
		delete[] Q;
		// A = L*C
		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m, n, m,
			      one, L, m, U, n, zero, A, n );
		delete[] L;
		delete[] U;
		
		for (size_t i = 0; i < m;++i)
			for (size_t j = 0; j < n; ++j)
				if (!F.areEqual( *(A+i*n+j), *(Abis+i*n+j))){
					ret=false;
				}
		delete[] A;
		delete[] Abis;
		if (!ret){
// 			write_field(F,std::cerr<<"A="<<endl,A,m,n,n);
// 			write_field(F,std::cerr<<"Abis="<<endl,Abis,m,n,n);
// 			write_field(F,std::cerr<<"B="<<endl,B,m,n,n);
// 			write_field(F,std::cerr<<"C="<<endl,C,m,n,n);
		}


	}

	commentator.stop(MSG_STATUS (ret), (const char *) 0, "testLQUP");
    
	return ret;
}

template <class Field>
static bool testMinPoly (const Field& F, size_t n, int iterations) 
{
	typedef typename Field::Element                  Element;
	typedef typename Field::RandIter                RandIter;
	typedef vector<Element>                       Polynomial;
	//Commentator commentator;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start (pretty("Testing minpoly"),"testMinPoly",iterations);
	Element tmp, one, zero,mone;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	F.init(one, 1UL);
	F.init(zero, 0UL);
	F.neg(mone, one);
	//F.neg( mone, one);
	bool ret = true;
			
	for (int k=0;k<iterations;++k) {
    
		commentator.progress(k);    
		
		Element * A = new Element[n*n];
		Element * X = new Element[n*(n+1)];
		size_t * Perm = new size_t[n];
			
		Polynomial P;
		// Test MinPoly(In) = X-1
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<n;++j)
				F.assign(*(A+j+i*n),zero);
			F.assign(*(A+i*(n+1)),one);
		}
		
		FFPACK::MinPoly( F, P, n, A, n, X, n, Perm );
		
		if ( P.size() !=2 )
			ret = false;
		if ( !F.areEqual(P[0], mone) )
			ret = false;
		if ( !F.areEqual(P[1], one) )
			ret = false;
		if(!ret) cerr<<"MinP(In)!=X-1"<<endl;
		// Test MinPoly(a*In) = X-a
		G.random(tmp);
		
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<n;++j)
				F.assign(*(A+j+i*n),zero);
			
			F.assign(*(A+i*(n+1)),tmp);
		}
		
		F.negin(tmp);
		
		for (size_t i=0; i<n ;++i) 
			Perm[i]=0;
		FFPACK::MinPoly( F, P, n, A, n, X, n, Perm );
	
		if ( P.size() !=2 )
			ret = false;
		if ( !F.areEqual(P[0], tmp) ){
			cerr<<"P[0]="<<P[0]<<"tmp="<<tmp<<endl;
			ret = false;
		}
		if ( !F.areEqual(P[1], one) ){
			cerr<<"P[1]="<<P[1]<<endl;

			ret = false;
		}
		if(!ret) cerr<<"MinP(aIn)!=X-a"<<endl;
		
// 		for (size_t i=0;i<n-1;++i){
// 			for (size_t j=0; j<n; ++j)
// 				F.assign(*(A+i*n+j),zero);
// 			F.assign(*(A+i*n+i+1),one);
// 		}
// 		for (size_t j=0;j<n;++j)
// 			F.assign(*(A+(n-1)*n+j),zero);
		
// 		for (size_t i=0; i<n ;++i) 
// 			Perm[i]=0;
// 		FFPACK::MinPoly( F, P, n, A, n, X, n, Perm );
	
// 		if ( P.size() !=n+1 )
// 			ret = false;
// 		for (size_t i=0; i<n;++i)
// 			if ( !F.areEqual(P[i], zero) )
// 				ret = false;
// 		if ( !F.areEqual(P[n], one) )
// 			ret = false;
// 		if(!ret) cerr<<"MinP(J)!=X^n"<<endl;
		delete[] A;
		delete[] X;
		delete[] Perm;
	}

	commentator.stop(MSG_STATUS (ret), (const char *) 0, "testMinPoly");
	
	return ret;
}

template <class Field>
static bool testCharPoly (const Field& F, size_t n, int iterations) 
{
	typedef typename Field::Element                  Element;
	typedef typename Field::RandIter                RandIter;
	typedef vector<Element>                       Polynomial;
	//Commentator commentator;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start (pretty("Testing charpoly"),"testCharPoly",iterations);
	Element tmp, one, zero,mone;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	F.init(one, 1UL);
	F.init(zero, 0UL);
	F.neg(mone, one);
	bool ret = true;
	
	for (int k=0;k<iterations;++k) {
    
		commentator.progress(k);    

		Element * A = new Element[n*n];
		
		std::list<Polynomial> P;
		// Test CharPoly(In) = (X-1)^n
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				F.assign(*(A+i*n+j),zero);
			F.assign(*(A+i*(n+1)),one);
			for (size_t j=i+1;j<n;++j)
				 F.assign(*(A+i*n+j),zero);
		}
		P.clear();

		FFPACK::CharPoly (F, P, n, A, n);
	

		typename list<Polynomial>::const_iterator P_it = P.begin();
		while (P_it != P.end()){
			if ( P_it->size() !=2 ){
				ret = false;
			}
			if ( !F.areEqual(P_it->operator[](0), mone) ){
				ret = false;
			}
			if ( !F.areEqual(P_it->operator[](1), one) ){
				ret = false;
			}
			
			P_it++;
		}
		
		// Test CharPoly(a*In) = X-a
		G.random(tmp);
		
		for (size_t i=0;i<n;++i)
			F.assign(*(A+i*(n+1)),tmp);
		F.negin(tmp);
		P.clear();

		FFPACK::CharPoly( F, P, n, A, n );
	
		P_it = P.begin();

		while (P_it != P.end()){
			if ( P_it->size() !=2 )
				ret = false;
			if ( !F.areEqual(P_it->operator[](0), tmp) )
				ret = false;
			if ( !F.areEqual(P_it->operator[](1), one) )
			ret = false;
			P_it++;
		}
		delete[] A;
	}

	commentator.stop(MSG_STATUS (ret), (const char *) 0, "testCharPoly");
	
	return ret;
}

template <class Field>
static bool testInv (const Field& F,size_t n, int iterations) 
{

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
  
	//Commentator commentator;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start (pretty("Testing inverse"),"testInv",iterations);
	
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element zero,one;
	F.init(one,1UL);
	F.init(zero,0UL);
	
	bool ret = true;
	
	Element * Id = new Element[n*n];
	for (size_t i=0;i<n;++i){
		for (size_t j=0;j<n;++j)
			F.assign(*(Id+j+i*n),zero);
		F.assign( *(Id+i*(n+1)),one);
	}
	for (int k=0;k<iterations;++k) {
    
		commentator.progress(k);
   

		Element * A = new Element[n*n];
		Element * L = new Element[n*n];
		Element * S = new Element[n*n];

		

		// create S as an upper triangular matrix of full rank 
		// with nonzero random diagonal's element 
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)      
				F.assign(*(S+i*n+j),zero);
			Gn.random(*(S+i*(n+1)));
			for (size_t j=i+1;j<n;++j)      
				G.random(*(S+i*n+j));
		}
     
		// create L as a lower triangular matrix 
		// with only 1's on diagonal 
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				G.random(*(L+i*n+j));
			F.assign(*(L+i*(n+1)),one);
			for (size_t j=i+1;j<n;++j)
				F.assign(*(L+i*n+j),zero);

		}
        
		//  compute A=LS

		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, L, n, S, n, zero, A, n );

		Element * Ab = new Element[n*n];
		Element * invA = new Element[n*n];
		for (size_t i = 0; i<n*n; ++i)
			F.assign( *(Ab+i), *(A+i) );
		// compute the inverse of A
		int nullity;
		FFPACK::Invert2 ( F, n, A, n, invA, n, nullity);
		
		// compute Ainv*A and A*Ainv

		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, Ab, n, invA, n, zero, L, n );

		FFLAS::fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, n, n, n,
			      one, invA, n, Ab, n, zero, S, n );
   
		for (size_t i=0;i<n*n;++i)
			if ( !F.areEqual(*(L+i),*(Id+i)) || !F.areEqual(*(S+i),*(Id+i)) ) 
				ret =false;
		if (!ret){
			// write_field (F, cerr<<"ID1="<<endl, L, n,n,n);
// 			write_field (F, cerr<<"ID2="<<endl, S, n,n,n);

		}
		delete[] L;
		delete[] S;
		delete[] A;
		delete[] Ab;
		delete[] invA;
	}
	delete[] Id;
	
	commentator.stop(MSG_STATUS (ret), (const char *) 0, "testInv");
    
	return ret;
}

template <class Field>
static bool testapplyP (const Field& F,size_t n, int iterations) 
{

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
  
	//Commentator commentator;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start (pretty("Testing applyP"),"testapplyP",iterations);
	
	RandIter G(F);
	Element zero,one,tmp,tmp2;
	F.init(one,1UL);
	F.init(zero,0UL);
	
	bool ret = true;
	Field Z2(2);
	RandIter G2(Z2);
	
	for (int k=0;k<iterations;++k) {
		
		commentator.progress(k);
		

		Element * A = new Element[n*n];
		Element * Ab = new Element[n*n];
		size_t * P =new size_t[n];
		
		
		for (size_t i=0; i<n; ++i){
			G.random(tmp);
			if ( Z2.isZero(G2.random(tmp2) ) )
				P[i] = i + ( (size_t) tmp % (n-i) ); 
			else				
				P[i] = i;
		}
		
		// create S as an upper triangular matrix of full rank 
		// with nonzero random diagonal's element 
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<n;++j)      
				G.random(*(A+i*n+j));
		}
		
		for (size_t i = 0; i<n*n; ++i)
			F.assign( *(Ab+i), *(A+i) );
		
		//  compute A=LS

		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasNoTrans, n, 0, n, A, n, P );
		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans, n, 0, n, A, n, P );
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans, n, 0, n, A, n, P );
		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasTrans, n, 0, n, A, n, P );
		
		for (size_t i=0;i<n*n;++i)
			if ( !F.areEqual(*(Ab+i),*(A+i)) ) 
				ret =false;
		delete[] A;
		delete[] Ab;
		delete[] P;
	}
	
	commentator.stop(MSG_STATUS (ret), (const char *) 0, "testApplyP");
    
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

	static size_t n = 130;
	static size_t m = 130;
	static integer q = 65521;
	static int iterations =1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to MxN.",       TYPE_INT,     &n },
		{ 'm', "-m M", "Set dimension of test matrices to MxN.",       TYPE_INT,     &m },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q }, 
		{ 'i', "-i I", "Perform each test for I iterations.",           TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);
    
	Field F (q);

	srand(time (NULL));
	
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.start("ffpack test suite", "ffpack");

 	if (!testLUdivine (F, m,n, iterations)) pass=false;
  	if (!testRank (F, n, iterations))   pass = false;   
  	if (!testDet (F, n, iterations))   pass = false;   
   	if (!testTURBO (F, n, iterations))   pass = false;   
   	if (!testapplyP  (F, n, iterations)) pass = false;
   	if (!testInv  (F, n, iterations)) pass = false;
   	if (!testMinPoly (F,n,iterations)) pass=false;
 	if (!testCharPoly (F,n,iterations)) pass=false;
	
	commentator.stop(MSG_STATUS(pass),"ffpack test suite");
    
	return pass ? 0 : -1;
}














/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
