
/* tests/test-blas-domain.C
 * Copyright (C) 2004 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 *
 */
 // where is this used?
#define __MINP_CONSTRUCT
#include "linbox/linbox-config.h"
#include <iostream>
#include <string>
#include <linbox/integer.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/field/modular.h>
#include <linbox/randiter/nonzero.h>
#include <linbox/util/commentator.h>
#include <linbox/algorithms/blas-domain.h>

#include <vector>

#include "test-common.h"

using namespace LinBox;

const int maxpretty = 35;

string blank;

const char* pretty(string a) {

	blank = "     " + a;
	int msgsize= maxpretty - blank.size();
	string dot(".");
	for (int i=0;i<msgsize ;++i)
		 blank+=dot;
	 return blank.c_str();
}

#define mycommentator commentator

template <class Field>
static bool testMulAdd (const Field& F, size_t n, int iterations) {

	typedef typename Field::Element     Element;
	typedef typename Field::RandIter   RandIter;
	typedef BlasMatrix<Element>          Matrix;

	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing muladd"),"testMulAdd",iterations);
	
	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	MatrixDomain<Field>      MD(F);
	VectorDomain<Field>      VD(F);
	
	for (int k=0;k<iterations; ++k) {
    
		mycommentator.progress(k);
		Matrix A(n,n),B(n,n),C(n,n),D(n,n),T(n,n),R(n,n);
		std::vector<Element> x(n),y(n),z(n),t(n);

		Element alpha, beta,malpha,tmp;
		

		// Create 3 random n*n matrices
		for (size_t i=0;i<n;++i)
			for (size_t j=0;j<n;++j){
				A.setEntry(i,j,G.random(tmp));
				B.setEntry(i,j,G.random(tmp));
				C.setEntry(i,j,G.random(tmp));
			}

		// Create 2 random vectors
		for (size_t i=0;i<n;++i) {
			G.random(x[i]);
			G.random(y[i]);
		}			

		// create 2 random element
		G.random(alpha);
		G.random(beta);
	
		F.neg(malpha,alpha);

		// compute D = -alpha.(A*C+B*C) + alpha.(A+B)*C
		
		BMD.mul(D,A,C);
		BMD.mul(T,B,C);
		MD.addin(D,T);
		
		MD.add(T,A,B);
		BMD.muladd(R,malpha,D,alpha,T,C);
		
		if (!MD.isZero(R))
			ret=false;

		// compute z = beta.y + alpha.A*x

		BMD.muladd(z,beta,y,alpha,A,x);

		MD.vectorMul(t,A,x);
		for (size_t i=0;i<n;++i){
		  F.mulin(t[i],alpha);
		  F.axpyin(t[i],beta,y[i]);
		}

		if (!VD.areEqual(t,z))
			ret=false;
	}
	
	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testMulAdd");
	
	return ret;
}



/*
 *  Testing the rank of dense matrices using BlasDomain
 *  construct a n*n matrices of rank r and compute the rank
 */
template <class Field>
static bool testRank (const Field& F,size_t n, int iterations) {

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing rank"),"testRank",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element tmp;
	unsigned int r;
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations; ++k) {
    
		mycommentator.progress(k);
		BlasMatrix<Element> A(n,n),S(n,n), L(n,n);
     
		r = rand() % n;
		// create S as an upper triangular matrix with r nonzero rows
		for (size_t i=0;i<r;++i){
			S.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<n;++j)     
				S.setEntry(i,j,G.random(tmp));
		}
     
		// create L as a lower triangular matrix with nonzero elements on the diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,Gn.random(tmp));
		}
     
		//  compute A=LS
		BMD.mul(A,L,S);
     
		// compute the rank of A
		unsigned int rank= BMD.rankin(A);
		if (rank!=r)
			ret=false;
	}
  	
	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testRank");
    
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
  
	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing determinant"),"testDet",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element tmp,One,d;
	F.init(One,1UL);

	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {
		
		mycommentator.progress(k);

		G.random(d);

		BlasMatrix<Element> A(n,n),S(n,n), L(n,n);

		// create S as an upper triangular matrix of full rank 
		// with diagonal's element equal to 1 except the first entry wich equals to d
		for (size_t i=0;i<n;++i){
			S.setEntry(i,i,One);
			for (size_t j=i+1;j<n;++j)      
				S.setEntry(i,j,G.random(tmp));
		}
		S.setEntry(0,0,d);
    
		// create L as a lower triangular matrix with only 1's on diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,One);
		}

    
		//  compute A=LS
		BMD.mul(A,L,S);
    
		// compute the determinant of A
		Element det= BMD.detin(A);

		if (!F.areEqual(det,d))
			ret=false;
	}
  
	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testDet");
    
	return ret;
}

/*
 *  Testing the inverse of dense matrices using BlasDomain
 *  construct a non-singular n*n matrices 
 */
template <class Field>
static bool testInv (const Field& F,size_t n, int iterations) {

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef  BlasMatrix<Element> Matrix;
  
	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing inverse"),"testInv",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element One,tmp;
	F.init(One,1UL);

	bool ret = true;
	MatrixDomain<Field> MD(F);
	BlasMatrixDomain<Field> BMD(F);

	Matrix Id(n,n);
	for (size_t i=0;i<n;++i)
		Id.setEntry(i,i,One);
  
	for (int k=0;k<iterations;++k) {
    
		mycommentator.progress(k);
   

		Matrix A(n,n),S(n,n), L(n,n), invA(n,n);

		// create S as an upper triangular matrix of full rank 
		// with nonzero random diagonal's element 
		for (size_t i=0;i<n;++i){
			S.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<n;++j)      
				S.setEntry(i,j,G.random(tmp));
		}
     
		// create L as a lower triangular matrix 
		// with only 1's on diagonal 
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,One);
		}
        
		//  compute A=LS
		BMD.mul(A,L,S);
      
		// compute the inverse of A
		BMD.inv(invA,A);
        
		// compute Ainv*A and A*Ainv
		BMD.mul(L,invA,A);
		BMD.mul(S,A,invA);
   
		if (!MD.areEqual(L,Id) || !MD.areEqual(S,Id))
			ret=false;
	}
  
	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testInv");
    
	return ret;
}


/*
 * Test resolution of linear system with a triangular matrix 
 */
template <class Field>
static bool testTriangularSolve (const Field& F, size_t m, size_t n, int iterations) {

	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Element>                       Matrix;
	typedef TriangularBlasMatrix<Element>   TriangularMatrix;
	typedef typename Field::RandIter                RandIter;

	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing triangular solver"),"testTriangularSolve",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element One,tmp;
	F.init(One,1UL);
  
	bool ret = true;
	MatrixDomain<Field> MD(F);
	VectorDomain<Field>  VD(F);
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {
    
		mycommentator.progress(k);    
    
		Matrix Al(m,m),Au(m,m);
		Matrix X(m,n), B(m,n), C(m,n);  

		std::vector<Element> b(m),x(m),c(m);

		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		//create random vector b
		for( size_t i=0;i<m;++i)
			F.init(b[i],G.random(tmp));
    
		// Create Au a random full rank upper triangular matrix
		for (size_t i=0;i<m;++i){
			Au.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<m;++j)
				Au.setEntry(i,j,G.random(tmp));
		}

		// Create Al a random full rank lower triangular matrix
		for (size_t i=0;i<m;++i){      
			for (size_t j=0;j<i;++j)
				Al.setEntry(i,j,G.random(tmp));
			Al.setEntry(i,i,Gn.random(tmp));	  
		}

		// Create 2 trinagular matrix as view of matrix
		TriangularMatrix TAl(Al,BlasTag::low,BlasTag::nonunit), TAu(Au,BlasTag::up,BlasTag::nonunit);
    
		// testing solver with matrix right hand side
		BMD.left_solve(X,TAl,B);    
		BMD.mul(C,Al,X);          
		if (!MD.areEqual(C,B))
			ret=false;
        
		BMD.left_solve(X,TAu,B);
		BMD.mul(C,Au,X);   
		if (!MD.areEqual(C,B))
			ret=false;
        
		// testing solver with matrix left hand side
		BMD.right_solve(X,TAl,B);
		BMD.mul(C,X,Al);    
		if (!MD.areEqual(C,B))
			ret=false;
        
		BMD.right_solve(X,TAu,B);
		BMD.mul(C,X,Au);    
		if (!MD.areEqual(C,B))
			ret=false;
       

		// testing solver with vector right hand side
		BMD.left_solve(x,TAl,b);    
		BMD.mul(c,Al,x);          
		if (!VD.areEqual(c,b))
			ret=false;
        
		BMD.left_solve(x,TAu,b);
		BMD.mul(c,Au,x);   
		if (!VD.areEqual(c,b))
			ret=false;
        
		// testing solver with vector left hand side
		BMD.right_solve(x,TAl,b);
		BMD.mul(c,x,Al);    
		if (!VD.areEqual(c,b))
			ret=false;
        
		BMD.right_solve(x,TAu,b);
		BMD.mul(c,x,Au);    
		if (!VD.areEqual(c,b))
			ret=false;

	}

	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testTriangularSolve");
    
	return ret;
}

/*
 * Test resolution of linear system with a matrix 
 */
template <class Field>
static bool testSolve (const Field& F, size_t m, size_t n, int iterations) {

	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Element>                       Matrix;
	typedef typename Field::RandIter                RandIter;

	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing solver"),"testSolve",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element One,tmp;
	F.init(One,1UL);
  
	bool ret = true;
	MatrixDomain<Field> MD(F);
	VectorDomain<Field>  VD(F);
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {
    
		mycommentator.progress(k);    
    
		Matrix A(m,m),L(m,m),S(m,m);
		Matrix X(m,n), B(m,n), C(m,n);  

		std::vector<Element> b(m),x(m),c(m);

		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<n;++j)
				B.setEntry(i,j,G.random(tmp));

		//create a random vector b
		for( size_t i=0;i<m;++i)
			F.init(b[i],G.random(tmp));
    
   
		// create S as an upper triangular matrix of full rank 
		// with nonzero random diagonal's element 
		for (size_t i=0;i<m;++i){
			S.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<n;++j)      
				S.setEntry(i,j,G.random(tmp));
		}
    
		// create L as a lower triangular matrix 
		// with only 1's on diagonal 
		for (size_t i=0;i<m;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,One);
		}
        
		//  compute A=LS
		BMD.mul(A,L,S);


		// testing solver with matrix right hand side
		BMD.left_solve(X,A,B);    
		BMD.mul(C,A,X);          
		if (!MD.areEqual(C,B))
			ret=false;
        
		BMD.left_solve(X,A,B);
		BMD.mul(C,A,X);   
		if (!MD.areEqual(C,B))
			ret=false;
        
		// testing solver with matrix left hand side
		BMD.right_solve(X,A,B);
		BMD.mul(C,X,A);    
		if (!MD.areEqual(C,B))
			ret=false;
        
		BMD.right_solve(X,A,B);
		BMD.mul(C,X,A);    
		if (!MD.areEqual(C,B))
			ret=false;
       

		// testing solver with vector right hand side
		BMD.left_solve(x,A,b);    
		BMD.mul(c,A,x);          
		if (!VD.areEqual(c,b))
			ret=false;
        
		BMD.left_solve(x,A,b);
		BMD.mul(c,A,x);   
		if (!VD.areEqual(c,b))
			ret=false;
        
		// testing solver with vector left hand side
		BMD.right_solve(x,A,b);
		BMD.mul(c,x,A);    
		if (!VD.areEqual(c,b))
			ret=false;
        
		BMD.right_solve(x,A,b);
		BMD.mul(c,x,A);    
		if (!VD.areEqual(c,b))
			ret=false;

	}

	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testSolve");
    
	return ret;
}

/*
 * Test of the BlasPermutations
 */

template <class Field>
static bool testPermutation (const Field& F, size_t m, int iterations) {

	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Element>                       Matrix;
	typedef typename Field::RandIter                RandIter;

	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing permutations"),"testPermutation",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element One,zero,tmp,tmp2;
	F.init(One,1UL);
	F.init(zero,0UL);
  
	bool ret = true;
	MatrixDomain<Field> MD(F);
	VectorDomain<Field>  VD(F);
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {
    
		mycommentator.progress(k);    
    
		std::vector<size_t> P(m);

		Field Z2(2);
		RandIter G2(Z2);

		for (size_t i=0; i<m; ++i){
			G.random(tmp);
			if ( Z2.isZero(G2.random(tmp2) ) )
				P[i] = i + ( (size_t) tmp % (m-i) ); 
			else				
				P[i] = i;
		}
		
		//std::cerr<<P<<std::endl;
		Matrix A(m,m), Abis(m,m), B(m,m), C(m,m), D(m,m);
		std::vector<Element> a(m),abis(m),b(m),c(m), d(m);
		BlasPermutation Perm(P);
		
		// Create A a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				A.setEntry(i,j,Gn.random(tmp));
		// Create a a random vector
		for (size_t i=0;i<m;++i)
			F.assign(a[i],Gn.random(tmp));

		/*
		 * Test A.P.P^t == A
		 */
		
		// B = A.P
		BMD.mul( B, A, Perm);
		// C = B.P^t 
		BMD.mul( C, B, TransposedBlasMatrix<BlasPermutation>(Perm) );
		// Test C==A
		if (!MD.areEqual(A,C))
			ret=false;
		/*
		 * Test A.P^t.P == A
		 */
		
		// B = A.P^t
		BMD.mul( B, A, TransposedBlasMatrix<BlasPermutation>(Perm));
		// C = B.P 
		BMD.mul( C, B, Perm );
		// Test C==A
		if (!MD.areEqual(A,C))
			ret=false;
		/*
		 * Test P.P^t.A == A
		 */
		
		// B = P.A
		BMD.mul( B, Perm, A);
		// C = P^t.B 
		BMD.mul( C, TransposedBlasMatrix<BlasPermutation>(Perm) , B);
		// Test C==A
		if (!MD.areEqual(A,C))
			ret=false;
		/*
		 * Test P^t.P.A == A
		 */
		
		// B = P^t.A
		BMD.mul( B, TransposedBlasMatrix<BlasPermutation>(Perm), A);
		// C = P.B 
		BMD.mul( C, Perm, B);
		// Test C==A
		if (!MD.areEqual(A,C))
			ret=false;

		/*
		 * Test a.P.P^t == a
		 */
		
		// b = a.P
		BMD.mul( b, a, Perm);
		// c = b.P^t 
		BMD.mul( c, b, TransposedBlasMatrix<BlasPermutation>(Perm) );
		// Test c==a
		if (!VD.areEqual(a,c))
			ret=false;

		/*
		 * Test a.P^t.P == a
		 */
		
		// b = a.P^t
		BMD.mul( b, a, TransposedBlasMatrix<BlasPermutation>(Perm));
		// c = B.P 
		BMD.mul( c, b, Perm );
		// Test c==a
		if (!VD.areEqual(a,c))
			ret=false;
		/*
		 * Test P.P^t.a == a
		 */
		
		// b = P.a
		BMD.mul( b, Perm, a);
		// c = P^t.b 
		BMD.mul( c, TransposedBlasMatrix<BlasPermutation>(Perm) , b);
		// Test c==a
		if (!VD.areEqual(a,c))
			ret=false;

		/*
		 * Test P^t.P.a == a
		 */
		
		// b = P^t.a
		BMD.mul( b, TransposedBlasMatrix<BlasPermutation>(Perm), a);
		// c = P.b 
		BMD.mul( c, Perm, b);
		// Test c==a
		if (!VD.areEqual(a,c))
			ret=false;


		/*
		 * Test P^t.A.(P.A)^-1.B == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));
		
		// Abis = P.A
		BMD.mul( Abis, Perm, A);
		// C = (P.A)^-1.B
		BMD.left_solve( C, Abis, B);
		// D = A.C (= P^-1.B)
		BMD.mul(D, A, C);
		// D = P.D
		BMD.mulin_right( Perm,D);
		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test A.P^t.(A.P)^-1.B == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));
		
		// Abis = A.P
		BMD.mul( Abis, A, Perm);
		// C = (A.P)^-1.B
		BMD.left_solve( C, Abis, B);
		// C = P.C
		BMD.mulin_right( Perm,C);
		// D = A.C (= P^-1.B)
		BMD.mul(D, A, C);

		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test B.P^t.A.(P.A)^-1 == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));
		
		// Abis = P.A
		BMD.mul( Abis, Perm, A);
		// C = B.(P.A)^-1
		BMD.right_solve( C, Abis, B);
		// C = C.P
		BMD.mulin_left( C,Perm);
		// D = C.A (=B)
		BMD.mul(D, C, A);
		if (!MD.areEqual(D,B))
		  ret=false;
		
		/*
		 * Test B.A.P^t.(A.P)^-1 == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));
		
		// Abis = A.P
		BMD.mul( Abis, A, Perm);
		// C = B.(A.P)^-1
		BMD.right_solve( C, Abis, B);
		// D = C.A (= B.P^t)
		BMD.mul(D, C, A);
		// C = C.P
		BMD.mulin_left( D, Perm);
			
		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test P.A.(P^t.A)^-1.B == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));
		
		// Abis = P^t.A
		BMD.mul( Abis, TransposedBlasMatrix<BlasPermutation>(Perm), A);
		// C = (P^t.A)^-1.B
		BMD.left_solve( C, Abis, B);
		// D = A.C (= P.B)
		BMD.mul(D, A, C);
		// D = P^t.D
		BMD.mulin_right( TransposedBlasMatrix<BlasPermutation>(Perm),D);
		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test A.P.(A.P^t)^-1.B == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));
		
		// Abis = A.P^t
		BMD.mul( Abis, A, TransposedBlasMatrix<BlasPermutation>(Perm));
		// C = (A.P^t)^-1.B
		BMD.left_solve( C, Abis, B);
		// C = P^t.C
		BMD.mulin_right( TransposedBlasMatrix<BlasPermutation>(Perm),C);
		// D = A.C (= P.B)
		BMD.mul(D, A, C);

		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test B.P.A.(P^t.A)^-1 == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));
		
		// Abis = P^t.A
		BMD.mul( Abis, TransposedBlasMatrix<BlasPermutation>(Perm), A);
		// C = B.(P^t.A)^-1
		BMD.right_solve( C, Abis, B);
		// C = C.P^t
		BMD.mulin_left( C,TransposedBlasMatrix<BlasPermutation>(Perm));
		// D = C.A (=B)
		BMD.mul(D, C, A);
		if (!MD.areEqual(D,B))
		  ret=false;
		
		/*
		 * Test B.A.P.(A.P^t)^-1 == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));
		
		// Abis = A.P^t
		BMD.mul( Abis, A, TransposedBlasMatrix<BlasPermutation>(Perm));
		// C = B.(A.P^t)^-1
		BMD.right_solve( C, Abis, B);
		// D = C.A (= B.P)
		BMD.mul(D, C, A);
		// C = C.P^t
		BMD.mulin_left( D, TransposedBlasMatrix<BlasPermutation>(Perm));
			
		if (!MD.areEqual(D,B))
			ret=false;
	}
	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testPermutation");
	
	return ret;
}

/*
 * Test of the LQUPMatrix class
 */
template <class Field>
static bool testLQUP (const Field& F, size_t m, size_t n, int iterations) {

	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Element>                       Matrix;
	typedef typename Field::RandIter                RandIter;

	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing LQUP factorization"),"testLQUP",iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	Element One,zero,tmp;
	F.init(One,1UL);
	F.init(zero,0UL);
  
	bool ret = true;
	MatrixDomain<Field> MD(F);
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {
    
		mycommentator.progress(k);    
    
		Matrix A(m,n), Abis(m,n), B(m,m), C(m,n);
		

		// Create B a random matrix of rank n/2
		for (size_t j=0;j<m;++j)
			if ( j % 2 )
				for (size_t i=0;i<m;++i)
				  B.setEntry(i,j,G.random(tmp));
			else
			  for (size_t i=0;i<m;++i)
			    B.setEntry(i,j,zero);
		// Create C a random matrix of rank n/2
		for (size_t i=0;i<m;++i)
			if ( i % 2 )
				for (size_t j=0;j<n;++j)
					C.setEntry(i,j,G.random(tmp));
			else
				for (size_t j=0;j<n;++j)
					C.setEntry(i,j,zero);
		
		// A = B*C
		BMD.mul(A, B, C);
		
		Abis = A;

		LQUPMatrix<Field> X(F,A);
		
		TriangularBlasMatrix<Element> L(m,m,BlasTag::low,BlasTag::unit);
		TriangularBlasMatrix<Element> U(m,n,BlasTag::up,BlasTag::nonunit);
		X.getL(L);
		X.getU(U);
		BlasPermutation P,Q;
		P=X.getP();

		Q=X.getQ();
		
		// C = U*P
		BMD.mul( C, U, P);
		// C = Q*C
		BMD.mulin_right( Q, C);
		// A = L*C
		BMD.mul( A, L, C);
		
		if (!MD.areEqual(A,Abis))
			ret=false;

		// Second pass
		// A = B*C
		BMD.mul(A, B, C);
		
		Abis = A;

		LQUPMatrix<Field> Y(F,A);
		
		TriangularBlasMatrix<Element> L2(m,m,BlasTag::low,BlasTag::unit);
		TriangularBlasMatrix<Element> U2(m,n,BlasTag::up,BlasTag::nonunit);
		Y.getL(L2);
		Y.getU(U2);
		P=Y.getP();

		Q=Y.getQ();
		
		// C = Q*U2
		BMD.mul( C,Q,U2);
		// C = Q*C
		BMD.mulin_left(  C,P);
		// A = L*C
		BMD.mul( A, L2, C);
		
		if (!MD.areEqual(A,Abis))
			ret=false;
	}

	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testLQUP");
    
	return ret;
}

template <class Field>
static bool testMinPoly (const Field& F, size_t n, int iterations) {
	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Element>                       Matrix;
	typedef typename Field::RandIter                RandIter;
	typedef vector<Element>                       Polynomial;

	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing minpoly"),"testMinPoly",iterations);
	Element tmp, one, zero,mone;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	F.init(one, 1UL);
	F.init(zero, 0UL);
	F.neg(mone, one);
	//F.neg( mone, one);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
		
	for (int k=0;k<iterations;++k) {
    
		mycommentator.progress(k);    

		Matrix A(n,n);
		Polynomial P;
		// Test MinPoly(In) = X-1
		for (size_t i=0;i<n;++i)
			A.setEntry(i,i,one);
		
		BMD.minpoly( P, A );
		
		if ( P.size() !=2 )
			ret = false;
		if ( !F.areEqual(P[0], mone) )
			ret = false;
		if ( !F.areEqual(P[1], one) )
			ret = false;

		// Test MinPoly(a*In) = X-a
		G.random(tmp);
		
		for (size_t i=0;i<n;++i)
			A.setEntry(i,i,tmp);
		F.negin(tmp);
		BMD.minpoly( P, A );
		if ( P.size() !=2 )
			ret = false;
		if ( !F.areEqual(P[0], tmp) )
			ret = false;
		if ( !F.areEqual(P[1], one) )
			ret = false;

		for (size_t i=0;i<n-1;++i){
			for (size_t j=0; j<i+1; ++j)
				A.setEntry(i,j,zero);
			A.setEntry(i,i+1,one);
			if (i<n-2)
				for (size_t j=i+2; j<n; ++j)
					A.setEntry(i,j,zero);
		}
		for (size_t j=0;j<n;++j)
			A.setEntry(n-1,j,zero);
		
		BMD.minpoly( P, A );
		if ( P.size() !=n+1 )
			ret = false;
		for (size_t i=0; i<n;++i)
			if ( !F.areEqual(P[i], zero) )
				ret = false;
		if ( !F.areEqual(P[n], one) )
			ret = false;

	}

	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testMinPoly");
	
	return ret;
}

template <class Field>
static bool testCharPoly (const Field& F, size_t n, int iterations) {
	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Element>                       Matrix;
	typedef typename Field::RandIter                RandIter;
	typedef vector<Element>                       Polynomial;

	//Commentator mycommentator (std::cout);
	//mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator.start (pretty("Testing charpoly"),"testCharPoly",iterations);
	Element tmp, one, zero,mone;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	F.init(one, 1UL);
	F.init(zero, 0UL);
	F.neg(mone, one);
	//F.neg( mone, one);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	
	for (int k=0;k<iterations;++k) {
    
		mycommentator.progress(k);    

		Matrix A(n,n);
		list<Polynomial> P;
		// Test CharPoly(In) = (X-1)^n
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				A.setEntry(i,j,zero);
			A.setEntry(i,i,one);
			for (size_t j=i+1;j<n;++j)
				A.setEntry(i,j,zero);
		}
		P.clear();
		BMD.charpoly( P, A );

		typename list<Polynomial>::const_iterator P_it = P.begin();
		while (P_it != P.end()){
			if ( P_it->size() !=2 )
				ret = false;
			if ( !F.areEqual(P_it->operator[](0), mone) )
				ret = false;
			if ( !F.areEqual(P_it->operator[](1), one) )
				ret = false;
			
			P_it++;
		}
		
		// Test CharPoly(a*In) = X-a
		G.random(tmp);
		
		for (size_t i=0;i<n;++i)
			A.setEntry(i,i,tmp);
		F.negin(tmp);
		P.clear();
		BMD.charpoly( P, A );
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
	}

	mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testCharPoly");
	
	return ret;
}

template<class T, template <class T> class Container>
std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
	for(typename Container<T>::const_iterator refs =  C.begin();
	    refs != C.end() ;
	    ++refs )
		o << (*refs) << " " ;
	return o << std::endl;
}

int main(int argc, char **argv) {
	commentator.setBriefReportStream (cout);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (2);
	commentator.getMessageClass (BRIEF_REPORT).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	static size_t n = 400;
	static integer q = 1000003U;
	static int iterations = 3;
	
    static Argument args[] = {
        { 'n', "-n N", "Set dimension of test matrices to NxN", TYPE_INT,     &n },
        { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1]", TYPE_INTEGER, &q },
        { 'i', "-i I", "Perform each test for I iterations",    TYPE_INT,     &iterations },
		{ '\0' }
    };

	parseArguments (argc, argv, args);

	commentator.setBriefReportStream (cout);
	typedef Modular<double> Field;
	//typedef Modular<int> Field;
	//typedef Modular<float> Field;
        
	Field F (q);
    
	bool pass = true;

	srand (time (NULL));

    
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (20);
	//commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_IMPORTANT);

	commentator.start("BlasMatrixDomain benchmark", "BlasMatrixDomain");


	if (!testMulAdd (F,n,iterations)) pass=false;
 	if (!testRank (F, n, iterations))   pass = false;   
 	if (!testDet  (F, n, iterations)) pass = false;
 	if (!testInv  (F, n, iterations)) pass = false;
 	if (!testTriangularSolve (F,n,n,iterations)) pass=false;
 	if (!testSolve (F,n,n,iterations)) pass=false;
 	if (!testPermutation (F,n,iterations)) pass=false;
 	if (!testLQUP (F,n,n,iterations)) pass=false;
 	if (!testMinPoly (F,n,iterations)) pass=false;
	if (!testCharPoly (F,n,iterations)) pass=false;
	
	commentator.stop("done", 0, "BlasMatrixDomain"); 
	return pass ? 0 : -1;
}
  
  
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
