/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#include "linbox-config.h"

#include <iostream>

#include <linbox/integer.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/field/givaro-zpz.h>
#include <linbox/field/modular.h>
#include <linbox/randiter/nonzero.h>
#include <linbox/util/commentator.h>
#include <linbox/algorithms/blas-domain.h>


#include <vector>

#include "test-common.h"

using namespace LinBox;


/*
 *  Testing the rank of dense matrices using BlasDomain
 *  construct a n*n matrices of rank r and compute the rank
 */
template <class Field>
static bool testRank (const Field& F,size_t n, int iterations) {

  typedef typename Field::Element Element;
  typedef typename Field::RandIter RandIter;

  Commentator mycommentator;
  mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
  mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
  mycommentator.start ("     Testing rank","",iterations);

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
  
  Commentator mycommentator;
  mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
  mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
  mycommentator.start ("     Testing determinant","",iterations);

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
  
  Commentator mycommentator;
  mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
  mycommentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
  mycommentator.start ("     Testing inverse","",iterations);

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
  
  mycommentator.stop(MSG_STATUS (ret), (const char *) 0, "testDet");
    
  return ret;
}


  int main(int argc, char **argv) {

    typedef Modular<double> Field;
    

    bool pass = true;

    static size_t n = 256;
    static integer q = 101U;
    static int iterations =10;
    
    static Argument args[] = {
      { 'n', "-n N", "Set dimension of test matrices to NxN (default 256)",       TYPE_INT,     &n },
      { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q }, 
      { 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
    };

    parseArguments (argc, argv, args);
    
    Field F (q);
    

    srand (time (NULL));

    
    commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
    commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
    commentator.start("BlasMatrixDomain Test suite");
    std::cerr<<endl<<endl;

    if (!testRank (F, n, iterations))   pass = false;   
    if (!testDet  (F, n, iterations)) pass = false;
    if (!testInv  (F, n, iterations)) pass = false;
    
	
    std::cerr<<"\nBlasMatrixDomain Test suite...";
    commentator.stop(MSG_STATUS(pass),"BlasMatrixDomain Test suite");
    
    return pass ? 0 : -1;
  }
  
  
