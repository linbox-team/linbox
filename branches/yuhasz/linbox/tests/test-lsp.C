/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-lsp.C
 * Copyright (C) 2003 Pascal Giorgi
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

#include "linbox/integer.h"
#include "linbox/matrix/dense.h"
#include "linbox/field/givaro-zpz.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/util/commentator.h"
#include "linbox/algorithms/lsp.h"
#include "linbox/fflas/fflas.h"

#include <vector>

#include "test-common.h"

using namespace LinBox;



template <class Field,class Matrix>
bool lsp_testing (integer p, int m,int n) {

  typedef typename Field::Element Element;  
  typedef typename Field::RandIter RandIter;

  Field F(p);

  Element One;
  F.init(One,1UL);
  RandIter G(F);

  LinBox::Timer timer; timer.clear();
  Matrix M(m,n);
  Matrix L(m,m);
  Matrix S(m,n);
  std::vector<int> P_lapackstyle;
  //  Permutation<std::vector<int> > Pbb(n);
  typename Matrix::RawIterator iter;
  timer.start();
  for (iter=M.rawBegin();iter!=M.rawEnd();iter++)
    G.random(*iter);
  timer.stop();
  std::cerr << "Random matrix  --->  " << timer << std::endl;

  MatrixDomain<Field> MD(F);


  lsp<Field> LSP(F);
  timer.start();
  LSP.compute( M,L,S,P_lapackstyle);
  timer.stop();
  std::cerr << "P: " << p << ", m: " << m << ", n: " << n << "  --->  " << timer << std::endl;

  //Matrix LL(LSP.get_L());
  //Matrix SS(LSP.get_S());  
  Matrix PP(n,n);
  std::vector<int> P(n);
  int tmp;
  for (int i=0;i<n;i++)
      P[i] = i;
  for (int i=0;i<n;i++)
      if ( P_lapackstyle[i] != i ) {
  	  tmp = P[i];
  	  P[i] = P[P_lapackstyle[i]];
  	  P[P_lapackstyle[i]] =  tmp;
      }
  for (int i=0;i<n;i++){
      PP.setEntry(i,P[i],One);
  }

  
  Matrix MM(m,n);

//  MD.mul(MM,L,S);
//x  MD.mulin(MM,PP);

  Element Zero;
  F.init(Zero,0);
  FFLAS::fgemm (F,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,m,One,L.FullIterator(),m,S.FullIterator(),n,Zero,MM.FullIterator(),n);
  FFLAS::fgemm (F,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,n,One,MM.FullIterator(),n,PP.FullIterator(),n,Zero,MM.FullIterator(),n);
  
  

  if (MD.areEqual(M,MM)){
	  return true;
  }
  else {
	  Matrix Diff(m,n);
	  MD.sub(Diff,M, MM);
	  Diff.write(cerr,F);
          std::cerr << P_lapackstyle << std::endl;
	  ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	  report<<"LSP is not correct for matrix of size ("<<m<<"*"<<n<<") over GF("<<p<<")"<<endl;
          return false;
  }
 
return true; 
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

  static Argument args[] = {};
  parseArguments (argc, argv, args);
  
  commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (-1);
  commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);


    typedef GivaroZpz<Std32> Field;
    typedef Field::Element Element;
    typedef DenseMatrixBase<Element> Matrix;
    bool pass =true;
    
    //integer prime [5] = {5,17,1009,30011,65521};   
    integer prime [2] = {19,65521};   
    //int size [9]     = {10,20,50,80,100,150,200,250,300};
    int size [4] = {50,100,400,800};
    // int size [4] = {100,200,400,600};
    LinBox::commentator.start("Testing LSP with Dense matrices");
    std::cerr<<endl;
    for (int i=0;(i<2)&pass;i++) 
      for (int j=0;(j<4)&pass;j++){
	      //LinBox::commentator.progress();
	if (! lsp_testing<Field,LinBox::DenseMatrixBase<Element> > (prime[i],size[j],size[j]))
	  pass=false;
      }

    LinBox::commentator.stop (MSG_STATUS (pass));

    return 0;
  }
  
