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

  Matrix M(m,n);
  typename Matrix::RawIterator iter;
  for (iter=M.rawBegin();iter!=M.rawEnd();iter++)
    G.random(*iter);

  MatrixDomain<Field> MD(F);


  lsp<Field,Matrix> LSP(F,M);
  LSP.compute();

  Matrix LL(LSP.get_L());
  Matrix SS(LSP.get_S());  
  Matrix PP(n,n);
  for (int i=0;i<n;i++)
    PP.setEntry(i,LSP.get_P()[i],One);

  
  Matrix MM(m,n);
  MD.mul(MM,LL,SS);
  MD.mulin(MM,PP);

  
  

  if (MD.areEqual(M,MM))
	  return true;
  else {
	  ostream &report = commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	  report<<"LSP is not correct for matrix of size ("<<m<<"*"<<n<<") over GF("<<p<<")"<<endl;
	  return false;
  }
  
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
    
    integer prime [5] = {5,17,1009,30011,65521};   
    int size [9]     = {10,20,50,80,100,150,200,250,300};
    
    LinBox::commentator.start("Testing LSP with Dense matrices","",50); 

    for (int i=0;(i<5)&pass;i++) 
      for (int j=0;(j<9)&pass;j++){
	LinBox::commentator.progress();
	if (! lsp_testing<Field,Matrix> (prime[i],size[j],size[j]))
	  pass=false;
      }

    LinBox::commentator.stop (MSG_STATUS (pass));

    return 0;
  }
  
