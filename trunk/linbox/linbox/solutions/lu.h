/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/LU.h
 *
 * writtend by zhendong wan <wan@mail.eecis.udel.edu>
 */

#ifndef LU_H
#define LU_H

#include <iostream>
#include <linbox/blackbox/dense.h>
#include <linbox/matrix/dense-submatrix.h>

namespace LinBox
{
  /**
   * M <-- the LU decomposition of M.
   *
   *@param M is a dense matrix on input.  
   * Every leading principal minor must be nonzero.
   * M is modified to represent
   * the LU decomposition of the original M.
   * L is unit lower triangular and occupies the strictly lower triangular
   * part of M.  The diagonal one's are implicit.
   * U is upper triangular and occupies the rest of M.
   */
  template<class Field>
    void LU(DenseMatrix<Field>& M);
  
  // M <-- LU decomp of M (just as for DenseMatrix argument).
  template<class Field>
    void LU(DenseSubmatrix<typename Field::Element>& M);
  
  // M <-- M L^{-1}, where L is unit lower triangular with implicit diagonal.
  template<class Field>
    void LL_MULIN(DenseSubmatrix<typename Field::Element>& M, const DenseSubmatrix<typename Field::Element>& L);
  
  // R <-- RU, where U is upper triangular.
  template<class Field>
    void RU_MULIN(DenseSubmatrix<typename Field::Element>& R, const DenseSubmatrix<typename Field::Element>& U);
  
  template<class Field>
    void  AXMYIN(DenseSubmatrix<typename Field::Element>&, const DenseSubmatrix<typename Field::Element>&, const DenseSubmatrix<typename Field::Element>&);
  
  template<class Field>
    void LU(DenseMatrix<Field>& M) 
    {
      linbox_check(M.rowdim()==M.coldim());
      if(M.rowdim()<=1)
	return;
      LU(DenseSubmatrix<typename Field::Element>(M,0,M.rowdim()-1,0,M.coldim()-1));      
    }
  
  template<class Field>
    void LU(DenseSubmatrix<typename Field::Element> M) 
    {
      int dim=M.rowdim();
      if(dim<=1)
	return;
      --dim;
      int HALF=dim/2;
      DenseSubmatrix<typename Field::Element> M00(M,0,HALF,0,HALF);
      DenseSubmatrix<typename Field::Element> M01(M,0,HALF,HALF+1,dim);
      DenseSubmatrix<typename Field::Element> M10(M,HALF+1,dim,0,HALF);
      DenseSubmatrix<typename Field::Element> M11(M,HALF+1,dim,HALF+1,dim);
      LU(M00);
      LL_MULIN(M01,M00);
      RU_MULIN(M10,M00);
      AXMYIN(M11,M10,M01);
      LU(M11);     
    }

  /*
   *@param L is a lower triangle matrix.
   * All diagonal entries in L are one.
   * $ M <- L^{-1} M$
   */
  template<class Field>
    void LL_MULIN(DenseSubmatrix<typename Field::Element>& M, const DenseSubmatrix<typename Field::Element>& L) 
    {
      typedef DenseSubmatrix<typename Field::Element> Matrix;
      typename Matrix::ColIterator colp;
      typename Matrix::ColIterator ep, epo;

      typename Matrix::ConstRowIterator crowp;
      typename Matrix::ConstRowIterator cep;
      
      typename Matrix::Element e;

      for(colp=M.rowOfColsBegin();colp!=M.rowOfColsEnd();++colp)
	{
	  crowp=L.colOfRowsBegin();
	  for(ep=colp->begin();ep!=colp->end();++ep,++crowp)
	    for(cep= crowp->begin(),epo=colp->begin();epo!=ep;++cep,++epo)
	      {
		M.field().mul(e,*epo,*cep);
		M.field().subin(*ep,e);
	      }
	}
    }
  
  /*
   *param U is an upper triangle matrix.
   *$M <- M U^{-1}$
   */
  template<class Field>
    void RU_MULIN(DenseSubmatrix<typename Field::Element>& M, const DenseSubmatrix<typename Field::Element>& U)
    {
      typedef DenseSubmatrix<typename Field::Element> Matrix;
      typename Matrix::ColOfRowsIterator rowp;
      typename Matrix::ConstRowOfColsIterator ccolp;
      typename Matrix::Element e;
      typename Matrix::RowIterator ep,epo;
      typename Matrix::ConstColIterator cep;
      
      for(rowp=M.colOfRowsBegin();rowp!=M.colOfRowsEnd();++rowp)
	{
	  ccolp=U.rowOfColsBegin();
	  for(ep=rowp->begin();ep!=rowp->end();++ep,++ccolp)
	    {
	      for(cep= ccolp->begin(),epo=rowp->begin();epo!=ep;++cep,++epo)
		{
		  M.field().mul(e,*epo,*cep);
		  M.field().subin(*ep,e);
		}
	      M.field().divin(*ep,*cep);
	    }
	}
    }

  /*
   *@M1, M2, dense matrix.
   *M<-M-M1*M2
   */
  template<class Field>
    void  AXMYIN(DenseSubmatrix<typename Field::Element>& M, const DenseSubmatrix<typename Field::Element>& M1, const DenseSubmatrix<typename Field::Element>& M2)
    {
      typedef DenseSubmatrix<typename Field::Element> Matrix;
      typename Matrix::RowOfColsIterator colp;
      typename Matrix::ConstColOfRowsIterator crowp1;
      typename Matrix::ConstRowOfColsIterator ccolp2; 
      typename Matrix::Element e;
      typename Matrix::ColIterator ep;
      typename Matrix::ConstRowIterator cep1;
      typename Matrix::ConstColIterator cep2;
      for(colp=M.rowOfColsBegin(),ccolp2=M2.rowOfColsBegin();colp!=M.rowOfColsEnd();++colp,++ccolp2)
	for(ep=colp->begin(),crowp1=M1.colOfRowsBegin();ep!=colp->end();++ep,++crowp1)
	  {
	    M.field().init(e,0);
	    for(cep1=crowp1->begin(),cep2=ccolp2->begin();cep1!=crowp1->end();++cep1,++cep2)
	      M.field().axpyin(e,*cep1,*cep2);
	    M.field().subin(*ep,e);
	  }
    }
}

#endif
