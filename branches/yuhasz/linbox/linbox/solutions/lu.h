/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/LU.h
 *
 * writtend by zhendong wan <wan@mail.eecis.udel.edu>
 */

#ifndef LU_H
#define LU_H

#include <iostream>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/densesubmatrix.h>

namespace LinBox
{
  /**
   * M <-- LU decomp of M.
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
    void LU(DenseSubMatrix<Field> M);
  
  // M <-- M L^{-1}, where L is unit lower triangular with implicit diagonal.
  template<class Field>
    void LL_MULIN(DenseSubMatrix<Field>& M, const DenseSubMatrix<Field>& L);
  
  // R <-- RU, where U is upper triangular.
  template<class Field>
    void RU_MULIN(DenseSubMatrix<Field>& R, const DenseSubMatrix<Field>& U);
  
  template<class Field>
    void  AXMYIN(DenseSubMatrix<Field>&, const DenseSubMatrix<Field>&, const DenseSubMatrix<Field>&);
  
  template<class Field>
    void LU(DenseMatrix<Field>& M) 
    {
      linbox_check(M.rowdim()==M.coldim());
      if(M.rowdim()<=1)
	return;
      LU(DenseSubMatrix<Field>(&M,0,M.rowdim()-1,0,M.coldim()-1));      
    }
  
  template<class Field>
    void LU(DenseSubMatrix<Field> M) 
    {
      int dim=M.rowdim();
      if(dim<=1)
	return;
      --dim;
      int HALF=dim/2;
      DenseSubMatrix<Field> M00(M,0,HALF,0,HALF);
      DenseSubMatrix<Field> M01(M,0,HALF,HALF+1,dim);
      DenseSubMatrix<Field> M10(M,HALF+1,dim,0,HALF);
      DenseSubMatrix<Field> M11(M,HALF+1,dim,HALF+1,dim);
      LU(M00);
      LL_MULIN(M01,M00);
      RU_MULIN(M10,M00);
      AXMYIN(M11,M10,M01);
      LU(M11);     
    }

  /**
   *@param L is a lower triangle matrix.
   * All diagonal entries in L are one.
   * M<-(inverse L)M;
   */
  template<class Field>
    void LL_MULIN(DenseSubMatrix<Field>& M, const DenseSubMatrix<Field>& L) 
    {
      DenseSubMatrix<Field>::RowOfColsIterator colp;
      DenseSubMatrix<Field>::ColIterator ep, epo;

      DenseSubMatrix<Field>::ConstColOfRowsIterator crowp;
      DenseSubMatrix<Field>::ConstRowIterator cep;
      
      DenseSubMatrix<Field>::Element e;

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
  
  /**
   *param U is a uuper triangle matrix
   *M<-M(inverse R)
   */
  template<class Field>
    void RU_MULIN(DenseSubMatrix<Field>& M, const DenseSubMatrix<Field>& U)
    {
      DenseSubMatrix<Field>::ColOfRowsIterator rowp;
      DenseSubMatrix<Field>::ConstRowOfColsIterator ccolp;
      DenseSubMatrix<Field>::Element e;
      DenseSubMatrix<Field>::RowIterator ep,epo;
      DenseSubMatrix<Field>::ConstColIterator cep;
      
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

  /**
   *@M1, M2, dense matrix.
   *M<-M-M1*M2
   */
  template<class Field>
    void  AXMYIN(DenseSubMatrix<Field>& M, const DenseSubMatrix<Field>& M1, const DenseSubMatrix<Field>& M2)
    {
      DenseSubMatrix<Field>::RowOfColsIterator colp;
      DenseSubMatrix<Field>::ConstColOfRowsIterator crowp1;
      DenseSubMatrix<Field>::ConstRowOfColsIterator ccolp2; 
      DenseSubMatrix<Field>::Element e;
      DenseSubMatrix<Field>::ColIterator ep;
      DenseSubMatrix<Field>::ConstRowIterator cep1;
      DenseSubMatrix<Field>::ConstColIterator cep2;
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
