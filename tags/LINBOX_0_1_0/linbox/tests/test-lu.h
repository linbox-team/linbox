#include <linbox/blackbox/dense-matrix1.h>
#include <linbox/util/debug.h>

namespace LinBox
{
  /** test if M==L*U;
   * where L is lower traingle matrix, whose diagonal entries are 1.
   * U is an upper triangle matrix.
   */
  template<class Field>
    bool LU_MUL_TEST(const DenseMatrix<Field>& M, const DenseMatrix<Field>& L,const DenseMatrix<Field>& U)
    {
      linbox_check((M.rowdim()==L.rowdim())&&(M.coldim()==U.coldim())&&
		   (L.rowdim()==L.coldim())&&(U.rowdim()==U.coldim())&&
		   (L.coldim()==U.rowdim()));
      
      DenseMatrix<Field>::ConstRowOfColsIterator ccolp;
      DenseMatrix<Field>::Element e;
      DenseMatrix<Field>::ConstColOfRowsIterator crowp1;   
      DenseMatrix<Field>::ConstRowOfColsIterator ccolp2;       
      DenseMatrix<Field>::ConstColIterator cep, cepo;
      DenseMatrix<Field>::ConstRowIterator cep1,cep1e;
      DenseMatrix<Field>::ConstColIterator cep2;
      int i, j;
      i=0;
      for(ccolp=M.rowOfColsBegin(),ccolp2=U.rowOfColsBegin();
	  ccolp!=M.rowOfColsEnd();
	  ++ccolp,++ccolp2,++i)
	{
	  j=0;
	  for(cep=ccolp->begin(),crowp1=L.colOfRowsBegin();
	      cep!=ccolp->end();
	      ++cep,++crowp1,++j)
	    {
	      if(j<=i)
		cep1e=crowp1->begin()+j;
	      else
		cep1e=crowp1->begin()+i+1;
	      
	      M.field().init(e,0);
	      for(cep1=crowp1->begin(),cep2=ccolp2->begin();
		  cep1!=cep1e;
		  ++cep1,++cep2)
		M.field().axpyin(e,*cep1,*cep2);
	      
	      if(j<=i)
		M.field().addin(e ,*cep2);
	      if(!M.field().areEqual(e,*cep))
		return false;
	    }     
	}
      return true;
    }


  template<class Field>
    void LU_MUL(DenseMatrix<Field>& M, const DenseMatrix<Field>& L,const DenseMatrix<Field>& U)
    {
      linbox_check((M.rowdim()==L.rowdim())&&(M.coldim()==U.coldim())&&
		   (L.rowdim()==L.coldim())&&(U.rowdim()==U.coldim())&&
		   (L.coldim()==U.rowdim()));

      DenseMatrix<Field>::RowOfColsIterator ccolp;
      DenseMatrix<Field>::Element e;
      DenseMatrix<Field>::ConstColOfRowsIterator crowp1;   
      DenseMatrix<Field>::ConstRowOfColsIterator ccolp2;       
      DenseMatrix<Field>::ColIterator cep;
      DenseMatrix<Field>::ConstRowIterator cep1,cep1e;
      DenseMatrix<Field>::ConstColIterator cep2;
      int i, j;
      i=0;
      for(ccolp=M.rowOfColsBegin(),ccolp2=U.rowOfColsBegin();
	  ccolp!=M.rowOfColsEnd();
	  ++ccolp,++ccolp2,++i)
	{
	  j=0;
	  for(cep=ccolp->begin(),crowp1=L.colOfRowsBegin();
	      cep!=ccolp->end();
	      ++cep,++crowp1,++j)
	    {
	      if(j<=i)
		cep1e=crowp1->begin()+j;
	      else
		cep1e=crowp1->begin()+i+1;
	      
	      M.field().init(*cep,0);
	      for(cep1=crowp1->begin(),cep2=ccolp2->begin();
		  cep1!=cep1e;
		  ++cep1,++cep2)
		M.field().axpyin(*cep,*cep1,*cep2);
	      
	      if(j<=i)
		M.field().addin(*cep ,*cep2);
	    }     
	}
    }
}
