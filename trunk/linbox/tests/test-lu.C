#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <linbox/field/ntl-zz_p.h>
#include <linbox/solutions/lu.h>
#include <linbox/blackbox/dense.h>
#include <linbox/util/commentator.h>
#include <vector>
#include <iomanip>

typedef LinBox::UnparametricField<NTL::zz_p> Field;

namespace LinBox
{
  /*- test if M==L*U;
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
}

bool test(int _SIZE)
{
  Field field;
  typedef LinBox::DenseMatrix<Field>  Matrix;
  Matrix M(field,_SIZE,_SIZE);
  for(int i=0;i<M.rowdim();++i)
      M.setEntry(i,i,NTL::to_zz_p(rand()%1000000+1000000));


  for(int i=0;i<M.rowdim();++i)
    for(int j=0;j<M.coldim();++j)
      M.setEntry(i,j,NTL::to_zz_p(rand()%1000000+1));
 
  LinBox::DenseMatrix<Field> M_C(M);
  LinBox::LU(M);
  return LinBox::LU_MUL_TEST(M_C,M,M);
}
  
int main()
{
  NTL::zz_p::init(1073741789);
  srand(time(0));
  int iteration=120;
  LinBox::commentator.start("Test LU over random matrix in ntl-zz_p field","",iteration-10);
  for(int i=10;i<iteration;++i)
    {
      LinBox::commentator.progress();
      if(!test(i))
	{
	  LinBox::commentator.stop("Failed","");
	  return 0;
	}
    }
  LinBox::commentator.stop("passed","");
  return 0;
}
