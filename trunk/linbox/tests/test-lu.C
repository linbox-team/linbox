#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <linbox/field/ntl-zz_p.h>
#include <linbox/blackbox/LU.h>
#include "wan_lu.h"
#include <vector>
#include <iomanip>

typedef LinBox::UnparametricField<NTL::zz_p> Field;

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
  int iteration=100;
  std::cout<<"finished ... 0%"<<flush;
  for(int i=10;i<iteration;++i)
    {
      std::cout<<"\b\b\b"<<setw(2)<<setiosflags(std::ios::fixed)<<(i*100)/iteration<<"%"<<flush;
      if(!test(i))
	{
	  std::cout<<"\b\b\bfailed\n"<<i;
	  return 0;
	}
    }
  std::cout<<"\b\b\bpassed\n";
  return 0;
  
}
