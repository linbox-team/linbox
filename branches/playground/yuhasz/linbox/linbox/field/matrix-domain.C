#ifndef MATRIX_DOMAIN_C
#define MATRIX_DOMAIN_C
#include <linbox/util/debug.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/vector/vector-domain.h>
#include <vector>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/dense.h>

namespace LinBox
{
  template<class Matrix1, class Matrix2>
  bool MatrixDomain::areEqual(const Matrix1& M1, const Matrix2& M2) const
  {
    if((M1.rowdim()!=M2.rowdim())||
       (M1.coldim()!=M2.coldim()))
      return false;
    typename Matrix1::RowIterator p1=M1.rowBegin();
    typename Matrix2::RowIterator p2=M2.rowBegin();
    
    for(;p1!=M1.rowEnd();++p1,++p2)
      if(!VectorDomain_gen(M1.field()).areEqual(*p1,*p2))
	return false;

    return true;
  }
      
  template<class Matrix1>
  bool MatrixDomain::isZero(const Matrix1& M1)
  {
    typename Matrix1::ConstRowIterator p1=M1.rowBegin();
    for(;p1!=M1.rowEnd();++p1)
      if(!VectorDomain_gen(M1.field()).isZero(*p1))
	return false;
    
    return true;
  }

  template<class Matrix1,class Matrix2,class Matrix3>
  Matrix1& MatrixDomain::mul(Matrix1& M1, const Matrix2& M2, const Matrix3& M3) const
  {
    linbox_check((M1.rowdim()==M2.rowdim())&&
		 (M1.coldim()==M3.coldim())&&
		 (M2.coldim()==M3.rowdim()));
    typename Matrix1::ColIterator p1=M1.colBegin();
    typename Matrix3::ConstColIterator p3=M3.colBegin();
    for(;p1!=M1.colEnd();++p1,++p3)
      M2.apply(*p1,*p3);

    return M1;
  }

  

  template<class Matrix1, class Matrix2>
  Matrix1& MatrixDomain::mulin_L(Matrix1& M1, const Matrix2& M2) const
  {
    linbox_check((M1.coldim()==M2.rowdim())&&
		 (M2.rowdim()==M2.coldim()));

    typename Matrix1::RowIterator p1;
    for(p1=M1.rowBegin();p1!=M1.rowEnd();++p1)
      M2.applyTransposeIn(*p1);

    return M1;
  }


 
  template<class Matrix1, class Matrix2>
  Matrix1& MatrixDomain::mulin_R(Matrix1& M1, const Matrix2& M2) const
  {
    linbox_check((M1.rowdim()==M2.coldim())&&
		 (M2.coldim()==M2.rowdim()));

    typename Matrix1::ColIterator p1;
    for(p1=M1.colBegin();p1!=M1.colEnd();++p1)
      M2.applyIn(*p1);

    return M1;
  }

  template<class Matrix1, class Matrix2, class Matrix3>
  Matrix1& MatrixDomain::add(Matrix1& M1, const Matrix2& M2, const Matrix3& M3) const
  {
    linbox_check((M1.rowdim()==M2.rowdim())&&(M1.rowdim()==M3.rowdim())&&
		 (M1.coldim()==M2.coldim())&&(M1.coldim()==M3.coldim()));

          
    typename Matrix1::ColIterator p1=M1.colBegin();
    typename Matrix2::ConstColIterator p2=M2.colBegin();
    typename Matrix3::ConstColIterator p3=M3.colBegin();
    typename Matrix1::Col::iterator pe1;
    typename Matrix2::ConstCol::iterator pe2;
    typename Matrix3::ConstCol::iterator pe3;
    
    for(;p1!=M1.colEnd();++p1,++p2,++p3)
      for(pe1=p1->begin(),pe2=p2->begin(),pe3=p3->begin();pe1!=p1->end();++pe1,++pe2,++pe3)
	M1.field().add(*pe1,*pe2,*pe3);
    
    return M1;    
  }
  
  template<class Matrix1, class Matrix2>
  Matrix1& MatrixDomain::addin(Matrix1& M1, const Matrix2& M2) const
  {
    linbox_check((M1.rowdim()==M2.rowdim())&&
		 (M1.coldim()==M2.coldim()));
    typename Matrix1::ColIterator p1=M1.colBegin();
    typename Matrix2::ConstColIterator p2=M2.colBegin();
    typename Matrix1::Col::iterator pe1;
    typename Matrix2::ConstCol::iterator pe2;
    for(;p1!=M1.colEnd();++p1,++p2)
      for(pe1=p1->begin(),pe2=p2->begin();pe1!=p1->end();++pe1,++pe2)
	M1.field().addin(*pe1,*pe2);

    return M1;
  }

  template<class Matrix1, class Matrix2, class Matrix3>
  Matrix1& MatrixDomain::sub(Matrix1& M1, const Matrix2& M2, const Matrix3& M3) const
  {
    linbox_check((M1.rowdim()==M2.rowdim())&&(M1.rowdim()==M3.rowdim())&&
		 (M1.coldim()==M2.coldim())&&(M1.coldim()==M3.coldim()));

          
    typename Matrix1::ColIterator p1=M1.colBegin();
    typename Matrix2::ConstColIterator p2=M2.colBegin();
    typename Matrix3::ConstColIterator p3=M3.colBegin();
    typename Matrix1::Col::iterator pe1;
    typename Matrix2::ConstCol::iterator pe2;
    typename Matrix3::ConstCol::iterator pe3;
    
    for(;p1!=M1.colEnd();++p1,++p2,++p3)
      for(pe1=p1->begin(),pe2=p2->begin(),pe3=p3->begin();pe1!=p1->end();++pe1,++pe2,++pe3)
	M1.field().sub(*pe1,*pe2,*pe3);
    
    return M1;    
  }
  
  template<class Matrix1,class Matrix2>
  Matrix1& MatrixDomain::subin(Matrix1& M1, const Matrix2& M2) const
  {
    linbox_check((M1.rowdim()==M2.rowdim())&&
		 (M1.coldim()==M2.coldim()));
    typename Matrix1::ColIterator p1=M1.colBegin();
    typename Matrix2::ConstColIterator p2=M2.colBegin();
    typename Matrix1::Col::iterator pe1;
    typename Matrix2::ConstCol::iterator pe2;
    for(;p1!=M1.colEnd();++p1,++p2)
      for(pe1=p1->begin(),pe2=p2->begin();pe1!=p1->end();++pe1,++pe2)
	M1.field().subin(*pe1,*pe2);

    return M1;
  }

   /*M1<-M2**k;
       */
  template<class Matrix1, class Matrix2>
  Matrix1& MatrixDomain::pow_apply(Matrix1& M1, const Matrix2& M2, unsigned long int k) const
  {
    linbox_check((M1.rowdim()==M1.coldim())&&
		 (M2.rowdim()==M2.coldim())&&
		 (M1.rowdim()==M2.rowdim()));

  
    typename Matrix1::RawIterator p=M1.rawBegin();
    for(;p!=M1.rawEnd();++p)
      M1.field().init(*p,0);
    for(p=M1.rawBegin();p<M1.rawEnd();)
      {
	M1.field().init(*p,1);
	p=p+M1.rowdim()+1;
      }
    
    
    for(int i=0;i<k;++i)
      mulin_R(M1,M2);
    
    return M1;
  }
    
  
  template<class Matrix1, class Matrix2>
  Matrix1& MatrixDomain::pow_horn(Matrix1& M1, const Matrix2& M2, unsigned long int k) const
  {
    linbox_check((M1.rowdim()==M1.coldim())&&
		 (M2.rowdim()==M2.coldim())&&
		 (M1.rowdim()==M2.rowdim()));
    
    if(k==0)
      {
	typename Matrix1::RawIterator p=M1.rawBegin();
	for(;p!=M1.rawEnd();++p)
	  M1.field().init(*p,0);
	for(p=M1.rawBegin();p<M1.rawEnd();)
	  {
	    M1.field().init(*p,1);
	    p+=M1.rowdim()+1;
	  }
	return M1;
      }
    
    typename Matrix1::RawIterator p1;
    typename Matrix2::ConstRawIterator p2;
    for(p1=M1.rawBegin(),p2=M2.rawBegin();p1!=M1.rawEnd();++p1,++p2)
      M1.field().assign(*p1,*p2);
  
    std::vector<bool> bit;
    bit.reserve(sizeof(unsigned long)*4);
    while(k>0)
      {
	bit.push_back(k%2);
	k/=2;
      };

    
    std::vector<bool>::reverse_iterator p=bit.rbegin();
    ++p;
    Matrix1 temp(M1);
    for(;p!=bit.rend();++p)
      {
	temp=M1;
	mulin_L(M1,temp);
	if(*p)
	  mulin_L(M1,M2);

      }
      
    return M1;     
  }
    
}
#endif
