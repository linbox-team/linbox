#ifndef MATRIX_DOMAIN_H
#define MATRIX_DOMAIN_H

namespace LinBox
{
  class MD
    {
    public:
      
      /* M1==M2
       */              
      template<class Matrix1, class Matrix2>
	bool areEqual(const Matrix1& M1, const Matrix2& M2) const;
      
      /*M1==0
       */
      template<class Matrix1>
	bool isZero(const Matrix1& M1);
            
      /* M1<-M2*M3
       */
      template<class Matrix1,class Matrix2,class Matrix3>
	Matrix1& mul(Matrix1& M1, const Matrix2& M2, const Matrix3& M3) const;

      /*M1<-M1*M2
       */
      template<class Matrix1,class Matrix2>
	Matrix1& mulin_L(Matrix1& M1, const Matrix2& M2) const;
      
      /*M1<-M2*M1
       */
      template<class Matrix1, class Matrix2>
	Matrix1& mulin_R(Matrix1& M1, const Matrix2& M2) const;

       /*M1<-M2+M3;
       */
      template<class Matrix1,class Matrix2,class Matrix3>
	Matrix1& add(Matrix1& M1, const Matrix2& M2, const Matrix3& M3) const;
    
      /*M1<-M1+M2;
       */
      template<class Matrix1, class Matrix2>
	Matrix1& addin(Matrix1& M1, const Matrix2& M2) const;

       /*M1<-M2-M3;
       */
      template<class Matrix1, class Matrix2, class Matrix3>
	Matrix1& sub(Matrix1& M1, const Matrix2& M2, const Matrix3& M3) const;

      /*M1<-M1-M2;
       */
      template<class Matrix1,class Matrix2>
	Matrix1& subin(Matrix1& M1, const Matrix2& M2) const;

      /*M1<-M2**k;
       */
      template<class Matrix1, class Matrix2>
	Matrix1& pow_apply(Matrix1& M1, const Matrix2& M2, unsigned long int k) const;
      

      template<class Matrix1, class Matrix2>
	Matrix1& pow_horn(Matrix1& M1, const Matrix2& M2, unsigned long int k) const;
    };

}
#include <linbox/field/matrix-domain.C>
#endif
