/*file directsum.h
 *Author: Zhendong Wan for LinBox group
 */
#ifndef _DIRECT_SUM_H
#define _DIRECT_SUM_H
#include <algorithm>
#include <stddef.h>

namespace LinBox
{
  template<class sparse_mat>
    class DirSum
    {
    public:
      // Constructor from two sparse matrix.
      DirSum(const sparse_mat& A, const sparse_mat& B):_A(A),_B(B){}
      
      //copy constructor.
      DirSum(const DirSum<sparse_mat>& M):_A(M._A),_B(M._B){}
      
      template<class Vector>
	Vector& apply(Vector& y, const Vector& x) const
	{
	  if(x.size()!=coldim())
	    return y;
	  
	  int col_a;
	  col_a=_A.coldim();
	  
	  Vector x1(x.begin(),x.begin()+col_a),
	    x2(x.begin()+col_a,x.end());
	  
	  Vector y1,y2;
	  _A.apply(y1,x1);
	  _B.apply(y2,x2);
	  
	  y.resize(y1.size()+y2.size());
	  std::copy(y1.begin(),y1.end(),y.begin());
	  std::copy(y2.begin(),y2.end(),y.begin()+y1.size());
	  return y;
	}
	  
      template<class Vector>
	Vector& applyTranspose(Vector& y, const Vector& x) const
	{
	  if(x.size()!=rowdim())
	    return y;
	  
	  int row_a;
	  row_a=_A.rowdim();
	  
	  Vector x1(x.begin(),x.begin()+row_a),
	    x2(x.begin()+row_a,x.end());
	  
	  Vector y1,y2;
	  _A.applyTranspose(y1,x1);
	  _B.applyTranspsoe(y2,x2);
	  
	  y.resize(y1.size()+y2.size());
	  std::copy(y1.begin(),y1.end(),y.begin());
	  std::copy(y2.begin(),y2.end(),y.begin()+y1.size());
	  return y;
	}

      size_t coldim(void) const
	{
	  return _A.coldim()+_B.coldim();}

      size_t rowdim(void) const
	{
	  return _A.rowdim()+_B.rowdim();
	}
      
    private:
      sparse_mat _A;
      sparse_mat _B;
    };
}
#endif
