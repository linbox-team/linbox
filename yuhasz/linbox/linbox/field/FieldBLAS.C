#ifndef FIELDBLAS_C
#define FIELDBLAS_C
#include <linbox/field/FieldBLAS.h>
#include <linbox/util/debug.h>

namespace LinBox
{
	template<class Field>
	template<class Matrix>
	bool FieldBLAS<Field>::isZero(const Matrix& M) const
	{
	 typename Matrix::ConstRawIterator p=M.rawBegin();
	 for(;p!=M.rawEnd();++p)
	    if(!_F.isZero(*p)) return false;
	 return true;
	}

	template<class Field>
	template<class Vector>
	Vector& FieldBLAS<Field>:: waxpby(Vector& w, const Element& a, const Vector& x, const Element& b, const Vector& y) const
	{
	 typename Vector::iterator pw=w.begin();
	 typename Vector::const_iterator px=x.begin(),
 					 py=y.begin();
	 linbox_check((w.size()==x.size())&&(w.size()==y.size()));
	 Element tmp;
	 for(;pw!=w.end();++pw,++px,++py)
	   {
	      _F.mul(tmp, a ,*px);
	      _F.axpy(*pw,b,*py,tmp);
	   }
	 return w;
	}

	template<class Field>
	template<class Matrix, class Vector>
	Vector& FieldBLAS<Field>::apply(Vector& res, const Matrix& M, const Vector& V) const
	{
	  return M.apply(res,V);
	}

	template<class Field>
	template<class Matrix>
	Matrix& FieldBLAS<Field>::mul(Matrix& res, const Matrix& M1, const Matrix& M2) const
	{
	   linbox_check((res.rowdim()==M1.rowdim())&&(res.coldim()==M2.coldim())&&(M1.coldim()==M2.rowdim()));
	   typename Matrix::RowOfColsIterator pres=res.rowOfColsBegin();
	   typename Matrix::ConstRowOfColsIterator pm2=M2.rowOfColsBegin();
	   for(;pres!=res.rowOfColsEnd();++pres,++pm2)
              M1.apply(*pres, *pm2);
	   return res;
	}
 	template<class Field>
	template<class Matrix>
	Matrix& FieldBLAS<Field>::mulin(Matrix& M1, const Matrix& M2) const
	{
	  linbox_check((M1.coldim()==M2.rowdim())&&(M2.coldim()==M2.rowdim()));
	  Matrix tmp(M1);
	  mul(M1,tmp,M2);
	  return M1;
	} 
}

#endif
