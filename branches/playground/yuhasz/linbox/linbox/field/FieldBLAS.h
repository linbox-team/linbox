#ifndef FIELDBLAS_H
#define FIELDBLAS_H

namespace LinBox
{
	template<class Field>
	class FieldBLAS
	{
	public:
	   typedef typename Field::Element Element;

	   FieldBLAS(const Field& F =Field()):_F(F){}

	   template<class Matrix>
	   bool isZero(const Matrix& m) const;

	   template<class Vector>
	   Vector& waxpby(Vector& w, const Element& a, const Vector& x, const Element& b, const Vector& y) const;

	   template<class Matrix, class Vector>
	   Vector& apply(Vector& res, const Matrix& M, const Vector& V) const;

	   template<class Matrix>
	   Matrix& mul(Matrix& res, const Matrix& M1, const Matrix& M2) const;

	   template<class Matrix>
	   Matrix& mulin(Matrix& M1, const Matrix& M2) const;

	protected:
	   Field _F;
	};
}

#include "FieldBLAS.C"
#endif
