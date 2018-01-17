#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixVectorMulKaratsuba_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixVectorMulKaratsuba_H

#include "linbox/matrix/slicedpolynomialmatrix/SlicedPolynomialMatrix.h"
#include "linbox/vector/slicedpolynomialvector/SlicedPolynomialVector.h"

namespace LinBox
{ 
	template< class Field, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixVectorMulKaratsuba
	{
	private:
		typedef Operand1::IntField IntField;
		typedef BlasMatrix<IntField> Matrix;
		typedef BlasVector<IntField> Vector;
		typedef std::vector<Matrix> vecM;
		typedef std::vector<Vector> vecV;
		typedef Operand1::polynomial polynomial;
		vec& karatsuba(IntField& F, vecV& C, vecM& A, vecV& B);
	public:
		Operand1 &operator() (const Field &GF, Operand1 &C, const Operand2 &A, const Operand3 &B) const;
	}; 
} /* end of namespace LinBox */

#include "SlicedPolynomialMatrixVectorMulKaratsuba.inl"

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
