#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_H

#include "SlicedPolynomialMatrix.h"
#include "SlicedPolynomialMatrixAddSub.h"

namespace LinBox
{ 
	template< class Field, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixMulKaratsuba
	{
	private:
		typedef Operand1::IntField IntField;
		typedef BlasMatrix<IntField> Matrix;
		typedef std::vector<BlasMatrix<Operand1::IntField>> vec;
		typedef Operand1::polynomial polynomial;
		vec& modulo(vec& C, int n, polynomial irreducible);
		vec& karatsuba(IntField& F, vec& A, vec& B, vec& C);
	public:
		Operand1 &operator() (const Field &GF, Operand1 &C, const Operand2 &A, const Operand3 &B) const;
	}; 
} /* end of namespace LinBox */

#include "SlicedPolynomialMatrixMulKaratsuba.inl"

#endif
// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
