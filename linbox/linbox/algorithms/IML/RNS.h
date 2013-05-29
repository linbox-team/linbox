#ifndef __LINBOX_algorithm_iml_rns_H
#define __LINBOX_algorithm_iml_rns_H

namespace LinBox { namespace IML {



	/*
	 *
	 * RNSOp.h includes rountines related to the operations in Residue Number
	 * System (RNS).
	 *
	 * Functions:
	 *   - basisExt: given a representation of a matrix/vector in some RNS,
	 *     extend to compute the representation in another positive integer
	 *
	 *   - basisExtPos: given a representation of a non-negative matrix/vector in
	 *     some RNS, extend to compute the representation in another positive
	 *     integer
	 *
	 *   - basisProd: compute the product of elements of a RNS basis
	 *
	 *   - ChineseRemainder: given a representation of an integer in some RNS, use
	 *     Chinese Remainder Algorithm to reconstruct the integer
	 *
	 *   - ChineseRemainderPos: given a representation of a non-negative integer
	 *     in some RNS, use Chinese Remainder Algorithm to reconstruct the integer
	 *
	 *   - combBasis: compute the special combination of RNS basis
	 *
	 *   - cumProd: compute the representation of the combination of elements of
	 *     one RNS basis in another RNS basis
	 *
	 *   - findRNS: find a RNS basis and its special combination
	 *
	 *   - maxInter: compute the maximum interval of positive and negative results
	 *      of a matrix-matrix or matrix-vector product
	 *
	 *   - repBound: compute the mix radix coefficients of a special integer in a
	 *     RNS basis
	 *
	 * Note: the modular operations in these functions are implemented by function
	 *   fmod, which requires the bit-length of operands(e.g., RNS basis) could
	 *   fit into mantissa of floating point numbers, i.e., less than 53 bits.
	 *
	 */

	template<class Container>
	class RNS {
	public:
		typedef typename Container::Field Field;
		typedef typename Field::Element Element;
		typedef typename BlasVector<Field> Vect;
		Field & F;
		size_t len ;
		BlasVector<Field> basis ;
		BlasVector<Field> cmbasis ;
		Element cumprod;
		BlasVector<Field> bdcoeff;
	}



// }; // RNS

} // IML

} // LinBox

#include "RNS.inl"

#endif // __LINBOX_algorithm_iml_rns_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
