#ifndef __LINBOX_algorithm_iml_rns_H
#define __LINBOX_algorithm_iml_rns_H

#include "linbox/field/unparametric.h"

namespace LinBox { namespace iml {



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

	template<class Container, class FiniteField>
	class RNS {
	public:
		typedef typename FiniteField::Element          ModElement;
		typedef          UnparametricField<ModElement> NoField;
		typedef typename Container::Field              Field;
		typedef typename Field::Element                Element;
		typedef          BlasVector<NoField>           ModVect;
		typedef typename Field::Element                ModField;


	private :
		size_t RNS_bound ;
		Integer          mp_maxInter ;
		size_t              basislen ;
		NoField                  unF ;

		ModVect     RNSbasis ;
		ModVect     RNScombi ;
		ModElement   cumprod ;
		ModVect      bdcoeff ;
		Integer      mp_prod ;

		public:

		RNS():
			RNS_bound(0)
			, mp_maxInter(-1)
			, unF()
			, RNSbasis(unF)
			, RNScombi(unF)
			, bdcoeff(unF)
			, mp_prod(1) //!@bug computed in findRNS
		{
		}

		void combBasis();

		void findRNS();

		void setRNSbound(size_t b)
		{
			RNS_bound = b ;
		}

		void setMaxInter(const Integer& m)
		{
			mp_maxInter = m;
		}

		void
		basisExt ( FiniteField                  & F
			   , std::vector<Container>     & R
			   , Container                  & RE
			   , bool                         pos=false);

		void
		ChineseRemainder ( BlasVector<Field> &Ac, Integer & mp_Ac
				  , bool pos = false);
		void
		cumProd (ModVect &cumprod,
			 const ModVect &extbasis);

		void
		maxInter (const Integer& mp_prod, const Integer& mp_alpha, const long n, Integer& mp_b);

		void
		maxExtInter (const Integer & mp_alpha, const long n, Integer & mp_b);

		void repBound ();

		void RNSbound (const long n);
	};


	/*
	 *
	 * Calling Sequence:
	 *   basisProd(len, basis, mp_prod)
	 *
	 * Summary:
	 *   Compute the product of elements of a RNS basis
	 *
	 * Description:
	 *   Let a RNS basis be 'basis'. The function computes the product of its
	 *   elements basis[0]*basis[1]*...*basis[len-1].
	 *
	 * Input:
	 *     len: long, dimension of RNS basis
	 *   basis: 1-dim typename _Field::Element array length len, RNS basis
	 *
	 * Output:
	 *   mp_prod: mpz_t, product of elements in 'basis'
	 *
	 */

	template<class Field>
	void
	basisProd (BlasVector<Field> & basis, Integer &  mp_prod)
	{
		long i;

		// mpz_set_ui(mp_prod, basis[0]);
		mp_prod = (Integer)basis[0];
		for (i = 1; i < basis.size() ; i++) {
			Integer::mulin(mp_prod,(long int)basis[i]);
		}
	}



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
