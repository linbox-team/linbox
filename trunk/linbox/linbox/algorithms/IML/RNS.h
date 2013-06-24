#ifndef __LINBOX_algorithm_iml_rns_H
#define __LINBOX_algorithm_iml_rns_H

#include "linbox/field/unparametric.h"

namespace LinBox { namespace iml {

	/*
	 *
	 * Calling Sequence:
	 *   basisProd(len, basis, _mp_prod)
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
	 *   _mp_prod: mpz_t, product of elements in 'basis'
	 *
	 */
	template<class Field>
	void
	basisProd (const BlasVector<Field> & basis, Integer &  _mp_prod)
	{
		size_t i;

		// mpz_set_ui(_mp_prod, basis[0]);
		_mp_prod = (Integer)basis[0];
		for (i = 1; i < basis.size() ; i++) {
			Integer::mulin(_mp_prod,(long int)basis[i]);
		}
	}


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

	template<class FiniteField>
	class RNS {
	public:
		typedef typename FiniteField::Element          ModElement;
		typedef          UnparametricField<ModElement> NoField;
		// typedef typename Container::Field              Field;
		// typedef typename Field::Element                Element;
		typedef          BlasVector<NoField>           ModVect;
		// typedef typename Field::Element                ModField;


	private :
		size_t             _RNS_bound ;
		Integer          _mp_maxInter ;
		size_t              _basislen ;
		NoField                  _unF ;

		ModVect     _RNSbasis ;
		ModVect     _RNScombi ;
		ModElement   _cumprod ;
		ModVect      _bdcoeff ;
		Integer      _mp_prod ;

		public:

		RNS():
			_RNS_bound(0)
			, _mp_maxInter(-1)
			, _unF()
			, _RNSbasis(_unF)
			, _RNScombi(_unF)
			, _bdcoeff(_unF)
			, _mp_prod(1) //!@bug computed in findRNS
		{
		}

		void setUp(size_t b, Integer&m)
		{
			_RNS_bound = b;
			_mp_maxInter = m ;
			findRNS();
			// Integer toto = 1;
			// basisProd<NoField>(_RNSbasis,toto);
			// linbox_check(toto == _mp_prod);
			// std::cout << "prod :" << _mp_prod << std::endl;
			repBound();
			// std::cout << "bound:" << _mp_prod << std::endl;

		}

		void
		combBasis (ModVect& RNScombi, const ModVect &RNSbasis);

		void
		combBasis ()
		{
			//!  resize ?
			combBasis(_RNScombi,_RNSbasis);
		}

		void findRNS();


		size_t size() { return _basislen ; }
		size_t basisLength()
		{
			return _basislen ;
		}

		Integer & basisProd() { return _mp_prod ; }


		// *   basisExt(len, n, p, _RNSbasis, _RNScombi, _cumprod, _bdcoeff, R, RE)
		template<class Container>
		void
		basisExt ( FiniteField                  & F
			   , std::vector<Container>     & R
			   , Container                  & RE
			   , bool                         pos=false);

		void
		ChineseRemainder ( BlasVector<NoField> &Ac, Integer & mp_Ac
				  , bool pos = false);
		void
		cumProd (ModVect &_cumprod,
			 const ModVect &extbasis);

		cumProd(const ModVect &extbasis)
		{
			cumProd(_cumprod,extbasis);
		}

		void
		maxInter (const Integer& _mp_prod, const Integer& mp_alpha, const size_t n, Integer& mp_b);

		void
		maxExtInter (const Integer & mp_alpha, const size_t n, Integer & mp_b);

		void repBound ();
		// void repBound (ModVect& bdcoeff, const ModVect &RNSbasis, const ModVect& RNScombi);

		static ModElement RNSbound (const size_t n);

		ModVect & refRNSbasis() { return _RNSbasis ; }
		ModVect & refRNScombi() { return _RNScombi ; }

		ModElement & prime(size_t i) { return _RNSbasis[i] ;}
	};


	template<class  Field>
	class RNSmatrix {
	public:
		RNS<Field>                      & basis;
		std::vector<BlasMatrix<Field> > & matRNS;

		RNSmatrix(RNS<Field>                    & myBasis,
			  std::vector<BlasMatrix<Field> & myARNS):
			basis(myBasis), matRNS(myARNS)
		{}
	};


	template<class Matrix>
	void magnitude(typename Matrix::Element alpha, const Matrix & A)
	{
		BlasMatrixDomain<Ring> BMD(A.field());
		BMD.Magnitude(alpha,A);
	}

	template<class Field>
	void magnitude(Integer & alpha, const RNSmatrix<Field> & A)
	{
		basisProd(A.basis.refBasis(), alpha);
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
