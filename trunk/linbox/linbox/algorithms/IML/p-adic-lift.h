#ifndef __LINBOX_algorithm_iml_p_adic_lift_H
#define __LINBOX_algorithm_iml_p_adic_lift_H

#include "linbox/field/unparametric.h"
#include "linbox/algorithms/IML/RNS.h"

namespace LinBox{ namespace iml{


/*
 * padiclift.h includes routines performing p-adic lifting operations
 *
 * Functions:
 *   - liftInit: perform initialization operations before lifting where the
 *     left hand side input matrix is a signed long matrix
 *
 *   - liftInitLlhs: perform initialization operations before lifting where
 *     the left hand side input matrix is a mpz_t matrix
 *
 *   - liftInitRNS: perform initialization operations before lifting where the
 *     left hand side input matrix is represented in a RNS
 *
 *   - lift: compute p-adic lifting coefficients of system of linear equations
 *
 */

	template<class Field>
	class pAdicLift {
	public :
		typedef typename Field::Element ModElement;
		typedef UnparametricField<ModElement> NoField;
	private :
		RNS<Field> &               _liftbasis ;
		RNS<Field>                  _extbasis ;
		BlasVector<NoField>      _liftbasisInv;
		std::vector<BlasVector<Field> > _AInv ;
		BlasVector<NoField>               ARNS;
	public :
		pAdicLift(RNS<Field> & liftbasis) :
			_liftbasis(liftbasis)
		{}

		template<class Matrix>
		int liftInit(const Matrix &A);

		template<class Matrix>
		int liftInitRNS(const Matrix &A);

		void invBasis () ;

		void iml_lift(LinBoxTag::Side s,
			      std::vector<std::vector<BlasMatrix<Field> > > & C,
			      const size_t k,
			       BlasMatrix<PID_integer> & mp_r
			     );

		size_t size()
		{
			return _liftbasis.size() ;
		}
}; // pAdicLift


	}// IML
}// LinBox

#include "p-adic-lift.inl"

#endif // __LINBOX_algorithm_iml_p_adic_lift_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


