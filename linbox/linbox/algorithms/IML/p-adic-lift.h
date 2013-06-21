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
	class LiftStep {
	private :
		size_t _len ;
		size_t _dim ; // # of lifting steps
		size_t _row ;
		size_t _col ;
		std::vector<std::vector<Field> >  _dataC ;
		std::vector<unsigned long> & _primes ;
	public :
		LiftStep(size_t l,size_t k, size_t m, size_t n, const std::vector<unsigned long> &primes) :
			_len(l),_dim(k),_row(m),_col(n),_primes(primes)
		{
			typedef std::vector<BlasMatrix<Field > > MV ;
			Field F(101);
			const BlasMatrix<Field> MZ(F,m,n);
			MV VMZ (_len,MZ);
			for (size_t i = 0 ; i < _len ; ++i) {
				Field Fi(_primes[i]);
				VMZ[i].changeField(Fi);
			}
			_dataC.resize(_dim,(const MV) VMZ);
		}
		size_t rowdim() { return _row ;}
		size_t coldim() { return _col ;}
		size_t steps()  { return _dim ;}
		size_t bdim()   { return _len ;}

	} ;

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
		std::vector<BlasVector<NoField> > ARNS;
	public :
		pAdicLift(RNS<Field> & liftbasis) :
			_liftbasis(liftbasis)
		{
		}

		template<class Matrix>
		int liftInit(const Matrix &A);

		template<class Matrix>
		int liftInitRNS(const Matrix &A);

		void invBasis () ;

		void iml_lift(LinBoxTag::Side s,
			      LiftStep<Field> & C,
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


