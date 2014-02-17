#ifndef __LINBOX_algorithm_iml_reconstruct_solution_H
#define __LINBOX_algorithm_iml_reconstruct_solution_H


/*
 * reconsolu.c includes routines performing rational reconstructions
 *
 * Functions:
 *   - sumliftCoeff: recursively sum up the p-adic lifting coefficients
 *
 *   - find2exp: compute floor(log[2](len))
 *
 *   - iratrecon: perform rational number reconstruction
 *
 *   - soluRecon: try reconstructing the rational solution using p-adic
 *     lifting coefficients
 *
 *   - findNumer: certify correctness of the input denominator and, upon
 *     success, compute the numerator
 *
 */

namespace LinBox { namespace iml {
template<class FiniteField>
class SolutionReconstruct {


public:
	void
	sumliftCoeff (const Integer             &mp_basisprod
		      , const size_t             k
		      , BlasVector<PID_integer> &C
		      , Integer                 &mp_sum);

	void
	sumCoeff_rec (size_t                      start
		      , size_t                    len
		      , BlasVector<PID_integer> & C
		      , BlasVector<PID_integer> & mp_pow
		      , size_t                    splflag
		      , size_t                    savflag
		      , size_t                  & idx
		      , BlasVector<PID_integer> & mp_left
		      , Integer                 & mp_right);

	size_t
	find2exp (const size_t len);

	int
	iratrecon (const Integer &mp_u
		   , const Integer &mp_m
		   , const Integer &mp_nb
		   , const Integer &mp_db
		   , Integer &mp_N
		   , Integer &mp_D);

	int
	soluRecon (const Tag::Side solupos
		   , const size_t              k
		   , RNS<FiniteField>        & rns
		   , pAdicLift<FiniteField>  & C
		   , Integer                 & mp_nb
		   , Integer                 & mp_db
		   , BlasVector<PID_integer> & mp_N
		   , Integer                 & mp_D);

	size_t
	findNumer (const Integer &mp_u
		   , const Integer &mp_m
		   , const Integer &mp_D
		   , const Integer &mp_nb
		   , Integer &mp_N);

}; // reconsolu

} // iml
} // LinBox

#include "reconstruct-solution.inl"

#endif // __LINBOX_algorithm_iml_reconstruct_solution_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


