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

class reconsolu {

void sumliftCoeff(const mpz_t mp_basisprod, const long k, mpz_t* C, \
		  mpz_t mp_sum);

void sumCoeff_rec(long start, long len, mpz_t *C, mpz_t *mp_pow, long splflag,\
		  long savflag, long *idx, mpz_t *mp_left, mpz_t mp_right);

long find2exp(const long len);

long iratrecon(const mpz_t mp_u, const mpz_t mp_m, const mpz_t mp_nb, \
	       const mpz_t mp_db, mpz_t mp_N, mpz_t mp_D);

long soluRecon(const enum SOLU_POS solupos, const long k, const long basislen,\
	       const long n, const long m, const mpz_t mp_basisprod, \
	       const FiniteField *basis, const FiniteField *cmbasis, \
	       Double ***C, mpz_t mp_nb, mpz_t mp_db, mpz_t *mp_N, mpz_t mp_D);

long findNumer(const mpz_t mp_u, const mpz_t mp_m, const mpz_t mp_D, \
	       const mpz_t mp_nb, mpz_t mp_N);



}; // reconsolu

#endif // __LINBOX_algorithm_iml_reconstruct_solution_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


