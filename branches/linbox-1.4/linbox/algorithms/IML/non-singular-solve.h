#ifndef __LINBOX_algorithm_iml_non_singular_solve_H
#define __LINBOX_algorithm_iml_non_singular_solve_H

/*
 * nonsysolve.h includes routines solving the nonsingular system of linear
 * equations.
 *
 * Functions:
 *
 *   - findliftbasisSmall: compute the p-adic lifting basis
 *
 *   - findliftbasisLarge: compute the p-adic lifting basis
 *
 *   - adBasis: adjust the lifting basis if some element in the lifting basis
 *     is bad
 *
 *   - liftbd: compute the initial number of lifting steps, initial solution
 *     bounds, maximum number of lifting steps and maximum solution bounds
 *
 */

class NonSingularSolve {

void adBasis(const long idx, const long basislen, FiniteField *liftbasis);

FiniteField *findLiftbasisSmall(const long n, const mpz_t mp_alpha, \
				long *basislen);

FiniteField *findLiftbasisLarge(const long n, const mpz_t mp_alpha, \
				long *basislen);

void liftbd(const mpz_t mp_basisprod, const long n, const mpz_t mp_alpha, const mpz_t mp_beta, long *maxk, mpz_t mp_maxnb, mpz_t maxdb, long *k, mpz_t mp_nb, mpz_t mp_db);


};

#endif // __LINBOX_algorithm_iml_non_singular_solve_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
