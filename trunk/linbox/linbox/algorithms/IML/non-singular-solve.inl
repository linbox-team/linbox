#ifndef __LINBOX_algorithm_iml_non_singular_solve_INL
#define __LINBOX_algorithm_iml_non_singular_solve_INL





/*
 * Calling Sequence:
 *   adBasis(idx, basislen, liftbasis)
 *
 * Summary:
 *   Adjust the lifting basis if some element in the lifting basis is bad.
 *
 * Description:
 *   If A^(-1) mod liftbasis[idx] does not exist, then use this function to
 *   delete liftbasis[idx] and add the previous prime of liftbasis[basislen-1]
 *   into the lifting basis. Meanwhile, move the elements in lifting basis to
 *   maintain the decreasing order.
 *
 * Input:
 *         idx: long, index in liftbasis where A^(-1) mod liftbasis[idx]
 *              does not exist
 *    basislen: long, dimension of basis
 *   liftbasis: 1-dim FiniteField array length basislen, lifting basis
 *
 * Note:
 *   liftbasis is update inplace.
 *
 */

void
adBasis (const long idx, const long basislen, FiniteField *liftbasis)
{
  long i;
  Integer mp_temp;
  mpz_init(mp_temp);

  /* move the elements forward */
  for (i = idx+1; i < basislen; i++) { liftbasis[i-1] = liftbasis[i]; }
  mpz_set_ui(mp_temp, liftbasis[basislen-1]);
  mpz_sub_ui(mp_temp, mp_temp, 1);
  while (mpz_probab_prime_p(mp_temp, 10) == 0)
    mpz_sub_ui(mp_temp, mp_temp, 1);
  liftbasis[basislen-1] = mpz_get_ui(mp_temp);

  mpz_clear(mp_temp);
  return;
}


/*
 * Calling Sequence:
 *   liftbd(mp_basisprod, n, mp_alpha, mp_beta, maxk, mp_maxnb, mp_maxdb, k,
 *          mp_nb, mp_db)
 *
 * Summary:
 *   Compute the initial number of lifting steps, initial solution bounds,
 *   maximum number of lifting steps and maximum solution bounds.
 *
 * Description:
 *   Let the system of linear equations be Ax = b. The function computes
 *   the worst bound(Hadamard's bound) of denominator and numerators and
 *   corresponding number of lifting steps maxk, which satisfies
 *   mp_basisprod^maxk >= 2*N*D and mp_basisprod^(maxk-1) < 2*N*D, where
 *   D = n^ceil(n/2)*mp_alpha^n, worst case bound for denominator,
 *   N = n^ceil(n/2)*mp_alpha^(n-1)*mp_beta, worst case bound for denominators.
 *
 *   The function also computes the initial number of lifting step k = 20 and
 *   the corresponding bound of denominator and numerators, such that
 *   mp_basisprod^k >= 2*N1*D1, mp_basisprod^(k-1) < 2*N1*D1, and N1 = D1,
 *   where N1 is the estimate bound for numerators and D1 is the estimate
 *   bound for denominator.
 *
 * Input:
 *   mp_basisprod: mpz_t, product of elements of lifting basis
 *              n: long, dimension of matrix A
 *       mp_alpha: mpz_t, maximum magnitude of A
 *        mp_beta: mpz_t, maximum magnitude of b
 *
 * Output:
 *       maxk: pointer to long int, maximum number of lifting steps
 *   mp_maxnb: mpz_t, the worst numerator bound N
 *   mp_maxdb: mpz_t, the worst denominator bound D
 *          k: pointer to long int, estimated number of lifting steps
 *      mp_nb: mpz_t, estimated numerator bound N1
 *      mp_db: mpz_t, estimated denominator bound D1
 *
 */

void
liftbd (const Integer   & mp_basisprod
	, const size_t    n
	, const Integer & mp_alpha
	, const Integer & mp_beta
	, long          & maxk
	, Integer       & mp_maxnb
	, Integer       & mp_maxdb
	, long          & k
	, Integer       & mp_nb
	, Integer       & mp_db)
{
  size_t n1;
  Integer mp_t1, mp_t2, mp_prod;

  if ((n % 2) == 0) { n1 = n/2; }
  else { n1 = n/2+1; }

  /* compute the worst case bound */
  // mpz_init(mp_t1);
  // mpz_init(mp_t2);
  Integer::pow(mp_t1,(long unsigned)n,(long unsigned)n1);
  // mpz_ui_pow_ui(mp_t1, n, n1);
  Integer::pow(mp_db, mp_alpha, (long unsigned)n);
  // mpz_pow_ui(mp_db, mp_alpha, n);
  mpz_mul(mp_maxdb, mp_db, mp_t1);
  mpz_pow_ui(mp_nb, mp_alpha, n-1);
  mpz_mul(mp_nb, mp_nb, mp_beta);
  mpz_mul(mp_maxnb, mp_nb, mp_t1);
  mpz_init_set(mp_prod, mp_maxdb);
  mpz_mul(mp_prod, mp_prod, mp_maxnb);
  mpz_mul_ui(mp_prod, mp_prod, 2);
  mpz_add_ui(mp_prod, mp_prod, 1);

  /* compute maxk */
  *maxk = 1;
  mpz_set(mp_t1, mp_basisprod);
  while (mpz_cmp(mp_t1, mp_prod) < 0)
    {
      mpz_mul(mp_t1, mp_t1, mp_basisprod);
      ++(*maxk);
    }

  /* compute k and estimate bound */
  *k = 20;
  mpz_pow_ui(mp_prod, mp_basisprod, *k);
  mpz_sub_ui(mp_prod, mp_prod, 1);
  mpz_divexact_ui(mp_prod, mp_prod, 2);
  mpz_sqrt(mp_nb, mp_prod);
  mpz_set(mp_db, mp_nb);
  if (*k >= *maxk) { mpz_set(mp_nb, mp_maxnb); mpz_set(mp_db, mp_maxdb); }
  mpz_clear(mp_t1);
  mpz_clear(mp_t2);
  mpz_clear(mp_prod);
  return;
}

#endif // __LINBOX_algorithm_iml_non_singular_solve_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
