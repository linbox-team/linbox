#ifndef __LINBOX_algorithm_iml_reconstruct_solution_INL
#define __LINBOX_algorithm_iml_reconstruct_solution_INL

namespace LinBox { namespace iml {

/*
 * Calling Sequence:
 *   sumliftCoeff(mp_basisprod, k, C, mp_sum)
 *
 * Summary:
 *   Recursively sum up the p-adic lifting coefficients
 *
 * Description:
 *   Given p-adic lifing coefficients C after k lifting steps, the function
 *   recursively compute
 *     A^(-1).B mod p^k =
 *     C[0]+C[1]*mp_basisprod+ ... +C[k-1]*mp_basisprod^(k-1)
 *
 * Input:
 *   mp_basisprod: mpz_t, product of lifting basis
 *              k: long, dimension of array C, number of lifting steps
 *              C: 1-dim mpz_t matrix length k, storing lifting coefficients
 *                 computed each lifting step
 *
 * Output:
 *   mp_sum: mpz_t, the sum of lifting coefficients
 *
 */

void
sumliftCoeff (const mpz_t mp_basisprod, const long k, mpz_t* C, mpz_t mp_sum)
{
  long i, t, splflag;
  long idx[1];
  mpz_t mp_right, *mp_pow, *mp_left;

  idx[0] = 0;

  /* precomputation mp_basisprod^log(k) */
  t = find2exp(k);
  mp_pow = XMALLOC(mpz_t, t+1);
  mpz_init_set(mp_pow[0], mp_basisprod);
  for (i = 1; i < t+1; i++)
    {
      mpz_init(mp_pow[i]);
      mpz_pow_ui(mp_pow[i], mp_pow[i-1], 2);
    }
  mpz_init(mp_right);
  if (t == 0)
    {
      mpz_set(mp_sum, C[0]);
      for (i = 0; i < t+1; i++) { mpz_clear(mp_pow[i]); } { XFREE(mp_pow); }
      mpz_clear(mp_right);
      return;
    }
  if ((1 << t) == k) { splflag = 1; }
  else { splflag = 0; }
  mp_left = XMALLOC(mpz_t, t);
  for (i = 0; i < t; i++) { mpz_init(mp_left[i]); }
  sumCoeff_rec(0, k, C, mp_pow, splflag, 0, idx, mp_left, mp_right);
  mpz_set(mp_sum, mp_left[0]);

  for (i = 0; i < t+1; i++) { mpz_clear(mp_pow[i]); } { XFREE(mp_pow); }
  for (i = 0; i < t; i++) { mpz_clear(mp_left[i]); } { XFREE(mp_left); }
  mpz_clear(mp_right);
  return;
}


/*
 * Calling Sequence:
 *   sumCoeff_rec(s, len, C, mp_pow, splflag, 0, idx, mp_left, mp_right)
 *
 * Summary:
 *   compute C[s]+C[s+1]*p+C[s+2]*p^2+...+C[s+len-1]*p^(len-1) recursively
 *
 * Input:
 *      start: long, recursion start index in array C, usually 0
 *        len: long, number of entries in C to use
 *          C: 1-dim mpz_t array, storing values for computation
 *     mp_pow: 1-dim mpz_t array length t+1, storing different powers of p,
 *             mp_pow[i] = p^(2^i), (i = 0..t), t = floor(log(len)/log(2))
 *    splflag: 1/0, input flag,
 *           - if len=2^i and i >= 0, then splflag = 1
 *           - else splflag = 0
 *    savflag: 1/0, input flag, initially use 0
 *           - if savflag = 0, then current subproblem is in the left subtree
 *             of recursion tree
 *           - else, current subproblem is in the right subtree of recursion
 *             tree
 *        idx: pointer to a long integer, the index of available entries in
 *             array mp_left, initially *idx = 0
 *    mp_left: 1-dim mpz_t array length t, storing the intermidiate result of
 *             left subtree
 *   mp_right: mpz_t, storing the result of right subtree
 *
 * Output:
 *   mp_left[0] will store the computation result
 *
 * Idea: the result of subproblem on the right subtree of recursion tree is
 *       overwritten into mp_right
 *       the result of subproblem on the left subtree of recursion tree is
 *       stored into mp_left[*idx], and make the *idx always be availabe
 *       index in mp_left, e.g., the space of mp_left used in right subtree can
 *       can be reused by its siblings and parent
 *
 */

void
sumCoeff_rec (long start, long len, mpz_t *C, mpz_t *mp_pow, long splflag,\
	      long savflag, long *idx, mpz_t *mp_left, mpz_t mp_right)
{
  long lpow, llen, tidx;

  /* basecase */
  if (len == 1)
    {
      if (savflag == 0)
	{
	  mpz_set(mp_left[*idx], C[start]);
	  ++*idx;
	}
      else { mpz_set(mp_right, C[start]); }
      return;
    }
  /* if len is power of 2, then directly divide len by 2 to
     obtain two subproblems*/
  if (splflag == 1)
    {
      lpow = find2exp(len)-1;
      llen = 1 << lpow;
      sumCoeff_rec(start, llen, C, mp_pow, 1, 0, idx, mp_left, mp_right);
      tidx = *idx-1;
      sumCoeff_rec(start+llen, len-llen, C, mp_pow, 1, 1, idx, mp_left, \
		   mp_right);
    }
  /* otherwise, find maximal lpow, 2^lpow < len, to make two subproblems */
  else
    {
      lpow = find2exp(len);
      llen = 1 << lpow;
      sumCoeff_rec(start, llen, C, mp_pow, 1, 0, idx, mp_left, mp_right);
      tidx = *idx-1;

      /* if len-2^lpow is power of 2, then not compute the second subproblem */
      if (llen == len)
	{
	  mpz_set(mp_right, mp_left[tidx]);
	  return;
	}
      sumCoeff_rec(start+llen, len-llen, C, mp_pow, 0, 1, idx, mp_left, \
		   mp_right);
    }

  /* combine subproblems by position of two subproblems in the recursion tree*/
  /* if subproblem lies in the left subtree, store result in a new mpz_t int */
  if (savflag == 0)
    mpz_addmul(mp_left[tidx], mp_right, mp_pow[lpow]);

  /* if in the right subtree, store result in the pre-allocated mp_right */
  else
    {
      mpz_mul(mp_right, mp_right, mp_pow[lpow]);
      mpz_add(mp_right, mp_left[tidx], mp_right);
    }
  *idx = tidx+1;
  return;
}


/*
 * Calling Sequence:
 *   i <-- find2exp(len)
 *
 * Summary:
 *   Compute floor(log[2](len))
 *
 * Description:
 *   Compute i >= 0 such that 2^i <= len and 2^(i+1) > len
 *
 */

long
find2exp (const long len)
{
  long i=0, l=1;
  while ((l <<= 1) <= len) { ++i; }
  return i;
}


/*
 * Calling Sequence:
 *   1/0 <-- iratrecon(mp_u, mp_m, mp_nb, mp_db, mp_N, mp_D)
 *
 * Summary:
 *   Perform rational number reconstruction
 *
 * Description:
 *   Given image mp_u, modulus mp_m, the function reconstruct numerator mp_N
 *   and denominator mp_D, such that
 *   mp_N/mp_D = mp_u (mod mp_m), and
 *   abs(mp_N) <= mp_nb, 0 < mp_D <= mp_db
 *
 * Input:
 *    mp_u: mpz_t, image
 *    mp_m: mpz_t, modulus
 *   mp_nb: mpz_t, numerator bound of the rational number
 *   mp_db: mpz_t, denominator bound of the rational number
 *
 * Output:
 *   mp_N: mpz_t, storing the reconstructed numerator if the reconstruction
 *         succeeds
 *   mp_D: mpz_t, storing the reconstructed denominator if the reconstruction
 *         fails
 *
 * Return:
 *   - 1 if reconstruction succeeds
 *   - 0 if reconstruction fails
 *
 */

long
iratrecon (const mpz_t mp_u, const mpz_t mp_m, const mpz_t mp_nb, \
	   const mpz_t mp_db, mpz_t mp_N, mpz_t mp_D)
{
  mpz_t mp_v3, mp_u2, mp_v2, mp_u3, mp_t2, mp_t3, mp_q;

  mpz_init(mp_v3);
  mpz_tdiv_r(mp_v3, mp_u, mp_m);
  mpz_init_set_ui(mp_u2, 0);
  mpz_init_set(mp_u3, mp_m);
  mpz_init_set_ui(mp_v2, 1);
  mpz_init(mp_q);
  mpz_init(mp_t2);
  mpz_init(mp_t3);
  while (mpz_cmpabs(mp_v3, mp_nb) > 0)
    {
      mpz_tdiv_qr(mp_q, mp_t3, mp_u3, mp_v3);
      mpz_mul(mp_t2, mp_q, mp_v2);
      mpz_sub(mp_t2, mp_u2, mp_t2);
      mpz_set(mp_u2, mp_v2);
      mpz_set(mp_u3, mp_v3);
      mpz_set(mp_v2, mp_t2);
      mpz_set(mp_v3, mp_t3);
    }
  if (mpz_cmpabs(mp_v2, mp_db) > 0)
    {
      mpz_clear(mp_v2); mpz_clear(mp_v3); mpz_clear(mp_u2); mpz_clear(mp_q);
      mpz_clear(mp_t2); mpz_clear(mp_t3); mpz_clear(mp_u3);
      return 0;
    }
  mpz_gcd(mp_u2, mp_v2, mp_v3);
  if (mpz_cmp_ui(mp_u2, 1) != 0)
    {
      mpz_clear(mp_v2); mpz_clear(mp_v3); mpz_clear(mp_u2); mpz_clear(mp_q);
      mpz_clear(mp_t2); mpz_clear(mp_t3); mpz_clear(mp_u3);
      return 0;
    }
  if (mpz_sgn(mp_v2) < 0) { mpz_neg(mp_v2, mp_v2);  mpz_neg(mp_v3, mp_v3); }
  { mpz_set(mp_D, mp_v2);  mpz_set(mp_N, mp_v3); }

  { mpz_clear(mp_v2); mpz_clear(mp_v3); mpz_clear(mp_u2); mpz_clear(mp_u3); }
  { mpz_clear(mp_q); mpz_clear(mp_t2);  mpz_clear(mp_t3); }
  return 1;
}


/*
 * Calling Sequence:
 *   1/0 <-- soluRecon(solupos, k, basislen, n, m, mp_basisprod, basis
 *                     cmbasis, C, mp_nb, mp_db, mp_N, mp_D)
 *
 * Summary:
 *   Try reconstructing the rational solution using p-adic lifting coefficients
 *
 * Description:
 *   Given p-adic lifting coefficients computed from function lift, this
 *   function try reconstructing the rational solution of a nonsingular
 *   system. If the trial succeeds, the function outputs the reconstructed
 *   numerator matrix and denominator and returns 1. If the system is not
 *   lifted high enough, usually the trial will fail in the first several
 *   solution entries. Then the function returns 0.
 *
 * Input:
 *        solupos: enumerate, flag to indicate whether to transpose A or not
 *               - solupos = LeftSolu: system be XA = B, B m x n matrix
 *               - solupos = RightSolu: system be AX = B, B n x m matrix
 *              k: long, first dimension of C, number of lifting steps
 *       basislen: long, dimension of lifting basis
 *              n: long, dimension of left hand side input matrix A
 *              m: long, row/column dimension of right hand side input matrix B
 *   mp_basisprod: mpz_t, product of lifting basis
 *          basis: 1-dim FiniteField array length basislen, lifting basis
 *        cmbasis: 1-dim FiniteField array length len, computed by function
 *                 combBasis, inverses of special combination of lifting basis
 *              C: 3-dim Double array, dimension k x liftbasislen x n*m,
 *                 computed by function lift.
 *          mp_nb: mpz_t, bound of absoluate value of numerator of
 *                 reconstructed rational number
 *          mp_db: mpz_t, bound of denominator of reconstructed rational number
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length n*m, representation of a n x m or m x n
 *         reconstructed numerator matrix of the solution
 *   mp_D: mpz_t, denominator of the solution
 *
 * Return:
 *   - 1 if reconstruction succeeds
 *   - 0 if reconstruction fails
 *
 */

long
soluRecon (const enum SOLU_POS solupos, const long k, const long basislen,\
	   const long n, const long m, const mpz_t mp_basisprod, \
	   const FiniteField *basis, const FiniteField *cmbasis, \
	   Double ***C, mpz_t mp_nb, mpz_t mp_db, mpz_t *mp_N, mpz_t mp_D)
{
  long i, j, s, t, h, l, r, ri, len=1;
  mpz_t mp_sum, mp_m, mp_temp;
  mpz_t *mp_Dt, *mp_C;
  long *idx;
  Double *C1;

  mpz_init(mp_sum);
  mpz_init(mp_m);
  mpz_pow_ui(mp_m, mp_basisprod, k);
  mp_Dt = XMALLOC(mpz_t, len);
  mpz_init(mp_Dt[0]);

  /* idx stores the index in mp_N where denominator changes */
  idx = XCALLOC(long, len);
  mp_C = XMALLOC(mpz_t, k);
  for (j = 0; j < k; j++) { mpz_init(mp_C[j]); }
  C1 = XMALLOC(Double, basislen);
  for (i = 0; i < k; i++)
    {
      for (j = 0; j < basislen; j++) { C1[j] = C[i][j][0]; }
      ChineseRemainderPos(basislen, basis, cmbasis, C1, mp_C[i]);
    }
  sumliftCoeff(mp_basisprod, k, mp_C, mp_sum);
  ri = iratrecon(mp_sum, mp_m, mp_nb, mp_db, mp_N[0], mp_Dt[0]);
  if (ri == 0)
    {
      for (i = 0; i < len; i++) { mpz_clear(mp_Dt[i]); } { XFREE(mp_Dt); }
      for (i = 0; i < k; i++) { mpz_clear(mp_C[i]); } { XFREE(mp_C); }
      { mpz_clear(mp_sum);  mpz_clear(mp_m); }
      { XFREE(C1); XFREE(idx); }
      return 0;
    }
  mpz_set(mp_D, mp_Dt[0]);

  /* reconstruct mp_N[i], i = 1..n*m */
  for (i = 1; i < n*m; i++)
    {
      /* transform i to the corresponding index h in C when solupos==LeftSolu*/
      { s = (long)(i/n); t = i-s*n; h = t*m+s; }
      for (j = 0; j < k; j++)
	{
	  if (solupos == RightSolu)
	    for (l = 0; l < basislen; l++) { C1[l] = C[j][l][i]; }
	  else if (solupos == LeftSolu)
	    for (l = 0; l < basislen; l++) { C1[l] = C[j][l][h]; }
	  ChineseRemainderPos(basislen, basis, cmbasis, C1, mp_C[j]);
	}
      sumliftCoeff(mp_basisprod, k, mp_C, mp_sum);

      /* try finding numerator using denominator obtained in last operation */
      r = findNumer(mp_sum, mp_m, mp_D, mp_nb, mp_N[i]);

      /* fail to find the numerator, apply rational reconstruction, mark
	 the additional multiple of denominator, and update mp_D */
      if (r == 0)
	{
	  ++len;
	  mp_Dt = XREALLOC(mpz_t, mp_Dt, len);
	  mpz_init(mp_Dt[len-1]);
	  idx = XREALLOC(long, idx, len);
	  idx[len-1] = i;
	  ri = iratrecon(mp_N[i], mp_m, mp_nb, mp_db, mp_N[i], mp_Dt[len-1]);
	  if (ri == 0)
	    {
	      for (j = 0; j < len; j++) { mpz_clear(mp_Dt[j]); }
	      { XFREE(mp_Dt);  XFREE(C1);  XFREE(idx); }
	      for (j = 0; j < k; j++) { mpz_clear(mp_C[j]); } { XFREE(mp_C); }
	      { mpz_clear(mp_sum);  mpz_clear(mp_m); }
	      return 0;
	    }
	  mpz_mul(mp_D, mp_D, mp_Dt[len-1]);
	}
    }

  /* normalize numerators */
  mpz_init_set(mp_temp, mp_Dt[len-1]);
  for (i = len-2; i >= 0; i--)
    {
      for (j = idx[i]; j < idx[i+1]; j++)
	mpz_mul(mp_N[j], mp_temp, mp_N[j]);
      mpz_mul(mp_temp, mp_temp, mp_Dt[i]);
    }

  for (i = 0; i < len; i++) { mpz_clear(mp_Dt[i]); } { XFREE(mp_Dt); }
  for (i = 0; i < k; i++) { mpz_clear(mp_C[i]); } { XFREE(mp_C); }
  { mpz_clear(mp_sum);  mpz_clear(mp_m);  mpz_clear(mp_temp); }
  { XFREE(C1); XFREE(idx); }
  return 1;
}


/*
 * Calling Sequence:
 *   1/0 <-- findNumer(mp_u, mp_m, mp_D, mp_nb, mp_N)
 *
 * Summary:
 *   Certify correctness of the input denominator and, upon success, compute
 *   the numerator
 *
 * Description:
 *   Given a possible denominator mp_D, numerator bound mp_nb, the function
 *   certifies the denominator by computing mp_N = mp_u*mp_D mod mp_m,
 *   if abs(mp_N) < mp_nb, then mp_D is the real denominator and mp_N is the
 *   real numerator. In this case, the function returns 1. Otherwise, the
 *   function returns 0.
 *
 * Input:
 *    mp_u: mpz_t, image
 *    mp_m: mpz_t, modulus
 *    mp_D: mpz_t, the possible input denominator
 *   mp_nb: mpz_t, numerator bound of rational number
 *
 * Output:
 *   mp_N: mpz_t, if the function returns 1, then mp_N stores the numerator
 *         of the reconstructed rational number
 *
 * Return:
 *   -1, if mp_N is the real numerator
 *   -0, if mp_N is not the real numerator
 *
 * Note:
 *   mp_N will be modified no matter it is really numerator or not.
 *
 */

long
findNumer (const mpz_t mp_u, const mpz_t mp_m, const mpz_t mp_D, \
	   const mpz_t mp_nb, mpz_t mp_N)
{
  mpz_mul(mp_N, mp_u, mp_D);
  mpz_mod(mp_N, mp_N, mp_m);
  if (mpz_cmpabs(mp_N, mp_nb) > 0) { return 0; }
  else { return 1; }
}


} // iml
} // LinBox

#endif // __LINBOX_algorithm_iml_reconstruct_solution_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

