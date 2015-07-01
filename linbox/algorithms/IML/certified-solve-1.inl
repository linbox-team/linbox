#ifndef __LINBOX_algorithm_iml_certified_solve_INL
#define __LINBOX_algorithm_iml_certified_solve_INL

#define NULLSPACE_COLUMN 10











/*
 * Calling Sequence:
 *   1/0 <-- specialminSolveRNS(certflag, redflag, nullcol, n, m,
 *                              basislen, mp_bdC, basis, CRNS, mp_b,
 *                              mp_N, mp_D, mp_NZ, mp_DZ)
 *
 * Summary:
 *   Compute the minimal denominator solution of a full row rank system,
 *   represented by mpz_t integers, and corresponding certificate
 *   vector(optional). The solution size could be reduced.
 *
 * Description:
 *   Let the full row rank system be Cx = b, where b is a vector, and C is
 *
 *               m
 *       [     |         ]
 *   C = [  A  |    B    ] n,   m > n
 *       [     |         ]
 *
 *   Then certificate vector z satisfies that z.C is an integer vector and
 *   z.b has the same denominator as the solution vector.
 *
 *   The function takes CRNS, representation of C in a RNS, as input. Based on
 *   function minSolnoncompressRNS, the function can deal with the case
 *   when m >> n. In that case, the function compresses the system, and then
 *   computes the solution and the certificate of the compressed system.
 *   Then the function verifies the certificate of the compressed system over
 *   the uncompressed system as well as converts the solution of the
 *   compressed system to that of uncompressed system. If the verification of
 *   certificate succeeds, the function returns 1. Otherwise, reture 0.
 *
 *   On the other hand, if the nullspace of post-conditioned is designated
 *   larger than m-n, then the system is uncompressed and m-n is used.
 *
 *   If redflag is specified as 1, then reduce the solution by the kernel
 *   basis with dimension min(max(nullcol, NULLSPACE_COLUMN), dimnullspace),
 *   where dimnullspace is the dimension of nullspace of C.
 *
 *
 * Input:
 *   certflag: 1/0, flag to indicate whether to compute the certificate vector
 *             or not. If certflag = 1, then compute the certificate vector.
 *             Otherwise, not compute the certificate.
 *    redflag: 1/0, flag to indicate whether to reduce the solution size or not
 *           - redflag = 1, reduce the solution size
 *           - redflag = 0, not reduce the solution size
 *    nullcol: long, dimension of nullspace and kernel basis of conditioned
 *             system,
 *             if nullcol < NULLSPACE_COLUMN, use NULLSPACE_COLUMN instead
 *          n: long, row dimension of the input system
 *          m: long, column dimension of input matrix C
 *   basislen: FiniteField, dimension of the RNS basis
 *     mp_bdC: mpz_t, the maximum sum of absolute values of entries of
 *             each row in matrix C
 *      basis: 1-dim FiniteField array length basislen, RNS basis
 *       CRNS: 2-dim Double array, dimension basislen x n*m, representation of
 *             C in the RNS. CRNS[i] = C mod basis[i], i = 0..basislen-1
 *       mp_b: 1-dim mpz_t array length n, right hand side vector b
 *
 * Return:
 *   1: verification of the certificate vector succeeds
 *   0: verification of the certificate vector fails, needs to restart
 *
 * Output (only correct when returned value is 1):
 *    mp_N: 1-dim mpz_t array length m, numerator vector of the solution
 *          with minimal denominator
 *    mp_D: mpz_t, minimal denominator of the solution
 *   mp_NZ: 1-dim mpz_t array, the first n entries store the numerator
 *          vector of the certificate z
 *   mp_DZ: mpz_t, the denominator of the certificate z
 *
 * Note:
 *   - The space of the solution (mp_N, mp_D) is needed to preallocated.
 *   - If certflag is specified to be 1, then also needs to preallocate space
 *     for (mp_NZ, mp_DZ). Otherwise, set mp_NZ = NULL, and mp_DZ = any int
 *
 * Precondition:
 *   any element p in RNS basis must satisfy 2*(p-1)^2 <= 2^53-1.
 *
 */

long
specialminSolveRNS (const long certflag, const long redflag, \
		    const long nullcol, const long n, const long m, \
		    const long basislen, const mpz_t mp_bdC, \
		    const FiniteField *basis, Double **CRNS, mpz_t *mp_b, \
		    mpz_t *mp_N, mpz_t mp_D, mpz_t *mp_NZ, mpz_t mp_DZ)
{
  long i, j, l, ver, k, len=1;
  FiniteField p, qh, d;
  mpz_t mp_maxInter, mp_basisprod;
  long *P, *rp;
  double *cumprod, *cumprod1;
  FiniteField *bdcoeff, *cmbasis, *extbdcoeff, **extbasiscmb;
  Double *P1, *P2, *dtemp, *CExtRNS, *AP, **ARNS, **BRNS;
  mpz_t *mp_Bb, *mp_Nt;
  double tt;

#if HAVE_TIME_H
  clock_t ti;
#endif

  /* only when we need to reduce the solution, we may reset number
     of columns */
  k = NULLSPACE_COLUMN;
  if ((redflag == 1) && (nullcol > NULLSPACE_COLUMN)) { k = nullcol; }

  p = RandPrime(15, 19);
  cmbasis = combBasis(basislen, basis);
  bdcoeff = repBound(basislen, basis, cmbasis);
  mpz_init(mp_basisprod);

  /* in case m-n is small, do not need to compress the system */
  if (m-n <= k)
    {
      dtemp = XMALLOC(Double, basislen);
      mp_Bb = XMALLOC(mpz_t, (m-n+1)*n);
      for (i = 0; i < (m-n+1)*n; i++) { mpz_init(mp_Bb[i]); }

      /* use CRT to reconstruct B */
      basisProd(basislen, basis, mp_basisprod);
      for (i = 0; i < n; i++)
	for (j = 0; j < m-n; j++)
	  {
	    for (l = 0; l < basislen; l++) { dtemp[l] = CRNS[l][i*m+(j+n)]; }
	    ChineseRemainder(basislen, mp_basisprod, basis, cmbasis, bdcoeff, \
			     dtemp, mp_Bb[i*(m-n+1)+j]);
	  }
      for (i = 0; i < n; i++) { mpz_set(mp_Bb[i*(m-n+1)+m-n], mp_b[i]); }
      mpz_clear(mp_basisprod);
      { XFREE(dtemp); XFREE(cmbasis); XFREE(bdcoeff); }

      /* copy n x n submatrices ARNS[i] from CRNS[i] */
      ARNS = XMALLOC(Double *, basislen);
      for (i = 0; i < basislen; i++)
	{
	  ARNS[i] = XMALLOC(Double, n*n);
	  for (j = 0; j < n; j++)
	    for (l = 0; l < n; l++)
	      ARNS[i][j*n+l] = CRNS[i][j*m+l];
	}
      minSolnoncompressRNS(certflag, redflag, n, m-n, basislen, basis, mp_Bb, \
			   ARNS, mp_N, mp_D, mp_NZ, mp_DZ);
      for (i = 0; i < (m-n+1)*n; i++) { mpz_clear(mp_Bb[i]); }
      XFREE(mp_Bb);
      for (i = 0; i < basislen; i++) { XFREE(ARNS[i]); } { XFREE(ARNS); }
      return 1;
    }

  /* compute maximum mangnitude bound of elements of extended RNS basis */
  qh = RNSbound(m);

  /* compute maximum intermidiate result during the compression */
  mpz_init_set_ui(mp_maxInter, 1);
  mpz_addmul_ui(mp_maxInter, mp_bdC, 2);

  /* compute extended RNS basis */
  extbasiscmb = findRNS(qh, mp_maxInter, &len);
  { mpz_clear(mp_maxInter); }
  cumprod = cumProd(basislen, basis, len, extbasiscmb[0]);
  ARNS = XMALLOC(Double *, len);
  BRNS = XMALLOC(Double *, len);
  CExtRNS = XMALLOC(Double, n*m);
  cumprod1 = cumProd(len, extbasiscmb[0], 1, &p);
  extbdcoeff = repBound(len, extbasiscmb[0], extbasiscmb[1]);
  AP = XMALLOC(Double, n*n);
  P = XMALLOC(long, n+1);
  rp = XMALLOC(long, n+1);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* compress the large matrix to a small one
     generate random 0, 1 matrices P1, P2 for compression
     note: matrix might be <I_(n+k), *> */
  while (1)
    {
      P1 = randomDb(m, n, 1);
      P2 = randomDb(m, k, 1);
      for (i = 0; i < len; i++)
	{
	  ARNS[i] = XMALLOC(Double, n*n);
	  BRNS[i] = XMALLOC(Double, n*k);

	  /* extend CRNS to the extended RNS basis */
	  basisExt(basislen, n*m, extbasiscmb[0][i], basis, cmbasis, \
		   cumprod[i], bdcoeff, CRNS, CExtRNS);
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, m, \
		      1.0, CExtRNS, m, P1, n, 0.0, ARNS[i], n);
	  Dmod(extbasiscmb[0][i], ARNS[i], n, n, n);
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, k, m, \
		      1.0, CExtRNS, m, P2, k, 0.0, BRNS[i], k);
	  Dmod(extbasiscmb[0][i], BRNS[i], n, k, k);
	}
      /* check whether the rank does not decrease after compression */
      for (i = 0; i < n+1; i++) { P[i] = i;  rp[i] = 0; }
      basisExt(len, n*n, p, extbasiscmb[0], extbasiscmb[1], *cumprod1, \
	       extbdcoeff, ARNS, AP);
      d=1;
      RowEchelonTransform(p, AP, n, n, 0, 0, 0, 0, P, rp, &d);
      if (rp[0] == n) { break; }
      else
	{
	  for (i = 0; i < len; i++) { XFREE(ARNS[i]); XFREE(BRNS[i]); }
	  { XFREE(P1); XFREE(P2); }
	}
    }
  { XFREE(bdcoeff); XFREE(cumprod); XFREE(cmbasis); XFREE(CExtRNS); }
  { XFREE(extbdcoeff); XFREE(cumprod1); XFREE(AP); XFREE(P); XFREE(rp); }

  /* use CRT to reconstruct mpz_t matrix B */
  mp_Bb = XMALLOC(mpz_t, n*(k+1));
  for (i = 0; i < n*(k+1); i++) { mpz_init(mp_Bb[i]); }
  basisProd(len, extbasiscmb[0], mp_basisprod);
  bdcoeff = repBound(len, extbasiscmb[0], extbasiscmb[1]);
  dtemp = XMALLOC(Double, len);
  for (i = 0; i < n; i++)
    for (j = 0; j < k; j++)
      {
	for (l = 0; l < len; l++) { dtemp[l] = BRNS[l][i*k+j]; }
	ChineseRemainder(len, mp_basisprod, extbasiscmb[0], extbasiscmb[1], \
			 bdcoeff, dtemp, mp_Bb[i*(k+1)+j]);
      }
  for (i = 0; i < n; i++) { mpz_set(mp_Bb[i*(k+1)+k], mp_b[i]); }
  mpz_clear(mp_basisprod);
  { XFREE(dtemp); XFREE(bdcoeff); }
  for (i = 0; i < len; i++) { XFREE(BRNS[i]); } { XFREE(BRNS); }
  mp_Nt = XMALLOC(mpz_t, n+k);
  for (i = 0; i < n+k; i++) { mpz_init(mp_Nt[i]); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("    Compression Time: %f\n", tt);
  ti = clock();
#endif

  /* compute minimal denominator solution of system APv=b */
  minSolnoncompressRNS(certflag, redflag, n, k, len, extbasiscmb[0], mp_Bb, \
		       ARNS, mp_Nt, mp_D, mp_NZ, mp_DZ);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("    Solution of Compressed System Finding Time: %f\n", tt);
#endif

  for (i = 0; i < n*(k+1); i++) { mpz_clear(mp_Bb[i]); } { XFREE(mp_Bb); }
  { XFREE(extbasiscmb[0]); XFREE(extbasiscmb[1]); XFREE(extbasiscmb); }
  for (i = 0; i < len; i++) { XFREE(ARNS[i]); } { XFREE(ARNS); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* verify certificate vector of compressed system in uncompressed system */
  if ((certflag == 1) && (mpz_cmp_ui(mp_D, 1) != 0))
    {
      ver = certVerify(basislen, n, m, basis, CRNS, mp_DZ, mp_NZ);

      /* check failed */
      if (ver == 0)
	{
	  XFREE(P1); XFREE(P2);
	  for (i = 0; i < n+k; i++) { mpz_clear(mp_Nt[i]); } { XFREE(mp_Nt); }
	  return 0;
	}
    }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("    Certificate Verification Time: %f\n", tt);
  fflush(stdout);
  ti = clock();
#endif

  /* compute mp_N by Pmp_Nt */
  for (i = 0; i < m; i++)
    {
      mpz_set_ui(mp_N[i], 0);
      for (j = 0; j < n; j++)
	if (P1[i*n+j] == 1) { mpz_add(mp_N[i], mp_N[i], mp_Nt[j]); }
      for (j = 0; j < k; j++)
	if (P2[i*k+j] == 1) { mpz_add(mp_N[i], mp_N[i], mp_Nt[n+j]); }
    }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("    Solution of Uncompressed System Recovering Time: %f\n", tt);
  fflush(stdout);
#endif

  { XFREE(P1); XFREE(P2); }
  for (i = 0; i < n+k; i++) { mpz_clear(mp_Nt[i]); } { XFREE(mp_Nt); }
  return 1;
}
