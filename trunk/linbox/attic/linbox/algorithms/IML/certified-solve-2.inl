#ifndef __LINBOX_algorithm_iml_certified_solve_INL
#define __LINBOX_algorithm_iml_certified_solve_INL

#define NULLSPACE_COLUMN 10











/*
 * Calling Sequence:
 *   1/0 <-- specialminSolveLong(certflag, redflag, nullcol, n, m, mp_bdC,
 *                               C, mp_b, mp_N, mp_D, mp_NZ, mp_DZ)
 *
 * Summary:
 *   Compute the minimal denominator solution of a full row rank system,
 *   represented by signed long integers, and corresponding certificate
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
 *   The certificate vector z satisfies that z.C is an integer vector and
 *   z.b has the same denominator as the solution vector.
 *
 *   The function takes the signed long array C, as input. Based on
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
 *             or not
 *           - certflag = 1, compute the certificate vector
 *           - certflag = 0, not compute the certificate
 *    redflag: 1/0, flag to indicate whether to reduce the solution size or not
 *           - redflag = 1, reduce the solution size
 *           - redflag = 0, not reduce the solution size
 *    nullcol: long, dimension of nullspace and kernel basis of conditioned
 *             system,
 *             if nullcol < NULLSPACE_COLUMN, use NULLSPACE_COLUMN instead
 *          n: long, row dimension of the input system
 *          m: long, column dimension of input matrix C
 *     mp_bdC: mpz_t, the maximum sum of absolute values of entries of
 *             each row in matrix C
 *      basis: 1-dim FiniteField array length basislen, RNS basis
 *          C: 1-dim long array length n*m, input matrix C
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
 *   - The space of the solution (mp_N, mp_D) is needed to be preallocated and
 *     the integer mp_D and the entries in mp_N are needed to be initiated as
 *     any integer values.
 *   - If certflag is specified to be 1, then also needs to preallocate space
 *     for (mp_NZ, mp_DZ). Otherwise, set mp_NZ = NULL, and mp_DZ = any int
 *
 */

long
specialminSolveLong (const long certflag, const long redflag, \
		     const long nullcol, const long n, \
		     const long m, const mpz_t mp_bdC, const long *C, \
		     mpz_t *mp_b,  mpz_t *mp_N, mpz_t mp_D, \
		     mpz_t *mp_NZ, mpz_t mp_DZ)
{
  long i, j, l, q, temp, k, len=1;
  FiniteField p, qh, d;
  mpz_t mp_maxInter, mp_basisprod, mp_temp, mp_temp1;
  long *A, *P, *rp;
  double *cumprod;
  FiniteField *bdcoeff, *basis, **basiscmb;
  Double temp1, *DA, *DB, *DC, *P1, *P2, *dtemp, *CRNS, *AP, **ARNS, **BRNS;
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

  /* in case m-n is small, do not need to compress the system */
  if (m-n <= k)
    {
      mp_Bb = XMALLOC(mpz_t, (m-n+1)*n);
      A = XMALLOC(long, n*n);
      for (i = 0; i < n; i++)
	{
	  for (j = 0; j < n; j++)
	    A[i*n+j] = C[i*m+j];
	  for (j = 0; j < m-n; j++)
	    mpz_init_set_si(mp_Bb[i*(m-n+1)+j], C[i*m+(j+n)]);
	  mpz_init_set(mp_Bb[i*(m-n+1)+m-n], mp_b[i]);
	}

      minSolnoncompressLong(certflag, redflag, n, m-n, mp_Bb, A, mp_N, mp_D, \
			    mp_NZ, mp_DZ);
      for (i = 0; i < (m-n+1)*n; i++){ mpz_clear(mp_Bb[i]); } { XFREE(mp_Bb); }
      XFREE(A);
      return 1;
    }
  mp_Nt = XMALLOC(mpz_t, n+k);
  for (i = 0; i < n+k; i++) { mpz_init(mp_Nt[i]); }
  AP = XMALLOC(Double, n*n);
  P = XMALLOC(long, n+1);
  rp = XMALLOC(long, n+1);
  mp_Bb = XMALLOC(mpz_t, n*(k+1));
  for (i = 0; i < n*(k+1); i++) { mpz_init(mp_Bb[i]); }

  /* if fit into long and mantissa of double, then not compute in RNS */
  if (mpz_cmp_ui(mp_bdC, LongRNSbound()) <= 0)
    {
      DA = XMALLOC(Double, n*n);
      DB = XMALLOC(Double, n*k);
      DC = XMALLOC(Double, n*m);
      for (i = 0; i < m*n; i++) { DC[i] = (Double)C[i]; }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      printf("    Compressed Matrix Represented in Long\n");
      fflush(stdout);
      ti = clock();
#endif

      while (1)
	{
	  P1 = randomDb(m, n, 1);
	  P2 = randomDb(m, k, 1);
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, m, \
		      1.0, DC, m, P1, n, 0.0, DA, n);
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, k, m, \
		      1.0, DC, m, P2, k, 0.0, DB, k);

	  /* check whether the rank does not decrease after compression */
	  for (i = 0; i < n+1; i++) { P[i] = i;  rp[i] = 0; }
	  for (i = 0; i < n*n; i++)
	    AP[i] = (temp1 =fmod(DA[i], p)) >= 0 ? temp1 : p+temp1;
	  d=1;
	  RowEchelonTransform(p, AP, n, n, 0, 0, 0, 0, P, rp, &d);
	  if (rp[0] == n) { break; }
	  else { XFREE(P1); XFREE(P2); }
	}
      A = XMALLOC(long, n*n);
      for (i = 0; i < n*n; i++) { A[i] = (long)DA[i]; }
      for (i = 0; i < n; i++)
	{
	  for (j = 0; j < k; j++)
	    mpz_set_d(mp_Bb[i*(k+1)+j], DB[i*k+j]);
  mpz_set(mp_Bb[i*(k+1)+k], mp_b[i]);
	}
      { XFREE(AP); XFREE(P); XFREE(rp); XFREE(DA); XFREE(DB); XFREE(DC); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
      printf("    Compression Time: %f\n", tt);
      fflush(stdout);
      ti = clock();
#endif

      minSolnoncompressLong(certflag, redflag, n, k, mp_Bb, A, mp_Nt, mp_D, \
			    mp_NZ, mp_DZ);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
      printf("    Solution of Compressed System Finding Time: %f\n", tt);
      fflush(stdout);
#endif
      XFREE(A);
      for (i = 0; i < n*(k+1); i++) { mpz_clear(mp_Bb[i]); } { XFREE(mp_Bb); }
    }
  else
    {
      /* compute maximum mangnitude bound of elements of extended RNS basis */
      qh = RNSbound(m);

      /* compute maximum intermidiate result during the compression */
      mpz_init_set_ui(mp_maxInter, 1);
      mpz_addmul_ui(mp_maxInter, mp_bdC, 2);

      /* compute RNS basis */
      basiscmb = findRNS(qh, mp_maxInter, &len);
      basis = basiscmb[0];
      mpz_clear(mp_maxInter);
      ARNS = XMALLOC(Double *, len);
      BRNS = XMALLOC(Double *, len);
      for (i = 0; i < len; i++)
	{
	  ARNS[i] = XMALLOC(Double, n*n);
	  BRNS[i] = XMALLOC(Double, n*k);
	}
      CRNS = XMALLOC(Double, n*m);
      cumprod = cumProd(len, basis, 1, &p);
      bdcoeff = repBound(len, basiscmb[0], basiscmb[1]);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      printf("    Compressed Matrix Represented in RNS\n");
      fflush(stdout);
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
	      /* extend CRNS to the extended RNS basis */
	      q = (long)basis[i];
	      for (j = 0; j < n*m; j++)
		CRNS[j] = (Double)((temp = (C[j] % q)) >= 0 ? temp : q+temp);
	      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, m, \
			  1.0, CRNS, m, P1, n, 0.0, ARNS[i], n);
	      Dmod(basis[i], ARNS[i], n, n, n);
	      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, k, m, \
			  1.0, CRNS, m, P2, k, 0.0, BRNS[i], k);
	      Dmod(basis[i], BRNS[i], n, k, k);
	    }
	  /* check whether the rank does not decrease after compression */
	  for (i = 0; i < n+1; i++) { P[i] = i;  rp[i] = 0; }
	  basisExt(len, n*n, p, basiscmb[0], basiscmb[1], *cumprod, \
		   bdcoeff, ARNS, AP);
	  RowEchelonTransform(p, AP, n, n, 0, 0, 0, 0, P, rp, &d);
	  if (rp[0] == n) { break; }
	  else { XFREE(P1); XFREE(P2); }
	}
      { XFREE(cumprod); XFREE(CRNS); XFREE(AP); XFREE(P); XFREE(rp); }

      /* use CRT to reconstruct B */
      mpz_init(mp_basisprod);
      basisProd(len, basis, mp_basisprod);
      dtemp = XMALLOC(Double, len);
      for (i = 0; i < n; i++)
	{
	  for (j = 0; j < k; j++)
	    {
	      for (l = 0; l < len; l++) { dtemp[l] = BRNS[l][i*k+j]; }
	      ChineseRemainder(len, mp_basisprod, basiscmb[0], basiscmb[1], \
			       bdcoeff, dtemp, mp_Bb[i*(k+1)+j]);
	    }
	  mpz_set(mp_Bb[i*(k+1)+k], mp_b[i]);
	}
      mpz_clear(mp_basisprod);
      { XFREE(dtemp); XFREE(bdcoeff); }
      for (i = 0; i < len; i++) { XFREE(BRNS[i]); } { XFREE(BRNS); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
      printf("    Compression Time: %f\n", tt);
      fflush(stdout);
      ti = clock();
#endif

      /* compute minimal denominator solution of system APv=b */
      minSolnoncompressRNS(certflag, redflag, n, k, len, basis, mp_Bb, ARNS,\
			   mp_Nt, mp_D, mp_NZ, mp_DZ);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
      printf("    Solution of Compressed System Finding Time: %f\n", tt);
      fflush(stdout);
#endif

      for (i = 0; i < n*(k+1); i++) { mpz_clear(mp_Bb[i]); } { XFREE(mp_Bb); }
      { XFREE(basiscmb[0]); XFREE(basiscmb[1]); XFREE(basiscmb); }
      for (i = 0; i < len; i++) { XFREE(ARNS[i]); } { XFREE(ARNS); }
    }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* verify certificate vector of compressed system in uncompressed system */
  if ((certflag == 1) && (mpz_cmp_ui(mp_D, 1) != 0))
    {
      { mpz_init(mp_temp); mpz_init(mp_temp1); }
      for (i = 0; i < m; i++)
	{
	  mpz_mul_si(mp_temp, mp_NZ[0], C[i]);
	  for (j = 1; j < n; j++)
	    {
	      mpz_set_si(mp_temp1, C[j*m+i]);
	      mpz_addmul(mp_temp, mp_NZ[j], mp_temp1);
	    }
	  if (!mpz_divisible_p(mp_temp, mp_DZ))
	    {
	      /* check failed */
	      for (l = 0; l < n+k; l++) { mpz_clear(mp_Nt[l]); }
	      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
	      { XFREE(P1); XFREE(P2); XFREE(mp_Nt); }
	      return 0;
	    }
	}
      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
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
