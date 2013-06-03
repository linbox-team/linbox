#ifndef __LINBOX_algorithm_iml_certified_solve_INL
#define __LINBOX_algorithm_iml_certified_solve_INL

#define NULLSPACE_COLUMN 10


/*
 *
 * Calling Sequence:
 *   1/2/3 <-- certSolveLong(certflag, n, m, A, mp_b, mp_N, mp_D,
 *                           mp_NZ, mp_DZ)
 *
 * Summary:
 *   Certified solve a system of linear equations without reducing the
 *   solution size, where the left hand side input matrix is represented
 *   by signed long integers
 *
 * Description:
 *   Let the system of linear equations be Av = b, where A is a n x m matrix,
 *   and b is a n x 1 vector. There are three possibilities:
 *
 *   1. The system has more than one rational solution
 *   2. The system has a unique rational solution
 *   3. The system has no solution
 *
 *   In the first case, there exist a solution vector v with minimal
 *   denominator and a rational certificate vector z to certify that the
 *   denominator of solution v is really the minimal denominator.
 *
 *   The 1 x n certificate vector z satisfies that z.A is an integer vector
 *   and z.b has the same denominator as the solution vector v.
 *   In this case, the function will output the solution with minimal
 *   denominator and optional certificate vector z (if certflag = 1).
 *
 *   Note: if choose not to compute the certificate vector z, the solution
 *     will not garantee, but with high probability, to be the minimal
 *     denominator solution, and the function will run faster.
 *
 *   In the second case, the function will only compute the unique solution
 *   and the contents in the space for certificate vector make no sense.
 *
 *   In the third case, there exists a certificate vector q to certify that
 *   the system has no solution. The 1 x n vector q satisfies q.A = 0 but
 *   q.b <> 0. In this case, the function will output this certificate vector
 *   q and store it into the same space for certificate z. The value of
 *   certflag also controls whether to output q or not.
 *
 *   Note: if the function returns 3, then the system determinately does not
 *     exist solution, no matter whether to output certificate q or not.
 *
 * Input:
 *   certflag: 1/0, flag to indicate whether or not to compute the certificate
 *             vector z or q.
 *           - If certflag = 1, compute the certificate.
 *           - If certflag = 0, not compute the certificate.
 *          n: long, row dimension of the system
 *          m: long, column dimension of the system
 *          A: 1-dim signed long array length n*m, representation of n x m
 *             matrix A
 *       mp_b: 1-dim mpz_t array length n, representation of n x 1 vector b
 *
 * Return:
 *   1: the first case, system has more than one solution
 *   2: the second case, system has a unique solution
 *   3: the third case, system has no solution
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length m,
 *       - numerator vector of the solution with minimal denominator
 *         in the first case
 *       - numerator vector of the unique solution in the second case
 *       - make no sense in the third case
 *   mp_D: mpz_t,
 *       - minimal denominator of the solutions in the first case
 *       - denominator of the unique solution in the second case
 *       - make no sense in the third case
 *
 * The following will only be computed when certflag = 1
 *   mp_NZ: 1-dim mpz_t array length n,
 *        - numerator vector of the certificate z in the first case
 *        - make no sense in the second case
 *        - numerator vector of the certificate q in the third case
 *   mp_DZ: mpz_t,
 *        - denominator of the certificate z if in the first case
 *        - make no sense in the second case
 *        - denominator of the certificate q in the third case
 *
 * Note:
 *   - The space of (mp_N, mp_D) is needed to be preallocated, and entries in
 *     mp_N and integer mp_D are needed to be initiated as any integer values.
 *   - If certflag is specified to be 1, then also needs to preallocate space
 *     for (mp_NZ, mp_DZ), and initiate integer mp_DZ and entries in mp_NZ to
 *     be any integer values.
 *     Otherwise, set mp_NZ = NULL, and mp_DZ = any integer
 *
 */

long
certSolveLong (const long certflag, const long n, const long m, \
	       const long *A, mpz_t *mp_b, mpz_t *mp_N, mpz_t mp_D, \
	       mpz_t *mp_NZ, mpz_t mp_DZ)
{
  long i, j, l, r, r1, t, idx, cmp, ver, alpha, checkflag, temp;
  FiniteField p, d=1;
  long *C, *P, *rp, *Pt, *rpt, *A11;
  Double *DA, *Db, *Dc;
  mpz_t mp_bd, mp_beta, mp_Du, mp_temp, mp_temp1;
  mpz_t *mp_b1, *mp_A21, *mp_Nu;
  double tt;

#if HAVE_TIME_H
  clock_t ti;
#endif

  alpha = maxMagnLong(A, n, m, m);
  if (!alpha)
    {
      /* in case A is a zero matrix, check vector b */
      mpz_init(mp_beta);
      maxMagnMP(mp_b, n, 1, 1, mp_beta);
      if (!mpz_sgn(mp_beta))
	{
	  /* if b is also a zero matrix, set N = 0, D = 1, NZ = 0, DZ = 1 */
	  for (i = 0; i < m; i++) { mpz_set_ui(mp_N[i], 0); }
	  mpz_set_ui(mp_D, 1);
	  if (certflag == 1)
	    {
	      for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	      mpz_set_ui(mp_DZ, 1);
	    }
	  mpz_clear(mp_beta);
	  return 1;
	}
      else
	{
	  /* compute certificate q, q.A = 0, q.b <> 0 in this case */
	  if (certflag == 1)
	    {
	      for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	      mpz_set_ui(mp_DZ, 1);
	      for (i = 0; i < n; i++)
		if (mpz_sgn(mp_b[i]) != 0)
		  { mpz_set_ui(mp_NZ[i], 1);  break; }
	    }
	  /* if b has non-zero entries, no solution */
	  mpz_clear(mp_beta);
	  return 3;
	}
    }
  while (1)
    {
      p = RandPrime(15, 19); /* choose prime p randomly, 2^15 < p < 2^19 */

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      ti = clock();
      printf("prime: %ld\n", p);
#endif

      /* call RowEchelonTransform to reduce DA = A mod p, and obtain rank r */
      DA = XMALLOC(Double, n*m);
      for (i = 0; i < n*m; i++)
	{ DA[i] = (Double)((temp = (A[i] % ((long) p))) >= 0 ? temp : ((long) p)+temp); }
      P = XMALLOC(long, n+1);
      for (i = 0; i < n+1; i++) { P[i] = i; }
      rp = XCALLOC(long, n+1);
      RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
      r = rp[0];

      /* if rank is 0, then p is a bad prime, restart */
      if (r == 0)
	{

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  printf("rank deficient in decomposition, repeat!\n");
#endif

	  { XFREE(DA); XFREE(P); XFREE(rp); }
	  continue;
	}

      /* compute Pt, rpt satisfying that row Pt[i], column rpt[j] of A will be
	 changed to row i, column j respectively after the decomposition */
      Pt = revseq(r, n, P);
      rpt = revseq(r, m, rp);

      /* decomposite <A|b> mod p */
      Db = XMALLOC(Double, n);
      for (i = 0; i < n; i++)
	Db[i] = (Double)mpz_fdiv_ui(mp_b[Pt[i]], p);
      Dc = XMALLOC(Double, n);
      for (i = r; i < n; i++) { Dc[i] = Db[i]; }

      /* compute first r rows of Dc by DA[1..r, 1..r].Db[1..r] */
      cblas_dgemv(CblasRowMajor, CblasNoTrans, r, r, 1.0, DA, m, Db, 1, 0.0, \
		  Dc, 1);

      /* compute last n-r rows of Dc by DA[r+1..n, 1..r].Db[1..r]+Db[r+1..n] */
      if (r < n)
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n-r, r, 1.0, DA+r*m, m, \
		    Db, 1, 1.0, Dc+r, 1);
      Dmod(p, Dc, n, 1, 1);
      idx = r-1;
      while ((++idx < n) && (Dc[idx] == 0)) ;
      { XFREE(DA); XFREE(Db); XFREE(Dc); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
      printf("Decomposition Time: %f\n", tt);
      fflush(stdout);
#endif

      /* rank of A mod p == rank of <A|b> mod p */
      if (idx == n)
	{

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  ti = clock();
#endif

	  /* compute mp_b1 = P.b */
	  mp_b1 = XMALLOC(mpz_t, n);
	  for (i = 0; i < n; i++) { mpz_init_set(mp_b1[i], mp_b[Pt[i]]); }

	  /* C = [A_11, A_12] */
	  C = XMALLOC(long, r*m);
	  for (i = 0; i < r; i++)
	    for (j = 0; j < m; j++)
	      C[i*m+j] = A[Pt[i]*m+rpt[j]];

	  /* full column rank, the solution is unique */
	  if (r == m)
	    nonsingSolvMM(RightSolu, r, 1, C, mp_b1, mp_N, mp_D);
	  else
	    {
	      /* not in full column rank, compute minimal denominator
		 solution */
	      mpz_init(mp_bd);
	      compressBoundLong(r, m, Pt, A, mp_bd);
	      ver = specialminSolveLong(certflag, 0, NULLSPACE_COLUMN, r, \
					m, mp_bd, C, mp_b1, mp_N, mp_D, \
					mp_NZ, mp_DZ);
	      mpz_clear(mp_bd);
	      if (ver == 0)
		{
		  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); }
		  { XFREE(C); XFREE(mp_b1); }
		  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		  printf("certificate checking fails, repeat!\n");
#endif

		  continue;
		}
	    }
	  XFREE(C);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
	  printf("Minimial Denominator Solution Solving Time: %f\n", tt);
	  fflush(stdout);
	  ti = clock();
#endif

	  /* not in full row rank, need to check whether <A_31|A_32>x = b_3 */
	  if (r < n)
	    {
	      { mpz_init(mp_temp); mpz_init(mp_temp1); }
	      checkflag = 0;
	      for (i = 0; i < n-r; i++)
		{
		  mpz_mul_si(mp_temp, mp_N[0], A[Pt[r+i]*m+rpt[0]]);
		  for (j = 1; j < m; j++)
		    {
		      mpz_set_si(mp_temp1, A[Pt[r+i]*m+rpt[j]]);
		      mpz_addmul(mp_temp, mp_N[j], mp_temp1);
		    }
		  mpz_mul(mp_temp1, mp_D, mp_b1[r+i]);
		  if (mpz_cmp(mp_temp, mp_temp1) != 0)
		    {
		      /* if compare fails, restart */
		      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
		      for (l = 0; l < n; l++) { mpz_clear(mp_b1[i]); }
		      { XFREE(P); XFREE(rp); XFREE(Pt); }
		      { XFREE(rpt); XFREE(mp_b1);}

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		      printf("checking lower n-r rows fails, repeat!\n");
#endif
		      checkflag = 1;
		      break;
		    }
		}
	      if (checkflag == 1) { continue; }
	      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
	    }
	  { XFREE(Pt); XFREE(rpt); }
	  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); } { XFREE(mp_b1); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
	  printf("Solution Checking Time: %f\n", tt);
	  fflush(stdout);
#endif

	  /* recover the solution vector and the certificate vector */
	  if (r < m)
	    {
	      /* compute Qy */
	      for (i = r; i > 0; i--)
		if (rp[i] != i) { mpz_swap(mp_N[i-1], mp_N[rp[i]-1]); }

	      /* set last n-r entries in z be 0 and compute zP */
	      if (certflag == 1)
		{
		  for (i = r; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
		  for (i = r; i > 0; i--)
		    if (P[i] != i) { mpz_swap(mp_NZ[i-1], mp_NZ[P[i]-1]); }
		}
	      { XFREE(P); XFREE(rp); }

	      /* system has more than one solution, the first case */
	      return 1;
	    }
	  else
	    {
	      /* system has a unique solution, the second case */
	      { XFREE(P); XFREE(rp); }
	      return 2;
	    }
	}
      else
	{
	  if ((r == m) && (certflag == 0))
	    {
	      { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }
	      return 3;
	    }
	  P[r+1] = idx+1;     /* update P for permutation of row r+1 */

	  /* compute u = A_21.A_11^(-1) */
	  mpz_init(mp_Du);
	  mp_Nu = XMALLOC(mpz_t, r);
	  for (i = 0; i < r; i++) { mpz_init(mp_Nu[i]); }
	  A11 = XMALLOC(long, r*r);
	  for (i = 0; i < r; i++)
	    for (j = 0; j < r; j++)
	      A11[i*r+j] = A[Pt[i]*m+rpt[j]];
	  mp_A21 = XMALLOC(mpz_t, r);
	  for (i = 0; i < r; i++)
	    mpz_init_set_si(mp_A21[i], A[Pt[idx]*m+rpt[i]]);
	  nonsingSolvMM(LeftSolu, r, 1, A11, mp_A21, mp_Nu, mp_Du);
	  XFREE(A11);
	  for (i = 0; i < r; i++) { mpz_clear(mp_A21[i]); } { XFREE(mp_A21); }
	  if (r < m)
	    {
	      { mpz_init(mp_temp); mpz_init(mp_temp1); }
	      checkflag = 0;

	      /* compare u.A_12 with A_22 */
	      for (i = 0; i < m-r; i++)
		{
		  mpz_mul_si(mp_temp, mp_Nu[0], A[Pt[0]*m+rpt[r+i]]);
		  for (j = 1; j < r; j++)
		    {
		      mpz_set_si(mp_temp1, A[Pt[j]*m+rpt[r+i]]);
		      mpz_addmul(mp_temp, mp_Nu[j], mp_temp1);
		    }
		  mpz_mul_si(mp_temp1, mp_Du, A[Pt[idx]*m+rpt[r+i]]);
		  if (mpz_cmp(mp_temp, mp_temp1) != 0)
		    {
		      /* comparison fails, restart */
		      mpz_clear(mp_Du);
		      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
		      for (i = 0; i < r; i++) { mpz_clear(mp_Nu[i]); }
		      { XFREE(P); XFREE(rp); XFREE(Pt); }
		      { XFREE(rpt); XFREE(mp_Nu); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		      printf("checking no solution case fails, repeat!\n");
#endif
		      checkflag = 1;
		      break;
		    }
		}
	      if (checkflag == 1) { continue; }
	      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
	    }
	  if (certflag == 1)
	    {
	      /* set q = (u, -1, 0, ..., 0), which is store in mp_NZ */
	      mpz_set(mp_DZ, mp_Du);
	      for (i = 0; i < r; i++) { mpz_set(mp_NZ[i], mp_Nu[i]); }
	      mpz_set(mp_NZ[r], mp_Du);
	      mpz_neg(mp_NZ[r], mp_NZ[r]);
	      for (i = r+1; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }

	      /* compute qP */
	      for (i = r+1; i > 0; i--)
		if (P[i] != i) { mpz_swap(mp_NZ[i-1], mp_NZ[P[i]-1]); }
	    }

	  mpz_clear(mp_Du);
	  for (i = 0; i < r; i++) { mpz_clear(mp_Nu[i]); } { XFREE(mp_Nu); }
	  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }

	  /* system has no solution, the third case */
	  return 3;
	}
    }
}




/*
 *
 * Calling Sequence:
 *   1/2/3 <-- certSolveRedLong(certflag, nullcol, n, m, A, mp_b, mp_N, mp_D,
 *                              mp_NZ, mp_DZ)
 *
 * Summary:
 *   Certified solve a system of linear equations and reduce the solution
 *   size, where the left hand side input matrix is represented by signed
 *   long integers
 *
 * Description:
 *   Let the system of linear equations be Av = b, where A is a n x m matrix,
 *   and b is a n x 1 vector. There are three possibilities:
 *
 *   1. The system has more than one rational solution
 *   2. The system has a unique rational solution
 *   3. The system has no solution
 *
 *   In the first case, there exist a solution vector v with minimal
 *   denominator and a rational certificate vector z to certify that the
 *   denominator of solution v is really the minimal denominator.
 *
 *   The 1 x n certificate vector z satisfies that z.A is an integer vector
 *   and z.b has the same denominator as the solution vector v.
 *   In this case, the function will output the solution with minimal
 *   denominator and optional certificate vector z (if certflag = 1).
 *
 *   Note: if choose not to compute the certificate vector z, the solution
 *     will not garantee, but with high probability, to be the minimal
 *     denominator solution, and the function will run faster.
 *
 *   Lattice reduction will be used to reduce the solution size. Parameter
 *   nullcol designates the dimension of kernal basis we use to reduce the
 *   solution size as well as the dimension of nullspace we use to compute
 *   the minimal denominator. The heuristic results show that the solution
 *   size will be reduced by factor 1/nullcol.
 *
 *   To find the minimum denominator as fast as possible, nullcol cannot be
 *   too small. We use NULLSPACE_COLUMN as the minimal value of nullcol. That
 *   is, if the input nullcol is less than NULLSPACE_COLUMN, NULLSPACE_COLUMN
 *   will be used instead. However, if the input nullcol becomes larger, the
 *   function will be slower. Meanwhile, it does not make sense to make
 *   nullcol greater than the dimension of nullspace of the input system.
 *
 *   As a result, the parameter nullcol will not take effect unless
 *   NULLSPACE_COLUMN < nullcol < dimnullspace is satisfied, where
 *   dimnullspace is the dimension of nullspace of the input system.  If the
 *   above condition is not satisfied, the boundary value NULLSPACE_COLUMN or
 *   dimnullspace will be used.
 *
 *   In the second case, the function will only compute the unique solution
 *   and the contents in the space for certificate vector make no sense.
 *
 *   In the third case, there exists a certificate vector q to certify that
 *   the system has no solution. The 1 x n vector q satisfies q.A = 0 but
 *   q.b <> 0. In this case, the function will output this certificate vector
 *   q and store it into the same space for certificate z. The value of
 *   certflag also controls whether to output q or not.
 *
 *   Note: if the function returns 3, then the system determinately does not
 *     exist solution, no matter whether to output certificate q or not.
 *
 * Input:
 *   certflag: 1/0, flag to indicate whether or not to compute the certificate
 *             vector z or q.
 *           - If certflag = 1, compute the certificate.
 *           - If certflag = 0, not compute the certificate.
 *    nullcol: long, dimension of nullspace and kernel basis of conditioned
 *             system,
 *             if nullcol < NULLSPACE_COLUMN, use NULLSPACE_COLUMN instead
 *          n: long, row dimension of the system
 *          m: long, column dimension of the system
 *          A: 1-dim signed long array length n*m, representation of n x m
 *             matrix A
 *       mp_b: 1-dim mpz_t array length n, representation of n x 1 vector b
 *
 * Return:
 *   1: the first case, system has more than one solution
 *   2: the second case, system has a unique solution
 *   3: the third case, system has no solution
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length m,
 *       - numerator vector of the solution with minimal denominator
 *         in the first case
 *       - numerator vector of the unique solution in the second case
 *       - make no sense in the third case
 *   mp_D: mpz_t,
 *       - minimal denominator of the solutions in the first case
 *       - denominator of the unique solution in the second case
 *       - make no sense in the third case
 *
 * The following will only be computed when certflag = 1
 *   mp_NZ: 1-dim mpz_t array length n,
 *        - numerator vector of the certificate z in the first case
 *        - make no sense in the second case
 *        - numerator vector of the certificate q in the third case
 *   mp_DZ: mpz_t,
 *        - denominator of the certificate z if in the first case
 *        - make no sense in the second case
 *        - denominator of the certificate q in the third case
 *
 * Note:
 *   - The space of (mp_N, mp_D) is needed to be preallocated, and entries in
 *     mp_N and integer mp_D are needed to be initiated as any integer values.
 *   - If certflag is specified to be 1, then also needs to preallocate space
 *     for (mp_NZ, mp_DZ), and initiate integer mp_DZ and entries in mp_NZ to
 *     be any integer values.
 *     Otherwise, set mp_NZ = NULL, and mp_DZ = any integer
 *
 */

long
certSolveRedLong (const long certflag, const long nullcol, const long n,
		  const long m, const long *A, mpz_t *mp_b, mpz_t *mp_N,
		  mpz_t mp_D, mpz_t *mp_NZ, mpz_t mp_DZ)
{
  long i, j, l, r, r1, t, idx, cmp, ver, alpha, checkflag, temp;
  FiniteField p, d=1;
  long *C, *P, *rp, *Pt, *rpt, *A11;
  Double *DA, *Db, *Dc;
  mpz_t mp_bd, mp_beta, mp_Du, mp_temp, mp_temp1;
  mpz_t *mp_b1, *mp_A21, *mp_Nu;
  double tt;

#if HAVE_TIME_H
  clock_t ti;
#endif

  alpha = maxMagnLong(A, n, m, m);
  if (!alpha)
    {
      /* in case A is a zero matrix, check vector b */
      mpz_init(mp_beta);
      maxMagnMP(mp_b, n, 1, 1, mp_beta);
      if (!mpz_sgn(mp_beta))
	{
	  /* if b is also a zero matrix, set N = 0, D = 1, NZ = 0, DZ = 1 */
	  for (i = 0; i < m; i++) { mpz_set_ui(mp_N[i], 0); }
	  mpz_set_ui(mp_D, 1);
	  if (certflag == 1)
	    {
	      for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	      mpz_set_ui(mp_DZ, 1);
	    }
	  mpz_clear(mp_beta);
	  return 1;
	}
      else
	{
	  /* compute certificate q, q.A = 0, q.b <> 0 in this case */
	  if (certflag == 1)
	    {
	      for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	      mpz_set_ui(mp_DZ, 1);
	      for (i = 0; i < n; i++)
		if (mpz_sgn(mp_b[i]) != 0)
		  { mpz_set_ui(mp_NZ[i], 1);  break; }
	    }
	  /* if b has non-zero entries, no solution */
	  mpz_clear(mp_beta);
	  return 3;
	}
    }
  while (1)
    {

      p = RandPrime(15, 19); /* choose prime p randomly, 2^15 < p < 2^19 */

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      ti = clock();
      printf("prime: %ld\n", p);
#endif

      /* call RowEchelonTransform to reduce DA = A mod p, and obtain rank r */
      DA = XMALLOC(Double, n*m);
      for (i = 0; i < n*m; i++)
	{ DA[i] = (Double)((temp = (A[i] % ((long) p))) >= 0 ? temp : ((long) p)+temp); }
      P = XMALLOC(long, n+1);
      for (i = 0; i < n+1; i++) { P[i] = i; }
      rp = XCALLOC(long, n+1);
      RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
      r = rp[0];

      /* if rank is 0, then p is a bad prime, restart */
      if (r == 0)
	{

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  printf("rank deficient in decomposition, repeat!\n");
#endif

	  { XFREE(DA); XFREE(P); XFREE(rp); }
	  continue;
	}

      /* compute Pt, rpt satisfying that row Pt[i], column rpt[j] of A will be
	 changed to row i, column j respectively after the decomposition */
      Pt = revseq(r, n, P);
      rpt = revseq(r, m, rp);

      /* decomposite <A|b> mod p */
      Db = XMALLOC(Double, n);
      for (i = 0; i < n; i++)
	Db[i] = (Double)mpz_fdiv_ui(mp_b[Pt[i]], p);
      Dc = XMALLOC(Double, n);
      for (i = r; i < n; i++) { Dc[i] = Db[i]; }

      /* compute first r rows of Dc by DA[1..r, 1..r].Db[1..r] */
      cblas_dgemv(CblasRowMajor, CblasNoTrans, r, r, 1.0, DA, m, Db, 1, 0.0, \
		  Dc, 1);

      /* compute last n-r rows of Dc by DA[r+1..n, 1..r].Db[1..r]+Db[r+1..n] */
      if (r < n)
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n-r, r, 1.0, DA+r*m, m, \
		    Db, 1, 1.0, Dc+r, 1);
      Dmod(p, Dc, n, 1, 1);
      idx = r-1;
      while ((++idx < n) && (Dc[idx] == 0)) ;
      { XFREE(DA); XFREE(Db); XFREE(Dc); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
      printf("Decomposition Time: %f\n", tt);
      fflush(stdout);
#endif

      /* rank of A mod p == rank of <A|b> mod p */
      if (idx == n)
	{

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  ti = clock();
#endif

	  /* compute mp_b1 = P.b */
	  mp_b1 = XMALLOC(mpz_t, n);
	  for (i = 0; i < n; i++) { mpz_init_set(mp_b1[i], mp_b[Pt[i]]); }

	  /* C = [A_11, A_12] */
	  C = XMALLOC(long, r*m);
	  for (i = 0; i < r; i++)
	    for (j = 0; j < m; j++)
	      C[i*m+j] = A[Pt[i]*m+rpt[j]];

	  /* full column rank, the solution is unique */
	  if (r == m)
	    nonsingSolvMM(RightSolu, r, 1, C, mp_b1, mp_N, mp_D);
	  else
	    {
	      /* not in full column rank, compute minimal denominator
		 solution */
	      mpz_init(mp_bd);
	      compressBoundLong(r, m, Pt, A, mp_bd);
	      ver = specialminSolveLong(certflag, 1, nullcol, r, m, mp_bd, \
					C, mp_b1, mp_N, mp_D, mp_NZ, mp_DZ);
	      mpz_clear(mp_bd);
	      if (ver == 0)
		{
		  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); }
		  { XFREE(C); XFREE(mp_b1); }
		  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		  printf("certificate checking fails, repeat!\n");
#endif

		  continue;
		}
	    }
	  XFREE(C);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
	  printf("Minimial Denominator Solution Solving Time: %f\n", tt);
	  fflush(stdout);
	  ti = clock();
#endif

	  /* not in full row rank, need to check whether <A_31|A_32>x = b_3 */
	  if (r < n)
	    {
	      { mpz_init(mp_temp); mpz_init(mp_temp1); }
	      checkflag = 0;
	      for (i = 0; i < n-r; i++)
		{
		  mpz_mul_si(mp_temp, mp_N[0], A[Pt[r+i]*m+rpt[0]]);
		  for (j = 1; j < m; j++)
		    {
		      mpz_set_si(mp_temp1, A[Pt[r+i]*m+rpt[j]]);
		      mpz_addmul(mp_temp, mp_N[j], mp_temp1);
		    }
		  mpz_mul(mp_temp1, mp_D, mp_b1[r+i]);
		  if (mpz_cmp(mp_temp, mp_temp1) != 0)
		    {
		      /* if compare fails, restart */
		      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
		      for (l = 0; l < n; l++) { mpz_clear(mp_b1[i]); }
		      { XFREE(P); XFREE(rp); XFREE(Pt); }
		      { XFREE(rpt); XFREE(mp_b1);}

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		      printf("checking lower n-r rows fails, repeat!\n");
#endif

		      checkflag = 1;
		      break;
		    }
		}
	      if (checkflag == 1) { continue; }
	      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
	    }
	  { XFREE(Pt); XFREE(rpt); }
	  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); } { XFREE(mp_b1); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
	  printf("Solution Checking Time: %f\n", tt);
	  fflush(stdout);
#endif

	  /* recover the solution vector and the certificate vector */
	  if (r < m)
	    {
	      /* compute Qy */
	      for (i = r; i > 0; i--)
		if (rp[i] != i) { mpz_swap(mp_N[i-1], mp_N[rp[i]-1]); }

	      /* set last n-r entries in z be 0 and compute zP */
	      if (certflag == 1)
		{
		  for (i = r; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
		  for (i = r; i > 0; i--)
		    if (P[i] != i) { mpz_swap(mp_NZ[i-1], mp_NZ[P[i]-1]); }
		}
	      { XFREE(P); XFREE(rp); }

	      /* system has more than one solution, the first case */
	      return 1;
	    }
	  else
	    {
	      /* system has a unique solution, the second case */
	      { XFREE(P); XFREE(rp); }
	      return 2;
	    }
	}
      else
	{
	  if ((r == m) && (certflag == 0))
	    {
	      { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }
	      return 3;
	    }
	  P[r+1] = idx+1;     /* update P for permutation of row r+1 */

	  /* compute u = A_21.A_11^(-1) */
	  mpz_init(mp_Du);
	  mp_Nu = XMALLOC(mpz_t, r);
	  for (i = 0; i < r; i++) { mpz_init(mp_Nu[i]); }
	  A11 = XMALLOC(long, r*r);
	  for (i = 0; i < r; i++)
	    for (j = 0; j < r; j++)
	      A11[i*r+j] = A[Pt[i]*m+rpt[j]];
	  mp_A21 = XMALLOC(mpz_t, r);
	  for (i = 0; i < r; i++)
	    mpz_init_set_si(mp_A21[i], A[Pt[idx]*m+rpt[i]]);
	  nonsingSolvMM(LeftSolu, r, 1, A11, mp_A21, mp_Nu, mp_Du);
	  XFREE(A11);
	  for (i = 0; i < r; i++) { mpz_clear(mp_A21[i]); } { XFREE(mp_A21); }
	  if (r < m)
	    {
	      { mpz_init(mp_temp); mpz_init(mp_temp1); }
	      checkflag = 0;

	      /* compare u.A_12 with A_22 */
	      for (i = 0; i < m-r; i++)
		{
		  mpz_mul_si(mp_temp, mp_Nu[0], A[Pt[0]*m+rpt[r+i]]);
		  for (j = 1; j < r; j++)
		    {
		      mpz_set_si(mp_temp1, A[Pt[j]*m+rpt[r+i]]);
		      mpz_addmul(mp_temp, mp_Nu[j], mp_temp1);
		    }
		  mpz_mul_si(mp_temp1, mp_Du, A[Pt[idx]*m+rpt[r+i]]);
		  if (mpz_cmp(mp_temp, mp_temp1) != 0)
		    {
		      /* comparison fails, restart */
		      mpz_clear(mp_Du);
		      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
		      for (i = 0; i < r; i++) { mpz_clear(mp_Nu[i]); }
		      { XFREE(P); XFREE(rp); XFREE(Pt); }
		      { XFREE(rpt); XFREE(mp_Nu); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		      printf("checking no solution case fails, repeat!\n");
#endif

		      checkflag = 1;
		      break;
		    }
		}
	      if (checkflag == 1) { continue; }
	      { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
	    }
	  if (certflag == 1)
	    {
	      /* set q = (u, -1, 0, ..., 0), which is store in mp_NZ */
	      mpz_set(mp_DZ, mp_Du);
	      for (i = 0; i < r; i++) { mpz_set(mp_NZ[i], mp_Nu[i]); }
	      mpz_set(mp_NZ[r], mp_Du);
	      mpz_neg(mp_NZ[r], mp_NZ[r]);
	      for (i = r+1; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }

	      /* compute qP */
	      for (i = r+1; i > 0; i--)
		if (P[i] != i) { mpz_swap(mp_NZ[i-1], mp_NZ[P[i]-1]); }
	    }

	  mpz_clear(mp_Du);
	  for (i = 0; i < r; i++) { mpz_clear(mp_Nu[i]); } { XFREE(mp_Nu); }
	  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }

	  /* system has no solution, the third case */
	  return 3;
	}
    }
}





/*
 *
 * Calling Sequence:
 *   1/2/3 <-- certSolveMP(certflag, n, m, mp_A, mp_b, mp_N, mp_D,
 *                         mp_NZ, mp_DZ)
 *
 * Summary:
 *   Certified solve a system of linear equations without reducing the
 *   solution size, where the left hand side input matrix is represented
 *   by mpz_t integers
 *
 * Description:
 *   Let the system of linear equations be Av = b, where A is a n x m matrix,
 *   and b is a n x 1 vector. There are three possibilities:
 *
 *   1. The system has more than one rational solution
 *   2. The system has a unique rational solution
 *   3. The system has no solution
 *
 *   In the first case, there exist a solution vector v with minimal
 *   denominator and a rational certificate vector z to certify that the
 *   denominator of solution v is really the minimal denominator.
 *
 *   The 1 x n certificate vector z satisfies that z.A is an integer vector
 *   and z.b has the same denominator as the solution vector v.
 *   In this case, the function will output the solution with minimal
 *   denominator and optional certificate vector z (if certflag = 1).
 *
 *   Note: if choose not to compute the certificate vector z, the solution
 *     will not garantee, but with high probability, to be the minimal
 *     denominator solution, and the function will run faster.
 *
 *   In the second case, the function will only compute the unique solution
 *   and the contents in the space for certificate vector make no sense.
 *
 *   In the third case, there exists a certificate vector q to certify that
 *   the system has no solution. The 1 x n vector q satisfies q.A = 0 but
 *   q.b <> 0. In this case, the function will output this certificate vector
 *   q and store it into the same space for certificate z. The value of
 *   certflag also controls whether to output q or not.
 *
 *   Note: if the function returns 3, then the system determinately does not
 *     exist solution, no matter whether to output certificate q or not.
 *   In the first case, there exist a solution vector v with minimal
 *   denominator and a rational certificate vector z to certify that the
 *   denominator of solution v is the minimal denominator.
 *
 * Input:
 *   certflag: 1/0, flag to indicate whether or not to compute the certificate
 *             vector z or q.
 *           - If certflag = 1, compute the certificate.
 *           - If certflag = 0, not compute the certificate.
 *          n: long, row dimension of the system
 *          m: long, column dimension of the system
 *       mp_A: 1-dim mpz_t array length n*m, representation of n x m matrix A
 *       mp_b: 1-dim mpz_t array length n, representation of n x 1 vector b
 *
 * Return:
 *   1: the first case, system has more than one solution
 *   2: the second case, system has a unique solution
 *   3: the third case, system has no solution
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length m,
 *       - numerator vector of the solution with minimal denominator
 *         in the first case
 *       - numerator vector of the unique solution in the second case
 *       - make no sense in the third case
 *   mp_D: mpz_t,
 *       - minimal denominator of the solutions in the first case
 *       - denominator of the unique solution in the second case
 *       - make no sense in the third case
 *
 * The following will only be computed when certflag = 1
 *   mp_NZ: 1-dim mpz_t array length n,
 *        - numerator vector of the certificate z in the first case
 *        - make no sense in the second case
 *        - numerator vector of the certificate q in the third case
 *   mp_DZ: mpz_t,
 *        - denominator of the certificate z if in the first case
 *        - make no sense in the second case
 *        - denominator of the certificate q in the third case
 *
 * Note:
 *   - The space of (mp_N, mp_D) is needed to be preallocated, and entries in
 *     mp_N and integer mp_D are needed to be initiated as any integer values.
 *   - If certflag is specified to be 1, then also needs to preallocate space
 *     for (mp_NZ, mp_DZ), and initiate integer mp_DZ and entries in mp_NZ to
 *     be any integer values.
 *     Otherwise, set mp_NZ = NULL, and mp_DZ = any integer
 *
 */

long
certSolveMP (const long certflag, const long n, const long m, \
		  mpz_t *mp_A, mpz_t *mp_b, mpz_t *mp_N, mpz_t mp_D, \
		  mpz_t *mp_NZ, mpz_t mp_DZ)
{
  long i, j, l, r, r1, t, idx, cmp, ver, basislen=1, checkflag;
  FiniteField qh, p, d=1;
  long *P, *rp, *Pt, *rpt;
  FiniteField *basis, **basiscmb;
  Double *DA, *Db, *Dc, **CRNS, **C1RNS, **A11RNS, **A12RNS;
  mpz_t mp_maxInter, mp_bd, mp_alpha, mp_beta, mp_Du;
  mpz_t *mp_b1, *mp_A21, *mp_A22, *mp_Nu;
  double tt;

#if HAVE_TIME_H
  clock_t ti;
#endif

  mpz_init(mp_alpha);
  maxMagnMP(mp_A, n, m, m, mp_alpha);
  if (!mpz_sgn(mp_alpha))
    {
      /* in case A is a zero matrix, check vector b */
      mpz_init(mp_beta);
      maxMagnMP(mp_b, n, 1, 1, mp_beta);
      if (!mpz_sgn(mp_beta))
	{
	  /* if b is also a zero matrix, set N = 0, D = 1, NZ = 0, DZ = 1 */
	  for (i = 0; i < m; i++) { mpz_set_ui(mp_N[i], 0); }
	  mpz_set_ui(mp_D, 1);
	  if (certflag == 1)
	    {
	      for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	      mpz_set_ui(mp_DZ, 1);
	    }
	  { mpz_clear(mp_alpha); mpz_clear(mp_beta); }
	  return 1;
	}
      else
	{
	  /* compute certificate q, q.A = 0, q.b <> 0 in this case */
	  if (certflag == 1)
	    {
	      for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	      mpz_set_ui(mp_DZ, 1);
	      for (i = 0; i < n; i++)
		if (mpz_sgn(mp_b[i]) != 0)
		  { mpz_set_ui(mp_NZ[i], 1);  break; }
	    }
	  /* if b has non-zero entries, no solution */
	  mpz_clear(mp_alpha); mpz_clear(mp_beta);
	  return 3;
	}
    }
  /* generate RNS basis */
  mpz_init_set_ui(mp_maxInter, 1);
  mpz_addmul_ui(mp_maxInter, mp_alpha, 2);
  qh = RNSbound(m);
  basiscmb = findRNS(qh, mp_maxInter, &basislen);
  basis = basiscmb[0];
  { mpz_clear(mp_maxInter); mpz_clear(mp_alpha); }

  while (1)
    {
      p = RandPrime(15, 19); /* choose p randomly, 2^15 < p < 2^19 */

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      ti = clock();
      printf("prime: %ld\n", p);
#endif

      /* call RowEchelonTransform to reduce DA = A mod p, and obtain rank r */
      DA = XMALLOC(Double, n*m);
      for (i = 0; i < n*m; i++) { DA[i] = (Double)mpz_fdiv_ui(mp_A[i], p); }
      P = XMALLOC(long, n+1);
      for (i = 0; i < n+1; i++) { P[i] = i; }
      rp = XCALLOC(long, n+1);
      RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
      r = rp[0];

      /* if rank is 0, then p is a bad prime, restart */
      if (r == 0)
	{

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  printf("rank deficient in decomposition, repeat!\n");
#endif

	  { XFREE(DA); XFREE(P); XFREE(rp); }
	  continue;
	}

      /* compute Pt, rpt satisfying that row Pt[i], column rpt[j] of A will be
	 changed to row i, column j respectively after the decomposition */
      Pt = revseq(r, n, P);
      rpt = revseq(r, m, rp);

      /* decomposite <A|b> mod p */
      Db = XMALLOC(Double, n);
      for (i = 0; i < n; i++)
	Db[i] = (Double)mpz_fdiv_ui(mp_b[Pt[i]], p);
      Dc = XMALLOC(Double, n);
      for (i = r; i < n; i++) { Dc[i] = Db[i]; }

      /* compute first r rows of Dc by DA[1..r, 1..r].Db[1..r] */
      cblas_dgemv(CblasRowMajor, CblasNoTrans, r, r, 1.0, DA, m, Db, 1, 0.0, \
		  Dc, 1);

      /* compute last n-r rows of Dc by DA[r+1..n, 1..r].Db[1..r]+Db[r+1..n] */
      if (r < n)
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n-r, r, 1.0, DA+r*m, m, \
		    Db, 1, 1.0, Dc+r, 1);
      Dmod(p, Dc, n, 1, 1);
      idx = r-1;
      while ((++idx < n) && (Dc[idx] == 0)) ;
      { XFREE(DA); XFREE(Db); XFREE(Dc); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
      printf("Decomposition Time: %f\n", tt);
#endif

      /* rank of A mod p == rank of <A|b> mod p */
      if (idx == n)
	{

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  ti = clock();
#endif

	  /* compute mp_b1 = P.b */
	  mp_b1 = XMALLOC(mpz_t, n);
	  for (i = 0; i < n; i++) { mpz_init_set(mp_b1[i], mp_b[Pt[i]]); }

	  /* CRNS[i] = [A_11, A_12] mod basis[i] */
	  CRNS = XMALLOC(Double *, basislen);
	  for (i = 0; i < basislen; i++)
	    {
	      CRNS[i] = XMALLOC(Double, r*m);
	      for (j = 0; j < r; j++)
		for (l = 0; l < m; l++)
		  CRNS[i][j*m+l] = (Double)mpz_fdiv_ui(mp_A[Pt[j]*m+rpt[l]], \
						       basis[i]);
	    }
	  /* full column rank, the solution is unique */
	  if (r == m)
	    nonsingSolvRNSMM(RightSolu, r, 1, basislen, basis, CRNS, mp_b1, \
			     mp_N, mp_D);
	  else
	    {
	      /* not in full column rank, compute minimal denominator solution
		 compute the bound after compression in specialminSolveRNS */
	      mpz_init(mp_bd);
	      compressBoundMP(r, m, Pt, mp_A, mp_bd);
	      ver = specialminSolveRNS(certflag, 0, NULLSPACE_COLUMN, r, m,\
				       basislen, mp_bd, basis, CRNS, mp_b1, \
				       mp_N, mp_D, mp_NZ, mp_DZ);
	      mpz_clear(mp_bd);
	      if (ver == 0)
		{
		  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); }
		  for (i = 0; i < basislen; i++) { XFREE(CRNS[i]); }
		  { XFREE(CRNS); XFREE(mp_b1); }
		  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		  printf("certificate checking fails, repeat!\n");
#endif

		  continue;
		}
	    }
	  for (i = 0; i < basislen; i++) { XFREE(CRNS[i]); } { XFREE(CRNS); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
	  printf("Minimial Denominator Solution Solving Time: %f\n", tt);
	  ti = clock();
#endif

	  /* not in full row rank, need to check whether <A_31|A_32>x = b_3 */
	  if (r < n)
	    {
	      C1RNS = XMALLOC(Double *, basislen);
	      for (i = 0; i < basislen; i++)
		{
		  C1RNS[i] = XMALLOC(Double, (n-r)*m);
		  for (j = 0; j < n-r; j++)
		    for (l = 0; l < m; l++)
		      C1RNS[i][j*m+l] = \
			(Double)mpz_fdiv_ui(mp_A[Pt[r+j]*m+rpt[l]], basis[i]);
		}

	      cmp = LVecSMatMulCmp(RightMul, basislen, n-r, m, basis, C1RNS, \
				   mp_D, mp_N, mp_b1+r);
	      for (i = 0; i < basislen; i++) { XFREE(C1RNS[i]); }
	      XFREE(C1RNS);

	      /* if compare fails, restart */
	      if (cmp == 0)
		{
		  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); }
		  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); XFREE(mp_b1); }
#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		  printf("checking lower n-r rows fails, repeat!\n");
#endif

		  continue;
		}
	    }
	  { XFREE(Pt); XFREE(rpt); }
	  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); } { XFREE(mp_b1); }
	  { XFREE(basiscmb[0]); XFREE(basiscmb[1]); XFREE(basiscmb); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
	  printf("Solution Checking Time: %f\n", tt);
#endif

	  /* recover the solution vector and the certificate vector */
	  if (r < m)
	    {
	      /* system has more than one solution, the first case */
	      /* compute Qy */
	      for (i = r; i > 0; i--)
		if (rp[i] != i) { mpz_swap(mp_N[i-1], mp_N[rp[i]-1]); }

	      /* set last n-r entries in z be 0 and compute zP */
	      if (certflag == 1)
		{
		  for (i = r; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
		  for (i = r; i > 0; i--)
		    if (P[i] != i) { mpz_swap(mp_NZ[i-1], mp_NZ[P[i]-1]); }
		}
	      { XFREE(P); XFREE(rp); }
	      return 1;
	    }
	  else
	    {
	      /* system has a unique solution, the second case */
	      { XFREE(P); XFREE(rp); }
	      return 2;
	    }
	}
      else
	{
	  if ((r == m) && (certflag == 0))
	    {
	      { XFREE(basiscmb[0]); XFREE(basiscmb[1]); XFREE(basiscmb); }
	      { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }
	      return 3;
	    }
	  P[r+1] = idx+1;      /* update P for permutation of row r+1 */

	  /* compute u = A_21.A_11^(-1) */
	  mpz_init(mp_Du);
	  mp_Nu = XMALLOC(mpz_t, r);
	  for (i = 0; i < r; i++) { mpz_init(mp_Nu[i]); }
	  A11RNS = XMALLOC(Double *, basislen);
	  for (i = 0; i < basislen; i++)
	    {
	      A11RNS[i] = XMALLOC(Double, r*r);
	      for (j = 0; j < r; j++)
		for (l = 0; l < r; l++)
		  A11RNS[i][j*r+l] = \
		    (Double)mpz_fdiv_ui(mp_A[Pt[j]*m+rpt[l]], basis[i]);
	    }
	  mp_A21 = XMALLOC(mpz_t, r);
	  for (i = 0; i < r; i++)
	    mpz_init_set(mp_A21[i], mp_A[Pt[idx]*m+rpt[i]]);
	  nonsingSolvRNSMM(LeftSolu, r, 1, basislen, basis, A11RNS, mp_A21, \
			   mp_Nu, mp_Du);
	  for (i = 0; i < basislen; i++)  { XFREE(A11RNS[i]); }
	  for (i = 0; i < r; i++) { mpz_clear(mp_A21[i]); }
	  { XFREE(mp_A21); XFREE(A11RNS); }
	  if (r < m)
	    {
	      /* compare u.A_12 with A_22 */
	      A12RNS = XMALLOC(Double *, basislen);
	      for (i = 0; i < basislen; i++)
		{
		  A12RNS[i] = XMALLOC(Double, r*(m-r));
		  for (j = 0; j < r; j++)
		    for (l = 0; l < m-r; l++)
		      A12RNS[i][j*(m-r)+l] = \
			(Double)mpz_fdiv_ui(mp_A[Pt[j]*m+rpt[r+l]], basis[i]);
		}
	      mp_A22 = XMALLOC(mpz_t, m-r);
	      for (i = 0; i < m-r; i++)
		mpz_init_set(mp_A22[i], mp_A[Pt[idx]*m+rpt[r+i]]);
	      cmp = LVecSMatMulCmp(LeftMul, basislen, r, m-r, basis, A12RNS, \
				   mp_Du, mp_Nu, mp_A22);
	      for (i = 0; i < basislen; i++) { XFREE(A12RNS[i]); }
	      for (i = 0; i < m-r; i++) { mpz_clear(mp_A22[i]); }
	      { XFREE(A12RNS); XFREE(mp_A22); }

	      /* comparison fails, restart */
	      if (cmp == 0)
		{
		  mpz_clear(mp_Du);
		  for (i = 0; i < r; i++) { mpz_clear(mp_Nu[i]); }
		  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); XFREE(mp_Nu); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		  printf("checking no solution case fails, repeat!\n");
#endif

		  continue;
		}
	    }
	  if (certflag == 1)
	    {
	      /* set q = (u, -1, 0, ..., 0), which is store in mp_NZ */
	      mpz_set(mp_DZ, mp_Du);
	      for (i = 0; i < r; i++) { mpz_set(mp_NZ[i], mp_Nu[i]); }
	      mpz_set(mp_NZ[r], mp_Du);
	      mpz_neg(mp_NZ[r], mp_NZ[r]);
	      for (i = r+1; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }

	      /* compute qP */
	      for (i = r+1; i > 0; i--)
		if (P[i] != i) { mpz_swap(mp_NZ[i-1], mp_NZ[P[i]-1]); }
	    }

	  mpz_clear(mp_Du);
	  { XFREE(basiscmb[0]); XFREE(basiscmb[1]); XFREE(basiscmb); }
	  for (i = 0; i < r; i++) { mpz_clear(mp_Nu[i]); } { XFREE(mp_Nu); }
	  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }

	  /* system has no solution, the third case */
	  return 3;
	}
    }
}



/*
 *
 * Calling Sequence:
 *   1/2/3 <-- certSolveRedMP(certflag, nullcol, n, m, mp_A, mp_b, mp_N, mp_D,
 *                            mp_NZ, mp_DZ)
 *
 * Summary:
 *   Certified solve a system of linear equations and reduce the solution
 *   size, where the left hand side input matrix is represented by signed
 *   mpz_t integers
 *
 * Description:
 *   Let the system of linear equations be Av = b, where A is a n x m matrix,
 *   and b is a n x 1 vector. There are three possibilities:
 *
 *   1. The system has more than one rational solution
 *   2. The system has a unique rational solution
 *   3. The system has no solution
 *
 *   In the first case, there exist a solution vector v with minimal
 *   denominator and a rational certificate vector z to certify that the
 *   denominator of solution v is really the minimal denominator.
 *
 *   The 1 x n certificate vector z satisfies that z.A is an integer vector
 *   and z.b has the same denominator as the solution vector v.
 *   In this case, the function will output the solution with minimal
 *   denominator and optional certificate vector z (if certflag = 1).
 *
 *   Note: if choose not to compute the certificate vector z, the solution
 *     will not garantee, but with high probability, to be the minimal
 *     denominator solution, and the function will run faster.
 *
 *   Lattice reduction will be used to reduce the solution size. Parameter
 *   nullcol designates the dimension of kernal basis we use to reduce the
 *   solution size as well as the dimension of nullspace we use to compute
 *   the minimal denominator. The heuristic results show that the solution
 *   size will be reduced by factor 1/nullcol.
 *
 *   To find the minimum denominator as fast as possible, nullcol cannot be
 *   too small. We use NULLSPACE_COLUMN as the minimal value of nullcol. That
 *   is, if the input nullcol is less than NULLSPACE_COLUMN, NULLSPACE_COLUMN
 *   will be used instead. However, if the input nullcol becomes larger, the
 *   function will be slower. Meanwhile, it does not make sense to make
 *   nullcol greater than the dimension of nullspace of the input system.
 *
 *   As a result, the parameter nullcol will not take effect unless
 *   NULLSPACE_COLUMN < nullcol < dimnullspace is satisfied, where
 *   dimnullspace is the dimension of nullspace of the input system.  If the
 *   above condition is not satisfied, the boundary value NULLSPACE_COLUMN or
 *   dimnullspace will be used.
 *
 *   In the second case, the function will only compute the unique solution
 *   and the contents in the space for certificate vector make no sense.
 *
 *   In the third case, there exists a certificate vector q to certify that
 *   the system has no solution. The 1 x n vector q satisfies q.A = 0 but
 *   q.b <> 0. In this case, the function will output this certificate vector
 *   q and store it into the same space for certificate z. The value of
 *   certflag also controls whether to output q or not.
 *
 *   Note: if the function returns 3, then the system determinately does not
 *     exist solution, no matter whether to output certificate q or not.
 *   In the first case, there exist a solution vector v with minimal
 *   denominator and a rational certificate vector z to certify that the
 *   denominator of solution v is the minimal denominator.
 *
 * Input:
 *   certflag: 1/0, flag to indicate whether or not to compute the certificate
 *             vector z or q.
 *           - If certflag = 1, compute the certificate.
 *           - If certflag = 0, not compute the certificate.
 *    nullcol: long, dimension of nullspace and kernel basis of conditioned
 *             system,
 *             if nullcol < NULLSPACE_COLUMN, use NULLSPACE_COLUMN instead
 *          n: long, row dimension of the system
 *          m: long, column dimension of the system
 *       mp_A: 1-dim mpz_t array length n*m, representation of n x m matrix A
 *       mp_b: 1-dim mpz_t array length n, representation of n x 1 vector b
 *
 * Return:
 *   1: the first case, system has more than one solution
 *   2: the second case, system has a unique solution
 *   3: the third case, system has no solution
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length m,
 *       - numerator vector of the solution with minimal denominator
 *         in the first case
 *       - numerator vector of the unique solution in the second case
 *       - make no sense in the third case
 *   mp_D: mpz_t,
 *       - minimal denominator of the solutions in the first case
 *       - denominator of the unique solution in the second case
 *       - make no sense in the third case
 *
 * The following will only be computed when certflag = 1
 *   mp_NZ: 1-dim mpz_t array length n,
 *        - numerator vector of the certificate z in the first case
 *        - make no sense in the second case
 *        - numerator vector of the certificate q in the third case
 *   mp_DZ: mpz_t,
 *        - denominator of the certificate z if in the first case
 *        - make no sense in the second case
 *        - denominator of the certificate q in the third case
 *
 * Note:
 *   - The space of (mp_N, mp_D) is needed to be preallocated, and entries in
 *     mp_N and integer mp_D are needed to be initiated as any integer values.
 *   - If certflag is specified to be 1, then also needs to preallocate space
 *     for (mp_NZ, mp_DZ), and initiate integer mp_DZ and entries in mp_NZ to
 *     be any integer values.
 *     Otherwise, set mp_NZ = NULL, and mp_DZ = any integer
 *
 */

long
certSolveRedMP (const long certflag, const long nullcol, const long n, \
		const long m, mpz_t *mp_A, mpz_t *mp_b, \
		mpz_t *mp_N, mpz_t mp_D, mpz_t *mp_NZ, mpz_t mp_DZ)
{
  long i, j, l, r, r1, t, idx, cmp, ver, basislen=1;
  FiniteField qh, p, d=1;
  long *P, *rp, *Pt, *rpt;
  FiniteField *basis, **basiscmb;
  Double *DA, *Db, *Dc, **CRNS, **C1RNS, **A11RNS, **A12RNS;
  mpz_t mp_maxInter, mp_bd, mp_alpha, mp_beta, mp_Du;
  mpz_t *mp_b1, *mp_A21, *mp_A22, *mp_Nu;
  double tt;

#if HAVE_TIME_H
  clock_t ti;
#endif

  mpz_init(mp_alpha);
  maxMagnMP(mp_A, n, m, m, mp_alpha);
  if (!mpz_sgn(mp_alpha))
    {
      /* in case A is a zero matrix, check vector b */
      mpz_init(mp_beta);
      maxMagnMP(mp_b, n, 1, 1, mp_beta);
      if (!mpz_sgn(mp_beta))
	{
	  /* if b is also a zero matrix, set N = 0, D = 1, NZ = 0, DZ = 1 */
	  for (i = 0; i < m; i++) { mpz_set_ui(mp_N[i], 0); }
	  mpz_set_ui(mp_D, 1);
	  if (certflag == 1)
	    {
	      for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	      mpz_set_ui(mp_DZ, 1);
	    }
	  { mpz_clear(mp_alpha); mpz_clear(mp_beta); }
	  return 1;
	}
      else
	{
	  /* compute certificate q, q.A = 0, q.b <> 0 in this case */
	  if (certflag == 1)
	    {
	      for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	      mpz_set_ui(mp_DZ, 1);
	      for (i = 0; i < n; i++)
		if (mpz_sgn(mp_b[i]) != 0)
		  { mpz_set_ui(mp_NZ[i], 1);  break; }
	    }
	  /* if b has non-zero entries, no solution */
	  mpz_clear(mp_alpha); mpz_clear(mp_beta);
	  return 3;
	}
    }
  /* generate RNS basis */
  mpz_init_set_ui(mp_maxInter, 1);
  mpz_addmul_ui(mp_maxInter, mp_alpha, 2);
  qh = RNSbound(m);
  basiscmb = findRNS(qh, mp_maxInter, &basislen);
  basis = basiscmb[0];
  { mpz_clear(mp_maxInter); mpz_clear(mp_alpha); }
  while (1)
    {
      p = RandPrime(15, 19); /* choose p randomly, 2^15 < p < 2^19 */

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      ti = clock();
      printf("prime: %ld\n", p);
#endif

      /* call RowEchelonTransform to reduce DA = A mod p, and obtain rank r */
      DA = XMALLOC(Double, n*m);
      for (i = 0; i < n*m; i++) { DA[i] = (Double)mpz_fdiv_ui(mp_A[i], p); }
      P = XMALLOC(long, n+1);
      for (i = 0; i < n+1; i++) { P[i] = i; }
      rp = XCALLOC(long, n+1);
      RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
      r = rp[0];

      /* if rank is 0, then p is a bad prime, restart */
      if (r == 0)
	{

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  printf("rank deficient in decomposition, repeat!\n");
#endif

	  { XFREE(DA); XFREE(P); XFREE(rp); }
	  continue;
	}

      /* compute Pt, rpt satisfying that row Pt[i], column rpt[j] of A will be
	 changed to row i, column j respectively after the decomposition */
      Pt = revseq(r, n, P);
      rpt = revseq(r, m, rp);

      /* decomposite <A|b> mod p */
      Db = XMALLOC(Double, n);
      for (i = 0; i < n; i++)
	Db[i] = (Double)mpz_fdiv_ui(mp_b[Pt[i]], p);
      Dc = XMALLOC(Double, n);
      for (i = r; i < n; i++) { Dc[i] = Db[i]; }

      /* compute first r rows of Dc by DA[1..r, 1..r].Db[1..r] */
      cblas_dgemv(CblasRowMajor, CblasNoTrans, r, r, 1.0, DA, m, Db, 1, 0.0, \
		  Dc, 1);

      /* compute last n-r rows of Dc by DA[r+1..n, 1..r].Db[1..r]+Db[r+1..n] */
      if (r < n)
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n-r, r, 1.0, DA+r*m, m, \
		    Db, 1, 1.0, Dc+r, 1);
      Dmod(p, Dc, n, 1, 1);
      idx = r-1;
      while ((++idx < n) && (Dc[idx] == 0)) ;
      { XFREE(DA); XFREE(Db); XFREE(Dc); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
      printf("Decomposition Time: %f\n", tt);
#endif

      /* rank of A mod p == rank of <A|b> mod p */
      if (idx == n)
	{

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  ti = clock();
#endif

	  /* compute mp_b1 = P.b */
	  mp_b1 = XMALLOC(mpz_t, n);
	  for (i = 0; i < n; i++) { mpz_init_set(mp_b1[i], mp_b[Pt[i]]); }

	  /* CRNS[i] = [A_11, A_12] mod basis[i] */
	  CRNS = XMALLOC(Double *, basislen);
	  for (i = 0; i < basislen; i++)
	    {
	      CRNS[i] = XMALLOC(Double, r*m);
	      for (j = 0; j < r; j++)
		for (l = 0; l < m; l++)
		  CRNS[i][j*m+l] = (Double)mpz_fdiv_ui(mp_A[Pt[j]*m+rpt[l]], \
						       basis[i]);
	    }
	  /* full column rank, the solution is unique */
	  if (r == m)
	    nonsingSolvRNSMM(RightSolu, r, 1, basislen, basis, CRNS, mp_b1, \
			     mp_N, mp_D);
	  else
	    {
	      /* not in full column rank, compute minimal denominator solution
		 compute the bound after compression in specialminSolveRNS */
	      mpz_init(mp_bd);
	      compressBoundMP(r, m, Pt, mp_A, mp_bd);
	      ver = specialminSolveRNS(certflag, 1, nullcol, r, m,\
				       basislen, mp_bd, basis, CRNS, mp_b1, \
				       mp_N, mp_D, mp_NZ, mp_DZ);
	      mpz_clear(mp_bd);
	      if (ver == 0)
		{
		  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); }
		  for (i = 0; i < basislen; i++) { XFREE(CRNS[i]); }
		  { XFREE(CRNS); XFREE(mp_b1); }
		  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		  printf("certificate checking fails, repeat!\n");
#endif

		  continue;
		}
	    }
	  for (i = 0; i < basislen; i++) { XFREE(CRNS[i]); } { XFREE(CRNS); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
	  printf("Minimial Denominator Solution Solving Time: %f\n", tt);
	  ti = clock();
#endif

	  /* not in full row rank, need to check whether <A_31|A_32>x = b_3 */
	  if (r < n)
	    {
	      C1RNS = XMALLOC(Double *, basislen);
	      for (i = 0; i < basislen; i++)
		{
		  C1RNS[i] = XMALLOC(Double, (n-r)*m);
		  for (j = 0; j < n-r; j++)
		    for (l = 0; l < m; l++)
		      C1RNS[i][j*m+l] = \
			(Double)mpz_fdiv_ui(mp_A[Pt[r+j]*m+rpt[l]], basis[i]);
		}

	      cmp = LVecSMatMulCmp(RightMul, basislen, n-r, m, basis, C1RNS, \
				   mp_D, mp_N, mp_b1+r);
	      for (i = 0; i < basislen; i++) { XFREE(C1RNS[i]); }
	      XFREE(C1RNS);

	      /* if compare fails, restart */
	      if (cmp == 0)
		{
		  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); }
		  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); XFREE(mp_b1); }
#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		  printf("checking lower n-r rows fails, repeat!\n");
#endif

		  continue;
		}
	    }
	  { XFREE(Pt); XFREE(rpt); }
	  for (i = 0; i < n; i++) { mpz_clear(mp_b1[i]); } { XFREE(mp_b1); }
	  { XFREE(basiscmb[0]); XFREE(basiscmb[1]); XFREE(basiscmb); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
	  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
	  printf("Solution Checking Time: %f\n", tt);
#endif

	  /* recover the solution vector and the certificate vector */
	  if (r < m)
	    {
	      /* system has more than one solution, the first case */
	      /* compute Qy */
	      for (i = r; i > 0; i--)
		if (rp[i] != i) { mpz_swap(mp_N[i-1], mp_N[rp[i]-1]); }

	      /* set last n-r entries in z be 0 and compute zP */
	      if (certflag == 1)
		{
		  for (i = r; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
		  for (i = r; i > 0; i--)
		    if (P[i] != i) { mpz_swap(mp_NZ[i-1], mp_NZ[P[i]-1]); }
		}
	      { XFREE(P); XFREE(rp); }
	      return 1;
	    }
	  else
	    {
	      /* system has a unique solution, the second case */
	      { XFREE(P); XFREE(rp); }
	      return 2;
	    }
	}
      else
	{
	  if ((r == m) && (certflag == 0))
	    {
	      { XFREE(basiscmb[0]); XFREE(basiscmb[1]); XFREE(basiscmb); }
	      { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }
	      return 3;
	    }
	  P[r+1] = idx+1;      /* update P for permutation of row r+1 */

	  /* compute u = A_21.A_11^(-1) */
	  mpz_init(mp_Du);
	  mp_Nu = XMALLOC(mpz_t, r);
	  for (i = 0; i < r; i++) { mpz_init(mp_Nu[i]); }
	  A11RNS = XMALLOC(Double *, basislen);
	  for (i = 0; i < basislen; i++)
	    {
	      A11RNS[i] = XMALLOC(Double, r*r);
	      for (j = 0; j < r; j++)
		for (l = 0; l < r; l++)
		  A11RNS[i][j*r+l] = \
		    (Double)mpz_fdiv_ui(mp_A[Pt[j]*m+rpt[l]], basis[i]);
	    }
	  mp_A21 = XMALLOC(mpz_t, r);
	  for (i = 0; i < r; i++)
	    mpz_init_set(mp_A21[i], mp_A[Pt[idx]*m+rpt[i]]);
	  nonsingSolvRNSMM(LeftSolu, r, 1, basislen, basis, A11RNS, mp_A21, \
			   mp_Nu, mp_Du);
	  for (i = 0; i < basislen; i++)  { XFREE(A11RNS[i]); }
	  for (i = 0; i < r; i++) { mpz_clear(mp_A21[i]); }
	  { XFREE(mp_A21); XFREE(A11RNS); }
	  if (r < m)
	    {
	      /* compare u.A_12 with A_22 */
	      A12RNS = XMALLOC(Double *, basislen);
	      for (i = 0; i < basislen; i++)
		{
		  A12RNS[i] = XMALLOC(Double, r*(m-r));
		  for (j = 0; j < r; j++)
		    for (l = 0; l < m-r; l++)
		      A12RNS[i][j*(m-r)+l] = \
			(Double)mpz_fdiv_ui(mp_A[Pt[j]*m+rpt[r+l]], basis[i]);
		}
	      mp_A22 = XMALLOC(mpz_t, m-r);
	      for (i = 0; i < m-r; i++)
		mpz_init_set(mp_A22[i], mp_A[Pt[idx]*m+rpt[r+i]]);
	      cmp = LVecSMatMulCmp(LeftMul, basislen, r, m-r, basis, A12RNS, \
				   mp_Du, mp_Nu, mp_A22);
	      for (i = 0; i < basislen; i++) { XFREE(A12RNS[i]); }
	      for (i = 0; i < m-r; i++) { mpz_clear(mp_A22[i]); }
	      { XFREE(A12RNS); XFREE(mp_A22); }

	      /* comparison fails, restart */
	      if (cmp == 0)
		{
		  mpz_clear(mp_Du);
		  for (i = 0; i < r; i++) { mpz_clear(mp_Nu[i]); }
		  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); XFREE(mp_Nu); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
		  printf("checking no solution case fails, repeat!\n");
#endif

		  continue;
		}
	    }
	  if (certflag == 1)
	    {
	      /* set q = (u, -1, 0, ..., 0), which is store in mp_NZ */
	      mpz_set(mp_DZ, mp_Du);
	      for (i = 0; i < r; i++) { mpz_set(mp_NZ[i], mp_Nu[i]); }
	      mpz_set(mp_NZ[r], mp_Du);
	      mpz_neg(mp_NZ[r], mp_NZ[r]);
	      for (i = r+1; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }

	      /* compute qP */
	      for (i = r+1; i > 0; i--)
		if (P[i] != i) { mpz_swap(mp_NZ[i-1], mp_NZ[P[i]-1]); }
	    }

	  mpz_clear(mp_Du);
	  { XFREE(basiscmb[0]); XFREE(basiscmb[1]); XFREE(basiscmb); }
	  for (i = 0; i < r; i++) { mpz_clear(mp_Nu[i]); } { XFREE(mp_Nu); }
	  { XFREE(P); XFREE(rp); XFREE(Pt); XFREE(rpt); }

	  /* system has no solution, the third case */
	  return 3;
	}
    }
}





/*
 *
 * Calling Sequence:
 *   At <-- revseq(r, m, A)
 *
 * Summary:
 *   Perform operations on the permutation vector
 *
 * Description:
 *   Let A be a vector length m+1 to record the permutations over a matrix M.
 *   A[i] represents the switch of row/column i-1 of M with row/column A[i]-1
 *   of M, i = 1..r. The permutation order is A[r].A[r-1]. ... .A[1].M.
 *
 *   The function at first generates a vector At length m, At[i] = i,
 *   i = 0..m-1. Then apply A[r]. ... .A[1] in order on At. The outputed At
 *   means that row/column At[i] of matrix M will be changed to row/column i
 *   of M after applying all the permutations to M in order.
 *
 * Input:
 *   r: long, permutation happens in the first r row/column of M
 *   m: long, length of A
 *   A: 1-dim long array length m+1, record of permutations on matrix M
 *
 * Return value:
 *   At: 1-dim long array length m, explained above
 *
 */

long *
revseq(const long r, const long m, const long *A)
{
  long i, t, *At;

  At = XMALLOC(long, m);
  for (i = 0; i < m; i++) { At[i] = i; }
  for (i = 1; i < r+1; i++)
    if (A[i] != i)
      { t = At[i-1];  At[i-1] = At[A[i]-1];  At[A[i]-1] = t; }

  return At;
}




/*
 * Calling Sequence:
 *   compressBoundLong(n, m, Pt, A, mp_bd)
 *
 * Summary:
 *   Compute the maximum magnitude of a compressed signed long submatrix
 *
 * Description:
 *   Let A be a n1 x m matrix, B be a n x m submatrix of A. Pt[i] represents
 *   that row i of B comes from row Pt[i] of A. The function computes the
 *   maximum magnitude of entries in the compressed matrix B.P, where P is an
 *   m x n+k 0-1 matrix.
 *
 * Input:
 *      n: long, row dimension of B
 *      m: long, column dimension of B
 *     Pt: 1-dim long array length n1+1, storing row relations between B and A
 *   mp_A: 1-dim long array length n1*m, n1 x m matrix A
 *
 * Output:
 *   mp_bd: mpz_t, maximum magnitude of entries in B.P
 *
 */

void
compressBoundLong (const long n, const long m, const long *Pt, const long *A,\
		   mpz_t mp_bd)
{
  long i, j, temp;
  mpz_t mp_temp;

  mpz_init(mp_temp);
  mpz_set_ui(mp_bd, 0);
  for (i = 0; i < n; i++)
    {
      mpz_set_ui(mp_temp, 0);
      for (j = 0; j < m; j++)
	{
	  temp = labs(A[Pt[i]*m+j]);
	  mpz_add_ui(mp_temp, mp_temp, temp);
	}
      if (mpz_cmp(mp_bd, mp_temp) < 0) { mpz_set(mp_bd, mp_temp); }
    }

  mpz_clear(mp_temp);
}



/*
 * Calling Sequence:
 *   compressBoundMP(n, m, Pt, mp_A, mp_bd)
 *
 * Summary:
 *   Compute the maximum magnitude of a compressed mpz_t submatrix
 *
 * Description:
 *   Let A be a n1 x m matrix, B be a n x m submatrix of A. Pt[i] represents
 *   that row i of B comes from row Pt[i] of A. The function computes the
 *   maximum magnitude of entries in the compressed matrix B.P, where P is an
 *   m x n+k 0-1 matrix.
 *
 * Input:
 *      n: long, row dimension of B
 *      m: long, column dimension of B
 *     Pt: 1-dim long array length n1+1, storing row relations between B and A
 *   mp_A: 1-dim mpz_t array length n1*m, n1 x m matrix A
 *
 * Output:
 *   mp_bd: mpz_t, maximum magnitude of entries in B.P
 *
 */

void
compressBoundMP (const long n, const long m, const long *Pt, mpz_t *mp_A, \
		 mpz_t mp_bd)
{
  long i, j;
  mpz_t mp_temp, mp_temp1;

  { mpz_init(mp_temp); mpz_init(mp_temp1); }
  mpz_set_ui(mp_bd, 0);
  for (i = 0; i < n; i++)
    {
      mpz_set_ui(mp_temp, 0);
      for (j = 0; j < m; j++)
	{
	  mpz_abs(mp_temp1, mp_A[Pt[i]*m+j]);
	  mpz_add(mp_temp, mp_temp, mp_temp1);
	}
      if (mpz_cmp(mp_bd, mp_temp) < 0) { mpz_set(mp_bd, mp_temp); }
    }

  { mpz_clear(mp_temp); mpz_clear(mp_temp1); }
}






/*
 * Calling Sequence:
 *   1/0 <-- LVecSMatMulCmp(mulpos, basislen, n, m, basis, ARNS, mp_s,
 *                          mp_V, mp_b)
 *
 * Summary:
 *   Verify the equality of a matrix-vector product and a scalar-vector
 *   product
 *
 * Description:
 *   The function verifies whether a matrix-vector product A.V or V.A is equal
 *   to a scalar-vector product s.b, where A is a n x m matrix, b and V is a
 *   vector, and s is a scalar. The flag mulpos is used to specify whether the
 *   matrix-vector product is A.V or V.A. The comparison result is output
 *   sensitive, i.e., the comparsion will stop as long as one failure occurs.
 *   The input matrix A is represented in a RNS.
 *
 * Input:
 *     mulpos: LeftMul/RightMul, flag to indicate whether A.V or V.A is
 *             performed.
 *             If mulpos = LeftMul ==> V.A. If mulpos = RightMul ==> A.V
 *   basislen: long, dimension of the RNS basis
 *          n: long, row dimension of A
 *          m: long, column dimension of A
 *      basis: 1-dim FiniteField array length basislen, RNS basis
 *       ARNS: 2-dim Double array, dimension basislen x n*m, representation of
 *             A in the RNS. ARNS[i] = A mod basis[i], i = 0..basislen-1
 *       mp_s: mpz_t, scalar s
 *       mp_V: 1-dim mpz_t array length n or m depending on mulpos, vector V
 *       mp_b: 1-dim mpz_t array length n or m depending on mulpos, vector b
 *
 * Return:
 *   1: comparison succeeds
 *   0: comparison fails
 *
 */

long
LVecSMatMulCmp (const enum MULT_POS mulpos, const long basislen, \
		const long n, const long m, const FiniteField *basis, \
		Double **ARNS, mpz_t mp_s, mpz_t *mp_V, mpz_t *mp_b)
{
  long i, j, l, sharednum, unsharednum;
  FiniteField *bdcoeff, *cmbasis;
  double temp;
  mpz_t mp_AV, mp_temp, mp_AV1;
  Double **U;

  if (mulpos == LeftMul) { sharednum = n;  unsharednum = m; }
  else if (mulpos == RightMul) { sharednum = m;  unsharednum = n; }

  /* allocate matrix U to store mix radix coefficients of one row/column
     of matrix A */
  U = XMALLOC(Double *, basislen);
  for (i = 0; i < basislen; i++)
    U[i] = XMALLOC(Double, sharednum);
  cmbasis = combBasis(basislen, basis);
  bdcoeff = repBound(basislen, basis, cmbasis);
  mpz_init(mp_AV);
  mpz_init(mp_AV1);
  mpz_init(mp_temp);
  for (i = 0; i < unsharednum; i++)
    {
      /* apply Garner's algorithm to compute U in positive representation */
      if (mulpos == LeftMul)
	cblas_dcopy(sharednum, ARNS[0]+i, unsharednum, U[0], 1);
      else if (mulpos == RightMul)
	cblas_dcopy(sharednum, ARNS[0]+i*sharednum, 1, U[0], 1);
      for (j = 1; j < basislen; j++)
	{
	  cblas_dcopy(sharednum, U[j-1], 1, U[j], 1);
	  for (l = j-2; l >= 0; l--)
	    {
	      cblas_dscal(sharednum, (double)(basis[l] % basis[j]), U[j], 1);
	      cblas_daxpy(sharednum, 1.0, U[l], 1, U[j], 1);
	      Dmod((double)basis[j], U[j], 1, sharednum, sharednum);
	    }
	  temp = (double)cmbasis[j]*(double)(basis[j]-1);
	  temp = fmod(temp, (double)basis[j]);
	  if (mulpos == LeftMul) {
	    cblas_dscal(sharednum, temp, U[j], 1);
	    cblas_daxpy(sharednum, (double)cmbasis[j], ARNS[j]+i, unsharednum, U[j], 1);
          } else if (mulpos == RightMul) {
	    cblas_dscal(sharednum, temp, U[j], 1);
	    cblas_daxpy(sharednum, (double)cmbasis[j], ARNS[j]+i*sharednum,  1, U[j], 1);
          }
	  Dmod((double)basis[j], U[j], 1, sharednum, sharednum);
	}
      /* change U to symmetric representation */
      for (j = 0; j < sharednum; j++)
	for (l = basislen-1; l >= 0; l--)
	  {
	    if (U[l][j] > bdcoeff[l])
	      {
		U[basislen-1][j] = U[basislen-1][j]-basis[basislen-1];
		break;
	      }
	    else if (U[l][j] < bdcoeff[l]) { break; }
	  }
      /* do vector-vector product and sum up the product */
      mpz_set_ui(mp_AV, 0);
      for (j = 0; j < sharednum; j++)
	{
	  mpz_mul_si(mp_temp, mp_V[j], U[basislen-1][j]);
	  mpz_add(mp_AV, mp_AV, mp_temp);
	}
      for (j = basislen-2; j >= 0; j--)
	{

	  mpz_mul_ui(mp_AV, mp_AV, basis[j]);
	  mpz_set_ui(mp_AV1, 0);
	  for (l = 0; l < sharednum; l++)
	    {
	      mpz_mul_ui(mp_temp, mp_V[l], U[j][l]);
	      mpz_add(mp_AV1, mp_AV1, mp_temp);
	    }
	  mpz_add(mp_AV, mp_AV, mp_AV1);
	}
      /* compare with mp_s*mp_b */
      mpz_mul(mp_temp, mp_s, mp_b[i]);
      if (mpz_cmp(mp_temp, mp_AV) != 0)
	{
	  { mpz_clear(mp_AV); mpz_clear(mp_AV1); mpz_clear(mp_temp); }
	  for (j = 0; j < basislen; j++) { XFREE(U[j]); } { XFREE(U); }
	  { XFREE(bdcoeff); XFREE(cmbasis); }
	  return 0;
	}
    }

  { mpz_clear(mp_AV); mpz_clear(mp_AV1); mpz_clear(mp_temp); }
  for (j = 0; j < basislen; j++) { XFREE(U[j]); } { XFREE(U); }
  { XFREE(bdcoeff); XFREE(cmbasis); }
  return 1;
}



/*
 *
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


/* return max( upper bound of signed long, 2^52-1 ) */
long
LongRNSbound(void)
{
  long a;
  mpz_t mp_a, mp_b;

  mpz_init(mp_a);
  mpz_ui_pow_ui(mp_a, 2, sizeof(long)*8-1);
  a = mpz_get_ui(mp_a);
  a = a-1;                    /* a := maximal number in Signed Long */
  mpz_clear(mp_a);
  mpz_init(mp_b);
  mpz_ui_pow_ui(mp_b, 2, 52);
  mpz_sub_ui(mp_b, mp_b, 1);  /* mp_b := 2^52-1 */
  if (mpz_cmp_ui(mp_b, a) >= 0) { mpz_clear(mp_b); return a; }
  else
    {
      unsigned long tmp = mpz_get_ui(mp_b);
      mpz_clear(mp_b);
      return tmp;
    }
}




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
#endif

  { XFREE(P1); XFREE(P2); }
  for (i = 0; i < n+k; i++) { mpz_clear(mp_Nt[i]); } { XFREE(mp_Nt); }
  return 1;
}



/*
 * Calling Sequences:
 *   1/0 <-- certVerify(basislen, n, m, basis, ARNS, mp_DZ, mp_NZ)
 *
 * Summary:
 *   Verify the correctness of a certificate vector for a system of linear
 *   equations
 *
 * Description:
 *   Let the system of linear equations be Ax = b. The function verifies the
 *   rational vector z by checking whether z.A is an integer vector. The
 *   n x m matrix A is represented in RNS. z is represented by an integer
 *   numerator vector and an integer denominator.  The function is output
 *   sensitive, i.e., the verification will stop as long as one failure
 *   occurs.
 *
 * Input:
 *   basislen: long, dimension of the RNS basis
 *          n: long, row dimension of A
 *          m: long, column dimension of A
 *      basis: 1-dim FiniteField array length basislen, RNS basis
 *       ARNS: 2-dim Double array, dimension basislen x n*m, representation of
 *             A in the RNS. ARNS[i] = A mod basis[i], i = 0..basislen-1
 *      mp_DZ: mpz_t, denominator of vector z
 *      mp_NZ: 1-dim mpz_t array length n, numerator array of vector z
 *
 * Return:
 *   1: z is the certificate of A
 *   0: z is not the certificate of A
 *
 */

long
certVerify (const long basislen, const long n, const long m, \
	    const FiniteField *basis, Double **ARNS, mpz_t mp_DZ, \
	    mpz_t *mp_NZ)
{
  long i, j, l;
  FiniteField *bdcoeff, *cmbasis;
  double temp;
  mpz_t mp_AZ, mp_temp, mp_AZ1;
  Double **U;

  /* allocate matrix U to store mix radix coefficients of one row/column
     of matrix A */
  U = XMALLOC(Double *, basislen);
  for (i = 0; i < basislen; i++)
    U[i] = XMALLOC(Double, n);
  cmbasis = combBasis(basislen, basis);
  bdcoeff = repBound(basislen, basis, cmbasis);
  mpz_init(mp_AZ);
  mpz_init(mp_AZ1);
  mpz_init(mp_temp);
  for (i = 0; i < m; i++)
    {
      /* apply Garner's algorithm to compute U in positive representation */
      cblas_dcopy(n, ARNS[0]+i, m, U[0], 1);
      for (j = 1; j < basislen; j++)
	{
	  cblas_dcopy(n, U[j-1], 1, U[j], 1);
	  for (l = j-2; l >= 0; l--)
	    {
	      cblas_dscal(n,  (double)(basis[l] % basis[j]), U[j], 1);
	      cblas_daxpy(n, 1.0, U[l], 1, U[j], 1);
	      Dmod((double)basis[j], U[j], 1, n, n);
	    }
	  temp = (double)cmbasis[j]*(double)(basis[j]-1);
	  temp = fmod(temp, (double)basis[j]);
	  cblas_dscal(n, temp, U[j], 1);
	  cblas_daxpy(n, (double)cmbasis[j], ARNS[j]+i, m, U[j], 1);
	  Dmod((double)basis[j], U[j], 1, n, n);
	}
      /* change U to symmetric representation */
      for (j = 0; j < n; j++)
	for (l = basislen-1; l >= 0; l--)
	  {
	    if (U[l][j] > bdcoeff[l])
	      {
		U[basislen-1][j] = U[basislen-1][j]-basis[basislen-1];
		break;
	      }
	    else if (U[l][j] < bdcoeff[l]) { break; }
	  }
      /* do vector-vector product and sum up the product */
      mpz_set_ui(mp_AZ, 0);
      for (j = 0; j < n; j++)
	{
	  mpz_mul_si(mp_temp, mp_NZ[j], U[basislen-1][j]);
	  mpz_add(mp_AZ, mp_AZ, mp_temp);
	}
      for (j = basislen-2; j >= 0; j--)
	{

	  mpz_mul_ui(mp_AZ, mp_AZ, basis[j]);
	  mpz_set_ui(mp_AZ1, 0);
	  for (l = 0; l < n; l++)
	    {
	      mpz_mul_ui(mp_temp, mp_NZ[l], U[j][l]);
	      mpz_add(mp_AZ1, mp_AZ1, mp_temp);
	    }
	  mpz_add(mp_AZ, mp_AZ, mp_AZ1);
	}
      /* check whether mp_DZ|mp_AZ */
      if (!mpz_divisible_p(mp_AZ, mp_DZ))
	{
	  { mpz_clear(mp_AZ); mpz_clear(mp_AZ1); mpz_clear(mp_temp); }
	  for (j = 0; j < basislen; j++) { XFREE(U[j]); } { XFREE(U); }
	  { XFREE(bdcoeff); XFREE(cmbasis); }
	  return 0;
	}
    }
  { mpz_clear(mp_AZ); mpz_clear(mp_AZ1); mpz_clear(mp_temp); }
  for (j = 0; j < basislen; j++) { XFREE(U[j]); } { XFREE(U); }
  { XFREE(bdcoeff); XFREE(cmbasis); }
  return 1;
}



/*
 * Calling Seqeuence:
 *   minSolnoncompressLong(certflag, redflag, n, k, mp_Bb, A, mp_N, mp_D,
 *                         mp_NZ, mp_DZ)
 *
 * Summary:
 *   Compute the minimal denominator solution of a full row rank system
 *   without compression, represented by signed long integers, and the
 *   corresponding certificate vector(optionalf). The solution size could be
 *   reduced.
 *
 * Description:
 *   Let the full row rank system be Cx = b, where b is a vector, and C is like
 *          n    k
 *       [     |   ]
 *   C = [  A  | B ] n
 *       [     |   ]
 *
 *   The function deals with the case when C is a signed long matrix. It
 *   takes A and Bb, represented by signed long and mpz_t respectively as
 *   input, where Bb is the matrix composed of B and b, Bb[1..-1, 1..k] = B
 *   and Bb[1..-1, k+1] = b.
 *
 *   The certificate vector z satisfy z.C is an integer vector and z.b has the
 *   same denominator as the solution vector.
 *
 *   If redflag is specified as 1, then reduce the solution by the kernel
 *   basis with dimension k.
 *
 * Input:
 *   certflag: 1/0, flag to indicate whether to compute the certificate vector
 *             or not. If certflag = 1, then compute the certificate vector.
 *             Otherwise, not compute the certificate.
 *    redflag: 1/0, flag to indicate whether to reduce the solution size or not
 *           - redflag = 1, reduce the solution size
 *           - redflag = 0, not reduce the solution size
 *          n: long, row dimension of the input system
 *          k: long, column dimension of B
 *      mp_Bb: 1-dim mpz_t array length (k+1)*n, matrix Bb consisting of
 *             submatrix B and vector b
 *          A: 1-dim long array length n*n, submatrix A of C
 *
 * Output:
 *    mp_N: 1-dim mpz_t array length n+k, numerator vector of the solution
 *          with minimal denominator
 *    mp_D: mpz_t, denominator of the solution
 *   mp_NZ: 1-dim mpz_t array, the first n entries store the numerator vector
 *          of the certificate z
 *   mp_DZ: mpz_t, the denominator of the certificate z
 *
 * Note:
 *   - The space of the solution (mp_N, mp_D) is needed to preallocated.
 *   - If certflag is specified to be 1, then also needs to preallocate space
 *     for (mp_NZ, mp_DZ). Otherwise, set mp_NZ = NULL, and mp_DZ = any int.
 *
 */

void
minSolnoncompressLong (const long certflag, const long redflag, const long n, \
		       const long k, mpz_t *mp_Bb, const long* A, \
		       mpz_t *mp_N, mpz_t mp_D, mpz_t *mp_NZ, mpz_t mp_DZ)
{
  long i, j, s, t, h;
  mpz_t mp_DABb, mp_s, mp_ns, mp_temp, mp_h;
  mpz_t *mp_NABb, *mp_M, *mp_Cert=NULL, *mp_T, *mp_Lat;
  double tt;

#if HAVE_TIME_H
  clock_t ti, ti1;
#endif

  mp_NABb = XMALLOC(mpz_t, (k+1)*n);
  for (i = 0; i < (k+1)*n; i++) { mpz_init(mp_NABb[i]); }
  mpz_init(mp_DABb);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* compute A^(-1)B and A^(-1)b at the same time */
  nonsingSolvMM(RightSolu, n, k+1, A, mp_Bb, mp_NABb, mp_DABb);

  /* compute s */
  mpz_init_set(mp_s, mp_DABb);
  mpz_init(mp_ns);
  mpz_neg(mp_ns, mp_s);

  /* construct matrix M = <s*<A^(-1)B, I_k>|s*<A^(-1)b, 0>, <0|s>> */
  mp_M = XMALLOC(mpz_t, (n+k+1)*(k+1));
  for (i = 0; i < (n+k+1)*(k+1); i++) { mpz_init(mp_M[i]); }
  mpz_init_set_ui(mp_temp, 1);
  scalCpMP(n, k+1, k+1, k+1, mp_temp, mp_NABb, mp_M);
  for (i = 0; i < k; i++) { mpz_set(mp_M[n*(k+1)+i*(k+1)+i], mp_ns); }
  mpz_set(mp_M[(n+k)*(k+1)+k], mp_s);
  { mpz_clear(mp_DABb); mpz_clear(mp_ns); }
  for (i = 0; i < (k+1)*n; i++) { mpz_clear(mp_NABb[i]); } { XFREE(mp_NABb); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("          Matix M Construction Time(A^(-1)B, A^(-1)b): %f\n", tt);
#endif

  /* construct matrix mp_T = s*<A^(-1)B, -I_k> for kernel basis N */
  mp_T = XMALLOC(mpz_t, (n+k)*k);
  for (i = 0; i < (n+k)*k; i++) { mpz_init(mp_T[i]); }
  scalCpMP(n+k, k, k+1, k, mp_temp, mp_M, mp_T);

  /* store s*<A^(-1)b, 0> into solution array mp_N */
  scalCpMP(n, 1, k+1, 1, mp_temp, mp_M+k, mp_N);
  for (i = n; i < n+k; i++) { mpz_set_ui(mp_N[i], 0); }
  mpz_clear(mp_temp);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* compute special hermite form of M and certificate(if certflag == true) */
  if (certflag == 1)
    {
      mp_Cert = XMALLOC(mpz_t, n+k+1);
      for (i = 0; i < n+k+1; i++) { mpz_init(mp_Cert[i]); }
    }
  specialHermite(certflag, n, k, 1, mp_M, mp_Cert);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("          Hermite Form Computing Time: %f\n", tt);
#endif

  /* compute minimal denominator */
  mpz_init_set(mp_h, mp_M[k*(k+1)+k]);
  mpz_divexact(mp_D, mp_s, mp_h);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* compute kernel basis N = mp_TH^(-1), which is inplaced into mp_T */
  kernelBasis(n, k, 1, mp_M, mp_T);

  /* initialize the lattice */
  if (redflag == 1)
    {
      mp_Lat = XMALLOC(mpz_t, (n+k)*(k+1));
      for (i = 0; i < k*(n+k); i++)
	{
	  s = (long)(i/(n+k));  t = i-s*(n+k);  h = t*k+s;
	  mpz_init_set(mp_Lat[i], mp_T[h]);
	}
    }

  /* compute mp_T.M[1..k, k+1] and inplace the result into first column
     of mp_T */
  for (i = 0; i < n+k; i++)
    {
      mpz_mul(mp_T[i*k], mp_T[i*k], mp_M[k]);
      for (j = 1; j < k; j++)
	mpz_addmul(mp_T[i*k], mp_T[i*k+j], mp_M[j*(k+1)+k]);

      /* compute numerator vector of the solution at the same time */
      mpz_sub(mp_N[i], mp_N[i], mp_T[i*k]);
      mpz_divexact(mp_N[i], mp_N[i], mp_h);
    }
  { mpz_clear(mp_s); mpz_clear(mp_h); }
  for (i = 0; i < (n+k)*k; i++) { mpz_clear(mp_T[i]); } { XFREE(mp_T); }
  for (i = 0; i < (n+k+1)*(k+1); i++) { mpz_clear(mp_M[i]); } { XFREE(mp_M); }

  /* reduct solution by lattice reduction */
  if (redflag == 1)
    {
      for (i = 0; i < n+k; i++) { mpz_init_set(mp_Lat[k*(n+k)+i], mp_N[i]); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      ti1 = clock();
#endif

      ired(mp_Lat, k+1, n+k, k);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti1)/CLOCKS_PER_SEC;
      printf("           lattice reduce solution time: %f\n", tt);
#endif

      for (i = 0; i < n+k; i++) { mpz_set(mp_N[i], mp_Lat[k*(n+k)+i]); }
      for (i = 0; i < (n+k)*(k+1); i++) { mpz_clear(mp_Lat[i]); }
      XFREE(mp_Lat);
    }

  /* compute certificate */
  if (certflag == 1)
    {
      if (mpz_cmp_ui(mp_D, 1) != 0)
	nonsingSolvMM(LeftSolu, n, 1, A, mp_Cert, mp_NZ, mp_DZ);
      else
	{
	  for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	  mpz_set_ui(mp_DZ, 1);
	}
      for (i = 0; i< n+k+1; i++) { mpz_clear(mp_Cert[i]); } { XFREE(mp_Cert); }
    }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("          Kernel Basis, Certificate, Solution Finding Time: %f\n", \
	 tt);
#endif

  return;
}



/*
 * Calling Sequence:
 *   minSolnoncompressRNS(certflag, redflag, n, k, basislen, basis, mp_Bb,
 *                        ARNS, mp_N, mp_D, mp_NZ, mp_DZ)
 *
 * Summary:
 *   Compute the minimal denominator solution of a full row rank system
 *   without compression, represented in RNS, and the corresponding
 *   certificate vector(optional). The solution size could be reduced.
 *
 * Description:
 *   Let the full row rank system be Cx = b, where b is a vector, and C is like
 *          n    k
 *       [     |   ]
 *   C = [  A  | B ] n
 *       [     |   ]
 *
 *   The function doesn't take C as input, but takes ARNS, the representation
 *   of A in a RNS, and mpz_t matrix Bb, the matrix composed of B and b as
 *   input, where Bb[1..-1, 1..k] = B and Bb[1..-1, k+1] = b.
 *
 *   The certificate vector z satisfy z.C is an integer vector and z.b has the
 *   same denominator as the solution vector.
 *
 *   If redflag is specified as 1, then reduce the solution by the kernel
 *   basis with dimension k.
 *
 *
 * Input:
 *   certflag: 1/0, flag to indicate whether to compute the certificate vector
 *             or not. If certflag = 1, then compute the certificate vector.
 *             Otherwise, not compute the certificate.
 *    redflag: 1/0, flag to indicate whether to reduce the solution size or not
 *           - redflag = 1, reduce the solution size
 *           - redflag = 0, not reduce the solution size
 *          n: long, row dimension of the input system
 *          k: long, column dimension of B
 *   basislen: long, dimension of the RNS basis
 *      basis: 1-dim FiniteField array length basislen, RNS basis
 *      mp_Bb: 1-dim mpz_t array length (k+1)*n, matrix Bb consisting of
 *             submatrix B and vector b
 *       ARNS: 2-dim Double array, dimension basislen x n*n, representation of
 *             A in the RNS. ARNS[i] = A mod basis[i], i = 0..basislen-1
 *
 * Output:
 *    mp_N: 1-dim mpz_t array length n+k, numerator vector of the solution
 *          with minimal denominator
 *    mp_D: mpz_t, denominator of the solution
 *   mp_NZ: 1-dim mpz_t array, the first n entries store the numerator
 *          vector of the certificate z
 *   mp_DZ: mpz_t, the denominator of the certificate z
 *
 * Note:
 *   - The space of the solution (mp_N, mp_D) is needed to preallocated.
 *   - If certflag is specified to be 1, then also needs to preallocate space
 *       for (mp_NZ, mp_DZ). Otherwise, set mp_NZ = NULL, and mp_DZ = any int.
 *
 * Precondition:
 *   any element p in RNS basis must satisfy 2*(p-1)^2 <= 2^53-1.
 *
 */

void
minSolnoncompressRNS (const long certflag, const long redflag, const long n, \
		      const long k, const long basislen, \
		      const FiniteField *basis, mpz_t *mp_Bb, Double **ARNS, \
		      mpz_t *mp_N, mpz_t mp_D, mpz_t *mp_NZ, mpz_t mp_DZ)
{
  long i, j, s, t, h;
  mpz_t mp_DABb, mp_s, mp_ns, mp_temp, mp_h;
  mpz_t *mp_NABb, *mp_M, *mp_Cert=NULL, *mp_T, *mp_Lat;
  double tt;

#if HAVE_TIME_H
  clock_t ti, ti1;
#endif

  mp_NABb = XMALLOC(mpz_t, (k+1)*n);
  for (i = 0; i < (k+1)*n; i++) { mpz_init(mp_NABb[i]); }
  mpz_init(mp_DABb);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* compute A^(-1)B and A^(-1)b at the same time */
  nonsingSolvRNSMM(RightSolu, n, k+1, basislen, basis, ARNS, mp_Bb, \
		   mp_NABb, mp_DABb);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("A^(-1)B, A^(-1)b solving time: %f\n", tt);
  fflush(stdout);
  ti = clock();
#endif

  /* compute s */
  mpz_init_set(mp_s, mp_DABb);
  mpz_init(mp_ns);
  mpz_neg(mp_ns, mp_s);

  /* construct matrix M = <s*<A^(-1)B, I_k>|s*<A^(-1)b, 0>, <0|s>> */
  mp_M = XMALLOC(mpz_t, (n+k+1)*(k+1));
  for (i = 0; i < (n+k+1)*(k+1); i++) { mpz_init(mp_M[i]); }
  mpz_init_set_ui(mp_temp, 1);
  scalCpMP(n, k+1, k+1, k+1, mp_temp, mp_NABb, mp_M);
  for (i = 0; i < k; i++) { mpz_set(mp_M[n*(k+1)+i*(k+1)+i], mp_ns); }
  mpz_set(mp_M[(n+k)*(k+1)+k], mp_s);
  { mpz_clear(mp_DABb); mpz_clear(mp_ns); }
  for (i = 0; i < (k+1)*n; i++) { mpz_clear(mp_NABb[i]); } { XFREE(mp_NABb); }

  /* construct matrix mp_T = s*<A^(-1)B, -I_k> for kernel basis N */
  mp_T = XMALLOC(mpz_t, (n+k)*k);
  for (i = 0; i < (n+k)*k; i++) { mpz_init(mp_T[i]); }
  scalCpMP(n+k, k, k+1, k, mp_temp, mp_M, mp_T);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("          Matix M, N Construction Time: %f\n", tt);
#endif

  /* store s*<A^(-1)b, 0> into solution array mp_N */
  scalCpMP(n, 1, k+1, 1, mp_temp, mp_M+k, mp_N);
  for (i = n; i < n+k; i++) { mpz_set_ui(mp_N[i], 0); }
  mpz_clear(mp_temp);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* compute special hermite form of M and certificate(if certflag == true) */
  if (certflag == 1)
    {
      mp_Cert = XMALLOC(mpz_t, n+k+1);
      for (i = 0; i < n+k+1; i++) { mpz_init(mp_Cert[i]); }
    }
  specialHermite(certflag, n, k, 1, mp_M, mp_Cert);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("          Hermite Form Computing Time: %f\n", tt);
  fflush(stdout);
#endif

  /* compute minimal denominator */
  mpz_init_set(mp_h, mp_M[k*(k+1)+k]);
  mpz_divexact(mp_D, mp_s, mp_h);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  /* compute kernel basis N = mp_TH^(-1), which is inplaced into mp_T */
  kernelBasis(n, k, 1, mp_M, mp_T);

  /* initialize the lattice */
  if (redflag == 1)
    {
      mp_Lat = XMALLOC(mpz_t, (n+k)*(k+1));
      for (i = 0; i < k*(n+k); i++)
	{
	  s = (long)(i/(n+k));  t = i-s*(n+k);  h = t*k+s;
	  mpz_init_set(mp_Lat[i], mp_T[h]);
	}
    }

  /* compute mp_T.M[1..k, k+1] and inplace the result into first column
     of mp_T */
  for (i = 0; i < n+k; i++)
    {
      mpz_mul(mp_T[i*k], mp_T[i*k], mp_M[k]);
      for (j = 1; j < k; j++)
	mpz_addmul(mp_T[i*k], mp_T[i*k+j], mp_M[j*(k+1)+k]);

      /* compute numerator vector of the solution at the same time */
      mpz_sub(mp_N[i], mp_N[i], mp_T[i*k]);
      mpz_divexact(mp_N[i], mp_N[i], mp_h);
    }
  { mpz_clear(mp_s); mpz_clear(mp_h); }
  for (i = 0; i < (n+k)*k; i++) { mpz_clear(mp_T[i]); } { XFREE(mp_T); }
  for (i = 0; i < (n+k+1)*(k+1); i++) { mpz_clear(mp_M[i]); } { XFREE(mp_M); }

  /* reduct solution by lattice reduction */
  if (redflag == 1)
    {
      for (i = 0; i < n+k; i++) { mpz_init_set(mp_Lat[k*(n+k)+i], mp_N[i]); }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      ti1 = clock();
#endif

      ired(mp_Lat, k+1, n+k, k);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
      tt = (double)(clock()-ti1)/CLOCKS_PER_SEC;
      printf("           lattice reduce solution time: %f\n", tt);
#endif

      for (i = 0; i < n+k; i++) { mpz_set(mp_N[i], mp_Lat[k*(n+k)+i]); }
      for (i = 0; i < (n+k)*(k+1); i++) { mpz_clear(mp_Lat[i]); }
      XFREE(mp_Lat);
    }

  /* compute certificate */
  if (certflag == 1)
    {
      if (mpz_cmp_ui(mp_D, 1) != 0)
	nonsingSolvRNSMM(LeftSolu, n, 1, basislen, basis, ARNS, mp_Cert, \
			 mp_NZ, mp_DZ);
      else
	{
	  for (i = 0; i < n; i++) { mpz_set_ui(mp_NZ[i], 0); }
	  mpz_set_ui(mp_DZ, 1);
	}
      for (i = 0; i< n+k+1; i++) { mpz_clear(mp_Cert[i]); } { XFREE(mp_Cert); }
    }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("          Kernel Basis, Certificate, Solution Finding Time: %f\n", \
	 tt);
  fflush(stdout);
#endif

  return;
}



/*
 * Calling Sequence:
 *   specialHermite(certflag, n, k, t, M, Cert)
 *
 * Summary:
 *   Perform unimodular transformation inplace in matrix M and return an
 *   certificate matrix Cert(optional) at the same time
 *
 * Description:
 *   Matrix M consists of:
 *
 *               k        t
 *         [          |        ]
 *         [          |        ]
 *         [  A^(-1)B | A^(-1)b] n
 *   M = s*[          |        ]
 *         [          |        ]
 *         [-------------------]
 *         [                   ]
 *         [                   ] k+t
 *         [     I_(k+t)       ]
 *         [                   ]
 *
 *   where A a n x n matrix, B a n x k matrix, b a n x t matrix, and I_(k+t) a
 *   k+t x k+t identity matrix.
 *
 *   The inplaced M has the properties:
 *    - upper triangular part of M[1..k, 1..k] is in hermite form
 *
 *          [g_1  *         * ]
 *          [    g_2   . .  * ]
 *      H = [          .    * ]
 *          [            .  * ]
 *          [              g_k]
 *
 *    - M[k+1, k+i] = s/d_i, where d_i is the minimum denominator of the
 *      system of linear equations [A | B]x = b[1..-1, i], i = 1..t
 *    - Vectors M[1..k, i] is reduced by doing modular M[k+1, i], i = k+1..k+t
 *
 *   The t x n+k+t certificate matrix Cert has the property:
 *    - For each row of Cert, fraction number computed by
 *      Cert[i, 1..-1].(A^(-1).b[1..-1, i]) has denominator d_i, i = 1..t.
 *      And Cert[i, 1..-1].(A^(-1).B) is an integer vector
 *
 * Input:
 *   certflag: 1/0, flag to indicate whether to compute certificate Cert
 *             or not. If certflag = 1, then compute the certificate.
 *             Otherwise, not compute the certificate.
 *          n: long, row dimension of A
 *          k: long, column dimension of B
 *          t: long, column dimension of b
 *          M: 1-dim mpz_t array length (n+k+t)*(k+t), n+k+t x k+t matrix M
 *
 * Output:
 *   inplaced mp_M
 *   Cert: 1-dim mpz_t array length t*(n+k+t), t x n+k+t certificate matrix
 *         Cert as above
 *
 * Note:
 *   the space of mp_Cert needs to be preallocated if certflag = 1.
 *   Otherwise, set mp_Cert = NULL
 *
 */

void
specialHermite (const long certflag, const long n, const long k, const long t,\
		mpz_t *M, mpz_t *Cert)
{
  long i, j, l;
  unsigned *C, *c;
  mpz_t s, g, p, q, p1, q1, tmpa, tmpb;
  mpz_t *P, *P1, *Q, *Q1, *tmpz;

  C = XMALLOC(unsigned, k*(n-1));
  c = XMALLOC(unsigned, n+t);
  mpz_init_set(s, _M(n+1,1));
  { mpz_init(g); mpz_init(p); mpz_init(q); mpz_init(p1); mpz_init(q1); }
  { mpz_init(tmpa); mpz_init(tmpb); }
  P = mpz_init_array(k);
  P1 = mpz_init_array(k);
  Q = mpz_init_array(k);
  Q1 = mpz_init_array(k);
  tmpz = mpz_init_array(n+t);

  for (i = 1; i <= k; i++)
    {
      /* compute one row of matrix C */
      for (j = i; j <= n+i-1; j++)
	mpz_set(_tmpz(j-i+1), _M(j,i));
      migcdex(s, _tmpz(1), &_tmpz(2), n-1, &_C(i,1));
      for (j = 1; j <= n-1; j++)
	{
	  for (l = i; l <= k+t; l++)
	    {
	      mpz_addmul_ui(_M(i,l), _M(i+j,l), _C(i,j));
	      mpz_mod(_M(i,l), _M(i,l), s);
	    }
	}
      /* compute gcdex of M(i,i) and M(n+i,i), apply the unimodular
	 transformation to M and store the transformation to P, Q, P1, Q1 */
      mpz_neg(tmpa, s);
      mpz_gcdext(g, p, q, _M(i,i), tmpa);
      mpz_divexact(p1, s, g);	mpz_divexact(q1, _M(i,i), g);
      mpz_set(_P(i), p);   mpz_set(_Q(i), q);
      mpz_set(_P1(i), p1); mpz_set(_Q1(i), q1);
      mpz_set(_M(i,i), g); mpz_set_ui(_M(n+i,i), 0);
      for (j = i+1; j <= k+t; j++)
	{
	  mpz_mul(tmpa, p, _M(i,j));
	  mpz_addmul(tmpa, q, _M(n+i,j));
	  mpz_mul(_M(n+i,j), p1, _M(i,j));
	  mpz_mod(_M(n+i,j), _M(n+i,j), s);
	  mpz_mod(_M(i,j), tmpa, s);
	}
      /* zero out the entries of column i of M below M(i,i) and store the
	 quotients q into M */
      for (j = i+1; j <= n+i-1; j++)
	{
	  mpz_divexact(q, _M(j,i), _M(i,i));
	  for (l = i+1; l <= k+t; l++)
	    {
	      mpz_submul(_M(j,l), q, _M(i,l));
	      mpz_mod(_M(j,l), _M(j,l), s);
	    }
	  mpz_neg(_M(j,i), q);
	}
    }
  for (i = 1; i <= t; i++)
    {
      for (j = k+1; j <= n+k+i-1; j++)
	mpz_set(_tmpz(j-k), _M(j,k+i));
      migcdex(s, _tmpz(1), &_tmpz(2), n+i-2, &_c(1));
      for (j = k+2; j <= n+k+i-1; j++)
	mpz_addmul_ui(_M(k+1,k+i), _M(j,k+i), _c(j-(k+1)));
      mpz_gcd(_M(k+1,k+i), _M(k+1,k+i), _M(n+k+i,k+i));

      if (certflag == 1)
	{
	  /* when gcd(M(k+1..n+k+i-1, k+i), s) = s then do not need to compute
	     Cert(i, 1..-1).  carry on to compute next row of Cert */
	  if (!mpz_cmp(_M(k+1,k+i), s)) { continue; }

	  /* initializtion of Cert */
	  mpz_set_ui(_Cert(i,k+1), 1);
	  for (j = k+2; j <= n+k+i-1; j++)
	    mpz_set_ui(_Cert(i,j), _c(j-(k+1)));

	  /* apply the stored transforms to Cert in inversed sequence */
	  for (j = k; j >= 1; j--)
	    {
	      for (l = j+1; l <= n+j-1; l++)
		{
		  mpz_addmul(_Cert(i,j), _Cert(i,l), _M(l,j));
		  mpz_mod(_Cert(i,j), _Cert(i,j), s);
		}
	      mpz_mul(tmpa, _P(j), _Cert(i,j));
	      mpz_addmul(tmpa, _P1(j), _Cert(i,n+j));
	      mpz_mul(tmpb, _Q(j), _Cert(i,j));
	      mpz_addmul(tmpb, _Q1(j), _Cert(i,n+j));
	      mpz_mod(_Cert(i,n+j), tmpb, s);
	      mpz_mod(_Cert(i,j), tmpa, s);
	      for (l = j+1; l <= n+j-1; l++)
		{
		  mpz_addmul_ui(_Cert(i,l), _Cert(i,j), _C(j,l-j));
		  mpz_mod(_Cert(i,l), _Cert(i,l), s);
		}
	    }
	}
    }
  /* carry on to reduce the upper triangular entries in the first k columns
     of M using diagonal entries */
  for (i = 2; i <= k; i++)
    {
      for (j = 1; j <= i-1; j++)
	{
	  mpz_tdiv_qr(q, _M(j,i), _M(j,i), _M(i,i));
	  for (l = i+1; l <= k+t; l++)
	    {
	      mpz_submul(_M(j,l), q, _M(i,l));
	      mpz_mod(_M(j,l), _M(j,l), s);
	    }
	}
    }
  /* reduce the last t columns of M using M(k+1,j), j=k+1..k+t */
  for (i = 1; i <= t; i++)
    {
      for (j = 1; j < k+1; j++)
	mpz_mod(_M(j,k+i), _M(j,k+i), _M(k+1,k+i));
    }

  { XFREE(C); XFREE(c); }
  { mpz_clear(s); mpz_clear(g); mpz_clear(tmpa); mpz_clear(tmpb); }
  { mpz_clear(p); mpz_clear(q); mpz_clear(p1); mpz_clear(q1); }
  for (i = 0; i < k; i++)
    { mpz_clear(P[i]); mpz_clear(P1[i]); mpz_clear(Q[i]); mpz_clear(Q1[i]); }
  { XFREE(P); XFREE(P1); XFREE(Q); XFREE(Q1); }
  for (i = 0; i < n+t; i++) { mpz_clear(tmpz[i]); } { XFREE(tmpz); }
  return;
}


/*
 * Calling Sequence:
 *   mpz_init_array(n)
 *
 * Summary:
 *   Allocate a mpz_t array length n dynamically and set all the entries in
 *   the array to be zero.
 *
 */

mpz_t *
mpz_init_array (const long n)
{
  long i;
  mpz_t *t = XMALLOC(mpz_t, n);
  for (i = 0; i < n; i++) { mpz_init(t[i]); }
  return t;
}



/*
 * Calling Sequence:
 *   migcdex(N, a, b, n, c)
 *
 * Summary:
 *   Compute solutions to the modulo N extended GCD problem
 *
 * Description:
 *   Given two mpz_t integers a and N, mpz_t array b, this function computes
 *   non-negative integers c[0], ... , c[n] such that
 *   gcd(a + c[0]b[0] + ... + c[n]b[n], N) = gcd(a, b[0], ..., b[n], N).
 *
 * Input:
 *   N: mpz_t, explained above
 *   a: mpz_t, explained above
 *   b: 1-dim mpz_t array length n, explained above
 *   n: long, dimension of array b and c
 *
 * Output:
 *   c: 1-dim unsigned array length n, explained above
 *
 */

void
migcdex (const mpz_t N, const mpz_t a, mpz_t *b, const long n, unsigned *c)
{
  long i, j;
  mpz_t gAN, gAbN, A;

  mpz_init(gAN);
  mpz_init(gAbN);
  mpz_init_set(A, a);
  mpz_gcd(gAbN, a, N);
  for (i = 0; i < n; i++)
    {
      mpz_gcd(gAbN, gAbN, b[i]);
      for (j = 0; ; j++)
	{
	  mpz_gcd(gAN, A, N);
	  if (!mpz_cmp(gAbN, gAN)) { break; }
	  mpz_add(A, A, b[i]);
	}
      c[i] = j;
    }

  { mpz_clear(gAN); mpz_clear(gAbN); mpz_clear(A); }
  return;
}



/*
 * Calling Sequence:
 *   kernelBasis(n, k, t, mp_M, mp_N)
 *
 * Summary:
 *   Compute a kernel basis of a full row rank matrix
 *
 * Description:
 *   The function computes inplace a right integer kernel basis N of full
 *   row rank matrix [A | B] taking M, the matrix computed by function
 *   specialHermite, and N as input.
 *
 *   The input matrix N looks like
 *
 *              k
 *         [ A^(-1)B ] n
 *   N = s*[---------]
 *         [  -I_k   ] k
 *
 *   where A a n x n matrix, B a n x k matrix, I_k a k x k identity matrix.
 *
 *   N is inplaced by
 *
 *         [ A^(-1)B ]
 *   N = s*[---------].H^(-1)
 *         [  -I_k   ]
 *
 *   where k x k matrix H is the hermite form of M computed by specialHermite.
 *
 * Input:
 *      n: long, row dimension of A
 *      k: long, column dimension of B
 *      t: long, same as parameter passed to function specialHermite (usually 1)
 *   mp_M: 1-dim mpz_t array length (n+k+t)*(k+t), computed by function
 *         specialHermite
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length (n+k)*k, representing a n+k x k kernel
 *         basis M, as shown above
 *
 * Note:
 *   the space of mp_N needs to be preallocated
 *
 */

void
kernelBasis (const long n, const long k, const long t,  mpz_t *mp_M, mpz_t *mp_N)
{
  long i, j, l;
  mpz_t mp_temp;

  mpz_init(mp_temp);
  for (i = 0; i < k; i++)
    {
      for (j = 0; j < i; j++)
	/* column operation N[1..-1, i] = N[1..-1, i]-N[1..-1, j]*M[j, i] */
	for (l = 0; l < n+k; l++)
	  {
	    mpz_mul(mp_temp, mp_N[l*k+j], mp_M[j*(k+t)+i]);
	    mpz_sub(mp_N[l*k+i], mp_N[l*k+i], mp_temp);
	  }

      /* column operation N[1..-1, i] = N[1..-1, i]/M[i, i] */
      for (l = 0; l < n+k; l++)
	{
	  mpz_divexact(mp_N[l*k+i], mp_N[l*k+i], mp_M[i*(k+t)+i]);
	}
    }

  mpz_clear(mp_temp);
}

#endif // __LINBOX_algorithm_iml_certified_solve_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


