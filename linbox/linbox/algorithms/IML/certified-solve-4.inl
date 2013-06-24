




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
		  mpz_t *mp_NZ, mpz_t mp_DZ,
		  bool red = 0, const long nullcol = NULLSPACE_COLUMN)
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
	      ver = specialminSolveRNS(certflag, red, nullcol, r, m,\
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
	  fflush(stdout);
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
	  fflush(stdout);
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






