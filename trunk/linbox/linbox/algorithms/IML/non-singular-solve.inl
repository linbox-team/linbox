#ifndef __LINBOX_algorithm_iml_non_singular_solve_INL
#define __LINBOX_algorithm_iml_non_singular_solve_INL


/*
 * Calling Sequence:
 *   nonsingSolvMM(solupos, n, m, A, mp_B, mp_N, mp_D)
 *
 * Summary:
 *   Solve nonsingular system of linear equations, where the left hand side
 *   input matrix is a signed long matrix.
 *
 * Description:
 *   Given a n x n nonsingular signed long matrix A, a n x m or m x n mpz_t
 *   matrix mp_B, this function will compute the solution of the system
 *   1. AX = mp_B
 *   2. XA = mp_B.
 *   The parameter solupos controls whether the system is in the type of 1
 *   or 2.
 *
 *   Since the unique solution X is a rational matrix, the function will
 *   output the numerator matrix mp_N and denominator mp_D respectively,
 *   such that A(mp_N) = mp_D*mp_B or (mp_N)A = mp_D*mp_B.
 *
 * Input:
 *   solupos: enumerate, flag to indicate the system to be solved
 *          - solupos = LeftSolu: solve XA = mp_B
 *          - solupos = RightSolu: solve AX = mp_B
 *         n: long, dimension of A
 *         m: long, column or row dimension of mp_B depending on solupos
 *         A: 1-dim signed long array length n*n, representing the n x n
 *            left hand side input matrix
 *      mp_B: 1-dim mpz_t array length n*m, representing the right hand side
 *            matrix of the system
 *          - solupos = LeftSolu: mp_B a m x n matrix
 *          - solupos = RightSolu: mp_B a n x m matrix
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length n*m, representing the numerator matrix
 *         of the solution
 *       - solupos = LeftSolu: mp_N a m x n matrix
 *       - solupos = RightSolu: mp_N a n x m matrix
 *   mp_D: mpz_t, denominator of the solution
 *
 * Precondition:
 *   A must be a nonsingular matrix.
 *
 * Note:
 *   - It is necessary to make sure the input parameters are correct,
 *     expecially the dimension, since there is no parameter checks in the
 *     function.
 *   - Input and output matrices are row majored and represented by
 *     one-dimension array.
 *   - It is needed to preallocate the memory space of mp_N and mp_D.
 *
 */

void
nonsingSolvMM (const enum SOLU_POS solupos, const long n, const long m, \
	       const long *A, mpz_t *mp_B, mpz_t *mp_N, mpz_t mp_D)
{
  long alpha, basislen, extbasislen, k=0, maxk=0, i, rt, kincr, ks=0, j, \
    l, minv;
  mpz_t mp_m, mp_alpha, mp_temp, mp_basisprod, mp_extbasisprod, mp_beta, \
    mp_nb, mp_db, mp_maxnb, mp_maxdb;
  mpz_t *mp_r;
  FiniteField *liftbasis, *extbdcoeff, *cmbasis, **extbasis;
  Double *liftbasisInv, **ARNS, **AInv, ***C, ***C1;
  double tt, tt1=0, tt2=0;
#if HAVE_TIME_H
  clock_t ti;
#endif

  { mpz_init(mp_nb); mpz_init(mp_db); mpz_init(mp_maxnb); mpz_init(mp_maxdb); }
  { mpz_init(mp_m); mpz_init(mp_beta); mpz_init(mp_alpha); mpz_init(mp_temp); }
  { mpz_init(mp_basisprod); mpz_init(mp_extbasisprod); }

  /* find lifting basis such that the product of the basis elements is
     approximately n*alpha */
  alpha = maxMagnLong(A, n, n, n);
  mpz_set_ui(mp_alpha, alpha);
  liftbasis = findLiftbasisSmall(n, mp_alpha, &basislen);

  /* initialization step before lifting */
  AInv = XMALLOC(Double *, basislen);
  for (i = 0; i < basislen; i++)
    AInv[i] = XMALLOC(Double, n*n);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  minv = liftInit(basislen, liftbasis, n, A, mp_basisprod, mp_extbasisprod, \
		  &extbasislen, &cmbasis, &extbdcoeff, &liftbasisInv, AInv,
		  &extbasis, &ARNS);

  /* if A^(-1) mod liftbasis[minv] does not exist, adjust liftbasis */
  while (minv != -1)
    {
      adBasis(minv, basislen, liftbasis);
      minv = liftInit(basislen, liftbasis, n, A, mp_basisprod, \
		      mp_extbasisprod, &extbasislen, &cmbasis, &extbdcoeff, \
		      &liftbasisInv, AInv, &extbasis, &ARNS);
    }
#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("              lifting initialization time: %f\n", tt);
  printf("              lifting basis length: %ld\n", basislen);
  for (i = 0; i < basislen; i++) { printf("   %ld", liftbasis[i]); }
  printf("\n");
  printf("              extended lifting basis length: %ld\n", extbasislen);
  for (i = 0; i < extbasislen; i++) { printf("    %ld", extbasis[0][i]); }
  printf("\n");
  fflush(stdout);
#endif

  mp_r = XMALLOC(mpz_t, m*n);
  for (i = 0; i < m*n; i++) { mpz_init(mp_r[i]); }
  mpz_set_ui(mp_D, 1);
  if (solupos == RightSolu)
    {
      maxMagnMP(mp_B, n, m, m, mp_beta);
      for (i = 0; i < m*n; i++) { mpz_set(mp_r[i], mp_B[i]); }
    }
  else if (solupos == LeftSolu)
    {
      maxMagnMP(mp_B, m, n, n, mp_beta);
      for (i = 0; i < m; i++)
	for (j = 0; j < n; j++)
	  mpz_set(mp_r[j*m+i], mp_B[i*n+j]);
    }

  /* set up k and bound of numerator and denominator */
  liftbd(mp_basisprod, n, mp_alpha, mp_beta, &maxk, mp_maxnb, mp_maxdb, &k, \
	 mp_nb, mp_db);
  kincr = k;
  ks = 0;
  C = NULL;
  do {

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    ti = clock();
#endif

    /* lifting kincr more steps */
    C1 = iml_lift(solupos, kincr, n, m, basislen, extbasislen, mp_basisprod, \
	      mp_extbasisprod, liftbasis, cmbasis, extbdcoeff, liftbasisInv, \
	      mp_r, extbasis, AInv, ARNS);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    tt1 += (double)(clock()-ti);
#endif

    /* update the lifting coefficients */
    C = XREALLOC(Double **, C, k);
    for (i = 0; i < kincr; i++) { C[i+ks] = C1[i]; }
    XFREE(C1);
    ks = k;

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    ti = clock();
#endif

    /* rational reconstruction */
    rt = soluRecon(solupos, k, basislen, n, m, mp_basisprod, liftbasis, \
		   cmbasis, C, mp_nb, mp_db, mp_N, mp_D);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    tt2 += (double)(clock()-ti);
#endif

    /* break the loop  when maximum step reached */
    if (k == maxk) { break; }

    /* rational reconstruction succeeds, check the result by magnitude bound */
    if (rt == 1)
      {
	if (solupos == RightSolu) { maxMagnMP(mp_N, n, m, m, mp_nb); }
	else if (solupos == LeftSolu) { maxMagnMP(mp_N, m, n, n, mp_nb); }
	mpz_mul_si(mp_nb, mp_nb, alpha);
	mpz_mul_ui(mp_nb, mp_nb, n);
	mpz_pow_ui(mp_m, mp_basisprod, k);
	mpz_sub_ui(mp_m, mp_m, 1);
	mpz_divexact_ui(mp_m, mp_m, 2);
	if (mpz_cmp(mp_nb, mp_m) < 0)
	  {
	    mpz_mul(mp_db, mp_D, mp_beta);
	    if (mpz_cmp(mp_db, mp_m) < 0)
	      break;
	  }
      }
    /* increase k and reset mp_nb and mp_db */
    kincr = (long)(0.1*k) > 20 ? (long)(0.1*k) : 20;
    if (k+kincr >= maxk)
      {
	kincr = maxk-k;
	k = maxk;
	mpz_set(mp_nb, mp_maxnb);
	mpz_set(mp_db, mp_maxdb);
	continue;
      }
    k += kincr;
    mpz_pow_ui(mp_m, mp_basisprod, k);
    mpz_sub_ui(mp_m, mp_m, 1);
    mpz_divexact_ui(mp_m, mp_m, 2);
    mpz_sqrt(mp_nb, mp_m);
    mpz_set(mp_db, mp_nb);
  } while (1);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  printf("              total lifting steps: %ld\n", k);
  tt1 = tt1/CLOCKS_PER_SEC;
  tt2 = tt2/CLOCKS_PER_SEC;
  printf("              lifting time: %f\n", tt1);
  printf("              rational reconstruction time: %f\n", tt2);
  fflush(stdout);
#endif


  for (i = 0; i < k; i++)
    {
      for (l = 0; l < basislen; l++) { XFREE(C[i][l]); }
      XFREE(C[i]);
    }
  XFREE(C);
  for (i = 0; i < m*n; i++) { mpz_clear(mp_r[i]); } { XFREE(mp_r); }
  for (i = 0; i < 2; i++) { XFREE(extbasis[i]); } { XFREE(extbasis); }
  for (i = 0; i < basislen; i++) { XFREE(AInv[i]); } { XFREE(AInv); }
  for (i = 0; i < extbasislen; i++) { XFREE(ARNS[i]); } { XFREE(ARNS); }
  { XFREE(liftbasis); XFREE(extbdcoeff); XFREE(cmbasis); XFREE(liftbasisInv); }
  { mpz_clear(mp_nb); mpz_clear(mp_db); mpz_clear(mp_maxnb); mpz_clear(mp_m); }
  { mpz_clear(mp_maxdb); mpz_clear(mp_beta); mpz_clear(mp_alpha); }
  { mpz_clear(mp_basisprod); mpz_clear(mp_extbasisprod); mpz_clear(mp_temp); }


  return;
}



/*
 * Calling Sequence:
 *   nonsingSolvLlhsMM(solupos, n, m, mp_A, mp_B, mp_N, mp_D)
 *
 * Summary:
 *   Solve nonsingular system of linear equations, where the left hand side
 *   input matrix is a mpz_t matrix.
 *
 * Description:
 *   Given a n x n nonsingular mpz_t matrix A, a n x m or m x n mpz_t
 *   matrix mp_B, this function will compute the solution of the system
 *   1. (mp_A)X = mp_B
 *   2. X(mp_A) = mp_B.
 *   The parameter solupos controls whether the system is in the type of 1
 *   or 2.
 *
 *   Since the unique solution X is a rational matrix, the function will
 *   output the numerator matrix mp_N and denominator mp_D respectively,
 *   such that mp_Amp_N = mp_D*mp_B or mp_Nmp_A = mp_D*mp_B.
 *
 * Input:
 *   solupos: enumerate, flag to indicate the system to be solved
 *          - solupos = LeftSolu: solve XA = mp_B
 *          - solupos = RightSolu: solve AX = mp_B
 *         n: long, dimension of A
 *         m: long, column or row dimension of mp_B depending on solupos
 *      mp_A: 1-dim mpz_t array length n*n, representing the n x n left hand
 *            side input matrix
 *      mp_B: 1-dim mpz_t array length n*m, representing the right hand side
 *            matrix of the system
 *          - solupos = LeftSolu: mp_B a m x n matrix
 *          - solupos = RightSolu: mp_B a n x m matrix
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length n*m, representing the numerator matrix
 *         of the solution
 *       - solupos = LeftSolu: mp_N a m x n matrix
 *       - solupos = RightSolu: mp_N a n x m matrix
 *   mp_D: mpz_t, denominator of the solution
 *
 * Precondition:
 *   mp_A must be a nonsingular matrix.
 *
 * Note:
 *   - It is necessary to make sure the input parameters are correct,
 *     expecially the dimension, since there is no parameter checks in the
 *     function.
 *   - Input and output matrices are row majored and represented by
 *     one-dimension array.
 *   - It is needed to preallocate the memory space of mp_N and mp_D.
 *
 */

void
nonsingSolvLlhsMM (const enum SOLU_POS solupos, const long n, \
		   const long m, mpz_t *mp_A, mpz_t *mp_B, mpz_t *mp_N, \
		   mpz_t mp_D)
{
  long basislen, extbasislen, k=0, maxk=0, i, rt, kincr, ks=0, j, l, minv;
  mpz_t mp_m, mp_alpha, mp_temp, mp_basisprod, mp_extbasisprod, mp_beta, \
    mp_nb, mp_db, mp_maxnb, mp_maxdb;
  mpz_t *mp_r;
  FiniteField *liftbasis, *extbdcoeff, *cmbasis, **extbasis;
  Double *liftbasisInv, **ARNS, **AInv, ***C, ***C1;

  { mpz_init(mp_nb); mpz_init(mp_db); mpz_init(mp_maxnb); mpz_init(mp_maxdb); }
  { mpz_init(mp_m); mpz_init(mp_beta); mpz_init(mp_alpha); mpz_init(mp_temp); }
  { mpz_init(mp_basisprod); mpz_init(mp_extbasisprod); }

  /* find lifting basis such that the product of the basis elements is
     approximately n*mp_alpha */
  maxMagnMP(mp_A, n, n, n, mp_alpha);
  liftbasis = findLiftbasisLarge(n, mp_alpha, &basislen);

  /* initialize for lifting */
  AInv = XMALLOC(Double *, basislen);
  for (i = 0; i < basislen; i++)
    AInv[i] = XMALLOC(Double, n*n);
  minv = liftInitLlhs(basislen, liftbasis, n, mp_A, mp_basisprod, \
		      mp_extbasisprod, &extbasislen, &cmbasis,  &extbdcoeff, \
		      &liftbasisInv, AInv, &extbasis, &ARNS);

  /* if A^(-1) mod liftbasis[minv] does not exist, adjust liftbasis */
  while (minv != -1)
    {
      adBasis(minv, basislen, liftbasis);
      minv = liftInitLlhs(basislen, liftbasis, n, mp_A, mp_basisprod, \
			  mp_extbasisprod,  &extbasislen, &cmbasis, \
			  &extbdcoeff, &liftbasisInv, AInv, &extbasis, &ARNS);
    }
  mp_r = XMALLOC(mpz_t, m*n);
  for (i = 0; i < m*n; i++) { mpz_init(mp_r[i]); }
  mpz_set_ui(mp_D, 1);
  if (solupos == RightSolu)
    {
      maxMagnMP(mp_B, n, m, m, mp_beta);
      for (i = 0; i < m*n; i++) { mpz_set(mp_r[i], mp_B[i]); }
    }
  else if (solupos == LeftSolu)
    {
      maxMagnMP(mp_B, m, n, n, mp_beta);
      for (i = 0; i < m; i++)
	for (j = 0; j < n; j++)
	  mpz_set(mp_r[j*m+i], mp_B[i*n+j]);
    }
  /* set up k and bound of numerator and denominator */
  liftbd(mp_basisprod, n, mp_alpha, mp_beta, &maxk, mp_maxnb, mp_maxdb, \
	 &k, mp_nb, mp_db);
  kincr = k;
  ks = 0;
  C = NULL;
  do {
    /* lifting kincr more steps */
    C1 = iml_lift(solupos, kincr, n, m, basislen, extbasislen, mp_basisprod, \
	      mp_extbasisprod, liftbasis, cmbasis, extbdcoeff, liftbasisInv, \
	      mp_r, extbasis, AInv, ARNS);

    /* update the lifting coefficients */
    C = XREALLOC(Double **, C, k);
    for (i = 0; i < kincr; i++) { C[i+ks] = C1[i]; }
    XFREE(C1);
    ks = k;

    /* rational reconstruction */
    rt = soluRecon(solupos, k, basislen, n, m, mp_basisprod, liftbasis, \
		   cmbasis, C, mp_nb, mp_db, mp_N, mp_D);

    /* break the loop  when maximum step reached */
    if (k == maxk) { break; }

    /* rational reconstruction succeeds, check the result by magnitude bound */
    if (rt == 1)
      {
	if (solupos == RightSolu) { maxMagnMP(mp_N, n, m, m, mp_nb); }
	else if (solupos == LeftSolu) { maxMagnMP(mp_N, m, n, n, mp_nb); }
	mpz_mul(mp_nb, mp_nb, mp_alpha);
	mpz_mul_ui(mp_nb, mp_nb, n);
	mpz_pow_ui(mp_m, mp_basisprod, k);
	mpz_sub_ui(mp_m, mp_m, 1);
	mpz_divexact_ui(mp_m, mp_m, 2);
	if (mpz_cmp(mp_nb, mp_m) < 0)
	  {
	    mpz_mul(mp_db, mp_D, mp_beta);
	    if (mpz_cmp(mp_db, mp_m) < 0)
	      break;
	  }
      }
    /* increase k and reset mp_nb and mp_db */
    kincr = (long)(0.1*k) > 20 ? (long)(0.1*k) : 20;
    if (k+kincr >= maxk)
      {
	kincr = maxk-k;
	k = maxk;
	mpz_set(mp_nb, mp_maxnb);
	mpz_set(mp_db, mp_maxdb);
	continue;
      }
    k += kincr;
    mpz_pow_ui(mp_m, mp_basisprod, k);
    mpz_sub_ui(mp_m, mp_m, 1);
    mpz_divexact_ui(mp_m, mp_m, 2);
    mpz_sqrt(mp_nb, mp_m);
    mpz_set(mp_db, mp_nb);
  } while (1);

  for (i = 0; i < k; i++)
    {
      for (l = 0; l < basislen; l++) { XFREE(C[i][l]); }
      XFREE(C[i]);
    }
  XFREE(C);
  for (i = 0; i < m*n; i++) { mpz_clear(mp_r[i]); } { XFREE(mp_r); }
  for (i = 0; i < 2; i++) { XFREE(extbasis[i]); } { XFREE(extbasis); }
  for (i = 0; i < basislen; i++) { XFREE(AInv[i]); } { XFREE(AInv); }
  for (i = 0; i < extbasislen; i++) { XFREE(ARNS[i]); } { XFREE(ARNS); }
  { XFREE(liftbasis); XFREE(extbdcoeff); XFREE(cmbasis); XFREE(liftbasisInv); }
  { mpz_clear(mp_nb); mpz_clear(mp_db); mpz_clear(mp_maxnb); mpz_clear(mp_m); }
  { mpz_clear(mp_maxdb); mpz_clear(mp_beta); mpz_clear(mp_alpha); }
  { mpz_clear(mp_basisprod); mpz_clear(mp_extbasisprod); mpz_clear(mp_temp); }

  return;
}



/*
 * Calling Sequence:
 *   nonsingSolvRNSMM(solupos, basislen, n, m, basis, ARNS, mp_B, mp_N, mp_D)
 *
 * Summary:
 *   Solve nonsingular system of linear equations, where the left hand side
 *   input matrix is represented in a RNS.
 *
 * Description:
 *   Given a n x n nonsingular matrix A represented in a RNS, a n x m or m x n
 *   mpz_t matrix mp_B, this function will compute the solution of the system
 *   1. AX = mp_B
 *   2. XA = mp_B.
 *   The parameter solupos controls whether the system is in the type of 1
 *   or 2.
 *
 *   Since the unique solution X is a rational matrix, the function will
 *   output the numerator matrix mp_N and denominator mp_D respectively,
 *   such that A(mp_N) = mp_D*mp_B or (mp_N)A = mp_D*mp_B.
 *
 * Input:
 *    solupos: enumerate, flag to indicate the system to be solved
 *           - solupos = LeftSolu: solve XA = mp_B
 *           - solupos = RightSolu: solve AX = mp_B
 *   basislen: long, dimension of RNS basis
 *          n: long, dimension of A
 *          m: long, column or row dimension of mp_B depending on solupos
 *      basis: 1-dim FiniteField array length basislen, RNS basis
 *       ARNS: 2-dim Double array, dimension basislen x n*n, representation of
 *             n x n input matrix A in RNS, where ARNS[i] = A mod basis[i]
 *       mp_B: 1-dim mpz_t array length n*m, representing the right hand side
 *             matrix of the system
 *           - solupos = LeftSolu: mp_B a m x n matrix
 *           - solupos = RightSolu: mp_B a n x m matrix
 *
 * Output:
 *   mp_N: 1-dim mpz_t array length n*m, representing the numerator matrix
 *         of the solution
 *       - solupos = LeftSolu: mp_N a m x n matrix
 *       - solupos = RightSolu: mp_N a n x m matrix
 *   mp_D: mpz_t, denominator of the solution
 *
 * Precondition:
 *   - A must be a nonsingular matrix.
 *   - Any element p in RNS basis must satisfy 2*(p-1)^2 <= 2^53-1.
 *
 * Note:
 *   - It is necessary to make sure the input parameters are correct,
 *     expecially the dimension, since there is no parameter checks in the
 *     function.
 *   - Input and output matrices are row majored and represented by
 *     one-dimension array.
 *   - It is needed to preallocate the memory space of mp_N and mp_D.
 *
 */

void
nonsingSolvRNSMM (const enum SOLU_POS solupos, const long n, const long m,\
		  const long basislen, const FiniteField *basis, \
		  Double **ARNS, mpz_t *mp_B, mpz_t *mp_N, mpz_t mp_D)
{
  long liftbasislen, extbasislen, k=0, maxk=0, i, rt, kincr, ks=0, j, l, minv;
  mpz_t mp_m, mp_alpha, mp_liftbasisprod, mp_extbasisprod, mp_beta, \
    mp_nb, mp_db, mp_maxnb, mp_maxdb;
  mpz_t *mp_r;
  FiniteField *liftbasis, *extbdcoeff, *cmliftbasis, **extbasis;
  Double *liftbasisInv, **AExtRNS, **AInv, ***C, ***C1;
  double tt, tt1=0, tt2=0;

#if HAVE_TIME_H
  clock_t ti;
#endif


  { mpz_init(mp_nb); mpz_init(mp_db); mpz_init(mp_maxnb); mpz_init(mp_maxdb); }
  { mpz_init(mp_m); mpz_init(mp_beta); mpz_init(mp_alpha); }
  { mpz_init(mp_liftbasisprod); mpz_init(mp_extbasisprod); }

  /* find lifting basis such that the product of the basis elements is
     approximately n*mp_alpha */
  basisProd(basislen, basis, mp_alpha);
  liftbasis = findLiftbasisLarge(n, mp_alpha, &liftbasislen);

  /* initialization step before lifting */
  AInv = XMALLOC(Double *, liftbasislen);
  for (i = 0; i < liftbasislen; i++)
    AInv[i] = XMALLOC(Double, n*n);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  minv = liftInitRNS(liftbasislen, liftbasis, basislen, basis, n, ARNS, \
		     mp_liftbasisprod, mp_extbasisprod, &extbasislen, \
		     &cmliftbasis, &extbdcoeff, &liftbasisInv, AInv, \
		     &extbasis, &AExtRNS);

  /* if A^(-1) mod liftbasis[minv] does not exist, adjust liftbasis */
  while (minv != -1)
    {
      adBasis(minv, liftbasislen, liftbasis);
      minv = liftInitRNS(liftbasislen, liftbasis, basislen, basis, n, ARNS, \
			 mp_liftbasisprod, mp_extbasisprod, &extbasislen, \
			 &cmliftbasis, &extbdcoeff, &liftbasisInv, AInv, \
			 &extbasis, &AExtRNS);
    }

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("              lifting initialization time: %f\n", tt);
  printf("              lifting basis length: %ld\n", liftbasislen);
  for (i = 0; i < liftbasislen; i++)
    printf("          %ld  ", liftbasis[i]);
  printf("\n");
  printf("              extended lifting basis length: %ld\n", extbasislen);
  for (i = 0; i < extbasislen; i++)
    printf("           %ld  ", extbasis[0][i]);
  printf("\n");
  fflush(stdout);
#endif

  mp_r = XMALLOC(mpz_t, m*n);
  for (i = 0; i < m*n; i++) { mpz_init(mp_r[i]); }
  mpz_set_ui(mp_D, 1);
  if (solupos == RightSolu)
    {
      maxMagnMP(mp_B, n, m, m, mp_beta);
      for (i = 0; i < m*n; i++) { mpz_set(mp_r[i], mp_B[i]); }
    }
  else if (solupos == LeftSolu)
    {
      maxMagnMP(mp_B, m, n, n, mp_beta);
      for (i = 0; i < m; i++)
	for (j = 0; j < n; j++)
	  mpz_set(mp_r[j*m+i], mp_B[i*n+j]);
    }
  /* set up k and bound of numerator and denominator */
  liftbd(mp_liftbasisprod, n, mp_alpha, mp_beta, &maxk, mp_maxnb, mp_maxdb, \
	 &k, mp_nb, mp_db);
  kincr = k;
  ks = 0;
  C = NULL;
  do {

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    ti = clock();
#endif

    /* lifting kincr more steps */
    C1 = iml_lift(solupos, kincr, n, m, liftbasislen, extbasislen, \
	      mp_liftbasisprod,  mp_extbasisprod, liftbasis, cmliftbasis, \
	      extbdcoeff, liftbasisInv, mp_r, extbasis, AInv, AExtRNS);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    tt1 += (double)(clock()-ti);
#endif

    /* update the lifting coefficients */
    C = XREALLOC(Double **, C, k);
    for (i = 0; i < kincr; i++) { C[i+ks] = C1[i]; }
    XFREE(C1);
    ks = k;

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    ti = clock();
#endif

    /* rational reconstruction */
    rt = soluRecon(solupos, k, liftbasislen, n, m, mp_liftbasisprod, \
		   liftbasis, cmliftbasis, C, mp_nb, mp_db, mp_N, mp_D);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    tt2 += (double)(clock()-ti);
#endif

    /* break the loop  when maximum lifting step reached */
    if (k == maxk) { break; }

    /* rational reconstruction succeed, check the result by magnitude bound */
    if (rt == 1)
      {
	if (solupos == RightSolu) { maxMagnMP(mp_N, n, m, m, mp_nb); }
	else if (solupos == LeftSolu) { maxMagnMP(mp_N, m, n, n, mp_nb); }
	mpz_mul(mp_nb, mp_nb, mp_alpha);
	mpz_mul_ui(mp_nb, mp_nb, n);
	mpz_pow_ui(mp_m, mp_liftbasisprod, k);
	mpz_sub_ui(mp_m, mp_m, 1);
	mpz_divexact_ui(mp_m, mp_m, 2);
	if (mpz_cmp(mp_nb, mp_m) < 0)
	  {
	    mpz_mul(mp_db, mp_D, mp_beta);
	    if (mpz_cmp(mp_db, mp_m) < 0)
	      break;
	  }
      }

    /* increase k and reset mp_nb and mp_db */
    kincr = (long)(0.1*k) > 20 ? (long)(0.1*k) : 20;
    if (k+kincr >= maxk)
      {
	kincr = maxk-k;
	k = maxk;
	mpz_set(mp_nb, mp_maxnb);
	mpz_set(mp_db, mp_maxdb);
	continue;
      }
    k += kincr;
    mpz_pow_ui(mp_m, mp_liftbasisprod, k);
    mpz_sub_ui(mp_m, mp_m, 1);
    mpz_divexact_ui(mp_m, mp_m, 2);
    mpz_sqrt(mp_nb, mp_m);
    mpz_set(mp_db, mp_nb);
  } while (1);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  printf("              total lifting steps: %ld\n", k);
  tt1 = tt1/CLOCKS_PER_SEC;
  tt2 = tt2/CLOCKS_PER_SEC;
  printf("              lifting time: %f\n", tt1);
  printf("              rational reconstruction time: %f\n", tt2);
  fflush(stdout);
#endif

  for (i = 0; i < k; i++)
    {
      for (l = 0; l < liftbasislen; l++) { XFREE(C[i][l]); }
      XFREE(C[i]);
    }
  XFREE(C);
  for (i = 0; i < m*n; i++) { mpz_clear(mp_r[i]); } { XFREE(mp_r); }
  for (i = 0; i < 2; i++) { XFREE(extbasis[i]); } { XFREE(extbasis); }
  for (i = 0; i < liftbasislen; i++) { XFREE(AInv[i]); } { XFREE(AInv); }
  for (i = 0; i < extbasislen; i++) { XFREE(AExtRNS[i]); } { XFREE(AExtRNS); }
  { XFREE(liftbasis); XFREE(extbdcoeff); XFREE(cmliftbasis); }
  XFREE(liftbasisInv);
  { mpz_clear(mp_nb); mpz_clear(mp_db); mpz_clear(mp_maxnb); mpz_clear(mp_m); }
  { mpz_clear(mp_maxdb); mpz_clear(mp_beta); mpz_clear(mp_alpha); }
  { mpz_clear(mp_liftbasisprod); mpz_clear(mp_extbasisprod); }
  return;
}



/*
 * Calling Sequence:
 *   liftbasis <-- findLiftbasisSmall(n, mp_alpha, basislen)
 *
 * Summary:
 *   Compute the p-adic lifting basis
 *
 * Description:
 *   The function computes p-adic lifting basis liftbasis, such that after
 *   every single lifting step, log[2](p) solution bits will be gained, where
 *   p is the product of the elements in the lifting basis.
 *
 *   Let the lifting basis be 'liftbasis' and dimension be m, then
 *   - elements in liftbasis are all primes
 *   - liftbasis[0] is the largest prime among all the primes less than
 *     RNSbound(n)
 *   - liftbasis[i+1] is the next prime smaller than liftbasis[i] (i = 0..m-2)
 *   - product of first m-2 elements in liftbasis <= mp_maxInter
 *   - product of first m-1 elements in liftbasis > mp_maxInter
 *
 * Input:
 *          n: long, dimension of A
 *   mp_alpha: mpz_t, maximum magnitude of A
 *   basislen: pointer to a long integer, the dimension of lifting basis
 *
 * Return:
 *   liftbasis: 1-dime FiniteField length *basislen, storing the lifting basis
 */

FiniteField *
findLiftbasisSmall (const long n, const mpz_t mp_alpha, long *basislen)
{
  long temp, len=0, count;
  FiniteField *liftbasis;
  mpz_t mp_temp, mp_bd, mp_prod;

  /* compute the upper bound of lifting basis */
  temp = RNSbound(n);
  mpz_init_set_ui(mp_temp, temp);

  /* compute n*mp_alpha */
  mpz_init_set_ui(mp_bd, n);
  mpz_mul(mp_bd, mp_bd, mp_alpha);
  mpz_init_set_ui(mp_prod, 1);
  liftbasis = NULL;
  while (mpz_cmp(mp_bd, mp_prod) > 0)
    {
      ++len;
      liftbasis = XREALLOC(FiniteField, liftbasis, len);
      while (mpz_probab_prime_p(mp_temp, 10) == 0)
	mpz_sub_ui(mp_temp, mp_temp, 1);
      liftbasis[len-1] = mpz_get_ui(mp_temp);
      mpz_sub_ui(mp_temp, mp_temp, 1);
      mpz_mul_ui(mp_prod, mp_prod, liftbasis[len-1]);
    }

  /* increase length to 2*len+4 */
  count = len+4;
  while (--count > 0)
    {
      ++len;
      liftbasis = XREALLOC(FiniteField, liftbasis, len);
      while (mpz_probab_prime_p(mp_temp, 10) == 0)
	mpz_sub_ui(mp_temp, mp_temp, 1);
      liftbasis[len-1] = mpz_get_ui(mp_temp);
      mpz_sub_ui(mp_temp, mp_temp, 1);
    }
  *basislen = len;
  { mpz_clear(mp_temp); mpz_clear(mp_bd); mpz_clear(mp_prod); }
  return liftbasis;
}



/*
 * Calling Sequence:
 *   liftbasis <-- findLiftbasisLarge(n, mp_alpha, basislen)
 *
 * Summary:
 *   Compute the p-adic lifting basis
 *
 * Description:
 *   The function computes p-adic lifting basis liftbasis, such that after
 *   every single lifting step, log[2](p) solution bits will be gained, where
 *   p is the product of the elements in the lifting basis.
 *
 *   Let the lifting basis be 'liftbasis' and dimension be m, then
 *   - elements in liftbasis are all primes
 *   - liftbasis[0] is the largest prime among all the primes less than
 *     RNSbound(n)
 *   - liftbasis[i+1] is the next prime smaller than liftbasis[i] (i = 0..m-2)
 *   - product of first m-2 elements in liftbasis <= mp_maxInter
 *   - product of first m-1 elements in liftbasis > mp_maxInter
 *
 * Input:
 *          n: long, dimension of A
 *   mp_alpha: mpz_t, maximum magnitude of A
 *   basislen: pointer to a long integer, the dimension of lifting basis
 *
 * Return:
 *   liftbasis: 1-dime FiniteField length *basislen, storing the lifting basis
 */

FiniteField *
findLiftbasisLarge (const long n, const mpz_t mp_alpha, long *basislen)
{
  long temp, len=0, count;
  FiniteField *liftbasis;
  mpz_t mp_temp, mp_bd, mp_prod;

  /* compute the upper bound of lifting basis */
  temp = RNSbound(n);
  mpz_init_set_ui(mp_temp, temp);

  /* compute n*mp_alpha */
  mpz_init_set_ui(mp_bd, n);
  mpz_mul(mp_bd, mp_bd, mp_alpha);
  mpz_init_set_ui(mp_prod, 1);
  liftbasis = NULL;
  while (mpz_cmp(mp_bd, mp_prod) > 0)
    {
      ++len;
      liftbasis = XREALLOC(FiniteField, liftbasis, len);
      while (mpz_probab_prime_p(mp_temp, 10) == 0)
	mpz_sub_ui(mp_temp, mp_temp, 1);
      liftbasis[len-1] = mpz_get_ui(mp_temp);
      mpz_sub_ui(mp_temp, mp_temp, 1);
      mpz_mul_ui(mp_prod, mp_prod, liftbasis[len-1]);
    }

  count = 3;
  while (--count > 0)
    {
      ++len;
      liftbasis = XREALLOC(FiniteField, liftbasis, len);
      while (mpz_probab_prime_p(mp_temp, 10) == 0)
	mpz_sub_ui(mp_temp, mp_temp, 1);
      liftbasis[len-1] = mpz_get_ui(mp_temp);
      mpz_sub_ui(mp_temp, mp_temp, 1);
      }


  *basislen = len;
  { mpz_clear(mp_temp); mpz_clear(mp_bd); mpz_clear(mp_prod); }
  return liftbasis;
}



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
  mpz_t mp_temp;
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
