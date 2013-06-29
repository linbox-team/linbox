#ifndef __LINBOX_algorithm_iml_non_singular_solve_INL
#define __LINBOX_algorithm_iml_non_singular_solve_INL

namespace LinBox { namespace iml {


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
		// mpz_init(mp_temp);

		/* move the elements forward */
		for (i = idx+1; i < basislen; i++) {
			liftbasis[i-1] = liftbasis[i];
		}
		mp_temp = (unsigned long) liftbasis[basislen-1];
		// mpz_set_ui(mp_temp, liftbasis[basislen-1]);
		Integer::subin(mp_temp, 1UL);
		// mpz_sub_ui(mp_temp, mp_temp, 1);
		Integer::prevprime(mp_temp);
		// while (mpz_probab_prime_p(mp_temp, 10) == 0)
		// mpz_sub_ui(mp_temp, mp_temp, 1);
		liftbasis[basislen-1] = (unsigned int) (mp_temp);
		// liftbasis[basislen-1] = mpz_get_ui(mp_temp);

		// mpz_clear(mp_temp);
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
		Integer::mul(mp_maxdb, mp_db, mp_t1);
		// mpz_mul(mp_maxdb, mp_db, mp_t1);
		Integer::pow(mp_nb, mp_alpha, (unsigned long) n-1);
		// mpz_pow_ui(mp_nb, mp_alpha, n-1);
		Integer::mulin(mp_nb, mp_beta);
		// mpz_mul(mp_nb, mp_nb, mp_beta);
		Integer::mul(mp_maxnb, mp_nb, mp_t1);
		// mpz_mul(mp_maxnb, mp_nb, mp_t1);
		// mp_prod = mp_maxdb;
		// mpz_init_set(mp_prod, mp_maxdb);
		Integer::mul(mp_prod, mp_maxdb, mp_maxnb);
		// mpz_mul(mp_prod, mp_prod, mp_maxnb);
		Integer::mulin(mp_prod, 2UL);
		// mpz_mul_ui(mp_prod, mp_prod, 2);
		Integer::addin(mp_prod,1UL);
		// mpz_add_ui(mp_prod, mp_prod, 1);

		/* compute maxk */
		maxk = 1;
		mp_t1 = mp_basisprod ;
		// mpz_set(mp_t1, mp_basisprod);
		while ( mp_t1  < mp_prod )
		{
			Integer::mulin(mp_t1,mp_basisprod);
			// mpz_mul(mp_t1, mp_t1, mp_basisprod);
			++ maxk;
			// ++(*maxk);
		}

		/* compute k and estimate bound */
		k = 20;
		Integer::pow(mp_prod, mp_basisprod, k);
		// mpz_pow_ui(mp_prod, mp_basisprod, *k);
		Integer::subin(mp_prod, 1UL);
		// mpz_sub_ui(mp_prod, mp_prod, 1);
		Integer::divexact(mp_prod, mp_prod, 2UL);
		// mpz_divexact_ui(mp_prod, mp_prod, 2);
		Integer::sqrt(mp_nb,mp_prod);
		// mpz_sqrt(mp_nb, mp_prod);
		mp_db = mp_nb ;
		// mpz_set(mp_db, mp_nb);
		if (k >= maxk) {
			mp_nb = mp_maxnb ;
			// mpz_set(mp_nb, mp_maxnb);
			mp_db = mp_maxdb ;
			// mpz_set(mp_db, mp_maxdb);
		}
		// mpz_clear(mp_t1);
		// mpz_clear(mp_t2);
		// mpz_clear(mp_prod);
		return;
	}

	template<class Field>
void
findLiftbasis(std::vector<unsigned long> liftbasis,
	      const long & n,
	      const Integer & mp_alpha,
	      // long & basislen,
	      bool small)
{
  long temp, len=0, count;
  Integer mp_temp, mp_bd, mp_prod;

  /* compute the upper bound of lifting basis */
  temp = RNS<Field>::RNSbound(n);
  mp_temp = temp;
  // mpz_init_set_ui(mp_temp, temp);

  /* compute n*mp_alpha */
  // mp_bd = n;
  // mpz_init_set_ui(mp_bd, n);
  Integer::mul(mp_bd, n, mp_alpha);
  // mpz_mul(mp_bd, mp_bd, mp_alpha);
  mp_prod =  1UL;
  // mpz_init_set_ui(mp_prod, 1);
  // liftbasis = NULL;
  while (mp_bd > mp_prod )
    {
      ++len;
      // liftbasis = XREALLOC(FiniteField, liftbasis, len);
      liftbasis.resise(len);
      // while (mpz_probab_prime_p(mp_temp, 10) == 0)
	// mpz_sub_ui(mp_temp, mp_temp, 1);
      Integer::prevprime(mp_temp,mp_temp);
      liftbasis[len-1] = (unsigned int)(mp_temp);
      // liftbasis[len-1] = mpz_get_ui(mp_temp);
      Integer::subin(mp_temp, 1UL);
      // mpz_sub_ui(mp_temp, mp_temp, 1);
      Integer::mulin(mp_prod, liftbasis[len-1]);
      // mpz_mul_ui(mp_prod, mp_prod, liftbasis[len-1]);
    }

  /* increase length to 2*len+4 */

  if (small)
	  count = len+4;
  else // large
	  count = 3 ;

  while (--count > 0)
    {
      ++len;
      // liftbasis = XREALLOC(FiniteField, liftbasis, len);
      liftbasis.resize(len);
      // while (mpz_probab_prime_p(mp_temp, 10) == 0)
	// mpz_sub_ui(mp_temp, mp_temp, 1);
      Integer::prevprime(mp_temp,mp_temp);
      liftbasis[len-1] = (unsigned int)(mp_temp);
      // liftbasis[len-1] = mpz_get_ui(mp_temp);
      Integer::subin(mp_temp, 1UL);
      // mpz_sub_ui(mp_temp, mp_temp, 1);
    }

  // *basislen = len;
  // { mpz_clear(mp_temp); mpz_clear(mp_bd); mpz_clear(mp_prod); }
  // return liftbasis;
  return ;
}

/*
 * Calling Sequence:
 *   nonsingSolvLlhsMM(solupos, n, m, mp_A, mp_B, mp_N, mp_D)
 *
 * Summary:
 *   Solve nonsingular system of linear equations, where the left hand side
 *   input matrix is a Integer matrix.
 *
 * Description:
 *   Given a n x n nonsingular Integer matrix A, a n x m or m x n mpz_t
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
 *      mp_A: 1-dim Integer array length n*n, representing the n x n left hand
 *            side input matrix
 *      mp_B: 1-dim Integer array length n*m, representing the right hand side
 *            matrix of the system
 *          - solupos = LeftSolu: mp_B a m x n matrix
 *          - solupos = RightSolu: mp_B a n x m matrix
 *
 * Output:
 *   mp_N: 1-dim Integer array length n*m, representing the numerator matrix
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

template<class Matrix>
void
nonsingSolv(const Tag::Side solupos, Matrix & A, BlasVector<PID_Integer> & mp_B, BlasVector<PID_Integer> & mp_N,  Integer & mp_D)
{

  double tt, tt1=0, tt2=0;
#if HAVE_TIME_H
  clock_t ti;
#endif


  /* find lifting basis such that the product of the basis elements is
     approximately n*mp_alpha */
  magnitude(mp_alpha,mp_A);

  bool large = isInteger<Element>();
  findLiftbasis(liftbasis, n, mp_alpha, &basislen, large);

  /* initialization step before lifting */
  // BlasMatrix<Field> zz(Field(), n, n );
  // std::vector<BlasMatrix<Field> > Ainv(basislen, zz);

  // AInv = XMALLOC(Double *, basislen);
  // for (i = 0; i < basislen; i++)
    // AInv[i] = XMALLOC(Double, n*n);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  ti = clock();
#endif

  pAdicLift PAL ;
  minv = PAL.liftInit();
  // minv = liftInitLlhs(basislen, liftbasis,
		  // n, mp_A, mp_basisprod, \ mp_extbasisprod,
		      // &extbasislen, &cmbasis,  &extbdcoeff,  &liftbasisInv, AInv,
		      // &extbasis, &ARNS);

  /* if A^(-1) mod liftbasis[minv] does not exist, adjust liftbasis */
  while (minv != -1)
    {
      adBasis(minv, basislen, liftbasis);
      minv = PAL.liftInit(mp_A);
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

  // mp_r = XMALLOC(mpz_t, m*n);
  BlasMatrix<PID_Integer> mpr(Z,n,m);
  // for (i = 0; i < m*n; i++) { mpz_init(mp_r[i]); }
  // mpz_set_ui(mp_D, 1);
  mp_D = 1 ;
  magnitude(mp_beta,mp_B);
  if (solupos == RightSolu)
    {
      // maxMagnMP(mp_B, n, m, m, mp_beta);
      // XXX copy
      // for (i = 0; i < m*n; i++) {
	      // mpz_set(mp_r[i], mp_B[i]);
      // }
      mp_r.copy(mp_B);
    }
  else if (solupos == LeftSolu)
    {
      // maxMagnMP(mp_B, m, n, n, mp_beta);
      for (i = 0; i < m; i++)
	for (j = 0; j < n; j++)
	  // mpz_set(mp_r[j*m+i], mp_B[i*n+j]);
	  mp_r.setEntry(j,i,mp_B.getEntry(i,j);
    }

  /* set up k and bound of numerator and denominator */
  liftbd(mp_basisprod, n, mp_alpha, mp_beta, &maxk, mp_maxnb, mp_maxdb,  &k, \
		  mp_nb, mp_db);
  kincr = k;
  ks = 0;
  C = NULL;
  do {

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    ti = clock();
#endif

    /* lifting kincr more steps */

    PAL.iml_lift(solupos, C1, kincr,mp_r);
#if HAVE_VERBOSE_MODE && HAVE_TIME_H
    tt1 += (double)(clock()-ti);
#endif

    /* update the lifting coefficients */
    C.resize(k);
    // C = XREALLOC(Double **, C, k);
    //! @todo C pointer array : no copy !
    for (i = 0; i < kincr; i++) { C[i+ks].copy(C1[i]); }
    // XFREE(C1);
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
	magnitude(mp_nb,mp_N);
	// if (solupos == RightSolu) { maxMagnMP(mp_N, n, m, m, mp_nb); }
	// else if (solupos == LeftSolu) { maxMagnMP(mp_N, m, n, n, mp_nb); }
	Integer::mulin(mp_nb, mp_alpha);
	// mpz_mul(mp_nb, mp_nb, mp_alpha);
	Integer::mulin(mp_nb, n);
	// mpz_mul_ui(mp_nb, mp_nb, n);
	Integer::pow(mp_m, mp_basisprod, k);
	// mpz_pow_ui(mp_m, mp_basisprod, k);
	Integer::subin(mp_m, 1UL);
	// mpz_sub_ui(mp_m, mp_m, 1);
	Integer::divexact_ui(mp_m, mp_m, 2UL);
	// mpz_divexact_ui(mp_m, mp_m, 2);
	if (mp_nb < mp_m )
	  {
		  Integer::mul(mp_db, mp_D, mp_beta);
	    // mpz_mul(mp_db, mp_D, mp_beta);
	    if (mp_db, < mp_m )
	      break;
	  }
      }
    /* increase k and reset mp_nb and mp_db */
    kincr = (long)(0.1*k) > 20 ? (long)(0.1*k) : 20;
    if (k+kincr >= maxk)
      {
	kincr = maxk-k;
	k = maxk;
	mp_nb =  mp_maxnb ;
	// mpz_set(mp_nb, mp_maxnb);
	mp_db = mp_maxdb ;
	// mpz_set(mp_db, mp_maxdb);
	continue;
      }
    k += kincr;
    Integer::pow(mp_m, mp_basisprod, k);
    // mpz_pow_ui(mp_m, mp_basisprod, k);
    Integer::subin(mp_m, 1UL);
    // mpz_sub_ui(mp_m, mp_m, 1);
    Integer::divexact(mp_m, mp_m, 2UL);
    // mpz_divexact_ui(mp_m, mp_m, 2);
    Integer::sqrt(mp_nb, mp_m);
    // mpz_sqrt(mp_nb, mp_m);
    mp_db = mp_nb ;
    // mpz_set(mp_db, mp_nb);
  } while (1);

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  printf("              total lifting steps: %ld\n", k);
  tt1 = tt1/CLOCKS_PER_SEC;
  tt2 = tt2/CLOCKS_PER_SEC;
  printf("              lifting time: %f\n", tt1);
  printf("              rational reconstruction time: %f\n", tt2);
  fflush(stdout);
#endif

   return;
}

} // iml
} // LinBox

#endif // __LINBOX_algorithm_iml_non_singular_solve_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
