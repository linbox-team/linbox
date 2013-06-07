#ifndef __LINBOX_algorithm_iml_certified_solve_INL
#define __LINBOX_algorithm_iml_certified_solve_INL

#define NULLSPACE_COLUMN 10

/*
 * Calling Sequence:
 *   compressBoundMP(n, m, Pt, mp_A, mp_bd)
 *
 * Summary:
 *   Compute the maximum magnitude of a compressed Integer submatrix
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
 *   mp_A: 1-dim Integer array length n1*m, n1 x m matrix A
 *
 * Output:
 *   mp_bd: mpz_t, maximum magnitude of entries in B.P
 *
 */

template<class PID>
void
compressBoundMP ( const std::vector<size_t> & Pt,
		  BlasMatrix<Ring>          & mp_A,
		  PID                       & R,
		  typename PID::Element     & mp_bd)
{
	// long i, j;
	typename Ring::Element temp;
	typename PID::Element mp_temp ;
	// Integer mp_temp;

	size_t  m = mp_A.coldim();
	size_t n = mp_A.rowdim();
	mp_bd = 0UL;
	// mpz_set_ui(mp_bd, 0);
	for (size_t i = 0; i < n; i++) {
		// mpz_set_ui(mp_temp, 0);
		mp_temp = 0UL;
		for (size_t j = 0; j < m; j++) {
			//! @bug I want to use either abs, labs or Integer::abs here
			//! so we need LinBox::abs
			R.abs(temp,mp_A[Pt[i]*m+j]);
			R.addin(mp_temp, mp_temp1);
		}
		if (R.absCompare(mp_bd, mp_temp) < 0) {
			mp_bd = mp_temp ;
			// mpz_set(mp_bd, mp_temp);
		}
	}
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
migcdex (const Integer &N
	 , const Integer &a
	 , std::vector<PID_integer> &b
	 // , const long n
	 , std::vector<size_t> & c)
{
  long i, j;
  Integer gAN, gAbN, A;

  // mpz_init(gAN);
  // mpz_init(gAbN);
  // mpz_init_set(A, a);
  A = a;
  Integer::gcd(gAbN, a, N);
  // mpz_gcd(gAbN, a, N);
  for (i = 0; i < n; i++)
    {
	    Integer::gcd(gAbN, gAbN, b[i]);
      // mpz_gcd(gAbN, gAbN, b[i]);
      for (j = 0; ; j++) {
	      Integer::gcd(gAN, A, N);
	  // mpz_gcd(gAN, A, N);
	  if (!Integer::compare(gAbN, gAN)) { break; }
	  Integer::addin(A,b[i]);
	  // mpz_add(A, A, b[i]);
	}
      c[i] = j;
    }

  // { mpz_clear(gAN); mpz_clear(gAbN); mpz_clear(A); }
  return;
}

#if 0
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


#endif

/* return min( upper bound of signed long, 2^52-1 ) */
template<class T>
T
LongRNSbound(void)
{
#if 0
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
#endif
  return std::min(std::numeric_limits<T>::max(),4503599627370495)
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
kernelBasis (const long t
	     ,  BlasMatrix<PID_integer> & mp_M
	     , BlasMatrix<PID_integer> & mp_N)
{
	size_t k = mp_N.coldim();
	size_t n = mp_N.rowdim()-k;
	long i, j, l;
	Integer mp_temp ;
	// mpz_t mp_temp;

	// mpz_init(mp_temp);
	for (i = 0; i < k; i++) {
		for (j = 0; j < i; j++)
			// XXX do blas op
			/* column operation N[1..-1, i] = N[1..-1, i]-N[1..-1, j]*M[j, i] */
			for (l = 0; l < n+k; l++) {
				Integer::mul(mp_temp,mp_N.getEntry(l,j),mp_M.getEntry(j,i));
				// mpz_mul(mp_temp, mp_N[l*k+j], mp_M[j*(k+t)+i]);
				Integer::subin(mp_N.getEntry(l,i),mp_temp);
				// mpz_sub(mp_N[l*k+i], mp_N[l*k+i], mp_temp);
			}

		/* column operation N[1..-1, i] = N[1..-1, i]/M[i, i] */
		for (l = 0; l < n+k; l++)
		{
			Integer::divexact(mp_N.refEntry(l,i), mp_N.refEntry(l,i), mp_M.refEntry(i,i));
			// mpz_divexact(mp_N[l*k+i], mp_N[l*k+i], mp_M[i*(k+t)+i]);
		}
	}

	// mpz_clear(mp_temp);
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

template<class Field>
int
certVerify (const RNS<Field> & basis
	    , RNS_rep<Field> & ARNS
	     // Double **ARNS
	     , Integer & mp_DZ,
	    std::vector<PID_integer> & mp_NZ)
{
	typedef typename Field::Element FiniteField ;
  long i, j, l;
  // FiniteField *bdcoeff, *cmbasis;
  double temp;
  Integer mp_AZ, mp_temp, mp_AZ1;
  // Double **U;
  MixRadix<Field> U(basislen,n);

  /* allocate matrix U to store mix radix coefficients of one row/column
     of matrix A */
  // U = XMALLOC(Double *, basislen);
  // for (i = 0; i < basislen; i++)
    // U[i] = XMALLOC(Double, n);
  // cmbasis = combBasis(basislen, basis);
  cmbasis
  // bdcoeff = repBound(basislen, basis, cmbasis);
  bdcoeff
  // mpz_init(mp_AZ);
  // mpz_init(mp_AZ1);
  // mpz_init(mp_temp);
  for (i = 0; i < m; i++) {
      /* apply Garner's algorithm to compute U in positive representation */
      cblas_dcopy(n, ARNS[0]+i, m, U[0], 1);
      for (j = 1; j < basislen; j++) {
	  cblas_dcopy(n, U[j-1], 1, U[j], 1);
	  for (l = j-2; l >= 0; l--) {
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
	for (l = basislen-1; l >= 0; l--) {
	    if (U[l][j] > bdcoeff[l]) {
		U[basislen-1][j] = U[basislen-1][j]-basis[basislen-1];
		break;
	      }
	    else if (U[l][j] < bdcoeff[l]) { break; }
	  }
      /* do vector-vector product and sum up the product */
      mpz_set_ui(mp_AZ, 0);
      for (j = 0; j < n; j++) {
	  mpz_mul_si(mp_temp, mp_NZ[j], U[basislen-1][j]);
	  mpz_add(mp_AZ, mp_AZ, mp_temp);
	}
      for (j = basislen-2; j >= 0; j--) {

	  mpz_mul_ui(mp_AZ, mp_AZ, basis[j]);
	  mpz_set_ui(mp_AZ1, 0);
	  for (l = 0; l < n; l++) {
	      mpz_mul_ui(mp_temp, mp_NZ[l], U[j][l]);
	      mpz_add(mp_AZ1, mp_AZ1, mp_temp);
	    }
	  mpz_add(mp_AZ, mp_AZ, mp_AZ1);
	}
      /* check whether mp_DZ|mp_AZ */
      if (!mpz_divisible_p(mp_AZ, mp_DZ)) {
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


#endif // __LINBOX_algorithm_iml_certified_solve_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


