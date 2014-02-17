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

template<class Matrix>
void
minSolnoncompress (const long certflag, const long redflag, const long n, \
		      const long k,
		       mpz_t *mp_Bb,
		       Matrix & AA,
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
  nonsingSolv(RightSolu, n, k+1, AA , mp_Bb, \
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

#if HAVE_VERBOSE_MODE && HAVE_TIME_H
  tt = (double)(clock()-ti)/CLOCKS_PER_SEC;
  printf("          Matix M, N Construction Time: %f\n", tt);
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
	nonsingSolvRNSMM(LeftSolu, n, 1, basislen, basis, ARNS, mp_Cert,  mp_NZ, mp_DZ);
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


#endif // __LINBOX_algorithm_iml_certified_solve_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


