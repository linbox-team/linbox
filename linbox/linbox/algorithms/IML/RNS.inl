#ifndef __LINBOX_algorithm_iml_rns_INL
#define __LINBOX_algorithm_iml_rns_INL


/*
 *
 * Calling Sequence:
 *   basisExt(len, n, p, basis, cmbasis, cumprod, bdcoeff, R, RE)
 *
 * Summary:
 *   Given a representation of a matrix/vector in some RNS, extend to compute
 *   the representation in another positive integer
 *
 * Description:
 *   Let R be the representation of a matrix/vector M in a residue basis
 *   'basis', i.e., R[i] = mod(M, basis[i]) (i = 0..len-1). The function
 *   computes the representation of M in another positive integer p,
 *   RE = mod(M, p) using Garner's algorithm. 'mod' represents positive modular
 *   operation.
 *
 *   Let q be product of all elements in the RNS basis. Every entry m in M
 *   satisfies -(q-1)/2 <= m <= (q-1)/2. That is, M has both positive entries
 *   and negative entries.
 *
 *   if Pos : Let q be product of all elements in the RNS basis. Every entry m in M
 *   satisfies 0 <= m <= q-1. That is, M only contains non-negative entries.

 *
 *   To avoid repeat computations, the function takes some precomputed
 *   informations as input, which are listed below.
 *
 * Input:
 *       len: long, dimension of RNS basis
 *         n: long, length of array RE
 *         p: typename _Field::Element, modulus
 *     basis: 1-dim typename _Field::Element array length len, RNS basis
 *   cmbasis: 1-dim typename _Field::Element array length len, computed by function
 *            combBasis, inverses of special combination of RNS basis
 *   cumprod: (only in non Pos) double, computed by function cumProd,
 *            (-basis[0]*basis[1]*...*basis[len-1]) mod p
 *   bdcoeff: 1-dim typename _Field::Element array length len, computed by function repBound
 *         R: 1-dim Double array length n*len, representation of a len x n
 *            matrix, R[i]=mod(M, basis[i]) (i=0..len-1)
 *
 * Output:
 *   RE: 1-dim Double array length n, RE = mod(M, p), the space of RE
 *       should be allocated before calling the function.
 *
 * Precondition:
 *   t <= 2^53-1, where t is the maximal intermidiate value arised in this
 *   function,
 *   t = max(2*(basis[len-1]-1)^2, (p-1)^2+basis[len-1]-1)

 * Precondition if Positive:
 *   t <= 2^53-1, where t is the maximal intermidiate value arised in this
 *   function, t = max(2*(basis[len-1]-1)^2, (p-1)^2+basis[len-1]-1)

 *
 */


template<class Container, class _Field>
void
basisExtPos (//const long                                    &len
	  // , const long  &n
	  // , const typename _Container::Field::Element   &p
	    _Field & Fp
	  // , const BlasVector<Container::Field>          &basis
	  // , const BlasVector<Container::Field>          &cmbasis
	  // , const typename _Container::Field::Element   &cumprod
	  // , const BlasVector<Container::Field>          &bdcoeff
	  , std::vector<Container>                      &R
	  , Container                                   &RE)
{
	bool pos = FieldTraits<_Field>::Rep == Positive;
	typedef typename Container::Field   Field;
	typedef typename Field::Element   Element;
	typedef typename BlasVector<Field>  Vect;

	long i, j;
	Element temp;
	Element p = Fp.characteristic();
	UnparametricField<Element> UF ;
	// Element q, qinv;
	// Double **U;
	// const typename _Field::Element *q, *qinv;

	const Vect & q    = basis; //!@bug no copy ?
	const Vect & qinv = cmbbasis;
	// q = basis;
	// qinv = cmbasis;

	/* if p = q[i] then just copy the corresponding column to RE */
	//! @todo order basis ?
	for (i = 0; i < len; i++) {
		if (p == q[i]) {
			Re = R[i] ;
			// cblas_dcopy(n, R[i], 1, RE, 1);
			return;
		}
	}

	const Container Z(F);
	std::vector<Container> U(len,Z);
	// U = XMALLOC(Double *, len);
	// for (i = 0; i < len; i++) { U[i] = XMALLOC(Double, n); }

	/* compute the coefficients of mix radix in positive representation by
	   inplacing modular matrix U */
	//!@bug copy method for blas things (vector/blas and sub)
	U[0].resize(R[0]);
	FFLAS::fcopy(UF,
		     U[0].getWritePointer(),1,
		     R[0].getPointer(),1);
	// cblas_dcopy(n, R[0], 1, U[0], 1);

	for (i = 1; i < len; i++) {
		Field F(q[i]);
		FFLAS::fcopy(UF,
			     U[i].getWritePointer(),1,
			     U[i-1].getPointer(),1);
		U[i].resize(R[0]);
		// cblas_dcopy(n, U[i-1], 1, U[i], 1);
		for (j = i-2; j >= 0; j--)
		{
			FFLAS::fscal(UF,n,(Element)(q[j] % q[i]),
				     U[i].getWritePointer(),1);
			// cblas_dscal(n, (double)(q[j] % q[i]), U[i], 1);
			FFLAS::faxpy(F,n,F.one,
				     U[j].getPointer(),1,U[i].getWritePointer(),1);
			// cblas_daxpy(n, 1.0, U[j], 1, U[i], 1);
			// Dmod((double)q[i], U[i], 1, n, n);
		}
		F.init(temp,(Element)qinv[i]*(Element)(q[i]-1));
		// temp = (Element)qinv[i]*(Element)(q[i]-1);
		// temp = fmod(temp, (Element)q[i]);
		FFLAS::fscal(UF,n,temp,U[i].getWritePointer(),1);
		// cblas_dscal(n, temp, U[i], 1);
		FFLAS::faxpy(F,n,qinv[i],R[i].getPointer(),1,
			     U[i].getWritePointer(),1);
		// cblas_daxpy(n, (double)qinv[i], R[i], 1, U[i], 1);
		// Dmod((double)q[i], U[i], 1, n, n);
	}

	/* compute mod(r, p) in positive representation and store into RE */
	FFLAS::fcopy(Fp,n,U[len-1].getPointer(),1,
		     RE.getWritePointer(),1);
	// cblas_dcopy(n, U[len-1], 1, RE, 1);
	// Dmod((double)p, RE, 1, n, n);
	for (i = len-2; i >= 0; i--) {
		FFLAS::fscal(UF,n,(Element)(q[i] % p),RE.getWritePointer(),1);
		// cblas_dscal(n, (double)(q[i] % p), RE, 1);
		FFLAS::faxpy(Fp,n,1,U[i].getPointer(),1,
			     RE.getWritePointer(),1);
		// cblas_daxpy(n, 1.0, U[i], 1,  RE, 1);
		// Dmod((double)p, RE, 1, n, n);
	}

	if (pos == false) {
	/* convert to symmetric representation */
	//! @todo use iterators here for i
	for (i = 0; i < n; i++)
		for (j = len-1; j >= 0; j--) {
			//!@bug does not work for every field (already balanced ?)
			if (U[j].getPointer()[i] > bdcoeff[j]) {
				RE.getPointer()[i] = fmod(RE[i]+cumprod, p);
				break;
			}
			else if (U[j].getPointer()[i] < bdcoeff[j]) {
				break;
			}
		}
	}
	// for (i = 0; i < len; i++) { XFREE(U[i]); } { XFREE(U); }

	return;
}


/*
 *
 * Calling Sequence:
 *   basisProd(len, basis, mp_prod)
 *
 * Summary:
 *   Compute the product of elements of a RNS basis
 *
 * Description:
 *   Let a RNS basis be 'basis'. The function computes the product of its
 *   elements basis[0]*basis[1]*...*basis[len-1].
 *
 * Input:
 *     len: long, dimension of RNS basis
 *   basis: 1-dim typename _Field::Element array length len, RNS basis
 *
 * Output:
 *   mp_prod: mpz_t, product of elements in 'basis'
 *
 */

void
basisProd (Integer &  mp_prod)
{
  long i;

  // mpz_set_ui(mp_prod, basis[0]);
  mp_prod = basis[0];
  for (i = 1; i < len; i++) { Integer::mulin(mp_prod,basis[i]); }
}



/*
 *
 * Calling Sequence:
 *   ChineseRemainder(len, mp_prod, basis, cmbasis, bdcoeff, Ac, mp_Ac)
 *
 * Summary:
 *   Given a representation of an integer in some RNS, use Chinese Remainder
 *   Algorithm to reconstruct the integer
 *
 * Description:
 *   Let A be an integer, and Ac contains the representation of A in a RNS
 *   basis 'basis', i.e. Ac[i] = mod(A, basis[i]), (i = 0..len). Here 'mod'
 *   is in positive representation. The function reconstructs the integer A
 *   given the RNS basis 'basis' and Ac.
 *
 *   To avoid repeat computations, the function takes some precomputed
 *   informations as input, which are listed below.
 *
 * Input:
 *       len: long, dimension of RNS basis
 *   mp_prod: mpz_t, computed by function basisProd, product of RNS basis
 *     basis: 1-dim typename _Field::Element array length len, RNS basis
 *   cmbasis: 1-dim typename _Field::Element array length len, computed by function
 *            combBasis, inverses of special combination of RNS basis
 *   bdcoeff: 1-dim typename _Field::Element array length len, computed by function repBound
 *        Ac: 1-dim Double array length n, representation of A in RNS
 *
 * Output:
 *   mp_Ac: mpz_t, reconstructed integer A
 *
 * Precondition:
 *   Let q be product of all elements in the RNS basis. Then A must satisfy
 *   -(q-1)/2 <= A <= (q-1)/2.
 *
 *
 *
 *
 * Calling Sequence:
 *   ChineseRemainderPos(len, basis, cmbasis, Ac, mp_Ac)
 *
 * Summary:
 *   Given a representation of a non-negative integer in some RNS, use Chinese
 *   Remainder Algorithm to reconstruct the integer
 *
 * Description:
 *   Let A be a non-negative integer, and Ac contains the representation of A
 *   in a RNS basis 'basis', i.e. Ac[i] = mod(A, basis[i]), (i = 0..len).
 *   Here 'mod' is in positive representation. The function reconstructs the
 *   integer A given the RNS basis 'basis' and Ac.
 *
 *   To avoid repeat computations, the function takes some precomputed
 *   informations as input, which are listed below.
 *
 * Input:
 *       len: long, dimension of RNS basis
 *     basis: 1-dim typename _Field::Element array length len, RNS basis
 *   cmbasis: 1-dim typename _Field::Element array length len, computed by function
 *            combBasis, inverses of special combination of RNS basis
 *        Ac: 1-dim Double array length n, representation of A in RNS
 *
 * Output:
 *   mp_Ac: mpz_t, reconstructed integer A
 *
 * Precondition:
 *   Let q be product of all elements in the RNS basis. Then A must satisfy
 *   0 <= A <= q-1.
 *
 */

void
ChineseRemainder (const Integer mp_prod,
		  BlasVector<Field> &Ac, Integer  mp_Ac)
{
	long i, j;
	Element temp, tempq, tempqinv;
	bool pos = FieldTraits<_Field>::Rep == Positive;

	BlasVector<Field> U(Ac.field(),len);
	// Double *U;
	// U = XMALLOC(Double, len);

	/* compute the coefficients of mix radix in positive representation by
	   inplacing modular matrix U */
	U[0] = Ac[0];
	for (i = 1; i < len; i++) {
		U[i] = U[i-1];
		tempq = (Element)basis[i];
		Field Fq(tempq);
		tempqinv = (Element)cmbasis[i];
		for (j = i-2; j >= 0; j--) {
			Fq.mulin(U[i],basis[j]);
			Fq.addin(U[i],U[j]);
			// U[i] = U[j] + U[i]*fmod((Element)basis[j], tempq);
			// U[i] = fmod(U[i], tempq);
		}
		F.init(temp,tempqinv*(Element)(basis[i]-1));
		// temp = fmod(tempqinv*(Element)(basis[i]-1), tempq);
		F.init(U[i],tempqinv*Ac[i]+temp*U[i]);
		// U[i] = fmod(tempqinv*Ac[i]+temp*U[i], tempq);
	}
	/* compute Ac in positive representation */
	mp_Ac = U[len-1];
	// mpz_set_d(mp_Ac, U[len-1]);
	for (j = len-2; j >= 0; j--)
	{
		Integer::mulin(mp_Ac,basis[j]);
		Integer::addin(mp_Ac, (Element)U[j]);
		// mpz_mul_ui(mp_Ac, mp_Ac, basis[j]);
		// mpz_add_ui(mp_Ac, mp_Ac, (Element)U[j]);
	}
	/* transfer from positive representation to symmetric representation */
	if (pos == false) {
		for (j = len-1; j >= 0; j--)
		{
			if (U[j] > bdcoeff[j])
			{
				Integer::subin(mp_Ac,mp_prod);
				// mpz_sub(mp_Ac, mp_Ac, mp_prod);
				break;
			}
			else if (U[j] < bdcoeff[j]) { break; }
		}
	}
	// XFREE(U);

	return;
}





/*
 *
 * Calling Sequence:
 *   cmbasis <-- combBasis(basislen, basis)
 *
 * Summary:
 *   Compute the special combination of a RNS basis
 *
 * Description:
 *   Let 'basis' be RNS basis. The function computes an array cmbasis
 *   satisfying
 *   cmbasis[0] = 0, cmbasis[i] = mod(1/(basis[0]*...*basis[i-1]), basis[i])
 *                   (i = 1..basislen-1)
 *
 * Input:
 *   basislen: long, dimension of RNS basis
 *      basis: 1-dim typename _Field::Element array length basislen, RNS basis
 *
 * Return:
 *   cmbasis: 1-dim typename _Field::Element array length basislen, shown as above
 *
 */

typename _Field::Element *
combBasis ()
{
  long i, j;
  Element dtemp;
  // mpz_t mp_prod, mp_q;
  // typename _Field::Element *cmbasis;

  // cmbasis = XMALLOC(typename _Field::Element, basislen);
  cmbasis.resize(basislen);
  cmbasis[0] = 0;
  // mpz_init(mp_prod);
  // mpz_init(mp_q);
  for (i = 1; i < basislen; i++)
    {
	    Field Fi(basis[i]);
	    F.init(dtemp,basis[0]);
      // dtemp = fmod((double)basis[0], (double)basis[i]);
      for (j = 1; j <= i-1; j++){
	dtemp = fmod(dtemp*(double)basis[j], (double)basis[i]);
	F.mulin(dtemp,basis[j]);
      }
      // Integer mp_q = basis[i];
      // mpz_set_ui(mp_q, basis[i]);
      // Integer mp_prod = dtemp ;
      // mpz_set_d(mp_prod, dtemp);
      // mpz_invert(mp_prod, mp_prod, mp_q);
      // cmbasis[i] = mpz_get_ui(mp_prod);
      F.inv(cmbasis[i],dtemp,basis[i]);
    }
  // mpz_clear(mp_prod);
  // mpz_clear(mp_q);

  return cmbasis;
}



/*
 *
 * Calling Sequence:
 *   cumprod <-- cumProd(basislen, basis, extbasislen, extbasis)
 *
 * Summary:
 *   Compute the representation of the combination of elements of one RNS basis
 *   in another RNS basis
 *
 * Description:
 *   Let 'basis' be one RNS basis with dimension basislen, and 'extbasis' be
 *   another RNS basis with dimension extbasislen. The function computes an
 *   array cumprod length extbasislen satisfying
 *   cumprod[i] = modp(-basis[0]*...*basis[basislen-1], extbasis[i]),
 *   i = 0..extbasislen-1
 *
 * Input:
 *      basislen: long, dimension of RNS basis 'basis'
 *         basis: 1-dim typename _Field::Element array length basislen, one RNS basis
 *   extbasislen: long, dimension of RNS basis 'extbasis'
 *      extbasis: 1-dim typename _Field::Element array length basislen, another RNS basis
 *
 * Return:
 *   cumprod: 1-dim double array length extbasislen, shown above
 *
 */

cumProd (BlasVector<Field> &cumprod,
	 const BlasVector<Field> &extbasis)
{
  long i, j;
  Element dtemp, dextbasis;
  // double *cumprod;

  // cumprod = XMALLOC(double, extbasislen);
  size_t extbasislen = extbasis.size();
  cumprod.resize(extbasislen);
  for (i = 0; i < extbasislen; i++) {
      Field F(extbasis[i]);
      // dextbasis = (Element)extbasis[i];
      F.init(cumprod[i],basis[0]);
      // cumprod[i] = fmod((double)basis[0], dextbasis);
      for (j = 1; j < basis.size(); j++) {
	      F.init(dtemp,basis[j]);
	  // dtemp = fmod((double)basis[j], dextbasis);
	  F.mulin(cumprod[i],dtemp); //!@bug 2 in one ?
	  // cumprod[i] = fmod(cumprod[i]*dtemp, dextbasis);
	}
      // cumprod[i] = dextbasis-cumprod[i];
      F.negin(cumprod[i]);
    }

  return ;
}


/*
 *
 * Calling Sequence:
 *   basiscmb <-- findRNS(RNS_bound, mp_maxInter, len)
 *
 * Summary:
 *   Find a RNS basis and its special combination
 *
 * Description:
 *   Given RNS_bound, the upper bound of the RNS basis, and mp_maxInter, the
 *   function finds a best RNS basis and a combination of that basis.
 *
 *   The RNS basis 'basis' has the property:
 *   - its elements are all primes
 *   - basis[0] is the largest prime among all the primes at most RNS_bound
 *   - basis[i+1] is the next prime smaller than basis[i] (i = 0..len-2)
 *   - basis[0]*basis[1]*...*basis[len-1] >= mp_maxInter
 *
 *   After finding 'basis', the functions also computes the combination of
 *   'basis' as the operations in function combBasis.
 *
 * Input:
 *     RNS_bound: typename _Field::Element, the upper bound of the RNS basis
 *   mp_maxInter: mpz_t, the lower bound for the product of elements of basis
 *
 * Return:
 *   basiscmb: 2-dim typename _Field::Element array, dimension 2 x len, where
 *           - basiscmb[0] represents the RNS basis
 *           - basiscmb[1] represents the special combination of basis
 *
 * Output:
 *   len: pointer to a long int, storing the dimension of the computed
 *        RNS basis
 *
 */

void
findRNS ( const Element RNS_bound, const Integer mp_maxInter,
	  long &length, BlasVector<Field> & RNSbasis, BlasVector<Field> & RNScomb)
{
  long i, j, len=0;
  double prod;
  // mpz_t mp_l, mp_prod, mp_q;
  // typename _Field::Element **qqinv;

  Integer mp_prod = 1;
  // mpz_init_set_ui(mp_prod, 1);
  Integer mp_l = RNS_bound ;
  // mpz_init_set_ui(mp_l, RNS_bound);
  // qqinv = XMALLOC(typename _Field::Element *, 2);
  // RNSbasis = NULL;
  while (mp_maxInter > mp_prod)
  // while (mpz_cmp(mp_maxInter, mp_prod) > 0)
    {
      ++len;
      RNSbasis.resize(len);
      // RNSbasis = XREALLOC(typename _Field::Element, RNSbasis, len);
      while (Integer::probab_prime_p(mp_l, 10) == 0) {
	      Integer::subin(mp_l,1); //! @bug why remove 1 ?
	      // mpz_sub_ui(mp_l, mp_l, 1);
      }
      RNSbasis[len-1] = (mp_l);
      Integer::subin(mp_l,1);
      // mpz_sub_ui(mp_l, mp_l, 1);
      Integer::mulin(mp_prod,RNSbasis[len-1]);
      // mpz_mul_ui(mp_prod, mp_prod, RNSbasis[len-1]);
    }
  // mpz_clear(mp_prod);
  // mpz_clear(mp_l);
  // RNScomb = XMALLOC(typename _Field::Element, len);
  RNScomb.resize(len);
  RNScomb[0] = 0;
  // mpz_init(mp_prod);
  // mpz_init(mp_q);
  for (i = 1; i < len; i++) {
      Element prod = (double)(RNSbasis[0] % RNSbasis[i]);
      for (j = 1; j <= i-1; j++) {
	// prod = fmod(prod*(double)RNSbasis[j], (double)RNSbasis[i]);
	F.mulin(prod,RNSbasis[j]);
      }
      // mpz_set_ui(mp_q, RNSbasis[i]);
      // mpz_set_d(mp_prod, prod);
      // mpz_invert(mp_prod, mp_prod, mp_q);
      // RNScomb[i] = mpz_get_ui(mp_prod);
      F.invert(RNScomb[i],prod;)
    }
  // mpz_clear(mp_prod);
  // mpz_clear(mp_q);
  length = len;

  // return qqinv;
}



/*
 *
 * Calling Sequence:
 *   maxInter(mp_prod, mp_alpha, n, mp_b)
 *
 * Summary:
 *   Compute the maximum interval of positive and negative results of a
 *   matrix-matrix or matrix-vector product
 *
 * Description:
 *   Let mp_alpha be the maximum magnitude of a m x n matrix A, mp_prod-1 be
 *   the maximum magnitude of a n x k matrix C. The function computes the
 *   maximum interval of positive and negative entries of A.C. That is, the
 *   function computes mp_b satisfying
 *   (mp_b-1)/2 = n*mp_alpha*(mp_prod-1)
 *
 * Input:
 *    mp_prod: mpz_t, mp_prod-1 be the maximum magnitude of matrix C
 *   mp_alpha: mpz_t, maximum magnitude of matrix A
 *          n: long, column dimension of A
 *
 * Output:
 *   mp_b: mpz_t, shown above
 *
 */

void
maxInter (const mpz_t mp_prod, const mpz_t mp_alpha, const long n, mpz_t mp_b)
{
  Integer mp_temp;
  // mpz_t mp_temp;

  // mpz_init(mp_temp);
  Integer::sub(mp_temp,mp_prod,1);
  // mpz_sub_ui(mp_temp, mp_prod, 1);
  Integer mp_b = mp_alpha ;
  // mpz_set(mp_b, mp_alpha);
  Integer::mulin(mp_b,n); //! @bug 2n here ?
  // mpz_mul_ui(mp_b, mp_b, n);
  Integer::mulin(mp_b,mp_temp);
  // mpz_mul(mp_b, mp_b, mp_temp);
  Integer::mulin(mp_b,2);
  // mpz_mul_ui(mp_b, mp_b, 2);
  Integer::addin(mp_b,1);
  // mpz_add_ui(mp_b, mp_b, 1);
  // mpz_clear(mp_temp);
}



/*
 *
 * Calling Sequence:
 *   maxExtInter(mp_alpha, n, mp_b)
 *
 * Summary:
 *   Compute the maximum interval of positive and negative results for
 *   lifting
 *
 * Description:
 *   Let mp_alpha be the maximum magnitude of a m x n matrix A,
 *   The function computes the mp_b satisfying
 *   (mp_b-1)/2 = n*mp_alpha+1
 *
 * Input:
 *   mp_alpha: mpz_t, maximum magnitude of matrix A
 *          n: long, column dimension of A
 *
 * Output:
 *   mp_b: mpz_t, shown above
 *
 */

void
maxExtInter (const Integer & mp_alpha, const long n, Integer & mp_b)
{
	Integer mp_b = 1;
  // mpz_set_ui(mp_b, 1);
	Integer::axpyin(mp_b,mp_alpha,n);
  // mpz_addmul_ui(mp_b, mp_alpha, n);
	Integer::mulin(mp_b,2);
  // mpz_mul_ui(mp_b, mp_b, 2);
  Integer::addin(mp_b,1);
  // mpz_add_ui(mp_b, mp_b, 1);
}



/*
 *
 * Calling Sequence:
 *   bdcoeff <-- repBound(len, basis, cmbasis)
 *
 * Summary:
 *   Compute the mix radix coefficients of a special integer in a RNS basis
 *
 * Description:
 *   Given a RNS basis, suppose the product of elements in the basis be q,
 *   then this RNS basis is able to represent integers lying in
 *   [-(q-1)/2, (q-1)/2] and [0, q-1] respectively with symmetric
 *   representation and positive representation. To transfer the result from
 *   positive representation to symmetric representation, the function
 *   computes the mix radix coefficients of the boundary value (q-1)/2 in the
 *   positive representation.
 *
 *   Let RNS basis be P. The function computes coefficient array U, such that
 * (q-1)/2 = U[0] + U[1]*P[0] + U[2]*P[0]*P[1] +...+ U[len-1]*P[0]*...*P[len-2]
 *
 * Input:
 *       len: long, dimension of RNS basis
 *     basis: 1-dim typename _Field::Element array length len, RNS basis
 *   cmbasis: 1-dim typename _Field::Element array length len, computed by function
 *            combBasis, inverses of special combination of RNS basis
 *
 * Output:
 *   bdcoeff: 1-dim typename _Field::Element array length len, the coefficient array U above
 *
 */

// typename _Field::Element *
void
repBound ()
{
  long i, j;
  Element dtemp;
  Integer mp_bd, mp_prod;
  // typename _Field::Element *bdcoeff;

  const Vect & q    = basis;
  const Vect & qinv = cmbasis;

  /* set the bound of transformation from positive to negative */
  // mpz_init(mp_prod);
  basisProd(q, mp_prod);
  // mpz_init(mp_bd);
  Integer::sub(mp_bd,mp_prod,1);
  // mpz_sub_ui(mp_bd, mp_prod, 1);
  Integer::divein(mp_bd,2);
  // mpz_divexact_ui(mp_bd, mp_bd, 2);

  /* compute the coeffcients of bound of mix radix and store in bdcoeff */
  // bdcoeff = XMALLOC(typename _Field::Element, len);
  // bdcoeff[0] = mpz_fdiv_ui(mp_bd, q[0]);
  Integer::fdivv(bdcoeff[0],mp_bd,q[0]);
  for (i = 1; i < len; i++) {
      dtemp = (Element)bdcoeff[i-1];
      Field Fj(q[i]);
      for (j = i-2; j >= 0; j--) {
	  F.mulin(dtemp,q[j]);
	  // dtemp = fmod(dtemp*q[j], (double)q[i]);
	  F.addin(dtemp,bdcoeff[j]);
	  // dtemp = fmod(dtemp+bdcoeff[j], (double)q[i]);
	}
      Integer::dive(bdcoeff[i],mp_bd,q[i]);
      // bdcoeff[i] = mpz_fdiv_ui(mp_bd, q[i]);
      F.init(dtemp,bdcoeff[i]-dtemp);
      // dtemp = fmod((double)bdcoeff[i]-dtemp, (double)q[i]);
      //! @bug can it be <0 ?
      // if (dtemp < 0) {
	      // dtemp = q[i]+dtemp;
      // }
      F.init(bdcoeff[i],dtemp*qinv[i]);
      // bdcoeff[i] = (typename _Field::Element)fmod(dtemp*qinv[i], (double)q[i]);
    }
  // mpz_clear(mp_bd);
  // mpz_clear(mp_prod);

  // return bdcoeff;
}


/*
 * Calling Sequence:
 *   bd <-- RNSbound(n)
 *
 * Summary:
 *   Compute the upper bound of a RNS basis
 *
 * Description:
 *   Given a m x n mod p matrix A, and a n x k mod p matrix B, the maximum
 *   magnitude of A.B is n*(p-1)^2. To use BLAS, it is needed that
 *   n*(p-1)^2 <= 2^53-1 to make the result of product correct.
 *
 *   The function computes an integer bd, such that
 *      n*(bd-1)^2 <= 2^53-1 and n*((bd+1)-1)^2 > 2^53-1
 *
 * Input:
 *   n: long, column dimension of matrix A
 *
 * Output:
 *   bd: typename _Field::Element, shown above
 *
 */

Element
RNSbound (const long n)
{
  Element bd;
  mpz_t mp_n, mp_d, mp_q;

  // mpz_init(mp_n);
  // mpz_init_set_ui(mp_d, n);
  // mpz_init(mp_q);
  Integer::pow(mp_n,2,53);
  // mpz_ui_pow_ui(mp_n, 2, 53);
  // mpz_sub_ui(mp_n, mp_n, 1);
  Integer::subin(mp_n,1);
  Integer::divq(mpa,mp_n,n);
  // mpz_fdiv_q(mp_q, mp_n, mp_d);
  Integer::sqrtin(mp_q)
  // mpz_sqrt(mp_q, mp_q);
  // bd = mpz_get_ui(mp_q)+1;
  bd = (Element)mp_q + 1 ;
  // mpz_clear(mp_n);
  // mpz_clear(mp_d);
  // mpz_clear(mp_q);

  return bd;
}

#endif // __LINBOX_algorithm_iml_rns_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
