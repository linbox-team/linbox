#ifndef __LINBOX_algorithm_iml_rns_INL
#define __LINBOX_algorithm_iml_rns_INL

namespace LinBox { namespace iml {

/*
 *
 * Calling Sequence:
 *   basisExt(len, n, p, RNSbasis, RNScombi, cumprod, bdcoeff, R, RE)
 *
 * Summary:
 *   Given a representation of a matrix/vector in some RNS, extend to compute
 *   the representation in another positive integer
 *
 * Description:
 *   Let R be the representation of a matrix/vector M in a residue RNSbasis
 *   'basis', i.e., R[i] = mod(M, RNSbasis[i]) (i = 0..len-1). The function
 *   computes the representation of M in another positive integer p,
 *   RE = mod(M, p) using Garner's algorithm. 'mod' represents positive modular
 *   operation.
 *
 *   Let q be product of all elements in the RNS RNSbasis. Every entry m in M
 *   satisfies -(q-1)/2 <= m <= (q-1)/2. That is, M has both positive entries
 *   and negative entries.
 *
 *   if Pos : Let q be product of all elements in the RNS RNSbasis. Every entry m in M
 *   satisfies 0 <= m <= q-1. That is, M only contains non-negative entries.

 *
 *   To avoid repeat computations, the function takes some precomputed
 *   informations as input, which are listed below.
 *
 * Input:
 *       RNSbasislen: long, dimension of RNS RNSbasis
 *         n: long, length of array RE
 *         p: typename _Field::Element, modulus
 *     RNSbasis: 1-dim typename _Field::Element array length basislen, RNS RNSbasis
 *   RNScombi: 1-dim typename _Field::Element array length basislen, computed by function
 *            combBasis, inverses of special combination of RNS RNSbasis
 *   cumprod: (only in non Pos) double, computed by function cumProd,
 *            (-RNSbasis[0]*RNSbasis[1]*...*RNSbasis[basislen-1]) mod p
 *   bdcoeff: 1-dim typename _Field::Element array length basislen, computed by function repBound
 *         R: 1-dim Double array length n*basislen, representation of a basislen x n
 *            matrix, R[i]=mod(M, RNSbasis[i]) (i=0..basislen-1)
 *
 * Output:
 *   RE: 1-dim Double array length n, RE = mod(M, p), the space of RE
 *       should be allocated before calling the function.
 *
 * Precondition:
 *   t <= 2^53-1, where t is the maximal intermidiate value arised in this
 *   function,
 *   t = max(2*(RNSbasis[basislen-1]-1)^2, (p-1)^2+RNSbasis[basislen-1]-1)

 * Precondition if Positive:
 *   t <= 2^53-1, where t is the maximal intermidiate value arised in this
 *   function, t = max(2*(RNSbasis[basislen-1]-1)^2, (p-1)^2+RNSbasis[basislen-1]-1)

 *
 */


template<class Container, class FiniteField>
void
RNS<Container,FiniteField>::basisExt ( FiniteField & Fp
				       , std::vector<Container>                      &R
				       , Container                                   &RE
				       , bool                                         pos )
{
	// bool pos = FieldTraits<_Field>::Rep == Positive;

	long i, j;
	long n = RE.size();
	Element temp;
	Element p = Fp.characteristic();
	// Element q, qinv;
	// Double **U;
	// const typename _Field::Element *q, *qinv;

	const ModVect & q    = RNSbasis; //!@bug no copy ?
	const ModVect & qinv = RNScombi;
	// q = RNSbasis;
	// qinv = RNScombi;

	/* if p = q[i] then just copy the corresponding column to RE */
	//! @todo order RNSbasis ?
	for (i = 0; i < basislen ; i++) {
		if (p == q[i]) {
		RE = R[i] ;
			// cblas_dcopy(n, R[i], 1, RE, 1);
			return;
		}
	}
	Field F = RE.field();

	const Container Z(F);
	std::vector<Container> U(basislen,Z);
	// U = XMALLOC(Double *, basislen);
	// for (i = 0; i < basislen; i++) { U[i] = XMALLOC(Double, n); }

	/* compute the coefficients of mix radix in positive representation by
	   inplacing modular matrix U */
	//!@bug copy method for blas things (vector/blas and sub)
	U[0].resize(R[0]);
	FFLAS::fcopy(unF,
		     U[0].getWritePointer(),1,
		     R[0].getPointer(),1);
	// cblas_dcopy(n, R[0], 1, U[0], 1);

	for (i = 1; i < basislen; i++) {
		FiniteField Fq(q[i]);
		FFLAS::fcopy(unF,
			     U[i].getWritePointer(),1,
			     U[i-1].getPointer(),1);
		U[i].resize(R[0]);
		// cblas_dcopy(n, U[i-1], 1, U[i], 1);
		for (j = i-2; j >= 0; j--)
		{
			FFLAS::fscal(unF,n,(Element)(q[j] % q[i]),
				     U[i].getWritePointer(),1);
			// cblas_dscal(n, (double)(q[j] % q[i]), U[i], 1);
			FFLAS::faxpy(Fq,n,Fq.one,
				     U[j].getPointer(),1,U[i].getWritePointer(),1);
			// cblas_daxpy(n, 1.0, U[j], 1, U[i], 1);
			// Dmod((double)q[i], U[i], 1, n, n);
		}
		Fq.init(temp,(Element)qinv[i]*(Element)(q[i]-1));
		// temp = (Element)qinv[i]*(Element)(q[i]-1);
		// temp = fmod(temp, (Element)q[i]);
		FFLAS::fscal(unF,n,temp,U[i].getWritePointer(),1);
		// cblas_dscal(n, temp, U[i], 1);
		FFLAS::faxpy(Fq,n,qinv[i],R[i].getPointer(),1,
			     U[i].getWritePointer(),1);
		// cblas_daxpy(n, (double)qinv[i], R[i], 1, U[i], 1);
		// Dmod((double)q[i], U[i], 1, n, n);
	}

	/* compute mod(r, p) in positive representation and store into RE */
	FFLAS::fcopy(Fp,n,U[basislen-1].getPointer(),1,
		     RE.getWritePointer(),1);
	// cblas_dcopy(n, U[basislen-1], 1, RE, 1);
	// Dmod((double)p, RE, 1, n, n);
	for (i = basislen-2; i >= 0; i--) {
		FFLAS::fscal(unF,n,(Element)(q[i] % p),RE.getWritePointer(),1);
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
		for (j = basislen-1; j >= 0; j--) {
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
	// for (i = 0; i < basislen; i++) { XFREE(U[i]); } { XFREE(U); }

	return;
}


/*
 *
 * Calling Sequence:
 *   ChineseRemainder(basislen, mp_prod, RNSbasis, RNScombi, bdcoeff, Ac, mp_Ac)
 *
 * Summary:
 *   Given a representation of an integer in some RNS, use Chinese Remainder
 *   Algorithm to reconstruct the integer
 *
 * Description:
 *   Let A be an integer, and Ac contains the representation of A in a RNS
 *   RNSbasis 'basis', i.e. Ac[i] = mod(A, RNSbasis[i]), (i = 0..basislen). Here 'mod'
 *   is in positive representation. The function reconstructs the integer A
 *   given the RNS RNSbasis 'basis' and Ac.
 *
 *   To avoid repeat computations, the function takes some precomputed
 *   informations as input, which are listed below.
 *
 * Input:
 *       basislen: long, dimension of RNS RNSbasis
 *   mp_prod: mpz_t, computed by function basisProd, product of RNS RNSbasis
 *     RNSbasis: 1-dim typename _Field::Element array length basislen, RNS RNSbasis
 *   RNScombi: 1-dim typename _Field::Element array length basislen, computed by function
 *            combBasis, inverses of special combination of RNS RNSbasis
 *   bdcoeff: 1-dim typename _Field::Element array length basilen, computed by function repBound
 *        Ac: 1-dim Double array length n, representation of A in RNS
 *
 * Output:
 *   mp_Ac: mpz_t, reconstructed integer A
 *
 * Precondition:
 *   Let q be product of all elements in the RNS RNSbasis. Then A must satisfy
 *   -(q-1)/2 <= A <= (q-1)/2.
 *
 *
 *
 *
 * Calling Sequence:
 *   ChineseRemainderPos(basislen, RNSbasis, RNScombi, Ac, mp_Ac)
 *
 * Summary:
 *   Given a representation of a non-negative integer in some RNS, use Chinese
 *   Remainder Algorithm to reconstruct the integer
 *
 * Description:
 *   Let A be a non-negative integer, and Ac contains the representation of A
 *   in a RNS RNSbasis 'basis', i.e. Ac[i] = mod(A, RNSbasis[i]), (i = 0..basislen).
 *   Here 'mod' is in positive representation. The function reconstructs the
 *   integer A given the RNS RNSbasis 'basis' and Ac.
 *
 *   To avoid repeat computations, the function takes some precomputed
 *   informations as input, which are listed below.
 *
 * Input:
 *       basislen: long, dimension of RNS RNSbasis
 *     RNSbasis: 1-dim typename _Field::Element array length basislen, RNS RNSbasis
 *   RNScombi: 1-dim typename _Field::Element array length basislen, computed by function
 *            combBasis, inverses of special combination of RNS RNSbasis
 *        Ac: 1-dim Double array length n, representation of A in RNS
 *
 * Output:
 *   mp_Ac: mpz_t, reconstructed integer A
 *
 * Precondition:
 *   Let q be product of all elements in the RNS RNSbasis. Then A must satisfy
 *   0 <= A <= q-1.
 *
 */

template<class Container, class FiniteField>
void
RNS<Container,FiniteField>::
ChineseRemainder ( BlasVector<Field> &Ac, Integer & mp_Ac
		  , bool pos )
{
	long i, j;
	Element temp, tempq, tempqinv;
	// bool pos = FieldTraits<_Field>::Rep == Positive;

	BlasVector<Field> U(Ac.field(),basislen);
	// Double *U;
	// U = XMALLOC(Double, basislen);

	/* compute the coefficients of mix radix in positive representation by
	   inplacing modular matrix U */
	U[0] = Ac[0];
	for (i = 1; i < basislen; i++) {
		U[i] = U[i-1];
		tempq = (Element)RNSbasis[i];
		Field Fq(tempq);
		tempqinv = (Element)RNScombi[i];
		for (j = i-2; j >= 0; j--) {
			Fq.mulin(U[i],RNSbasis[j]);
			Fq.addin(U[i],U[j]);
			// U[i] = U[j] + U[i]*fmod((Element)RNSbasis[j], tempq);
			// U[i] = fmod(U[i], tempq);
		}
		Fq.init(temp,tempqinv*(Element)(RNSbasis[i]-1));
		// temp = fmod(tempqinv*(Element)(RNSbasis[i]-1), tempq);
		Fq.init(U[i],tempqinv*Ac[i]+temp*U[i]);
		// U[i] = fmod(tempqinv*Ac[i]+temp*U[i], tempq);
	}
	/* compute Ac in positive representation */
	mp_Ac = U[basislen-1];
	// mpz_set_d(mp_Ac, U[basislen-1]);
	for (j = basislen-2; j >= 0; j--)
	{
		Integer::mulin(mp_Ac,RNSbasis[j]);
		Integer::addin(mp_Ac, (Element)U[j]);
		// mpz_mul_ui(mp_Ac, mp_Ac, RNSbasis[j]);
		// mpz_add_ui(mp_Ac, mp_Ac, (Element)U[j]);
	}
	/* transfer from positive representation to symmetric representation */
	if (pos == false) {
		for (j = basislen-1; j >= 0; j--)
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
 *   RNScombi <-- combBasis(basislen, RNSbasis)
 *
 * Summary:
 *   Compute the special combination of a RNS RNSbasis
 *
 * Description:
 *   Let 'basis' be RNS RNSbasis. The function computes an array RNScombi
 *   satisfying
 *   RNScombi[0] = 0, RNScombi[i] = mod(1/(RNSbasis[0]*...*RNSbasis[i-1]), RNSbasis[i])
 *                   (i = 1..basislen-1)
 *
 * Input:
 *   basislen: long, dimension of RNS RNSbasis
 *      RNSbasis: 1-dim typename _Field::Element array length basislen, RNS RNSbasis
 *
 * Return:
 *   RNScombi: 1-dim typename _Field::Element array length basislen, shown as above
 *
 */

template<class Container, class FiniteField>
void
RNS<Container,FiniteField>::
combBasis ()
{
	long i, j;
	Element prod;
	// mpz_t mp_prod, mp_q;
	// typename _Field::Element *RNScombi;

	// RNScombi = XMALLOC(typename _Field::Element, basislen);
	RNScombi.resize(basislen);
	RNScombi[0] = 0;
	// mpz_init(mp_prod);
	// mpz_init(mp_q);
	for (i = 1; i < basislen; i++)
	{
		FiniteField Fi(RNSbasis[i]);
		Fi.init(prod,RNSbasis[0]);
		// prod = fmod((double)RNSbasis[0], (double)RNSbasis[i]);
		for (j = 1; j <= i-1; j++){
			prod = fmod(prod*(double)RNSbasis[j], (double)RNSbasis[i]);
			Fi.mulin(prod,RNSbasis[j]);
		}
		// Integer mp_q = RNSbasis[i];
		// mpz_set_ui(mp_q, RNSbasis[i]);
		// Integer mp_prod = prod ;
		// mpz_set_d(mp_prod, prod);
		// mpz_invert(mp_prod, mp_prod, mp_q);
		// RNScombi[i] = mpz_get_ui(mp_prod);
		Fi.inv(RNScombi[i],prod);
	}
	// mpz_clear(mp_prod);
	// mpz_clear(mp_q);
	return;

}



/*
 *
 * Calling Sequence:
 *   cumprod <-- cumProd(basislen, RNSbasis, extbasislen, extbasis)
 *
 * Summary:
 *   Compute the representation of the combination of elements of one RNS RNSbasis
 *   in another RNS RNSbasis
 *
 * Description:
 *   Let 'basis' be one RNS RNSbasis with dimension basislen, and 'extbasis' be
 *   another RNS RNSbasis with dimension extbasislen. The function computes an
 *   array cumprod length extbasislen satisfying
 *   cumprod[i] = modp(-RNSbasis[0]*...*RNSbasis[basislen-1], extbasis[i]),
 *   i = 0..extbasislen-1
 *
 * Input:
 *      basislen: long, dimension of RNS RNSbasis 'basis'
 *         RNSbasis: 1-dim typename _Field::Element array length basislen, one RNS RNSbasis
 *   extbasislen: long, dimension of RNS RNSbasis 'extbasis'
 *      extbasis: 1-dim typename _Field::Element array length basislen, another RNS RNSbasis
 *
 * Return:
 *   cumprod: 1-dim double array length extbasislen, shown above
 *
 */

template<class Container, class FiniteField>
void
RNS<Container,FiniteField>::
cumProd (ModVect &cumprod_v,
	 const ModVect &extbasis)
{
  long i, j;
  Element dtemp, dextbasis;
  // double *cumprod_v;

  // cumprod_v = XMALLOC(double, extbasislen);
  size_t extbasislen = extbasis.size();
  cumprod_v.resize(extbasislen);
  for (i = 0; i < (long)extbasislen; i++) {
      FiniteField Fq(extbasis[i]);
      // dextbasis = (Element)extbasis[i];
      Fq.init(cumprod_v[i],RNSbasis[0]);
      // cumprod_v[i] = fmod((double)RNSbasis[0], dextbasis);
      for (j = 1; j < RNSbasis.size(); j++) {
	      Fq.init(dtemp,RNSbasis[j]);
	  // dtemp = fmod((double)RNSbasis[j], dextbasis);
	  Fq.mulin(cumprod_v[i],dtemp); //!@bug 2 in one ?
	  // cumprod_v[i] = fmod(cumprod_v[i]*dtemp, dextbasis);
	}
      // cumprod_v[i] = dextbasis-cumprod_v[i];
      Fq.negin(cumprod_v[i]);
    }

  return ;
}


/*
 *
 * Calling Sequence:
 *   basiscmb <-- findRNS(RNS_bound, mp_maxInter, len)
 *
 * Summary:
 *   Find a RNS RNSbasis and its special combination
 *
 * Description:
 *   Given RNS_bound, the upper bound of the RNS RNSbasis, and mp_maxInter, the
 *   function finds a best RNS RNSbasis and a combination of that RNSbasis.
 *
 *   The RNS RNSbasis 'basis' has the property:
 *   - its elements are all primes
 *   - RNSbasis[0] is the largest prime among all the primes at most RNS_bound
 *   - RNSbasis[i+1] is the next prime smaller than RNSbasis[i] (i = 0..len-2)
 *   - RNSbasis[0]*basis[1]*...*basis[len-1] >= mp_maxInter
 *
 *   After finding 'basis', the functions also computes the combination of
 *   'basis' as the operations in function combBasis.
 *
 * Input:
 *     RNS_bound: typename _Field::Element, the upper bound of the RNS RNSbasis
 *   mp_maxInter: mpz_t, the lower bound for the product of elements of RNSbasis
 *
 * Return:
 *   basiscmb: 2-dim typename _Field::Element array, dimension 2 x basislen, where
 *           - basiscmb[0] represents the RNS RNSbasis
 *           - basiscmb[1] represents the special combination of RNSbasis
 *
 * Output:
 *   len: pointer to a long int, storing the dimension of the computed
 *        RNS RNSbasis
 *
 */

template<class Container, class FiniteField>
void
RNS<Container,FiniteField>::
findRNS ()
{
  long i, j, len=0;
  double prod;
  // mpz_t mp_l, mp_prod, mp_q;
  // typename _Field::Element **qqinv;

  // Integer mp_prod = 1;
  // mpz_init_set_ui(mp_prod, 1);
  linbox_check(this->RNS_bound != 0);
  Integer mp_l = this->RNS_bound ;
  // mpz_init_set_ui(mp_l, RNS_bound);
  // qqinv = XMALLOC(typename _Field::Element *, 2);
  // RNSbasis = NULL;
  while (mp_maxInter > mp_prod)
  // while (mpz_cmp(mp_maxInter, mp_prod) > 0)
    {
      ++len;
      this->RNSbasis.resize(len);
      // RNSbasis = XREALLOC(typename _Field::Element, RNSbasis, len);
      // while (Integer::probab_prime_p(mp_l, 10) == 0) {
	      // Integer::subin(mp_l,1); //! @bug why remove 1 ?
	      // mpz_sub_ui(mp_l, mp_l, 1);
      // }
      Givaro::prevprime(mp_l,mp_l);
      this->RNSbasis[len-1] = (mp_l);
      Integer::subin(mp_l,1L);
      // mpz_sub_ui(mp_l, mp_l, 1);
      Integer::mulin(mp_prod,this->RNSbasis[len-1]);
      // mpz_mul_ui(mp_prod, mp_prod, RNSbasis[len-1]);
    }
  // mpz_clear(mp_prod);
  // mpz_clear(mp_l);
  // RNScombi = XMALLOC(typename _Field::Element, len);
  //
  this->basislen = len;
  this->combBasis();
  return;


}



/*
 *
 * Calling Sequence:
 *   maxInter(mp_mag, mp_alpha, n, mp_b)
 *
 * Summary:
 *   Compute the maximum interval of positive and negative results of a
 *   matrix-matrix or matrix-vector product
 *
 * Description:
 *   Let mp_alpha be the maximum magnitude of a m x n matrix A, mp_mag-1 be
 *   the maximum magnitude of a n x k matrix C. The function computes the
 *   maximum interval of positive and negative entries of A.C. That is, the
 *   function computes mp_b satisfying
 *   (mp_b-1)/2 = n*mp_alpha*(mp_mag-1)
 *
 * Input:
 *    mp_mag: mpz_t, mp_mag-1 be the maximum magnitude of matrix C
 *   mp_alpha: mpz_t, maximum magnitude of matrix A
 *          n: long, column dimension of A
 *
 * Output:
 *   mp_b: mpz_t, shown above
 *
 */

template<class Container, class FiniteField>
void
RNS<Container,FiniteField>::
maxInter (const Integer& mp_mag, const Integer& mp_alpha, const long n, Integer& mp_b)
{
  Integer mp_temp;
  // mpz_t mp_temp;

  // mpz_init(mp_temp);
  Integer::sub(mp_temp,mp_mag,1L);
  // mpz_sub_ui(mp_temp, mp_mag, 1);
  mp_b = mp_alpha ;
  // mpz_set(mp_b, mp_alpha);
  Integer::mulin(mp_b,n); //! @bug 2n here ?
  // mpz_mul_ui(mp_b, mp_b, n);
  Integer::mulin(mp_b,mp_temp);
  // mpz_mul(mp_b, mp_b, mp_temp);
  Integer::mulin(mp_b,2L);
  // mpz_mul_ui(mp_b, mp_b, 2);
  Integer::addin(mp_b,1L);
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

template<class Container, class FiniteField>
void
RNS<Container,FiniteField>::
maxExtInter (const Integer & mp_alpha, const long n, Integer & mp_b)
{
	mp_b = 1L;
  // mpz_set_ui(mp_b, 1);
	Integer::axpyin(mp_b,mp_alpha,(long)n);
  // mpz_addmul_ui(mp_b, mp_alpha, n);
	Integer::mulin(mp_b,2L);
  // mpz_mul_ui(mp_b, mp_b, 2);
  Integer::addin(mp_b,1L);
  // mpz_add_ui(mp_b, mp_b, 1);
}



/*
 *
 * Calling Sequence:
 *   bdcoeff <-- repBound(len, RNSbasis, RNScombi)
 *
 * Summary:
 *   Compute the mix radix coefficients of a special integer in a RNS RNSbasis
 *
 * Description:
 *   Given a RNS RNSbasis, suppose the product of elements in the RNSbasis be q,
 *   then this RNS RNSbasis is able to represent integers lying in
 *   [-(q-1)/2, (q-1)/2] and [0, q-1] respectively with symmetric
 *   representation and positive representation. To transfer the result from
 *   positive representation to symmetric representation, the function
 *   computes the mix radix coefficients of the boundary value (q-1)/2 in the
 *   positive representation.
 *
 *   Let RNS RNSbasis be P. The function computes coefficient array U, such that
 * (q-1)/2 = U[0] + U[1]*P[0] + U[2]*P[0]*P[1] +...+ U[len-1]*P[0]*...*P[len-2]
 *
 * Input:
 *       len: long, dimension of RNS RNSbasis
 *     RNSbasis: 1-dim typename _Field::Element array length basislen, RNS RNSbasis
 *   RNScombi: 1-dim typename _Field::Element array length basislen, computed by function
 *            combBasis, inverses of special combination of RNS RNSbasis
 *
 * Output:
 *   bdcoeff: 1-dim typename _Field::Element array length basislen, the coefficient array U above
 *
 */

template<class Container, class FiniteField>
void
RNS<Container,FiniteField>::
repBound ()
{
  long i, j;
  Element dtemp;
  Integer mp_bd, mp_prod;
  // typename _Field::Element *bdcoeff;

  const Vect & q    = RNSbasis;
  const Vect & qinv = RNScombi;

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
 *   Compute the upper bound of a RNS RNSbasis
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

template<class Container, class FiniteField>
RNS<Container,FiniteField>::
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

} // iml
} // LinBox

#endif // __LINBOX_algorithm_iml_rns_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
