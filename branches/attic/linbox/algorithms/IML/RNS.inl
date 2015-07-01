#ifndef __LINBOX_algorithm_iml_rns_INL
#define __LINBOX_algorithm_iml_rns_INL

namespace LinBox { namespace iml {

	/*
	 *
	 * Calling Sequence:
	 *   basisExt(len, n, p, _RNSbasis, _RNScombi, _cumprod, _bdcoeff, R, RE)
	 *
	 * Summary:
	 *   Given a representation of a matrix/vector in some RNS, extend to compute
	 *   the representation in another positive integer
	 *
	 * Description:
	 *   Let R be the representation of a matrix/vector M in a residue _RNSbasis
	 *   'basis', i.e., R[i] = mod(M, _RNSbasis[i]) (i = 0..len-1). The function
	 *   computes the representation of M in another positive integer p,
	 *   RE = mod(M, p) using Garner's algorithm. 'mod' represents positive modular
	 *   operation.
	 *
	 *   Let q be product of all elements in the RNS _RNSbasis. Every entry m in M
	 *   satisfies -(q-1)/2 <= m <= (q-1)/2. That is, M has both positive entries
	 *   and negative entries.
	 *
	 *   if Pos : Let q be product of all elements in the RNS _RNSbasis. Every entry m in M
	 *   satisfies 0 <= m <= q-1. That is, M only contains non-negative entries.

	 *
	 *   To avoid repeat computations, the function takes some precomputed
	 *   informations as input, which are listed below.
	 *
	 * Input:
	 *       RNS_basislen: long, dimension of RNS _RNSbasis
	 *         n: long, length of array RE
	 *         p: FiniteField, modulus
	 *     _RNSbasis: 1-dim FiniteField array length _basislen, RNS _RNSbasis
	 *   _RNScombi: 1-dim FiniteField array length _basislen, computed by function
	 *            combBasis, inverses of special combination of RNS _RNSbasis
	 *   _cumprod: (only in non Pos) double, computed by function cumProd,
	 *            (-_RNSbasis[0]*_RNSbasis[1]*...*_RNSbasis[_basislen-1]) mod p
	 *   _bdcoeff: 1-dim FiniteField array length _basislen, computed by function repBound
	 *         R: 1-dim Double array length n*_basislen, representation of a _basislen x n
	 *            matrix, R[i]=mod(M, _RNSbasis[i]) (i=0.._basislen-1)
	 *
	 * Output:
	 *   RE: 1-dim Double array length n, RE = mod(M, p), the space of RE
	 *       should be allocated before calling the function.
	 *
	 * Precondition:
	 *   t <= 2^53-1, where t is the maximal intermidiate value arised in this
	 *   function,
	 *   t = max(2*(_RNSbasis[_basislen-1]-1)^2, (p-1)^2+_RNSbasis[_basislen-1]-1)

	 * Precondition if Positive:
	 *   t <= 2^53-1, where t is the maximal intermidiate value arised in this
	 *   function, t = max(2*(_RNSbasis[_basislen-1]-1)^2, (p-1)^2+_RNSbasis[_basislen-1]-1)

	 *
	 */


	template<class FiniteField>
	template<class Container>
	void
	RNS<FiniteField>::basisExt ( FiniteField                       & Fp
				     , std::vector<Container>          & R
				     , Container                       & RE
				     , bool                            pos )
	{
		// bool pos = FieldTraits<_Field>::Rep == Positive;

		size_t i, j;
		size_t n = RE.size();
		ModElement temp;
		ModElement p = Fp.characteristic();
		// ModElement q, qinv;
		// Double **U;
		// const FiniteField *q, *qinv;

		const ModVect & q    = _RNSbasis; //!@bug no copy ?
		const ModVect & qinv = _RNScombi;
		// q = _RNSbasis;
		// qinv = _RNScombi;

		/* if p = q[i] then just copy the corresponding column to RE */
		//! @todo order _RNSbasis ?
		for (i = 0; i < _basislen ; i++) {
			if (p == q[i]) {
				RE = R[i] ;
				// cblas_dcopy(n, R[i], 1, RE, 1);
				return;
			}
		}
		typename Container::Field  F = RE.field();

		const Container Z(F);
		std::vector<Container> U(_basislen,Z);
		// U = XMALLOC(Double *, _basislen);
		// for (i = 0; i < _basislen; i++) { U[i] = XMALLOC(Double, n); }

		/* compute the coefficients of mix radix in positive representation by
		   inplacing modular matrix U */
		//!@bug copy method for blas things (vector/blas and sub)
		U[0].resize(R[0]);
		FFLAS::fcopy(_unF,
			     R[0].getPointer(),1,
			     U[0].getWritePointer(),1);
		// cblas_dcopy(n, R[0], 1, U[0], 1);

		for (i = 1; i < _basislen; i++) {
			FiniteField Fq(q[i]);
			FFLAS::fcopy(_unF,
				     U[i-1].getPointer(),1,
				     U[i].getWritePointer(),1);
			U[i].resize(R[0]);
			// cblas_dcopy(n, U[i-1], 1, U[i], 1);
			// for (j = i-2; j >= 0; j--) {}
			for (j = i-1; j-- ; )
			{
				FFLAS::fscal(_unF,n,(ModElement)(q[j] % q[i]),
					     U[i].getWritePointer(),1);
				// cblas_dscal(n, (double)(q[j] % q[i]), U[i], 1);
				FFLAS::faxpy(Fq,n,Fq.one,
					     U[j].getPointer(),1,U[i].getWritePointer(),1);
				// cblas_daxpy(n, 1.0, U[j], 1, U[i], 1);
				// Dmod((double)q[i], U[i], 1, n, n);
			}
			Fq.init(temp,(ModElement)qinv[i]*(ModElement)(q[i]-1));
			// temp = (ModElement)qinv[i]*(ModElement)(q[i]-1);
			// temp = fmod(temp, (ModElement)q[i]);
			FFLAS::fscal(_unF,n,temp,U[i].getWritePointer(),1);
			// cblas_dscal(n, temp, U[i], 1);
			FFLAS::faxpy(Fq,n,qinv[i],R[i].getPointer(),1,
				     U[i].getWritePointer(),1);
			// cblas_daxpy(n, (double)qinv[i], R[i], 1, U[i], 1);
			// Dmod((double)q[i], U[i], 1, n, n);
		}

		/* compute mod(r, p) in positive representation and store into RE */
		FFLAS::fcopy(Fp,n,
			     RE.getWritePointer(),1,U[_basislen-1].getPointer(),1);
		// cblas_dcopy(n, U[_basislen-1], 1, RE, 1);
		// Dmod((double)p, RE, 1, n, n);
		ModElement qmp ;
		// for (i = _basislen-2; i >= 0; i--) {}
		for (i = _basislen-1; i-- ; ) {
			FFLAS::fscal(_unF,n,(ModElement)Fp.init(qmp,q[i]),RE.getWritePointer(),1);
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
				// for (j = _basislen-1; j >= 0; j--) {}
				for (j = _basislen; j-- ;) {
					//!@bug does not work for every field (already balanced ?)
					if (U[j].getPointer()[i] > _bdcoeff[j]) {
						Fp.init(RE.getPointer()[i],RE[i]+_cumprod);
						break;
					}
					else if (U[j].getPointer()[i] < _bdcoeff[j]) {
						break;
					}
				}
		}
		// for (i = 0; i < _basislen; i++) { XFREE(U[i]); } { XFREE(U); }

		return;
	}


	/*
	 *
	 * Calling Sequence:
	 *   ChineseRemainder(_basislen, _mp_prod, _RNSbasis, _RNScombi, _bdcoeff, Ac, mp_Ac)
	 *
	 * Summary:
	 *   Given a representation of an integer in some RNS, use Chinese Remainder
	 *   Algorithm to reconstruct the integer
	 *
	 * Description:
	 *   Let A be an integer, and Ac contains the representation of A in a RNS
	 *   _RNSbasis 'basis', i.e. Ac[i] = mod(A, _RNSbasis[i]), (i = 0.._basislen). Here 'mod'
	 *   is in positive representation. The function reconstructs the integer A
	 *   given the RNS _RNSbasis 'basis' and Ac.
	 *
	 *   To avoid repeat computations, the function takes some precomputed
	 *   informations as input, which are listed below.
	 *
	 * Input:
	 *       _basislen: long, dimension of RNS _RNSbasis
	 *   _mp_prod: mpz_t, computed by function basisProd, product of RNS _RNSbasis
	 *     _RNSbasis: 1-dim FiniteField array length _basislen, RNS _RNSbasis
	 *   _RNScombi: 1-dim FiniteField array length _basislen, computed by function
	 *            combBasis, inverses of special combination of RNS _RNSbasis
	 *   _bdcoeff: 1-dim FiniteField array length basilen, computed by function repBound
	 *        Ac: 1-dim Double array length n, representation of A in RNS
	 *
	 * Output:
	 *   mp_Ac: mpz_t, reconstructed integer A
	 *
	 * Precondition:
	 *   Let q be product of all elements in the RNS _RNSbasis. Then A must satisfy
	 *   -(q-1)/2 <= A <= (q-1)/2.
	 *
	 *
	 *
	 *
	 * Calling Sequence:
	 *   ChineseRemainderPos(_basislen, _RNSbasis, _RNScombi, Ac, mp_Ac)
	 *
	 * Summary:
	 *   Given a representation of a non-negative integer in some RNS, use Chinese
	 *   Remainder Algorithm to reconstruct the integer
	 *
	 * Description:
	 *   Let A be a non-negative integer, and Ac contains the representation of A
	 *   in a RNS _RNSbasis 'basis', i.e. Ac[i] = mod(A, _RNSbasis[i]), (i = 0.._basislen).
	 *   Here 'mod' is in positive representation. The function reconstructs the
	 *   integer A given the RNS _RNSbasis 'basis' and Ac.
	 *
	 *   To avoid repeat computations, the function takes some precomputed
	 *   informations as input, which are listed below.
	 *
	 * Input:
	 *       _basislen: long, dimension of RNS _RNSbasis
	 *     _RNSbasis: 1-dim FiniteField array length _basislen, RNS _RNSbasis
	 *   _RNScombi: 1-dim FiniteField array length _basislen, computed by function
	 *            combBasis, inverses of special combination of RNS _RNSbasis
	 *        Ac: 1-dim Double array length n, representation of A in RNS
	 *
	 * Output:
	 *   mp_Ac: mpz_t, reconstructed integer A
	 *
	 * Precondition:
	 *   Let q be product of all elements in the RNS _RNSbasis. Then A must satisfy
	 *   0 <= A <= q-1.
	 *
	 */

	template<class FiniteField>
	void
	RNS<FiniteField>::
	ChineseRemainder ( BlasVector<NoField> &Ac, Integer & mp_Ac
			   , bool pos )
	{
		// for (size_t z = 0 ; z < _basislen ; ++z) {
			// std::cout << (long)Ac[z] << " (mod " << (long)_RNSbasis[z] << ") ," ;
		// }
		// std::cout << std::endl;
		// for (size_t z = 0 ; z < _basislen ; ++z) {
			// std::cout << (long)_RNScombi[z] << " ( " << (long)_bdcoeff[z] << ") ," ;
		// }
		// std::cout << std::endl;


		size_t i, j;
		ModElement temp, tempqinv;

		BlasVector<NoField> U(Ac.field(),_basislen);

		/* compute the coefficients of mix radix in positive representation by
		   inplacing modular matrix U */
		U[0] = Ac[0];
		for (i = 1; i < _basislen; i++) {
			U[i] = U[i-1];
			FiniteField Fq(_RNSbasis[i]);
			tempqinv = (ModElement)_RNScombi[i];
			for (j = i-1;  j-- ;) {
				Fq.mulin(U[i],_RNSbasis[j]);
				Fq.addin(U[i],U[j]);
			}
			Fq.init(temp,tempqinv*(ModElement)(_RNSbasis[i]-1));
			Fq.init(U[i],tempqinv*Ac[i]+temp*U[i]);
		}
		/* compute Ac in positive representation */
		mp_Ac = U[_basislen-1];
		for (j = _basislen-1; j--; ) {
			Integer::mulin(mp_Ac, (long int)_RNSbasis[j]);
			Integer::addin(mp_Ac, (long int)U[j]);
		}
		/* transfer from positive representation to symmetric representation */
		if (pos == false) {
			for (j = _basislen; j-- ; ) {
				if (U[j] > _bdcoeff[j]) {
					Integer::subin(mp_Ac,_mp_prod);
					break;
				}
				else if (U[j] < _bdcoeff[j]) { break; }
			}
		}

		return;
	}

	/*
	 *
	 * Calling Sequence:
	 *   _RNScombi <-- combBasis(_basislen, _RNSbasis)
	 *
	 * Summary:
	 *   Compute the special combination of a RNS _RNSbasis
	 *
	 * Description:
	 *   Let 'basis' be RNS _RNSbasis. The function computes an array _RNScombi
	 *   satisfying
	 *   _RNScombi[0] = 0, _RNScombi[i] = mod(1/(_RNSbasis[0]*...*_RNSbasis[i-1]), _RNSbasis[i])
	 *                   (i = 1.._basislen-1)
	 *
	 * Input:
	 *   _basislen: long, dimension of RNS _RNSbasis
	 *      _RNSbasis: 1-dim FiniteField array length _basislen, RNS _RNSbasis
	 *
	 * Return:
	 *   _RNScombi: 1-dim FiniteField array length _basislen, shown as above
	 *
	 */

	template<class FiniteField>
	void
	RNS<FiniteField>::
	combBasis (BlasVector<NoField>& RNScombi, const BlasVector<NoField> &RNSbasis)
	{
		size_t i, j;
		ModElement prod;
		size_t basislen = RNSbasis.size();

		RNScombi.resize(basislen);
		RNScombi[0] = 0;
		for (i = 1; i < basislen; i++)
		{
			FiniteField Fi(RNSbasis[i]);
			Fi.init(prod,RNSbasis[0]);
			for (j = 1; j <= i-1; j++){
				Fi.mulin(prod,RNSbasis[j]);
			}
			Fi.inv(RNScombi[i],prod);
		}
		return;
	}


	/*
	 *
	 * Calling Sequence:
	 *   _cumprod <-- cumProd(_basislen, _RNSbasis, ext_basislen, extbasis)
	 *
	 * Summary:
	 *   Compute the representation of the combination of elements of one RNS _RNSbasis
	 *   in another RNS _RNSbasis
	 *
	 * Description:
	 *   Let 'basis' be one RNS _RNSbasis with dimension _basislen, and 'extbasis' be
	 *   another RNS _RNSbasis with dimension ext_basislen. The function computes an
	 *   array _cumprod length ext_basislen satisfying
	 *   _cumprod[i] = modp(-_RNSbasis[0]*...*_RNSbasis[_basislen-1], extbasis[i]),
	 *   i = 0..ext_basislen-1
	 *
	 * Input:
	 *      _basislen: long, dimension of RNS _RNSbasis 'basis'
	 *         _RNSbasis: 1-dim FiniteField array length _basislen, one RNS _RNSbasis
	 *   ext_basislen: long, dimension of RNS _RNSbasis 'extbasis'
	 *      extbasis: 1-dim FiniteField array length _basislen, another RNS _RNSbasis
	 *
	 * Return:
	 *   _cumprod: 1-dim double array length ext_basislen, shown above
	 *
	 */

	template<class FiniteField>
	void
	RNS<FiniteField>::
	cumProd (ModVect &_cumprod_v,
		 const ModVect &extbasis)
	{
		size_t i, j;
		ModElement dtemp, dextbasis;
		// double *_cumprod_v;

		// _cumprod_v = XMALLOC(double, ext_basislen);
		size_t ext_basislen = extbasis.size();
		_cumprod_v.resize(ext_basislen);
		for (i = 0; i < ext_basislen; i++) {
			FiniteField Fq(extbasis[i]);
			// dextbasis = (ModElement)extbasis[i];
			Fq.init(_cumprod_v[i],_RNSbasis[0]);
			// _cumprod_v[i] = fmod((double)_RNSbasis[0], dextbasis);
			for (j = 1; j < _RNSbasis.size(); j++) {
				Fq.init(dtemp,_RNSbasis[j]);
				// dtemp = fmod((double)_RNSbasis[j], dextbasis);
				Fq.mulin(_cumprod_v[i],dtemp); //!@bug 2 in one ?
				// _cumprod_v[i] = fmod(_cumprod_v[i]*dtemp, dextbasis);
			}
			// _cumprod_v[i] = dextbasis-_cumprod_v[i];
			Fq.negin(_cumprod_v[i]);
		}

		return ;
	}


	/*
	 *
	 * Calling Sequence:
	 *   basiscmb <-- findRNS(RNS_bound, _mp_maxInter, len)
	 *
	 * Summary:
	 *   Find a RNS _RNSbasis and its special combination
	 *
	 * Description:
	 *   Given RNS_bound, the upper bound of the RNS _RNSbasis, and _mp_maxInter, the
	 *   function finds a best RNS _RNSbasis and a combination of that _RNSbasis.
	 *
	 *   The RNS _RNSbasis 'basis' has the property:
	 *   - its elements are all primes
	 *   - _RNSbasis[0] is the largest prime among all the primes at most RNS_bound
	 *   - _RNSbasis[i+1] is the next prime smaller than _RNSbasis[i] (i = 0..len-2)
	 *   - _RNSbasis[0]*basis[1]*...*basis[len-1] >= _mp_maxInter
	 *
	 *   After finding 'basis', the functions also computes the combination of
	 *   'basis' as the operations in function combBasis.
	 *
	 * Input:
	 *     RNS_bound: FiniteField, the upper bound of the RNS _RNSbasis
	 *   _mp_maxInter: mpz_t, the lower bound for the product of elements of _RNSbasis
	 *
	 * Return:
	 *   basiscmb: 2-dim FiniteField array, dimension 2 x _basislen, where
	 *           - basiscmb[0] represents the RNS _RNSbasis
	 *           - basiscmb[1] represents the special combination of _RNSbasis
	 *
	 * Output:
	 *   len: pointer to a long int, storing the dimension of the computed
	 *        RNS _RNSbasis
	 *
	 */

	template<class FiniteField>
	void
	RNS<FiniteField>::
	findRNS ()
	{
		size_t i, j, len=0;
		double prod;

		linbox_check(_RNS_bound != 0);
		Integer mp_l = _RNS_bound ;
		while (_mp_maxInter > _mp_prod) {
			++len;
			_RNSbasis.resize(len);
			Givaro::prevprime(mp_l,mp_l);
			_RNSbasis[len-1] = (mp_l);
			Integer::subin(mp_l,1L);
			Integer::mulin(_mp_prod,(long)_RNSbasis[len-1]);
		}
		_basislen = len;
		combBasis(_RNScombi, _RNSbasis);
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

	template<class FiniteField>
	void
	RNS<FiniteField>::
	maxInter (const Integer& mp_mag, const Integer& mp_alpha, const size_t n, Integer& mp_b)
	{
		Integer mp_temp;
		// mpz_t mp_temp;

		// mpz_init(mp_temp);
		Integer::sub(mp_temp,mp_mag,1L);
		// mpz_sub_ui(mp_temp, mp_mag, 1);
		mp_b = mp_alpha ;
		// mpz_set(mp_b, mp_alpha);
		Integer::mulin(mp_b,(long unsigned)n); //! @bug 2n here ?
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

	template<class FiniteField>
	void
	RNS<FiniteField>::
	maxExtInter (const Integer & mp_alpha, const size_t n, Integer & mp_b)
	{
		mp_b = 1L;
		// mpz_set_ui(mp_b, 1);
		Integer::axpyin(mp_b,mp_alpha,(unsigned long)n);
		// mpz_addmul_ui(mp_b, mp_alpha, n);
		Integer::mulin(mp_b,2L);
		// mpz_mul_ui(mp_b, mp_b, 2);
		Integer::addin(mp_b,1L);
		// mpz_add_ui(mp_b, mp_b, 1);
	}



	/*
	 *
	 * Calling Sequence:
	 *   _bdcoeff <-- repBound(len, _RNSbasis, _RNScombi)
	 *
	 * Summary:
	 *   Compute the mix radix coefficients of a special integer in a RNS _RNSbasis
	 *
	 * Description:
	 *   Given a RNS _RNSbasis, suppose the product of elements in the _RNSbasis be q,
	 *   then this RNS _RNSbasis is able to represent integers lying in
	 *   [-(q-1)/2, (q-1)/2] and [0, q-1] respectively with symmetric
	 *   representation and positive representation. To transfer the result from
	 *   positive representation to symmetric representation, the function
	 *   computes the mix radix coefficients of the boundary value (q-1)/2 in the
	 *   positive representation.
	 *
	 *   Let RNS _RNSbasis be P. The function computes coefficient array U, such that
	 * (q-1)/2 = U[0] + U[1]*P[0] + U[2]*P[0]*P[1] +...+ U[len-1]*P[0]*...*P[len-2]
	 *
	 * Input:
	 *       len: long, dimension of RNS _RNSbasis
	 *     _RNSbasis: 1-dim FiniteField array length _basislen, RNS _RNSbasis
	 *   _RNScombi: 1-dim FiniteField array length _basislen, computed by function
	 *            combBasis, inverses of special combination of RNS _RNSbasis
	 *
	 * Output:
	 *   _bdcoeff: 1-dim FiniteField array length _basislen, the coefficient array U above
	 *
	 */

	template<class FiniteField>
	void
	RNS<FiniteField>::
	repBound ()
	{
		size_t i, j;
		ModElement dtemp;
		Integer mp_bd, loc__mp_prod;

		const ModVect & q    = _RNSbasis;
		const ModVect & qinv = _RNScombi;

		/* set the bound of transformation from positive to negative */
		loc__mp_prod = _mp_prod ;
		Integer::sub(mp_bd,loc__mp_prod,1L);
		Integer::divexact(mp_bd,mp_bd,2L);

		/* compute the coeffcients of bound of mix radix and store in _bdcoeff */
		_bdcoeff.resize(_basislen);
		_bdcoeff[0] = (ModElement) Integer::frem(mp_bd,(long unsigned)q[0]);
		for (i = 1; i < _basislen; i++) {
			dtemp = (ModElement)_bdcoeff[i-1];
			FiniteField Fj(q[i]);
			for (j = i-1; j-- ; ) {
				Fj.mulin(dtemp,q[j]);
				Fj.addin(dtemp,_bdcoeff[j]);
			}
			_bdcoeff[i] = (ModElement)Integer::frem(mp_bd,(long unsigned)q[i]);
			Fj.init(dtemp,_bdcoeff[i]-dtemp);
			//! @bug can it be <0 ?
#if 0
			if (dtemp < 0) {
				dtemp = q[i]+dtemp;
			}
#endif
			Fj.init(_bdcoeff[i],dtemp*qinv[i]);
		}
		return;
	}


	/*
	 * Calling Sequence:
	 *   bd <-- RNSbound(n)
	 *
	 * Summary:
	 *   Compute the upper bound of a RNS _RNSbasis
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
	 *   bd: FiniteField, shown above
	 *
	 */

	template<class FiniteField>
	typename RNS<FiniteField>::ModElement
	RNS<FiniteField>::
	RNSbound (const size_t n)
	{
		ModElement bd;
		Integer mp_n, mp_d, mp_q;

		// mpz_init(mp_n);
		// mpz_init_set_ui(mp_d, n);
		// mpz_init(mp_q);
		Givaro::pow(mp_n,2UL,53UL);
		// mpz_ui_pow_ui(mp_n, 2, 53);
		// mpz_sub_ui(mp_n, mp_n, 1);
		Integer::subin(mp_n,1L);
		Integer::floor(mp_q,mp_n,(unsigned long)n);
		// mpz_fdiv_q(mp_q, mp_n, mp_d);
		Givaro::sqrt(mp_q,mp_q);
		// mpz_sqrt(mp_q, mp_q);
		// bd = mpz_get_ui(mp_q)+1;
		bd = (ModElement)mp_q + 1 ;
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
