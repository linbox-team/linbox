#ifndef __LINBOX_algorithm_iml_reconstruct_solution_INL
#define __LINBOX_algorithm_iml_reconstruct_solution_INL

namespace LinBox { namespace iml {

	/*
	 * Calling Sequence:
	 *   sumliftCoeff(mp_basisprod, k, C, mp_sum)
	 *
	 * Summary:
	 *   Recursively sum up the p-adic lifting coefficients
	 *
	 * Description:
	 *   Given p-adic lifing coefficients C after k lifting steps, the function
	 *   recursively compute
	 *     A^(-1).B mod p^k =
	 *     C[0]+C[1]*mp_basisprod+ ... +C[k-1]*mp_basisprod^(k-1)
	 *
	 * Input:
	 *   mp_basisprod: Integer, product of lifting basis
	 *              k: size_t, dimension of array C, number of lifting steps
	 *              C: 1-dim Integer matrix length k, storing lifting coefficients
	 *                 computed each lifting step
	 *
	 * Output:
	 *   mp_sum: Integer, the sum of lifting coefficients
	 *
	 */

template<class FiniteField>
	void
	SolutionReconstruct<FiniteField>::
	sumliftCoeff (const Integer             &mp_basisprod
		      , const size_t             k
		      , BlasVector<PID_integer> &C
		      , Integer                 &mp_sum)
	{
		size_t i, splflag;
		size_t idx = 0 ;
		size_t t = find2exp(k);
		Integer mp_right;
		PID_integer Z;
		/* precomputation mp_basisprod^log(k) */
		BlasVector<PID_integer> mp_pow  (Z,t+1);
		// mp_pow = XMALLOC(Integer, t+1);
		mp_pow[0] = mp_basisprod ;
		// mpz_init_set(mp_pow[0], mp_basisprod);
		for (i = 1; i < t+1; i++)
		{
			// mpz_init(mp_pow[i]);
			Givaro::pow(mp_pow[i],mp_pow[i-1],2L);
			// mpz_pow_ui(mp_pow[i], mp_pow[i-1], 2);
		}
		// mpz_init(mp_right);
		if (t == 0)
		{
			mp_sum = C[0];
			// mpz_set(mp_sum, C[0]);
			// for (i = 0; i < t+1; i++) { mpz_clear(mp_pow[i]); } { XFREE(mp_pow); }
			// mpz_clear(mp_right);
			return;
		}
		if ((1 << t) == k) {
			splflag = 1;
		}
		else {
			splflag = 0;
		}
		BlasVector<PID_integer> mp_left (Z,t);
		// mp_left = XMALLOC(Integer, t);
		// for (i = 0; i < t; i++) { mpz_init(mp_left[i]); }
		sumCoeff_rec(0, k, C, mp_pow, splflag, 0, idx, mp_left, mp_right);
		mp_sum= mp_left[0];
		// mpz_set(mp_sum, mp_left[0]);

		// for (i = 0; i < t+1; i++) { mpz_clear(mp_pow[i]); } { XFREE(mp_pow); }
		// for (i = 0; i < t; i++) { mpz_clear(mp_left[i]); } { XFREE(mp_left); }
		// mpz_clear(mp_right);
		return;
	}


	/*
	 * Calling Sequence:
	 *   sumCoeff_rec(s, len, C, mp_pow, splflag, 0, idx, mp_left, mp_right)
	 *
	 * Summary:
	 *   compute C[s]+C[s+1]*p+C[s+2]*p^2+...+C[s+len-1]*p^(len-1) recursively
	 *
	 * Input:
	 *      start: size_t, recursion start index in array C, usually 0
	 *        len: size_t, number of entries in C to use
	 *          C: 1-dim Integer array, storing values for computation
	 *     mp_pow: 1-dim Integer array length t+1, storing different powers of p,
	 *             mp_pow[i] = p^(2^i), (i = 0..t), t = floor(log(len)/log(2))
	 *    splflag: 1/0, input flag,
	 *           - if len=2^i and i >= 0, then splflag = 1
	 *           - else splflag = 0
	 *    savflag: 1/0, input flag, initially use 0
	 *           - if savflag = 0, then current subproblem is in the left subtree
	 *             of recursion tree
	 *           - else, current subproblem is in the right subtree of recursion
	 *             tree
	 *        idx: pointer to a size_t integer, the index of available entries in
	 *             array mp_left, initially *idx = 0
	 *    mp_left: 1-dim Integer array length t, storing the intermidiate result of
	 *             left subtree
	 *   mp_right: Integer, storing the result of right subtree
	 *
	 * Output:
	 *   mp_left[0] will store the computation result
	 *
	 * Idea: the result of subproblem on the right subtree of recursion tree is
	 *       overwritten into mp_right
	 *       the result of subproblem on the left subtree of recursion tree is
	 *       stored into mp_left[*idx], and make the *idx always be availabe
	 *       index in mp_left, e.g., the space of mp_left used in right subtree can
	 *       can be reused by its siblings and parent
	 *
	 */

template<class FiniteField>
	void
	SolutionReconstruct<FiniteField>::
	sumCoeff_rec (size_t                      start
		      , size_t                    len
		      , BlasVector<PID_integer> & C
		      , BlasVector<PID_integer> & mp_pow
		      , size_t                    splflag
		      , size_t                    savflag
		      , size_t                  & idx
		      , BlasVector<PID_integer> & mp_left
		      , Integer                 & mp_right)
	{
		size_t lpow, llen, tidx;

		/* basecase */
		if (len == 1) {
			if (savflag == 0) {
				mp_left[idx] = C[start];
				// mpz_set(mp_left[*idx], C[start]);
				++idx;
			}
			else {
				mp_right = C[start];
				// mpz_set(mp_right, C[start]);
			}
			return;
		}
		/* if len is power of 2, then directly divide len by 2 to
		   obtain two subproblems*/
		if (splflag == 1) {
			lpow = find2exp(len)-1;
			llen = 1 << lpow;
			sumCoeff_rec(start, llen, C, mp_pow, 1, 0, idx, mp_left, mp_right);
			tidx = idx-1;
			sumCoeff_rec(start+llen, len-llen, C, mp_pow, 1, 1, idx, mp_left, \
				     mp_right);
		}
		/* otherwise, find maximal lpow, 2^lpow < len, to make two subproblems */
		else {
			lpow = find2exp(len);
			llen = 1 << lpow;
			sumCoeff_rec(start, llen, C, mp_pow, 1, 0, idx, mp_left, mp_right);
			tidx = idx-1;

			/* if len-2^lpow is power of 2, then not compute the second subproblem */
			if (llen == len)
			{
				mp_right = mp_left[tidx];
				// mpz_set(mp_right, mp_left[tidx]);
				return;
			}
			sumCoeff_rec(start+llen, len-llen, C, mp_pow, 0, 1, idx, mp_left, \
				     mp_right);
		}

		/* combine subproblems by position of two subproblems in the recursion tree*/
		/* if subproblem lies in the left subtree, store result in a new Integer int */
		if (savflag == 0) {
			Integer::axpyin(mp_left[idx],mp_right, mp_pow[lpow]);
			// mpz_addmul(mp_left[tidx], mp_right, mp_pow[lpow]);
		}

		/* if in the right subtree, store result in the pre-allocated mp_right */
		else {
			Integer::mulin(mp_right,mp_pow[lpow]);
			// mpz_mul(mp_right, mp_right, mp_pow[lpow]);
			Integer::addin(mp_right,mp_left[tidx]);
			// mpz_add(mp_right, mp_left[tidx], mp_right);
		}
		idx = tidx+1;
		return;
	}


	/*
	 * Calling Sequence:
	 *   i <-- find2exp(len)
	 *
	 * Summary:
	 *   Compute floor(log[2](len))
	 *
	 * Description:
	 *   Compute i >= 0 such that 2^i <= len and 2^(i+1) > len
	 *
	 */

template<class FiniteField>
	size_t
	SolutionReconstruct<FiniteField>::
	find2exp (const size_t len)
	{
		return std::floor(log2(len));
		// size_t i=0, l=1;
		// while ((l <<= 1) <= len) { ++i; }
		// return i;
	}


	/*
	 * Calling Sequence:
	 *   1/0 <-- iratrecon(mp_u, mp_m, mp_nb, mp_db, mp_N, mp_D)
	 *
	 * Summary:
	 *   Perform rational number reconstruction
	 *
	 * Description:
	 *   Given image mp_u, modulus mp_m, the function reconstruct numerator mp_N
	 *   and denominator mp_D, such that
	 *   mp_N/mp_D = mp_u (mod mp_m), and
	 *   abs(mp_N) <= mp_nb, 0 < mp_D <= mp_db
	 *
	 * Input:
	 *    mp_u: Integer, image
	 *    mp_m: Integer, modulus
	 *   mp_nb: Integer, numerator bound of the rational number
	 *   mp_db: Integer, denominator bound of the rational number
	 *
	 * Output:
	 *   mp_N: Integer, storing the reconstructed numerator if the reconstruction
	 *         succeeds
	 *   mp_D: Integer, storing the reconstructed denominator if the reconstruction
	 *         fails
	 *
	 * Return:
	 *   - 1 if reconstruction succeeds
	 *   - 0 if reconstruction fails
	 *
	 */

template<class FiniteField>
	int
	SolutionReconstruct<FiniteField>::
	iratrecon (const Integer &mp_u
		   , const Integer &mp_m
		   , const Integer &mp_nb
		   , const Integer &mp_db
		   , Integer &mp_N
		   , Integer &mp_D)
	{
		Integer mp_v3, mp_u2, mp_v2, mp_u3, mp_t2, mp_t3, mp_q;

		// mpz_init(mp_v3);
		mp_v3 = mp_u % mp_m ;
		// mpz_tdiv_r(mp_v3, mp_u, mp_m);
		mp_u2 = 0UL;
		// mpz_init_set_ui(mp_u2, 0);
		mp_u3 = mp_m ;
		// mpz_init_set(mp_u3, mp_m);
		mp_v2 = 1UL ;
		// mpz_init_set_ui(mp_v2, 1);
		// mpz_init(mp_q);
		// mpz_init(mp_t2);
		// mpz_init(mp_t3);
		while (absCompare(mp_v3, mp_nb) > 0)
		{
			Integer::divmod(mp_q,mp_t3,mp_u3,mp_v3);
			// linbox_check(a >=0 || r==0); // XXX otherwise divmod != mpz_tdiv_qr
			// mpz_tdiv_qr(mp_q, mp_t3, mp_u3, mp_v3);
			Integer::mul(mp_t2,mp_q,mp_v2);
			// mpz_mul(mp_t2, mp_q, mp_v2);
			Integer::sub(mp_t2,mp_u2,mp_t2);
			// mpz_sub(mp_t2, mp_u2, mp_t2);
			mp_u2 = mp_v2;
			// mpz_set(mp_u2, mp_v2);
			mp_u3 = mp_v3;
			// mpz_set(mp_u3, mp_v3);
			mp_v2 = mp_t2;
			// mpz_set(mp_v2, mp_t2);
			mp_v3 = mp_t3;
			// mpz_set(mp_v3, mp_t3);
		}
		if (absCompare(mp_v2, mp_db) > 0)
		{
			// mpz_clear(mp_v2); mpz_clear(mp_v3); mpz_clear(mp_u2); mpz_clear(mp_q);
			// mpz_clear(mp_t2); mpz_clear(mp_t3); mpz_clear(mp_u3);
			return 0;
		}
		gcd(mp_u2,mp_v2,mp_v3);
		linbox_check(mp_u2>=0); // XXX otherwise, different
		// mpz_gcd(mp_u2, mp_v2, mp_v3);
		if (absCompare(mp_u2, 1UL) != 0)
		{
			// mpz_clear(mp_v2); mpz_clear(mp_v3); mpz_clear(mp_u2); mpz_clear(mp_q);
			// mpz_clear(mp_t2); mpz_clear(mp_t3); mpz_clear(mp_u3);
			return 0;
		}
		if ( mp_v2.priv_sign() < 0) {
			Integer::negin(mp_v2);
			// mpz_neg(mp_v2, mp_v2);
			Integer::negin(mp_v3);
			// mpz_neg(mp_v3, mp_v3);
		}
		{
			mp_D = mp_v2;
			// mpz_set(mp_D, mp_v2);
			mp_N = mp_v3;
			// mpz_set(mp_N, mp_v3);
		}

		// { mpz_clear(mp_v2); mpz_clear(mp_v3); mpz_clear(mp_u2); mpz_clear(mp_u3); }
		// { mpz_clear(mp_q); mpz_clear(mp_t2);  mpz_clear(mp_t3); }
		return 1;
	}


	/*
	 * Calling Sequence:
	 *   1/0 <-- soluRecon(solupos, k, basislen, n, m, mp_basisprod, basis
	 *                     cmbasis, C, mp_nb, mp_db, mp_N, mp_D)
	 *
	 * Summary:
	 *   Try reconstructing the rational solution using p-adic lifting coefficients
	 *
	 * Description:
	 *   Given p-adic lifting coefficients computed from function lift, this
	 *   function try reconstructing the rational solution of a nonsingular
	 *   system. If the trial succeeds, the function outputs the reconstructed
	 *   numerator matrix and denominator and returns 1. If the system is not
	 *   lifted high enough, usually the trial will fail in the first several
	 *   solution entries. Then the function returns 0.
	 *
	 * Input:
	 *        solupos: enumerate, flag to indicate whether to transpose A or not
	 *               - solupos = LeftSolu: system be XA = B, B m x n matrix
	 *               - solupos = RightSolu: system be AX = B, B n x m matrix
	 *              k: size_t, first dimension of C, number of lifting steps
	 *       basislen: size_t, dimension of lifting basis
	 *              n: size_t, dimension of left hand side input matrix A
	 *              m: size_t, row/column dimension of right hand side input matrix B
	 *   mp_basisprod: Integer, product of lifting basis
	 *          basis: 1-dim FiniteField array length basislen, lifting basis
	 *        cmbasis: 1-dim FiniteField array length len, computed by function
	 *                 combBasis, inverses of special combination of lifting basis
	 *              C: 3-dim Double array, dimension k x liftbasislen x n*m,
	 *                 computed by function lift.
	 *          mp_nb: Integer, bound of absoluate value of numerator of
	 *                 reconstructed rational number
	 *          mp_db: Integer, bound of denominator of reconstructed rational number
	 *
	 * Output:
	 *   mp_N: 1-dim Integer array length n*m, representation of a n x m or m x n
	 *         reconstructed numerator matrix of the solution
	 *   mp_D: Integer, denominator of the solution
	 *
	 * Return:
	 *   - 1 if reconstruction succeeds
	 *   - 0 if reconstruction fails
	 *
	 */

template<class FiniteField>
	int
	SolutionReconstruct<FiniteField>::
	soluRecon (const LinBoxTag::Side solupos
		   , const size_t              k
		   , RNS<FiniteField>        & rns
		   , pAdicLift<FiniteField>  & C
		   , Integer                 & mp_nb
		   , Integer                 & mp_db
		   , BlasVector<PID_integer> & mp_N
		   , Integer                 & mp_D)
	{
		size_t i, j, s, t, h, l, r, ri, len=1;
		Integer mp_sum, mp_m, mp_temp;
		PID_integer Z;
		// Integer *mp_Dt, *mp_C;
		// size_t *idx ;
		// Double *C1;

		// mpz_init(mp_sum);
		// mpz_init(mp_m);
		Givaro::pow(mp_m,rns.basisProd(),(unsigned long)k);
		// mpz_pow_ui(mp_m, mp_basisprod, k);
		BlasVector<PID_integer> mp_Dt(Z,len);
		// mp_Dt = XMALLOC(Integer, len);
		mp_Dt[0] = Integer(0) ; // XXX try Z.zero
		// mpz_init(mp_Dt[0]);

		/* idx stores the index in mp_N where denominator changes */
		std::vector<size_t> idx(len);
		// idx = XCALLOC(size_t, len);
		BlasVector<PID_integer> mp_C(Z,k);
		// mp_C = XMALLOC(Integer, k);
		// for (j = 0; j < k; j++) { mpz_init(mp_C[j]); }
		BlasVector<PID_integer> C1(Z,rns.size());
		// C1 = XMALLOC(Double, basislen);
		for (i = 0; i < k; i++) {
			for (j = 0; j < rns.size(); j++) {
				C1[j] = C.getEntry(i,j,0);
			}
			rns.ChineseRemainder(C1, mp_C[i],LinBoxTag::Positive); //positive
		}
		sumliftCoeff(rns.basisProd(), k, mp_C, mp_sum);
		ri = iratrecon(mp_sum, mp_m, mp_nb, mp_db, mp_N[0], mp_Dt[0]);
		if (ri == 0) {
			// for (i = 0; i < len; i++) { mpz_clear(mp_Dt[i]); } { XFREE(mp_Dt); }
			// for (i = 0; i < k; i++) { mpz_clear(mp_C[i]); } { XFREE(mp_C); }
			// { mpz_clear(mp_sum);  mpz_clear(mp_m); }
			// { XFREE(C1); XFREE(idx); }
			return 0;
		}
		mp_D = mp_Dt[0];
		// mpz_set(mp_D, mp_Dt[0]);

		/* reconstruct mp_N[i], i = 1..n*m */
		size_t n = C.rowdim();
		size_t m = C.coldim();
		for (i = 1; i < m * n ; i++) {
			/* transform i to the corresponding index h in C when solupos==LeftSolu*/
			{ s = (size_t)(i/n); t = i-s*n; h = t*m+s; }
			for (j = 0; j < k; j++) {
				if (solupos == RightSolu)
					for (l = 0; l < basislen; l++) {
						C1[l] = C.getEntry(j,l,i);
						// C1[l] = C[j][l][i];
					}
				else if (solupos == LeftSolu)
					for (l = 0; l < basislen; l++) {
						C1[l] = C.getEntry(j,l,h);
						// C1[l] = C[j][l][h];
					}
				ChineseRemainder(basislen, basis, cmbasis, C1, mp_C[j],LinBoxTag::Positive);
			}

			sumliftCoeff(mp_basisprod, k, mp_C, mp_sum);

			/* try finding numerator using denominator obtained in last operation */
			r = findNumer(mp_sum, mp_m, mp_D, mp_nb, mp_N[i]);

			/* fail to find the numerator, apply rational reconstruction, mark
			   the additional multiple of denominator, and update mp_D */
			if (r == 0) {
				++len;
				mp_Dt.resize(len);
				// mp_Dt = XREALLOC(Integer, mp_Dt, len);
				// mpz_init(mp_Dt[len-1]);
				idx.resize(len);
				// idx = XREALLOC(size_t, idx, len);
				idx[len-1] = i;
				ri = iratrecon(mp_N[i], mp_m, mp_nb, mp_db, mp_N[i], mp_Dt[len-1]);
				if (ri == 0) {
					// for (j = 0; j < len; j++) { mpz_clear(mp_Dt[j]); }
					// { XFREE(mp_Dt);  XFREE(C1);  XFREE(idx); }
					// for (j = 0; j < k; j++) { mpz_clear(mp_C[j]); } { XFREE(mp_C); }
					// { mpz_clear(mp_sum);  mpz_clear(mp_m); }
					return 0;
				}
				Integer::mulin(mp_D,mp_Dt[len-1]);
				// mpz_mul(mp_D, mp_D, mp_Dt[len-1]);
			}
		}

		/* normalize numerators */
		// mpz_init_set(mp_temp, mp_Dt[len-1]);
		mp_temp = mp_Dt[len-1];
		for (i = len-2; i >= 0; i--) {
			for (j = idx[i]; j < idx[i+1]; j++) {
				Integer::mul(mp_N[j],mp_temp, mp_N[j]);
				// mpz_mul(mp_N[j], mp_temp, mp_N[j]);
				}
			Integer::mulin(mp_temp,mp_Dt[i]);
			// mpz_mul(mp_temp, mp_temp, mp_Dt[i]);
		}

		// for (i = 0; i < len; i++) { mpz_clear(mp_Dt[i]); } { XFREE(mp_Dt); }
		// for (i = 0; i < k; i++) { mpz_clear(mp_C[i]); } { XFREE(mp_C); }
		// { mpz_clear(mp_sum);  mpz_clear(mp_m);  mpz_clear(mp_temp); }
		// { XFREE(C1); XFREE(idx); }
		return 1;
	}


	/*
	 * Calling Sequence:
	 *   1/0 <-- findNumer(mp_u, mp_m, mp_D, mp_nb, mp_N)
	 *
	 * Summary:
	 *   Certify correctness of the input denominator and, upon success, compute
	 *   the numerator
	 *
	 * Description:
	 *   Given a possible denominator mp_D, numerator bound mp_nb, the function
	 *   certifies the denominator by computing mp_N = mp_u*mp_D mod mp_m,
	 *   if abs(mp_N) < mp_nb, then mp_D is the real denominator and mp_N is the
	 *   real numerator. In this case, the function returns 1. Otherwise, the
	 *   function returns 0.
	 *
	 * Input:
	 *    mp_u: Integer, image
	 *    mp_m: Integer, modulus
	 *    mp_D: Integer, the possible input denominator
	 *   mp_nb: Integer, numerator bound of rational number
	 *
	 * Output:
	 *   mp_N: Integer, if the function returns 1, then mp_N stores the numerator
	 *         of the reconstructed rational number
	 *
	 * Return:
	 *   -1, if mp_N is the real numerator
	 *   -0, if mp_N is not the real numerator
	 *
	 * Note:
	 *   mp_N will be modified no matter it is really numerator or not.
	 *
	 */

template<class FiniteField>
	size_t
	SolutionReconstruct<FiniteField>::
	findNumer (const Integer &mp_u
		   , const Integer &mp_m
		   , const Integer &mp_D
		   , const Integer &mp_nb
		   , Integer &mp_N)
	{
		Integer::mul(mp_N,mp_U,mp_D);
		// mpz_mul(mp_N, mp_u, mp_D);
		Ingeger::modin(mp_N,mp_m);
		// mpz_mod(mp_N, mp_N, mp_m);
		if ( absCompare(mp_N,mp_nb)>0 ) {
			return 0;
		}
		else {
			return 1;
		}
	}


} // iml
} // LinBox

#endif // __LINBOX_algorithm_iml_reconstruct_solution_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

