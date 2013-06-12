#ifndef __LINBOX_algorithm_iml_p_adic_lift_INL
#define __LINBOX_algorithm_iml_p_adic_lift_INL

namespace LinBox{ namespace iml{



	/*
	 * Calling Sequence:
	 *   -1/i <-- liftInit(liftbasislen, liftbasis, n, A, mp_basisprod,
	 *            mp_extbasisprod, extbasislen, cmbasis, extbdcoeff,
	 *            liftbasisInv, AInv, extbasis, ARNS)
	 *
	 * Summary:
	 *   Perform initialization operations before lifting, where the left hand
	 *   side input matrix is a signed size_t matrix
	 *
	 * Description:
	 *   Initialization step before lifting computes necessory data and stores
	 *   them to avoid recomputations if continuous lifting is needed. Seperating
	 *   initialization step and lifting step makes it possible to lift adaptively.
	 *   Pointers are passed into the calling sequence to store the outputs of
	 *   initialization.
	 *
	 * Input:
	 *   liftbasislen: size_t, dimension of lifting basis
	 *      liftbasis: 1-dim ModElement array length liftbasislen, lifting basis
	 *              n: size_t, dimension of input matrix A
	 *              A: 1-dim size_t array length n*n, representation of n x n input
	 *                 matrix
	 *
	 * Return:
	 *   The first index i such that A^(-1) mod liftbasis[i] doesn't exist, where
	 *   i starts from 0. Otherwise, return -1 if A^(-1) mod liftbasis[i] exists
	 *   for any i
	 *
	 * Output:
	 *     mp_basisiprod: mpz_t, product of lifting basis
	 *   mp_extbasisprod: mpz_t, product of extended RNS basis
	 *       extbasislen: pointer to size_t int, storing the dimension of extended
	 *                    RNS basis. Extended basis is used to compute the
	 *                    matrix-vector product AC_i during lifting
	 *           cmbasis: pointer to a 1-dim ModElement array, storing the
	 *                    special combination of lifting basis computed by
	 *                    function combBasis
	 *        extbdcoeff: pointer to a 1-dim ModElement array, storing the mix
	 *                    radix coefficients of a special integer, computed by
	 *                    function repBound, in extended RNS basis
	 *      liftbasisInv: pointer to a 1-dim Double array, storing
	 *                    (1/mp_basisprod) mod extbasis[i]
	 *              AInv: pointer to a 1-dim Double array, storing the modp matrix
	 *                    A^(-1) mod liftbasis[i]
	 *          extbasis: pointer to a 2-dim ModElement array, where
	 *                  - (*extbasis)[0] = extended RNS basis
	 *                  - (*extbasis)[1] = the special combination of extended RNS
	 *                    basis computed by function combBasis
	 *              ARNS: pointer to a 2-dim Double array, where the last
	 *                    dimension (*ARNS)[i] stores A mod ith element in
	 *                    extended RNS basis
	 */



	//! @todo liftbasis/cmbasis in a RNS ?
	int
	liftInit (//const size_t liftbasislen
		  // , const BlasVector<NoField> & liftbasis
		  // , const size_t n
		  , const BlasMatrix<Ring> &A
		  // , Integer & mp_basisprod
		  // , Integer & mp_extbasisprod
		  // , size_t &extbasislen
		  // , BlasVector<NoField> &cmbasis
		  // , BlasVector<NoField> &extbdcoeff
		  // , BlasVector<NoField> &liftbasisInv
		  // , std::vector<BlasVector<Field> > &AInv
		  // , BlasVector<NoField> &RNSbasis
		  // , BlasVector<NoField> &RNScombi
		  // , RNS<Field> & extbasis
		  // , BlasVector<NoField> &ARNS
		  )
	{
		size_t n = A.coldim();
		size_t liftbasislen = liftbasis.size();

		size_t i, j, minv, p, temp, len=0;
		typename Ring::Element alpha;
		mpz_t mp_maxInter, mp_alpha, mp_temp;
		ModElement *q, *qinv;

		for (i = 0; i < liftbasislen; i++) {
			// p = (size_t)liftbasis[i];
			Field Fp(liftbasis[i]);
			for (j = 0; j < n*n; j++)
				Fp.init(AInv[i].getWritePointer()[j],A.getPointer()[j]);
				// AInv[i][j] = (double)((temp = (A[j] % p)) >= 0 ? temp : p+temp);
			AInv[i].changeField(Fp);
			minv = mInverseIn(AInv[i]);
			// minv = mInverse(liftbasis[i], AInv[i], n);

			/* if fail to find inverse of A mod liftbasis[i] */
			if (minv == 0) { return i; }
		}
		extbasis.combBasis(cmbasis,liftbasis);
		// *cmbasis = combBasis(liftbasislen, liftbasis);
		basisProd(liftbasis,mp_basisprod);
		// basisProd(liftbasislen, liftbasis, mp_basisprod);

		/* compute maximum intermediate result mp_maxInter */
		BlasMatrixDomain BMD(A.field());
		BMD.Magnitude(alpha,A);
		// alpha = maxMagnLong(A, n, n, n);

		Integer mp_alpha = alpha;
		// mpz_init_set_ui(mp_alpha, alpha);
		// mpz_init(mp_maxInter);
		extbasis.maxExtInter(mp_alpha, n, mp_maxInter);
		// mpz_clear(mp_alpha);

		extbasis.setUp(liftbasis[liftbasislen-1]-1, mp_maxInter);
		// *extbasis = findRNS(liftbasis[liftbasislen-1]-1, mp_maxInter, &len);
		// mpz_clear(mp_maxInter);
		// *extbasislen = len;
		typename RNS<Field>::ModVect & q    = extbasis.refRNSbasis();
		typename RNS<Field>::ModVect & qinv = extbasis.refRNScombi();
		// q = *(*extbasis);
		// qinv = *((*extbasis)+1);
		invBasis(liftbasisInv, q, mp_basisprod);
		// done in setUp
		// basisProd(len, q, mp_extbasisprod);
		// *extbdcoeff = repBound(len, q, qinv);
		// *ARNS = XMALLOC(Double *, len);
		const BlasMatrix<typename RNS<Field>::NoField> Z(extbasis.unF,n,n);
		ARS.resize(len,Z);
		for (i = 0; i < len; i++) {
			FiniteField Fp(q[i]);
			// p = (size_t)q[i];
			// (*ARNS)[i] = XMALLOC(Double, n*n);
			for (j = 0; j < n*n; j++)
				Fq.init(ARNS[i].getWritePointer()[j],A.getPointer()[j]);
			// (*ARNS)[i][j] = (double)((temp = (A[j] % p)) >= 0 ?  temp : p+temp);
		}
		return -1;
	}


#if 0 /* redundant */
	/*
	 * Calling Sequence:
	 *   -1/i <-- liftInitLlhs(liftbasislen, liftbasis, n, mp_A, mp_basisprod,
	 *                         mp_extbasisprod, extbasislen, cmbasis, extbdcoeff,
	 *                         liftbasisInv, AInv, extbasis, ARNS)
	 *
	 * Summary:
	 *   Perform initialization operations before lifting, where the left hand
	 *   side input matrix is a mpz_t matrix
	 *
	 * Description:
	 *   Initialization step before lifting computes necessory data and stores
	 *   them to avoid recomputations if continuous lifting is needed. Seperating
	 *   initialization step and lifting step makes it possible to lift adaptively.
	 *   Pointers are passed into the calling sequence to store the outputs of
	 *   initialization.
	 *
	 * Input:
	 *   liftbasislen: size_t, dimension of lifting basis
	 *      liftbasis: 1-dim ModElement array length liftbasislen, lifting basis
	 *              n: size_t, dimension of input matrix A
	 *           mp_A: 1-dim mpz_t array length n*n, representation of n x n input
	 *                 matrix
	 *
	 * Return:
	 *   The first index i such that A^(-1) mod liftbasis[i] doesn't exist, where
	 *   i starts from 0. Otherwise, return -1 if A^(-1) mod liftbasis[i] exists
	 *   for any i
	 *
	 * Output:
	 *     mp_basisiprod: mpz_t, product of lifting basis
	 *   mp_extbasisprod: mpz_t, product of extended RNS basis
	 *       extbasislen: pointer to size_t int, storing the dimension of extended
	 *                    RNS basis. Extended basis is used to compute the
	 *                    matrix-vector product AC_i during lifting
	 *           cmbasis: pointer to a 1-dim ModElement array, storing the
	 *                    special combination of lifting basis computed by
	 *                    function combBasis
	 *        extbdcoeff: pointer to a 1-dim ModElement array, storing the mix
	 *                    radix coefficients of a special integer, computed by
	 *                    function repBound, in extended RNS basis
	 *      liftbasisInv: pointer to a 1-dim Double array, storing
	 *                    (1/mp_basisprod) mod extbasis[i]
	 *              AInv: pointer to a 1-dim Double array, storing the modp matrix
	 *                    A^(-1) mod liftbasis[i]
	 *          extbasis: pointer to a 2-dim ModElement array, where
	 *                  - (*extbasis)[0] = extended RNS basis
	 *                  - (*extbasis)[1] = the special combination of extended RNS
	 *                    basis computed by function combBasis
	 *              ARNS: pointer to a 2-dim Double array, where the last
	 *                    dimension (*ARNS)[i] stores A mod ith element in
	 *                    extended RNS basis
	 */

	size_t
	liftInitLlhs (//const size_t liftbasislen
		      , const ModElement *liftbasis
		      // , const size_t n
		      , const BlasMatrix<PID_integer> &mp_A
		      , mpz_t mp_basisprod
		      // , mpz_t mp_extbasisprod
		      // , size_t *extbasislen
		      // , ModElement **cmbasis
		      , BlasVector<NoField> &cmbasis
		      // , ModElement **extbdcoeff
		      , BlasVector<NoField> &liftbasisInv
		      // , Double **liftbasisInv
		      , std::vector<BlasVector<Field> > &AInv
		      // , Double **AInv
		      // , ModElement ***extbasis
		      , RNS<Field> & extbasis
		      // , Double ***ARNS
		      , BlasVector<NoField> &ARNS)
	{
		size_t i, j, minv, len=0;
		double dtemp;
		Integer mp_maxInter, mp_alpha;
		// ModElement *q, *qinv;

		for (i = 0; i < liftbasislen; i++) {
			FiniteField Fp(liftbasis[i]);
			for (j = 0; j < n*n; j++) {
				Fp.init(AInv[i].getWritePointer()[j],mp_A.getPointer()[j]);
				// AInv[i][j] = (Double)mpz_fdiv_ui(mp_A[j], liftbasis[i]);
			}
			AInv[i].changeField(Fp);
			minv = mInverse(AInv[i]);

			/* if fail to find inverse of A mod liftbasis[i] */
			if (minv == 0) { return i; }
		}
		extbasis.combBasis(cmbasis,liftbasis);
		// *cmbasis = combBasis(liftbasislen, liftbasis);
		basisProd(liftbasis,mp_basisprod);
		// basisProd(liftbasislen, liftbasis, mp_basisprod);

		/* compute maximum intermediate result mp_maxInter */
		BlasMatrixDomain BMD(A.field());
		BMD.Magnitude(mp_alpha,A);
		// mpz_init(mp_alpha);
		// maxMagnMP(mp_A, n, n, n, mp_alpha);
		extbasis.maxExtInter(mp_alpha, n, mp_maxInter);
		// mpz_init(mp_maxInter);
		// maxExtInter(mp_alpha, n, mp_maxInter);
		// mpz_clear(mp_alpha);

		extbasis.setUp(liftbasis[liftbasislen-1]-1, mp_maxInter);
		// *extbasis = findRNS(liftbasis[liftbasislen-1]-1, mp_maxInter, &len);
		// mpz_clear(mp_maxInter);
		// *extbasislen = len;
		// q = *(*extbasis);
		qinv = *((*extbasis)+1);
		*liftbasisInv = invBasis(len, q, mp_basisprod);
		basisProd(len, q, mp_extbasisprod);
		*extbdcoeff = repBound(len, q, qinv);
		*ARNS = XMALLOC(Double *, len);
		for (i = 0; i < len; i++)
		{
			(*ARNS)[i] = XMALLOC(Double, n*n);
			for (j = 0; j < n*n; j++)
				(*ARNS)[i][j] = (Double)mpz_fdiv_ui(mp_A[j], q[i]);
		}

		return -1;
	}
#endif

#if 0 /* later */
	/*
	 * Calling Sequence:
	 * -1/i <-- liftInitRNS(liftbasislen, liftbasis, basislen, basis, n, ARNS,
	 *          mp_liftbasisprod, mp_extbasisprod, extbasislen, cmliftbasis,
	 *	    liftbasisInv, extbdcoeff, AInv, extbasis, AExtRNS)
	 *
	 * Summary:
	 *   Perform initialization operations before lifting where the left hand
	 *   side input matrix is represented in a RNS
	 *
	 * Description:
	 *   Initialization step before lifting computes necessory data and stores
	 *   them to avoid recomputations if continuous lifting is needed. Seperating
	 *   initialization step and lifting step makes it possible to lift adaptively.
	 *   Pointers are passed into the calling sequence to store the outputs of
	 *   initialization.
	 *
	 * Input:
	 *   liftbasislen: size_t, dimension of lifting basis
	 *      liftbasis: 1-dim ModElement array length liftbasislen, lifting basis
	 *       basislen: size_t, dimension of RNS basis used to represent A
	 *          basis: 1-dim ModElement array length basislen, RNS basis used to
	 *                 represent A
	 *              n: size_t, dimension of A
	 *           ARNS: 2-dim Double array, dimension basislen x n^2,
	 *                 representation of A in RNS, ARNS[i] = A mod basis[i]
	 *
	 * Return:
	 *   The first index i such that A^(-1) mod liftbasis[i] doesn't exist, where
	 *   i starts from 0. Otherwise, return -1 if A^(-1) mod liftbasis[i] exists
	 *   for any i
	 *
	 * Output:
	 *   mp_liftbasisiprod: mpz_t, product of lifting basis
	 *     mp_extbasisprod: mpz_t, product of extended RNS basis
	 *         extbasislen: pointer to size_t int, storing the dimension of extended
	 *                      RNS basis. Extended basis is used to compute the
	 *                      matrix-vector product AC_i during lifting
	 *             cmbasis: pointer to a 1-dim ModElement array, storing the
	 *                      special combination of lifting basis computed by
	 *                      function combBasis
	 *          extbdcoeff: pointer to a 1-dim ModElement array, storing the mix
	 *                      radix coefficients of a special integer, computed by
	 *                      function repBound, in extended RNS basis
	 *      liftbasisInv: pointer to a 1-dim Double array, storing
	 *                    (1/mp_basisprod) mod extbasis[i]
	 *                AInv: pointer to a 1-dim Double array, storing the modp
	 *                      matrix A^(-1) mod liftbasis[i]
	 *            extbasis: pointer to a 2-dim ModElement array, where
	 *                    - (*extbasis)[0] = extended RNS basis
	 *                    - (*extbasis)[1] = the special combination of extended
	 *                      RNS basis computed by function combBasis
	 *                ARNS: pointer to a 2-dim Double array, where the last
	 *                      dimension (*ARNS)[i] stores A mod ith element in
	 *                      extended RNS basis
	 */

	size_t
	liftInitRNS (const size_t liftbasislen
		     , const ModElement *liftbasis
		     , const size_t basislen
		     , const ModElement *basis
		     , const size_t n
		     , Double **ARNS
		     , mpz_t mp_liftbasisprod
		     , mpz_t mp_extbasisprod
		     , size_t *extbasislen
		     , ModElement **cmliftbasis
		     , ModElement **extbdcoeff
		     , Double **liftbasisInv
		     , Double **AInv
		     , ModElement ***extbasis
		     , Double ***AExtRNS
		     )
	{
		size_t i, j, alpha, minv, len=0;
		double dtemp, *cumprodRNS;
		mpz_t mp_maxInter, mp_alpha;
		ModElement *q, *qinv, *cmbasis, *bdcoeff;

		cmbasis = combBasis(basislen, basis);
		bdcoeff = repBound(basislen, basis, cmbasis);
		cumprodRNS = cumProd(basislen, basis, liftbasislen, liftbasis);
		for (i = 0; i < liftbasislen; i++)
		{
			/* compute A mod liftbasis[i] from ARNS */
			basisExt(basislen, n*n, liftbasis[i], basis, cmbasis, cumprodRNS[i], \
				 bdcoeff, ARNS, AInv[i]);
			minv = mInverse(liftbasis[i], AInv[i], n);

			/* if fail to find inverse of A mod basis[i] */
			if (minv == 0)
			{ XFREE(bdcoeff); XFREE(cmbasis); XFREE(cumprodRNS); return i; }
		}
		XFREE(cumprodRNS);
		*cmliftbasis = combBasis(liftbasislen, liftbasis);
		basisProd(liftbasislen, liftbasis, mp_liftbasisprod);

		/* compute maximum intermediate result mp_maxInter */
		mpz_init(mp_alpha);
		basisProd(basislen, basis, mp_alpha);
		mpz_init(mp_maxInter);
		maxExtInter(mp_alpha, n, mp_maxInter);
		mpz_clear(mp_alpha);

		*extbasis = findRNS(liftbasis[liftbasislen-1]-1, mp_maxInter, &len);
		mpz_clear(mp_maxInter);
		*extbasislen = len;
		q = *(*extbasis);
		qinv = *((*extbasis)+1);
		*liftbasisInv = invBasis(len, q, mp_liftbasisprod);
		basisProd(len, q, mp_extbasisprod);
		*extbdcoeff = repBound(len, q, qinv);
		*AExtRNS = XMALLOC(Double *, len);
		cumprodRNS = cumProd(basislen, basis, len, q);
		for (i = 0; i < len; i++)
		{
			(*AExtRNS)[i] = XMALLOC(Double, n*n);

			/* compute A mod extbasis[i] from ARNS */
			basisExt(basislen, n*n, q[i], basis, cmbasis, cumprodRNS[i], \
				 bdcoeff, ARNS, (*AExtRNS)[i]);
		}

		{ XFREE(bdcoeff); XFREE(cmbasis); XFREE(cumprodRNS); }

		return -1;
	}
#endif



	/*
	 * Calling Sequence:
	 *   C <-- lift(solupos, k, n, m, liftbasislen, extbasislen, mp_basisprod,
	 *              mp_extbasisprod, liftbasis, cmbasis, extbdcoeff,
	 *              liftbasiInv, mp_r, extbasis, AInv, ARNS)
	 *
	 * Summary:
	 *   Compute p-adic lifting coefficients of system of linear equations
	 *
	 * Description:
	 *   Given a system of linear equations AX = mp_r or Transpose(A)X = mp_r,
	 *   where A is a n x n nonsingular matrix and mp_r is a n x m matrix, the
	 *   function computes and stores lifting coefficients upto k lifting steps.
	 *
	 *   The data computed from initialization function are reused each time this
	 *   function is called. The right hand side matrix mp_r is updated such that
	 *   we could continue to call this function to perform lifting using updated
	 *   mp_r if we do not lift high enough.
	 *
	 * Input:
	 *           solupos: enumerate, flag to indicate whether to transpose A or not
	 *                  - solupos = LeftSolu: system be Transpose(A)X = mp_r
	 *                  - solupos = RightSolu: system be AX = mp_r
	 *                 k: size_t, number of lifting steps
	 *                 n: size_t, dimension of A
	 *                 m: size_t, column dimension of right hand side matrix mp_r
	 *      liftbasislen: size_t, dimension of lifting basis
	 *       extbasislen: size_t, dimension of extended RNS basis
	 *     mp_basisiprod: mpz_t, product of lifting basis
	 *   mp_extbasisprod: mpz_t, product of extended lifting basis
	 *         liftbasis: 1-dim ModElement array length liftbasislen, storing
	 *                    lifting basis
	 *           cmbasis: 1-dim ModElement array length liftbasislen, storing the
	 *                    special combination of lifting basis computed in
	 *                    initialization step
	 *        extbdcoeff: 1-dim ModElement array length liftbasislen, storing the
	 *                    mix radix coefficients of a special integer in extended
	 *                    RNS basis computed in initialization step
	 *      liftbasisInv: a 1-dim Double array, storing
	 *                    (1/mp_basisprod) mod extbasis[i]
	 *              mp_r: 1-dim mpz_t array length n*m, representation of n x m
	 *                    right hand side lifting matrix
	 *          extbasis: 2-dim ModElement array, dimension 2 x extbasislen,
	 *                    computed in initialization step
	 *              AInv: 1-dim Double array length liftbasislen, storing the modp
	 *                    matrix A^(-1) mod liftbasis[i]
	 *              ARNS: 2-dim Double array, dimension extbasislen x n^2, where
	 *                    the second dimension (*ARNS)[i] stores A mod ith element
	 *                    in extended RNS basis
	 *
	 * Output:
	 *   C: 3-dim Double array, dimension k x liftbasislen x n*m
	 *    - If solupos = RightSolu, then C[i][j] represents the n x m coefficient
	 *      matrix computed by A^(-1)mp_r mod liftbasis[j] at the ith lifting step.
	 *    - If solupos = LeftSolu, then C[i][j] represents the n x m coefficient
	 *      matrix computed by Transpose(A)^(-1)mp_r mod liftbasis[j] at the ith
	 *      lifting step.
	 *
	 * Precondition:
	 *   Any element p in array liftbasis must satisfy n*(p-1)^2 <= 2^53-1.
	 */

	void
	iml_lift (const enum SOLU_POS solupos
		  , std::vector<std::vector<BlasMatrix<Field> > > & C
		  , const size_t k
		  // , const size_t n
		  // , const size_t m
		  // , const size_t liftbasislen
		  // , const size_t extbasislen
		  // , const Integer & mp_basisprod
		  // , const Integer & mp_extbasisprod
		  // , const ModElement *liftbasis
		  // , const ModElement *cmbasis
		  // , const ModElement *extbdcoeff
		  // , const BlasVector<NoField> & liftbasisInv
		  , BlasMatrix<PID_integer> & mp_r
		  // , ModElement **extbasis
		  // , RNS<Field> & liftbasis
		  // , RNS<Field> & extbasis
		  // , std::vector<BlasMatrix<Field> > & AInv
		  // , std::vector<BlasMatrix<Field> > & ARNS
		  )
	{
		size_t i, j, l;
		mpz_t mp_r1;
		ModElement *q, *qinv;
		Double *dtemp, *dtemp1, **Ac, ***C;

		/* initialize lifting coefficient matrix C[k][liftbasislen][n] */
		// C = XMALLOC(Double **, k);
		const BlasMatrix<Field> Z(F,n,m);
		const std::vector<BlasMatrix<Field> > ZZ(liftbasislen,Z);
		C.resize(k,ZZ);
		// for (i = 0; i < k; i++)
		// {
			// C[i] = XMALLOC(Double *, liftbasislen);
			// for (j = 0; j < liftbasislen; j++)
				// C[i][j] = XMALLOC(Double, m*n);
		// }
		ModVect & q    = extbasis.refRNSbasis() ;
		ModVect & qinv = extbasis.refRNScombi() ;
		// mpz_init(mp_r1);
		std::vector<BlasMatrix<Field> > Ac(extbasislen,Z);
		// Ac = XCALLOC(Double *, extbasislen);
		// for (i = 0; i < extbasislen; i++)
			// Ac[i] = XCALLOC(Double, m*n);
		BlasMatrix<NoField> dtemp(unF,n,m);
		BlasVector<NoField> dtemp1(unF,extbasislen);
		// dtemp = XMALLOC(Double, m*n);
		// dtemp1 = XMALLOC(Double, extbasislen);

		/* start lifting */
		for (i = 0; i < k; i++)
		{
			/* compute coefficients of p-adic lifting C[i][l] */
			for (l = 0; l < liftbasislen; l++)
			{
				Field Fq(liftbasis[l]);
				/* mod(mp_r, liftbasis[j]) */
				for (j = 0; j < m*n; j++) {
					//!@todo use hom/rebind
					Fq.init(dtemp.getWritePointer()[j],mp_r.getPointer()[j]);
					// dtemp[j] = (Double)mpz_fdiv_ui(mp_r[j], liftbasis[l]);
				}

				/* compute the coefficients of p-adic lifting */
				FFLAS::fgemm(Fq,(solupos==LinBoxTag::Left)?(LinBoxTag::Trans):(LinBoxTag::NoTrans)
					     ,LinBoxTag::NoTrans,n,m,n
					     ,Fq.one,AInv[l].getPointer(),n,dtemp.getPointer(),m
					     ,F.zero,C[i][l].getWritePointer(),m);
#if 0
				if (solupos == LeftSolu)
				{
					if (m == 1)
						cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0, AInv[l], \
							    n, dtemp, 1, 0.0, C[i][l], 1);
					else
						cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, \
							    n, 1.0, AInv[l], n, dtemp, m, 0.0, C[i][l], m);
				}
				else if (solupos == RightSolu)
				{
					if (m == 1)
						cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, AInv[l],\
							    n, dtemp, 1, 0.0, C[i][l], 1);
					else
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, \
							    n, 1.0, AInv[l], n, dtemp, m, 0.0, C[i][l], m);
				}
#endif
				// Dmod((double)liftbasis[l], C[i][l], n, m, m);
			}

			/* compute Ac mod extbasis[j] */
			for (j = 0; j < extbasislen; j++)
			{
				liftbasis.basisExtPos(liftbasislen, m*n, q[j], liftbasis, cmbasis, C[i], dtemp);

				FFLAS::fgemm(Fq,(solupos==LinBoxTag::Left)?(LinBoxTag::Trans):(LinBoxTag::NoTrans)
					     ,LinBoxTag::NoTrans,n,m,n
					     ,Fq.one,ARNS[j].getPointer(),n,dtemp.getPointer(),m
					     ,F.zero,Ac[j].getWritePointer(),m);

#if 0
				if (solupos == LeftSolu)
				{
					if (m == 1)
						cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0, ARNS[j],\
							    n, dtemp, 1, 0.0, Ac[j], 1);
					else
						cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, n, \
							    1.0, ARNS[j], n, dtemp, m, 0.0, Ac[j], m);
				}
				else if (solupos == RightSolu)
				{
					if (m == 1)
						cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, ARNS[j],\
							    n, dtemp, 1, 0.0, Ac[j], 1);
					else
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m,\
							    n, 1.0, ARNS[j], n, dtemp, m, 0.0, Ac[j], m);
				}
				Dmod((double)q[j], Ac[j], n, m, m);
#endif
			}

			/* compute r_quo_p+(r mod p-Ac)/p */
			for (j = 0; j < m*n; j++)
			{
				/* mp_r[j] := Quo(mp_r[j], p), mp_r1 := Mod(mp_r[j], p) */
				Integer::divmod(mp_r[j],mp_r1,mp_r[j],p);
				// mpz_fdiv_qr(mp_r[j], mp_r1, mp_r[j], mp_basisprod);

				/* compute ((r mod p) mod q[l] - Ac mod q[l])(1/p mod q[l]) */
				for (l = 0; l < extbasislen; l++)
				{
					dtemp1[l]= Integer::frem(mp_r1,(long)q[l]);
					// dtemp1[l] = (Double)mpz_fdiv_ui(mp_r1, q[l]);
					Fp.addin(dtemp1[l],(q[l]-1)*Ac[l][j]);
					// dtemp1[l] = fmod(dtemp1[l]+(q[l]-1)*Ac[l][j], q[l]);
					Fp.mulin(dtemp1[l],liftbasisInv[l]);
					// dtemp1[l] = fmod(dtemp1[l]*liftbasisInv[l], q[l]);
				}

				/* compute (r mod p-Ac)(1/p) by CRT */
				extbasis.ChineseRemainder(dtemp1,mp_r1);
				// ChineseRemainder(extbasislen, mp_extbasisprod, q, qinv, extbdcoeff,
				// dtemp1, mp_r1);
				// mpz_add(mp_r[j], mp_r[j], mp_r1);
				Integer::addin(mp_r[j],mp_r1);
			}
		}

		// mpz_clear(mp_r1);
		// { XFREE(dtemp); XFREE(dtemp1); }
		// for (i = 0; i < extbasislen; i++) { XFREE(Ac[i]); } { XFREE(Ac); }

		return C;
	}


	/* return (1/mp_basisprod) mod basis[i] */

	void
	invBasis( )
	{
		size_t basislen = _extbasis.size();
		size_t i;
		// mpz_t mp_temp, mp_basis;
		// Double *inv;

		// { mpz_init(mp_temp); mpz_init(mp_basis); }
		// inv = XMALLOC(Double, basislen);
		_liftbasisInv .resize(basislen);
		for (i = 0; i < basislen; i++) {
			Field Fi(_extbasis.prime(i));
			F.inv(_liftbasisInv[i],_extbasis.basisProd());
			// mpz_set_ui(mp_basis, basis[i]);
			// mpz_invert(mp_temp, mp_basisprod, mp_basis);
			// inv[i] = mpz_get_d(mp_temp);
		}
		// { mpz_clear(mp_temp); mpz_clear(mp_basis); }
		return;
		// return inv;

	}

} // iml
} // LinBox

#endif // __LINBOX_algorithm_iml_p_adic_lift_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


