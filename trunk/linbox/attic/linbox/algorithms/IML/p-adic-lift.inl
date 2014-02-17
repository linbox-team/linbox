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
	template<class Field>
	template<class Matrix>
	int
	pAdicLift<Field>::liftInit (//const size_t liftbasislen
		  // , const BlasVector<NoField> & liftbasis
		  // , const size_t n
		  const Matrix &A
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
		size_t liftbasislen = _liftbasis.size();

		typedef typename Matrix::Ring Ring ;
		size_t i, j, minv, p, temp, len=0;
		typename Ring::Element alpha;
		mpz_t mp_maxInter,  mp_temp;
		// ModElement *q, *qinv;

		for (i = 0; i < liftbasislen; i++) {
			// p = (size_t)liftbasis[i];
			Field Fp(_liftbasis[i]);
			for (j = 0; j < n*n; j++)
				Fp.init(_AInv[i].getWritePointer()[j],A.getPointer()[j]);
				// AInv[i][j] = (double)((temp = (A[j] % p)) >= 0 ? temp : p+temp);
			_AInv[i].changeField(Fp);
			minv = mInverseIn(_AInv[i]);
			// minv = mInverse(liftbasis[i], AInv[i], n);

			/* if fail to find inverse of A mod liftbasis[i] */
			if (minv == 0) { return (int)i; }
		}
		_liftbasis.combBasis();
		// *cmbasis = combBasis(liftbasislen, liftbasis);
		Integer & mp_basisprod = _liftbasis.basisProd();
		basisProd(_liftbasis,mp_basisprod);
		// basisProd(liftbasislen, liftbasis, mp_basisprod);

		/* compute maximum intermediate result mp_maxInter */
		magnitude(alpha,A);

		// alpha = maxMagnLong(A, n, n, n);

		Integer mp_alpha = alpha;
		// mpz_init_set_ui(mp_alpha, alpha);
		// mpz_init(mp_maxInter);
		_extbasis.maxExtInter(mp_alpha, n, mp_maxInter);
		// mpz_clear(mp_alpha);

		_extbasis.setUp(_liftbasis[liftbasislen-1]-1, mp_maxInter);
		_extbasis.findRNS();
		// *extbasis = findRNS(liftbasis[liftbasislen-1]-1, mp_maxInter, &len);
		// mpz_clear(mp_maxInter);
		// *extbasislen = len;
		typename RNS<Field>::ModVect & q    = _extbasis.refRNSbasis();
		typename RNS<Field>::ModVect & qinv = _extbasis.refRNScombi();
		// q = *(*extbasis);
		// qinv = *((*extbasis)+1);
		invBasis(_liftbasisInv, q, mp_basisprod);
		// done in setUp
		// basisProd(len, q, mp_extbasisprod);
		// *extbdcoeff = repBound(len, q, qinv);
		// *ARNS = XMALLOC(Double *, len);
		const BlasMatrix<typename RNS<Field>::NoField> Z(_extbasis.unF,n,n);
		ARNS.resize(len,Z);
		for (i = 0; i < len; i++) {
			Field Fq(q[i]);
			// p = (size_t)q[i];
			// (*ARNS)[i] = XMALLOC(Double, n*n);
			for (j = 0; j < n*n; j++)
				Fq.init(ARNS[i].getWritePointer()[j],A.getPointer()[j]);
			// (*ARNS)[i][j] = (double)((temp = (A[j] % p)) >= 0 ?  temp : p+temp);
		}
		return -1;
	}


	template<class Field>
	int
	pAdicLift<Field>::liftInit(const RNSMatrix<Field> & AinRNS)
	{
		// long i, j, alpha, minv, len=0;
		// double dtemp, *cumprodRNS;
		// mpz_t mp_maxInter, mp_alpha;
		// FiniteField *q, *qinv, *cmbasis, *bdcoeff;

		RNS<Field> & AbaseRNS = AinRNS.basis ;
		std::vector<BlasMatrix<Field> > & AmatRNS = AinRNS.matRNS ;

		AbaseRNS.combBasis();
		// cmbasis = combBasis(basislen, basis);
		AbaseRNS.repBound();
		// bdcoeff = repBound(basislen, basis, cmbasis);
		AbaseRNS.cumProd(_liftbasis);
		// cumprodRNS = cumProd(basislen, basis, liftbasislen, liftbasis);
		for (i = 0; i < size() ; i++)
		{
			/* compute A mod liftbasis[i] from AbaseRNS */
			Field Fi(_liftbasis[i])
			AInv[i].changeField(Fi); // XXX needed ?
			AbaseRNS.basisExt(Fi,AbaseRNS,AInv[i]);
			// basisExt(basislen, n*n, liftbasis[i], basis, cmbasis, cumprodRNS[i],
				 // bdcoeff, AbaseRNS, AInv[i]);

			int minv = mInversein(AInv[i], n);
			// minv = mInverse(liftbasis[i], AInv[i], n);

			/* if fail to find inverse of A mod basis[i] */
			if (minv == 0)
			{
				// XFREE(bdcoeff); XFREE(cmbasis); XFREE(cumprodRNS);
				return i;
			}
		}

		// XFREE(cumprodRNS);
		_liftbasis.combBasis();
		// *cmliftbasis = combBasis(liftbasislen, liftbasis);
		_liftbasis.basisProd();
		// basisProd(liftbasislen, liftbasis, mp_liftbasisprod);

		/* compute maximum intermediate result mp_maxInter */
		// mpz_init(mp_alpha);
		//! @bug need to set it
		magnitude(mp_alpha,AbaseRNS);
		// basisProd(basislen, basis, mp_alpha);
		// maxExtInter(mp_alpha, n, mp_maxInter);
		// mpz_clear(mp_alpha);


		// *extbasis = findRNS(liftbasis[liftbasislen-1]-1, mp_maxInter, &len);
		_extbasis.setUp(liftbasis[liftbasislen-1]-1, mp_maxInter);
		_extbasis = _liftbasis.findRNS();
		// mpz_clear(mp_maxInter);
		// *extbasislen = len;
		typename RNS<Field>::ModVect & q    = _extbasis.refRNSbasis();
		typename RNS<Field>::ModVect & qinv = _extbasis.refRNScombi();

		// q = *(*extbasis);
		// qinv = *((*extbasis)+1);
		invBasis(_liftbasisInv, q, mp_basisprod);
		// *liftbasisInv = invBasis(len, q, mp_liftbasisprod);
		// basisProd(len, q, mp_extbasisprod);
		// *extbdcoeff = repBound(len, q, qinv);

		const BlasMatrix<typename RNS<Field>::NoField> Z(_extbasis.unF,n,n);
		ARNS.resize(len,Z);

		// *AExtRNS = XMALLOC(Double *, len);
		// cumprodRNS = cumProd(basislen, basis, len, q);
		for (i = 0; i < len; i++)
		{
			// (*AExtRNS)[i] = XMALLOC(Double, n*n);

			Field Fq(q[i]);
			/* compute A mod extbasis[i] from ARNS */
			// basisExt(basislen, n*n, q[i], basis, cmbasis, cumprodRNS[i], \
				 // bdcoeff, ARNS, (*AExtRNS)[i]);
			AbaseRNS.basisExt(Fq,AinRNS,ARNS[i]);
		}

		// { XFREE(bdcoeff); XFREE(cmbasis); XFREE(cumprodRNS); }

		return -1;
	}



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

	template<class Field>
	void
	pAdicLift<Field>::iml_lift (const Tag::Side solupos
		  // , std::vector<std::vector<BlasMatrix<Field> > > & C
				    , LiftStep<Field> & C
		  // , const size_t k
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
		// mpz_t mp_r1;
		// ModElement *q, *qinv;
		// Double *dtemp, *dtemp1, **Ac, ***C;
		typedef BlasVector<NoField> ModVect;

		/* initialize lifting coefficient matrix C[k][liftbasislen][n] */
		// C = XMALLOC(Double **, k);
		size_t n = C.rowdim() ;
		size_t m = C.coldim() ;
		size_t k = C.steps();
		size_t liftbasislen = _liftbasis.size();
		size_t extbasislen  = _extbasis. size();
		Field F(); // XXX
		const BlasMatrix<Field> Z(F,n,m);
		// const std::vector<BlasMatrix<Field> > ZZ(liftbasislen,Z);
		// C.resize(k,ZZ);
		// for (i = 0; i < k; i++)
		// {
			// C[i] = XMALLOC(Double *, liftbasislen);
			// for (j = 0; j < liftbasislen; j++)
				// C[i][j] = XMALLOC(Double, m*n);
		// }
		ModVect & q    = _extbasis.refRNSbasis() ;
		ModVect & qinv = _extbasis.refRNScombi() ;
		// mpz_init(mp_r1);
		std::vector<BlasMatrix<Field> > Ac(extbasislen,Z);
		for (size_t fi = 0 ; fi < extbasislen ; ++fi) {
			Field Fp(_extbasis[fi]);
			Ac[fi].changeField(Fp);
		}
		// Ac = XCALLOC(Double *, extbasislen);
		// for (i = 0; i < extbasislen; i++)
			// Ac[i] = XCALLOC(Double, m*n);
		NoField unF;
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
				Field Fq(_liftbasis[l]);
				/* mod(mp_r, liftbasis[j]) */
				for (j = 0; j < m*n; j++) {
					//!@todo use hom/rebind
					Fq.init(dtemp.getWritePointer()[j],mp_r.getPointer()[j]);
					// dtemp[j] = (Double)mpz_fdiv_ui(mp_r[j], liftbasis[l]);
				}

				/* compute the coefficients of p-adic lifting */
				FFLAS::fgemm(Fq,(solupos==Tag::Left)?(Tag::Trans):(Tag::NoTrans)
					     ,Tag::NoTrans,n,m,n
					     ,Fq.one,_AInv[l].getPointer(),n,dtemp.getPointer(),m
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
				Field Fq(q[j]);
				_liftbasis.basisExtPos(Fq, C[i], dtemp);
				// _liftbasis.basisExtPos(liftbasislen, m*n, q[j], _liftbasis, cmbasis, C[i], dtemp);

				FFLAS::fgemm(Fq,(solupos==Tag::Left)?(Tag::Trans):(Tag::NoTrans)
					     ,Tag::NoTrans,n,m,n
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
			for (j = 0; j < m*n; j++) {
				Integer mp_r1;
				/* mp_r[j] := Quo(mp_r[j], p), mp_r1 := Mod(mp_r[j], p) */
				Integer::divmod(mp_r[j],mp_r1,mp_r[j],_liftbasis.basisProd());
				// mpz_fdiv_qr(mp_r[j], mp_r1, mp_r[j], mp_basisprod);

				/* compute ((r mod p) mod q[l] - Ac mod q[l])(1/p mod q[l]) */
				for (l = 0; l < extbasislen; l++) {
					Field Fp(q[l]);
					dtemp1[l]= Integer::frem(mp_r1,(long)q[l]);
					// dtemp1[l] = (Double)mpz_fdiv_ui(mp_r1, q[l]);
					Fp.addin(dtemp1[l],(q[l]-1)*Ac[l][j]);
					// dtemp1[l] = fmod(dtemp1[l]+(q[l]-1)*Ac[l][j], q[l]);
					Fp.mulin(dtemp1[l],_liftbasisInv[l]);
					// dtemp1[l] = fmod(dtemp1[l]*liftbasisInv[l], q[l]);
				}

				/* compute (r mod p-Ac)(1/p) by CRT */
				_extbasis.ChineseRemainder(dtemp1,mp_r1);
				// ChineseRemainder(extbasislen, mp_extbasisprod, q, qinv, extbdcoeff,
				// dtemp1, mp_r1);
				// mpz_add(mp_r[j], mp_r[j], mp_r1);
				Integer::addin(mp_r.getPointer()[j],mp_r1);
			}
		}

		// mpz_clear(mp_r1);
		// { XFREE(dtemp); XFREE(dtemp1); }
		// for (i = 0; i < extbasislen; i++) { XFREE(Ac[i]); } { XFREE(Ac); }

		return C;
	}


	/* return (1/mp_basisprod) mod basis[i] */

	template<class Field>
	void
	pAdicLift<Field>::invBasis( )
	{
		size_t basislen = _extbasis.size();
		size_t i;
		// mpz_t mp_temp, mp_basis;
		// Double *inv;

		// { mpz_init(mp_temp); mpz_init(mp_basis); }
		// inv = XMALLOC(Double, basislen);
		_liftbasisInv .resize(basislen);
		for (i = 0; i < basislen; i++) {
			Field Fq(_extbasis.prime(i));
			Fq.inv(_liftbasisInv[i],_extbasis.basisProd());
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


