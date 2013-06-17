/* ---------------------------------------------------------------------
 *
 * -- Integer Matrix Library (IML)
 *    (C) Copyright 2004, 2006 All Rights Reserved
 *
 * -- IML routines -- Version 1.0.1 -- November, 2004
 *
 * Author         : Zhuliang Chen
 * Contributor(s) : Arne Storjohann
 * University of Waterloo -- School of Computer Science
 * Waterloo, Ontario, N2L3G1 Canada
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the IML group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef __LINBOX_algorithms_iml_nullspace_INL
#define __LINBOX_algorithms_iml_nullspace_INL


#include "linbox/matrix/factorized-matrix.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/algorithms/rational-solver.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/integer.h"
#include "linbox/util/debug.h"
// #include "linbox/algorithms/integer-tools.h"



namespace LinBox { namespace iml {

	template<class Ring, class Field>
	void
	Nullspace<Ring,Field>::
	justTry(  )
	{
		komp = true ;
		size_t n = A.rowdim();
		size_t m = A.coldim();
		bool ok = false ;

		std::vector<size_t> P(n+1);
		// std::vector<size_t> rp(n+1);
		rp.resize(n+1);
		// P = XCALLOC(long, n + 1);
		// rp = XCALLOC(long, n + 1);
		// size_t r,s ;
		{
			RandomPrimeIter PrimeGen(19) ;
			unsigned long p =  PrimeGen.random_between(15);
			Field F(p);
			typedef typename Field::Element Element;
			// p = RandPrime(15, 19);

			BlasMatrix<Field> DA(F,n,m);
			// DA = XCALLOC(double, n * m);
			for (size_t i = 0; i < n * m; i++)
				F.init(DA.getWritePointer()[i],A[i]); //!@todo hom
			// DA[i] = (double) ((temp = (A[i] % ((long) p))) >= 0 ? temp : ((long) p) + temp);
			for (size_t i = 0; i < n + 1; i++) {
				P[i] = i;
				rp[i] = 0;
			}
			Element d = F.one;
			RowEchelonTransform<Field> RET(DA);
			RET.reduce(1,1,0,0,P,rp,d);
			// RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
			// XFREE(DA);
			rank = rp[0];
			size_t s = m - rank;
			if ((s == 0) || (rank==0)) {
				mp_N.resize(0,0);
			}
#if 0
			else if (r == 0) {
				BlasMatrixDomain BMD(A.field());
				bool zero = BMD.isZero(A);
				// for (i = 0; i < n * m; i++)
				// if (A[i] != 0)
				// flag = 0;
				if (zero == false) {
					continue;
					//! @todo add a safe guard here, lb_max_repet
				}
				mp_N.resize(m,m);
				// mp_N = XCALLOC(mpz_t, m * m);
				BMD.setIdentity(mp_N);
				// for (i = 0; i < m; i++) {
				// for (j = 0; j < m; j++)
				// mpz_init_set_ui(mp_N[i * m + j], 0);
				// mpz_init_set_ui(mp_N[i * m + i], 1);
				// }
				// *mp_N_pass = mp_N;
				return (ok = true) ;
			}
#endif
			else {		/* r>0 and s>0 */

				// std::vector<size_t> Pt ;
				revseq(Pt,rank,P);
				// Pt = revseq(r, n, P);
				std::vector<size_t> rpt(m,0) ;
				rp.resize(m); //! @bug why ?
				revseq(rpt,rank,rp);
				// rpt = revseq(r, m, rp);

				BlasMatrix<Field> C(A.field(),rank,rank);
				// C = XCALLOC(long, r * r);
				for (size_t i = 0; i < rank; i++)
					for (size_t j = 0; j < rank; j++){
						C.setEntry(i, j, A.getEntry(Pt[i], rpt[j]));
						// C[i * r + j] = A[Pt[i] * m + rpt[j]];
					}

				// mp_B = XCALLOC(mpz_t, r * s);
				//! @todo put in the previous loop ?
				BlasMatrix<Ring> mp_B(R,rank,s);
				for (size_t i = 0; i < rank; i++)
					for (size_t j = 0; j < s; j++) {
						mp_B.setEntry(i,j, A.getEntry(Pt[i], rpt[rank+j]));
						// mpz_init_set_si(mp_B[i * s + j],
						// A[Pt[i] * m + rpt[r + j]]);
					}

				// Integer mp_D;
				// mpz_init(mp_D);
				// BlasMatrix<PID_integer> mp_N(R,s,s);
				// XXX not doing last m-s rows
				mp_N.resize(s,s);
				// mp_N = XCALLOC(mpz_t, m * s);
				// for (i = 0; i < m * s; i++)
				// mpz_init(mp_N[i]);

				nonsingSolve(LinBoxTag::Right, C, mp_B, mp_N, mp_D);

				Integer::negin(mp_D);

			}

		}
	}

	template<class Ring, class Field>
	size_t
	Nullspace<Ring,Field>::
	getNullspaceCompressed( int verfied )
	{
		justTry() ;
		bool ok = verify_compressed(verfied);
		while (!ok) {
			justTry();
			ok = verify_compressed(verfied) ;
		}
		komp = true ;
		ver1 = verfied ;
		return size;
	}


	template<class Ring, class Field>
	bool
	Nullspace<Ring,Field>::
	verify (unsigned verified)
	{
		linbox_check(verified<=2);
		linbox_check(komp= false);

		if (verified == 0)
			return true;

		if (ver1 && verified<=ver1)
			return true; // verified while compressed

		if (verified == 1) {
			std::cout << "not implemented" << std::endl;
			return true;
		}

		if (verified == 2) {
			std::cout << "not implemented" << std::endl;
			return true;
		}

		return false;
	}

	template<class Ring, class Field>
	bool
	Nullspace<Ring,Field>::
	verify_compressed (unsigned verified)
	{
		linbox_check(komp==true);

		ver1 = false ;
		if (verified == 0)
			return true ;
		linbox_check(verified<=2);

		if (verified ==1) {
			std::cout << "not implemented" << std::endl;
			return true;
		}

		if (verified ==2) {
			std::cout << "not implemented" << std::endl;
			return true;
		}
		return false;
	}

	template<class Ring, class Field>
	size_t
	Nullspace<Ring,Field>::
	getNullspace( int verified )
	{
		if (!komp)
			getNullspaceCompressed(verified);
		reconstruct();
		bool ok = verify(verified);
		while (!ok) {
			getNullspaceCompressed(verified);
			reconstruct();
			ok = verify(verified);
		}
		komp = false ;
		ver1=verified;
		return size ;
	}


	template<class Ring, class Field>
	void
	Nullspace<Ring,Field>::
	reconstruct()
	{
		komp = false ;
		mp_N.resize(A.coldim(),size);
		if (size == A.coldim()) {
			BlasMatrixDomain<Ring> BMD(R);
			BMD.setIdentity(mp_N);
			return;
		}
		for (size_t i = 0; i < size; i++)
			mp_N.setEntry((rank + i), i, mp_D);

		//!@todo use blas permutation
		for (size_t i = rank; i >= 1; i--) {
			//!@bug can skip some
			FFLAS::fswap(R,size,mp_N.getWritePointer()+((i-1)*size),1,
				     mp_N.getWritePointer()+((rp[i]-1)*size),1);
			// for (size_t j = 0; j < size; j++)
			// mpz_swap(mp_N[(i - 1) * s + j],mp_N[(rp[i] - 1) * s + j]);
		}
		return;

	}

} // IML
} // LinBox

#endif // __LINBOX_algorithms_iml_nullspace_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
