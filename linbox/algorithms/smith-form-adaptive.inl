/* Copyright (C)  LinBox
 *
 * authors: bds and zw
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#ifndef __LINBOX_smith_form_adaptive_INL
#define __LINBOX_smith_form_adaptive_INL

#include <cmath>
#include <vector>
#include <givaro/modular-integral.h>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/util/debug.h"
#include "linbox/ring/pir-modular-int32.h"
#include "linbox/ring/local2_32.h"
#include "linbox/ring/local-pir-modular.h"
#include "linbox/algorithms/smith-form-iliopoulos.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/algorithms/rational-solver-adaptive.h"
#include "linbox/algorithms/last-invariant-factor.h"
#include "linbox/algorithms/one-invariant-factor.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/blackbox/random-matrix.h"
#include "linbox/blackbox/scompose.h"
#include <fflas-ffpack/ffpack/ffpack.h>
#include "linbox/algorithms/smith-form-binary.h"
#include "linbox/algorithms/smith-form-adaptive.inl"
#include "linbox/solutions/valence.h"

namespace LinBox
{

	/* Compute the local smith form at prime p, when modular (p^e) fits in long
	*/
	template <class Matrix>
	void SmithFormAdaptive::compute_local_long (BlasVector<Givaro::ZRing<Integer> >& s, const Matrix& A, int64_t p, int64_t e)
	{
		//std::ostream& report(std::cout);
		std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);

		int order = (int)(A. rowdim() < A. coldim() ? A. rowdim() : A. coldim());
		linbox_check ((s. size() >= (size_t)order) && (p > 0) && ( e >= 0));
		if (e == 0) return;

		if (p == 2) {
			report << "      Compute local smith at 2^32 using special Local2_32\n";
			Local2_32 R;
			std::list <Local2_32::Element> l;
			SmithFormLocal<Local2_32> SF;
			BlasMatrix <Local2_32> A_local(R, A.rowdim(),A.coldim());
			MatrixHom::map (A_local, A);
			SF (l, A_local, R);
			std::list <Local2_32::Element>::iterator l_p;
			BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
			for (s_p = s. begin(), l_p = l. begin(); s_p != s. begin() +(ptrdiff_t) order; ++ s_p, ++ l_p) {
				*s_p = *l_p;
				/* deb */report << '[' << *s_p << ',' << *l_p << ']';
			}
			report << "     Done\n";
		}
		else if (e == 1) {
			report << "      Compute local smith at prime " << p << ", by rank.\n";
#if 0 /*Meet trouble to call ffpack routine*/
			Givaro::Modular<double> F (p); Givaro::Modular<double>::Element elt;
			int n = A. rowdim(); int m = A. coldim();
			Givaro::Modular<double>::Element* A_local = new Givaro::Modular<double>::Element [n * m];
			typename Matrix::ConstIndexedIterator rawi_p;
			typename Matrix::ConstIterator raw_p;
			Givaro::Modular<double>::Element* A_local_p;
			for (A_local_p = A_local; A_local_p != A_local + (n*m); ++ A_local_p)
				F. assign (*A_local_p, F.zero);
			integer tmp;
			for (rawi_p = A. IndexedBegin(), raw_p = A. Begin(), A_local_p = A_local; rawi_p != A. IndexedEnd(); ++ rawi_p, ++ raw_p, ++ A_local_p) {
				//F. init (*A_local_p, *raw_p);
				A. field(). convert (tmp, *raw_p);
				F. init (elt, tmp);
				F. assign (*(A_local + (int(rawi_p. rowIndex()) * m + int(rawi_p. colIndex()))), elt);
			}

			std::cout << "Initialize matrix done\n";
			for (int i = 0; i < n * m; ++ i)
				std::cout << *(A_local + i) << '\t';
			std::cout << "\nbegin to call ffpack:\n";

			unsigned int rank = FPACK::Rank ((Givaro::Modular<double>::Father_t)(F, n, m, A_local, m);
			std::cout << "Call of ffpack is done\n";
			delete[] A_local;
#endif
			size_t rank; integral_rank(rank, A, Method::Auto());

			BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
			for (s_p = s. begin(); s_p != s. begin() + (long) rank; ++ s_p)
				*s_p = 1;
			for (; s_p != s. begin() +(ptrdiff_t) order; ++ s_p)
				*s_p = 0;
			report << "      Done\n";
		}
		else {
			report << "      Compute local smith at " << p <<'^' << e << " using PIRModular<int32_t>\n";
			int64_t m = 1;
			int i = 0;
			for (i = 0; i < e; ++ i)
				m *= p;
			typedef PIRModular<int32_t> PIR;
			PIR R((int32_t)m);
			BlasMatrix <PIR> A_local(R, A.rowdim(), A.coldim());
			SmithFormLocal <PIR> SF;
			std::list <PIR::Element> l;
			MatrixHom::map (A_local, A);
			SF (l, A_local, R);
			std::list <PIR::Element>::iterator l_p;
			BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
			for (s_p = s. begin(), l_p = l. begin(); s_p != s. begin() +(ptrdiff_t) order; ++ s_p, ++ l_p)
				*s_p = *l_p;
			report <<  "      Done\n";
		}

	}

	/* Compute the local smith form at prime p, when modular (p^e) doesnot fit in long
	*/
	template <class Matrix>
	void SmithFormAdaptive::compute_local_big (BlasVector<Givaro::ZRing<Integer> >& s, const Matrix& A, int64_t p, int64_t e)
	{

		//std::ostream& report(std::cout);
		std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		int order = (int)(A. rowdim() < A. coldim() ? A. rowdim() : A. coldim());
		linbox_check ((s. size() >= (size_t) order) && (p > 0) && ( e >= 0));
		integer T; T = order; T <<= 20; T = pow (T, (int) sqrt((double)order));
		Givaro::Integer m(Integer::one);
        for (int i = 0; i < e; ++ i) m *= p;
        report << "      Compute local Smith at " << p << '^' << e << " over ZZ_p\n";
        typedef LocalPIRModular<Givaro::Integer> GZZ_p;
        GZZ_p R(m);
        BlasMatrix <GZZ_p> A_local(R, A.rowdim(), A.coldim());
        SmithFormLocal <GZZ_p> SF;
        std::list <GZZ_p::Element> l;
        MatrixHom::map (A_local, A);
        SF (l, A_local, R);
        std::list <GZZ_p::Element>::iterator l_p;
        BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
        for (s_p = s. begin(), l_p = l. begin(); s_p != s. begin() +(ptrdiff_t) order; ++ s_p, ++ l_p)
            R. convert(*s_p, *l_p);
        report << "      Done \n";
		return;
	}

	/* Compute the local smith form at prime p
	*/
	template <class Matrix>
	void SmithFormAdaptive::compute_local (BlasVector<Givaro::ZRing<Integer> >& s, const Matrix& A, int64_t p, int64_t e)
	{

		linbox_check ((p > 0) && ( e >= 0));
		integer m = 1; int i = 0; for ( i = 0; i < e; ++ i) m *= p;
		if (((p == 2) && (e <= 32)) || (m <= FieldTraits<PIRModular<int32_t> >::maxModulus()))
			compute_local_long (s, A, p, e);
		else
			compute_local_big (s, A, p, e);

		// normalize the answer
		for (BlasVector<Givaro::ZRing<Integer> >::iterator p_it = s. begin(); p_it != s. end(); ++ p_it)
			*p_it = gcd (*p_it, m);
	}

	/* Compute the k-smooth part of the invariant factor, where k = 100.
	 * @param sev is the exponent part ...
	 * By local smith form and rank computation
	 * r >= 2;
	 */
	template <class Matrix>
	void SmithFormAdaptive::smithFormSmooth (BlasVector<Givaro::ZRing<Integer> >& s, const Matrix& A, long r, const std::vector<int64_t>& sev)
	{
		Givaro::ZRing<Integer> Z;
		//std::ostream& report(std::cout);
		std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation the k-smooth part of the invariant factors starts(via local and rank):" << std::endl;
		int order = (int)(A. rowdim() < A. coldim() ? A. rowdim() : A. coldim());
		linbox_check (s. size() >= (size_t)order);
		std::vector<int64_t>::const_iterator sev_p; const int64_t* prime_p; BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
		BlasVector<Givaro::ZRing<Integer> > local(Z,(size_t)order);
		BlasVector<Givaro::ZRing<Integer> >::iterator local_p;

		for (s_p = s. begin(); s_p != s. begin() +(ptrdiff_t) r; ++ s_p)
			*s_p = 1;
		for (; s_p != s. end(); ++ s_p)
			*s_p = 0;
		if (r == 0) return;

		for (sev_p = sev. begin(), prime_p = prime; sev_p != sev. begin() +(ptrdiff_t) NPrime; ++ sev_p, ++ prime_p) {
			int extra = 1;
			do {

				if ((*prime_p == 2) && (*sev_p < 32))
					extra =  32 -(int) *sev_p;
				// put in a warning if over 2^32
				integer m = 1;
				for (int i = 0; i < *sev_p + extra; ++ i) m *= * prime_p;
				report << "   Compute the local smith form mod " << *prime_p <<"^" << *sev_p + extra << std::endl;
				compute_local (local, A, *prime_p, *sev_p + extra);
				//check
				report << "   Check if it agrees with the rank: ";
				if ((local[(size_t)r-1] % m != 0 ) && ((r == order) ||(local[(size_t)r] % m == 0))) {report << "yes.\n"; break;}
				report << "no. \n";
				extra *= 2;
			} while (true);
			for (s_p = s. begin(), local_p = local. begin(); s_p != s. begin() +(ptrdiff_t) order; ++ s_p, ++ local_p)
				*s_p *= *local_p;
		}
		report << "Computation of the smooth part is done.\n";

	}


	/* Compute the k-rough part of the invariant factor, where k = 100.
	 * By EGV+ algorithm or Iliopoulos' algorithm for Smith form.
	 */
	template <class Matrix>
	void SmithFormAdaptive::smithFormRough  (BlasVector<Givaro::ZRing<Integer> >& s, const Matrix& A, integer m)
	{

		//std::ostream& report(std::cout);
		std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Compuation of the k-rough part f the invariant factors starts(via EGV+ or Iliopolous):\n";
		int order = (int)(A. rowdim() < A. coldim() ? A. rowdim() : A. coldim());
		integer T; T = order; T <<= 20; T = pow (T, (int) sqrt((double)order));
		linbox_check ((s. size() >= (size_t)order) && (m > 0));
		if (m == 1)
			report << "   Not rough part." << std::endl;
		else if ( m <=  FieldTraits< PIRModular<int32_t> >::maxModulus() ) {
			report << "    Elimination starts:\n";
			PIRModular<int32_t> R (m);
			BlasMatrix<PIRModular<int32_t> > A_ilio(R, A.rowdim(), A.coldim());
			MatrixHom::map (A_ilio, A);
			SmithFormIliopoulos::smithFormIn (A_ilio);
			int i; BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
			for (i = 0, s_p = s. begin(); s_p != s. begin() +(ptrdiff_t) order; ++ i, ++ s_p)
				R. convert(*s_p, A_ilio[(size_t)i][(size_t)i]);
			report << "    Elimination ends.\n";
		}
		// else if bisection possible
		else if (m > T)  {
			report << "   Big rough part, bisection starts:\n";
			typedef typename Matrix::Field Ring;
			typedef RationalSolverAdaptive Solver;
			typedef LastInvariantFactor<Ring, Solver> LIF;
			typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
			SmithFormBinary<Ring, OIF> sf;;
			sf. setOIFThreshold (2);
			sf. setLIFThreshold (2);
			std::vector<int64_t> primeL (prime, prime + NPrime);
			std::vector<typename Ring::Element> out ((size_t)order);
			sf. smithForm (out, A, primeL);
			typename std::vector<typename Ring::Element>::iterator out_p;
			BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
			for (s_p = s. begin(), out_p = out. begin(); out_p != out. end(); ++ out_p, ++ s_p)
				A. field(). convert (*s_p, *out_p);
			report << "   Big rough part, bisection ends.\n";
		}
		else {
			report << "    Elimination start:\n";
			typedef LocalPIRModular<Givaro::Integer> GZZ_p;
            GZZ_p R (m);
			BlasMatrix<GZZ_p> A_ilio(R, A.rowdim(), A.coldim());
			MatrixHom::map (A_ilio, A);
			SmithFormIliopoulos::smithFormIn (A_ilio);
			int i; BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
			for (i = 0, s_p = s. begin(); s_p != s. begin() +(ptrdiff_t) order; ++ i, ++ s_p)
				R. convert(*s_p, A_ilio[(size_t)i][(size_t)i]);
			report << "    Elimination ends.\n";
		}
		report << "Compuation of the k-rough part of the invariant factors finishes.\n";
	}

	/* Compute the Smith form via valence algorithms
	 * Compute the local Smtih form at each possible prime
	 * r >= 2;
	 */
	template <class Matrix>
	void SmithFormAdaptive::smithFormVal (BlasVector<Givaro::ZRing<Integer> >&s, const Matrix& A, long r, const std::vector<int64_t>& sev)
	{
		//....
		Givaro::ZRing<Integer> Z;
		//std::ostream& report(std::cout);
		std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation the local smith form at each possible prime:\n";
		int order = (int)(A. rowdim() < A. coldim() ? A. rowdim() : A. coldim());
		linbox_check (s. size() >= (size_t)order);

		std::vector<int64_t>::const_iterator sev_p;
		const int64_t* prime_p;

		BlasVector<Givaro::ZRing<Integer> >::iterator s_p;
		BlasVector<Givaro::ZRing<Integer> > local(Z,(size_t)order);

		BlasVector<Givaro::ZRing<Integer> >::iterator local_p;

		for (s_p = s. begin(); s_p != s. begin() +(ptrdiff_t) r; ++ s_p)
			*s_p = 1;
		for (; s_p != s. end(); ++ s_p)
			*s_p = 0;
		if (r == 0) return;

		for (sev_p = sev. begin(), prime_p = prime; sev_p != sev. begin() +(ptrdiff_t) NPrime; ++ sev_p, ++ prime_p) {
			if (*sev_p <= 0) continue;
			//only compute the local Smith form at each possible prime
			int extra = 2;
			if (*prime_p == 2) extra = 32;
			else {
				// cheating here, try to use the max word size modular
				double log_max_mod = log((double) FieldTraits<PIRModular<int32_t> >:: maxModulus() - 1) ;
				extra = (int)(floor(log_max_mod / log (double(*prime_p))));
			}
			do {
				integer m = 1;
				for (int i = 0; i < extra; ++ i) m *= * prime_p;
				report << "   Compute the local smith form mod " << *prime_p <<"^" << extra << std::endl;
				compute_local (local, A, *prime_p, extra);
				//check
				report << "   Check if it agrees with the rank: ";
				if ((local[(size_t)r-1] % m != 0 ) && ((r == order) ||(local[(size_t)r] % m == 0))) {report << "yes.\n"; break;}
				report << "no. \n";
				extra *= 2;
			} while (true);
			for (s_p = s. begin(), local_p = local. begin(); s_p != s. begin() +(ptrdiff_t) order; ++ s_p, ++ local_p)
				*s_p *= *local_p;
		}
		report << "Computation of the smith form done.\n";

	}

	/* Compute the Smith form of a dense matrix
	 * By adaptive algorithm.
	 * Compute the valence possible, or valence is rough,
	 * Otherwise, compute the largest invariant factor,
	 * then based on that, compute the rough and smooth part, seperately.
	 */
	template <class Matrix>
	void SmithFormAdaptive::smithForm (BlasVector<Givaro::ZRing<Integer> >& s, const Matrix& A)
	{
		//commentator().start ("Smith Form starts", "Smithform");

		//std::ostream& report(std::cout);
		std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation of the invariant factors starts (via an adaptive alg):" << std::endl;

		// compute the rank over a random prime field.
		int order = (A. rowdim() < A. coldim()) ? (int)A. rowdim() : (int)A. coldim();
		report << "Computation of the rank starts:\n";
		typedef typename Matrix::Field Ring;
		size_t r; integral_rank(r, A, Method::Auto());
		report << "   Matrix rank over a random prime field: " << r << '\n';
		report << "Computation of the rank finished.\n";
		const int64_t* prime_p;
		std::vector<int64_t> e(NPrime); std::vector<int64_t>::iterator e_p;

		report <<"   Compute the degree of min poly of AA^T: \n";
		typedef Givaro::Modular<int32_t> Field;
		integer Val; Field::Element v; size_t degree;
        PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<Field>::bestBitSize(A.coldim()));
		Field F (*genprime);
		typename MatrixHomTrait<Matrix, Field>::value_type Ap(F, A.rowdim(), A.coldim());
		MatrixHom::map (Ap, A);
		Valence::one_valence (v, degree, Ap);
		report <<"   Degree of minimal polynomial of AA^T = " << degree << '\n';
		// if degree is small
		if (degree < sqrt(double(order))) {
			report << "   Computation of the valence starts:\n";
			Valence::valence (Val, degree, A);
			report << "      Valence = " << Val << std::endl;
			report << "   Computation of the valence ends.\n";
			Val = abs (Val);
			//Factor the k-smooth part of Val
			for (prime_p = prime, e_p = e. begin(); e_p != e. end(); ++ prime_p, ++ e_p) {
				*e_p = 0;
				while (Val % *prime_p == 0) {
					++ *e_p;
					Val = Val / *prime_p;
				}
			}
			if (Val == 1) {
				smithFormVal (s, A, r, e);
				report << "Computation of the invariant factors ends." << std::endl;
				return;
			}
			else
				report << "   Valence is rough.\n";
		}

		report << "Computation of the largest invariant factor with bonus starts:\n";
		typedef RationalSolverAdaptive Solver;
		typedef LastInvariantFactor<Ring, Solver> LIF;
		typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
		OIF oif; oif. setThreshold  (4); oif.getLastInvariantFactor().setThreshold (6);
		typename Ring::Element _lif, _bonus; integer lif, bonus;
		//Chnage A to BlasMatrix
		Ring R(A. field());
		BlasMatrix<Ring> DA(R,A.rowdim(),A.coldim());
		MatrixHom::map (DA, A);
		do {
			oif. oneInvariantFactor_Bonus (_lif, _bonus, DA, (int)r);

			A. field(). convert (lif, _lif); A. field(). convert (bonus, _bonus);
			//oif. oneInvariantFactor (bonus, A, (int)r);
			report << "   The largest invariant factor: " << lif << std::endl;
			report << "   Bonus (previous one): " << bonus << std::endl;
		} while (lif == 0);
		report << "Computation of the largest invariant factor with bonus finished.\n";
		// bonus = smooth * rough;
		integer r_mod; r_mod = lif;
		for (prime_p = prime, e_p = e. begin(); e_p != e. end(); ++ prime_p, ++ e_p) {
			*e_p = 0;
			while (r_mod % *prime_p == 0) {
				++ *e_p;
				r_mod = r_mod / *prime_p;
			}
		}
		// bonus assigns to its rough part
		bonus = gcd (bonus, r_mod);
		Givaro::ZRing<Integer> Z;
		BlasVector<Givaro::ZRing<Integer> > smooth (Z,(size_t)order), rough (Z,(size_t)order);
		smithFormRough (rough, DA, bonus);
		smithFormSmooth (smooth, A, r, e);
		//fixed the rough largest invariant factor
		if (r > 0) rough[r-1] = r_mod;

		BlasVector<Givaro::ZRing<Integer> >::iterator s_p, rough_p, smooth_p;

		/*
		   report << "Smooth part\n";
		   for (smooth_p = smooth. begin(); smooth_p != smooth. end(); ++ smooth_p)
		   report<< *smooth_p << ' ';
		   report<< '\n';
		   report<<"Rough part\n";
		   for (rough_p = rough. begin(); rough_p != rough. begin() +(ptrdiff_t) r; ++ rough_p)
		   report<< *rough_p << ' ';
		   report<< '\n';
		   */

		for (rough_p = rough. begin(); rough_p != rough. begin() +(ptrdiff_t) r; ++ rough_p)
			if (* rough_p == 0) *rough_p = bonus;

		for (s_p = s. begin(), smooth_p = smooth. begin (), rough_p = rough. begin(); s_p != s. begin() +(ptrdiff_t) order; ++s_p, ++ smooth_p, ++ rough_p)
			*s_p = *smooth_p * *rough_p;

		report << "Computation of the invariant factors ends." << std::endl;
		//commentator().stop ("done", NULL, "Smithform");
	}

	/* Compute the Smith form of a dense matrix
	 * By adaptive algorithm.
	 * Compute the valence possible, or valence is rough,
	 * Otherwise, compute the largest invariant factor,
	 * then based on that, compute the rough and smooth part, seperately.
	 */
	template <class IRing, class _Rep>
	void SmithFormAdaptive::smithForm (BlasVector<Givaro::ZRing<Integer> >& s, const BlasMatrix<IRing, _Rep>& A)
	{
		//commentator().start ("Smith Form starts", "Smithform");
		Givaro::ZRing<Integer> Z;

		//std::ostream& report(std::cout);
		std::ostream& report = commentator().report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation of the invariant factors starts (via an adaptive alg):" << std::endl;

		// compute the rank over a random prime field.
		const size_t order = (A. rowdim() < A. coldim() ? A. rowdim() : A. coldim());

		report << "Computation of the rank starts:" << std::endl;
		typedef typename BlasMatrix<IRing,_Rep>::Field Ring;
		size_t r; integral_rank(r, A, Method::Auto());
		report << "   Matrix rank over a random prime field: " << r << std::endl;
		report << "Computation of the rank finished.\n";
		// a hack
		if (r == 0) { for (size_t i = 0; i <  order; ++i) s[i]=0; return; }
		const int64_t* prime_p;
		std::vector<int64_t> e(NPrime); std::vector<int64_t>::iterator e_p;

		report <<"   Compute the degree of min poly of AA^T: \n";
		typedef Givaro::Modular<int32_t> Field;
		integer Val; Field::Element v; size_t degree;
		//integer Val; Field::Element v; size_t degree;
        PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<Field>::bestBitSize(A.coldim()));
		Field F (*genprime);
		typename MatrixHomTrait<BlasMatrix <IRing, _Rep>, Field>::value_type Ap(F,A.rowdim(),A.coldim());
		MatrixHom::map (Ap, A);
		Valence::one_valence (v, degree, Ap);
		//valence (v, Ap);
		report <<"   Degree of one_valence minpol of AA^T = " << degree << '\n';
		report <<"   value of one_valence minpol of AA^T = " << v << '\n';
		// if degree is small
		if (degree < sqrt(double(order))) {
			report << "   Computation of the valence starts:\n";
			//valence (Val, A);
			Valence::valence (Val, degree, A);
			report << "      Valence = " << Val << std::endl;
		report <<"   degree of minpol of AA^T = " << degree << '\n';
			if (Val == 0) valence (Val, A);
			report << " Corrected valence = " << Val << std::endl;
			report << "   Computation of the valence ends.\n";
			Val = abs (Val);
			//Factor the k-smooth part of Val
			for (prime_p = prime, e_p = e. begin(); e_p != e. end(); ++ prime_p, ++ e_p) {
				*e_p = 0;
				while (Val % *prime_p == 0) {
					++ *e_p;
					Val = Val / *prime_p;
				}
			}
			if (Val == 1) {
				smithFormVal (s, A, (long)r, e);
				report << "Computation of the invariant factors ends." << std::endl;
				return;
			}
			else
				report << "   Valence is rough.\n";
		}

		report << "Computation of the largest invariant factor with bonus starts:\n";
		typedef RationalSolverAdaptive Solver;
		typedef LastInvariantFactor<Ring, Solver> LIF;
		typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
		OIF oif; oif. setThreshold  (10); oif.getLastInvariantFactor().setThreshold (6);
		typename Ring::Element _lif, _bonus; integer lif, bonus;
		//Chnage A to BlasMatrix
		do {
			oif. oneInvariantFactor_Bonus (_lif, _bonus, A, (int)r);
			A. field(). convert (lif, _lif); A. field(). convert (bonus, _bonus);
			//oif. oneInvariantFactor (bonus, A, (int)r);
			report << "   The largest invariant factor: " << lif << std::endl;
			report << "   Bonus (previous one): " << bonus << std::endl;
		} while (lif == 0);
		report << "Computation of the largest invariant factor with bonus finished.\n";
		// bonus = smooth * rough;
		integer r_mod; r_mod = lif;
		for (prime_p = prime, e_p = e. begin(); e_p != e. end(); ++ prime_p, ++ e_p) {
			*e_p = 0;
			while (r_mod % *prime_p == 0) {
				++ *e_p;
				r_mod = r_mod / *prime_p;
			}
		}
		// bonus assigns to its rough part
		bonus = gcd (bonus, r_mod);
		BlasVector<Givaro::ZRing<Integer> > smooth (Z,order), rough (Z,order);
		report << "Computation of smooth part begins.\n";
		smithFormSmooth (smooth, A, (long)r, e);
		report << "Computation of rough part begins.\n";
		smithFormRough (rough, A, bonus);
		report << "Computation of rough/smooth parts finished.\n";
		// fixed the rough largest invariant factor
		if (r > 0) rough[r-1] = r_mod;

		BlasVector<Givaro::ZRing<Integer> >::iterator s_p, rough_p, smooth_p;

		/*
		   report << "Smooth part\n";
		   for (smooth_p = smooth. begin(); smooth_p != smooth. end(); ++ smooth_p)
		   report<< *smooth_p << ' ';
		   report<< '\n';
		   report<<"Rough part\n";
		   for (rough_p = rough. begin(); rough_p != rough. begin() +(ptrdiff_t) r; ++ rough_p)
		   report<< *rough_p << ' ';
		   report<< '\n';
		   */

		for (rough_p = rough. begin(); rough_p != rough. begin() + (ptrdiff_t)r; ++ rough_p)
			if (* rough_p == 0) *rough_p = bonus;

		for (s_p = s. begin(), smooth_p = smooth. begin (), rough_p = rough. begin(); s_p != s. begin() + (ptrdiff_t)order; ++s_p, ++ smooth_p, ++ rough_p)
			*s_p = *smooth_p * *rough_p;

		report << "Computation of the invariant factors ends." << std::endl;
		//commentator().stop ("done", NULL, "Smithform");
	}

}

#endif //__LINBOX_smith_form_adaptive_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
