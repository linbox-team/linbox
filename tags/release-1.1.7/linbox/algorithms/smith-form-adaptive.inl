/* Copyright (C)  LinBox
 *
 * authors: bds and zw
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __LINBOX_smith_form_adaptive_INL
#define __LINBOX_smith_form_adaptive_INL

#include <math.h>
#include <vector>
#include <linbox/linbox-config.h>
#include <linbox/integer.h>
#include <linbox/util/debug.h>
#include <linbox/field/PIR-modular-int32.h>
#include <linbox/field/local2_32.h>
#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/smith-form-iliopoulos.h>
#include <linbox/algorithms/smith-form-local.h>
#include <linbox/algorithms/rational-solver-adaptive.h>
#include <linbox/algorithms/last-invariant-factor.h>
#include <linbox/algorithms/one-invariant-factor.h>
#include <linbox/algorithms/matrix-rank.h>
#include <linbox/algorithms/matrix-hom.h>
#include <linbox/blackbox/random-matrix.h>
#include <linbox/blackbox/scompose.h>
#include <linbox/ffpack/ffpack.h>
#include <linbox/algorithms/smith-form-binary.h>
#include <linbox/algorithms/smith-form-adaptive.inl>
#include <linbox/solutions/valence.h>



#ifdef __LINBOX_HAVE_NTL
#include <linbox/field/PIR-ntl-ZZ_p.h>
#endif

namespace LinBox 
{ 

	/* Compute the local smith form at prime p, when modular (p^e) fits in long
	*/
	template <class Matrix>
	void SmithFormAdaptive::compute_local_long (std::vector <integer>& s, const Matrix& A, long p, long e) {
		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		linbox_check ((s. size() >= (unsigned long)order) && (p > 0) && ( e >= 0));
		if (e == 0) return;

		if (p == 2) {
			 report << "      Compute local smith at 2^32 using special Local2_32\n";
			 Local2_32 R;
			 std::list <Local2_32::Element> l;
			 SmithFormLocal<Local2_32> SF;
			 DenseMatrix <Local2_32> A_local(R, A.rowdim(),A.coldim()); 
			 MatrixHom::map (A_local, A, R);
			 SF (l, A_local, R);
			 std::list <Local2_32::Element>::iterator l_p;
			 std::vector <integer>::iterator s_p;
			 for (s_p = s. begin(), l_p = l. begin(); s_p != s. begin() + order; ++ s_p, ++ l_p) 
			 	*s_p = *l_p;
			report << "     Done\n";
		}
		else if (e == 1) {
			report << "      Compute local smith at prime " << p << ", by rank.\n";
			/** Meet trouble to call ffpack routine
			Modular<double> F (p); Modular<double>::Element elt;
			int n = A. rowdim(); int m = A. coldim();
			Modular<double>::Element* A_local = new Modular<double>::Element [n * m];
			typename Matrix::ConstRawIndexedIterator rawi_p;
			typename Matrix::ConstRawIterator raw_p;
			Modular<double>::Element* A_local_p;
			for (A_local_p = A_local; A_local_p != A_local + (n*m); ++ A_local_p)
				F. init (*A_local_p, 0);
			integer tmp;
			for (rawi_p = A. rawIndexedBegin(), raw_p = A. rawBegin(), A_local_p = A_local; rawi_p != A. rawIndexedEnd(); ++ rawi_p, ++ raw_p, ++ A_local_p) {
				 //F. init (*A_local_p, *raw_p);
				 A. field(). convert (tmp, *raw_p); 
				 F. init (elt, tmp); 
				 F. assign (*(A_local + (int(rawi_p. rowIndex()) * m + int(rawi_p. colIndex()))), elt);
			}

			std::cout << "Initialize matrix done\n";
			for (int i = 0; i < n * m; ++ i)
				std::cout << *(A_local + i) << '\t';
			std::cout << "\nbegin to call ffpack:\n";

			unsigned int rank = FFPACK::Rank (F, n, m, A_local, m);
			std::cout << "Call of ffpack is done\n";
			delete[] A_local;
			*/
			typedef Modular<int32> Field;
			typedef DenseMatrix<Field> FMatrix;
			MatrixRank<typename Matrix::Field, Field> MR;
			Field F(p); FMatrix A_local(A, F);
			long rank = MR. rankIn (A_local);

			std::vector <integer>::iterator s_p;
			for (s_p = s. begin(); s_p != s. begin() + (long) rank; ++ s_p)
				*s_p = 1;
			for (; s_p != s. begin() + order; ++ s_p)
				*s_p = 0;
			report << "      Done\n";
		}
		else {
			report << "      Compute local smith at " << p <<'^' << e << " using PIRModular<int32>\n";
			long m = 1; int i = 0; for (i = 0; i < e; ++ i) m *= p;
			typedef PIRModular<int32> PIR;
			PIR R(m);
			DenseMatrix <PIR> A_local(R, A.rowdim(), A.coldim()); 
			SmithFormLocal <PIR> SF;
			std::list <PIR::Element> l; 
			MatrixHom::map (A_local, A, R);
			SF (l, A_local, R);
			std::list <PIR::Element>::iterator l_p;
			std::vector <integer>::iterator s_p;
			for (s_p = s. begin(), l_p = l. begin(); s_p != s. begin() + order; ++ s_p, ++ l_p)
				*s_p = *l_p;
			report <<  "      Done\n";
		}
	
	}
	
#ifdef __LINBOX_HAVE_NTL
	/* Compute the local smith form at prime p, when modular (p^e) doesnot fit in long
	*/
	template <class Matrix>
	void SmithFormAdaptive::compute_local_big (std::vector<integer>& s, const Matrix& A, long p, long e) {
		
		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		linbox_check ((s. size() >= (unsigned long) order) && (p > 0) && ( e >= 0));
		integer T; T = order; T <<= 20; T = pow (T, (int) sqrt((double)order));
		NTL::ZZ m;  NTL::conv(m, 1); int i = 0; for (i = 0; i < e; ++ i) m *= p;
		//if (m < T) {
		if (1) {
			report << "      Compute local Smith at " << p << '^' << e << " over PIR-ntl-ZZ_p\n";
			PIR_ntl_ZZ_p R(m);
			DenseMatrix <PIR_ntl_ZZ_p> A_local(R, A.rowdim(), A.coldim()); 
			SmithFormLocal <PIR_ntl_ZZ_p> SF;
			std::list <PIR_ntl_ZZ_p::Element> l; 
			MatrixHom::map (A_local, A, R);
			SF (l, A_local, R);
			std::list <PIR_ntl_ZZ_p::Element>::iterator l_p;
			std::vector <integer>::iterator s_p;
			for (s_p = s. begin(), l_p = l. begin(); s_p != s. begin() + order; ++ s_p, ++ l_p)
				R. convert(*s_p, *l_p);
			report << "      Done \n";
		}
		else {
			std::cout << "Compute local smith at " << p << '^' << e << std::endl;
			std::cerr << "Not implemented yet.\n";
		}
		return;
	}
#endif

	/* Compute the local smith form at prime p
	*/
	template <class Matrix>
	void SmithFormAdaptive::compute_local (std::vector<integer>& s, const Matrix& A, long p, long e) {

		linbox_check ((p > 0) && ( e >= 0));
		integer m = 1; int i = 0; for ( i = 0; i < e; ++ i) m *= p;
		if (((p == 2) && (e <= 32)) || (m <= PIRModular<int32>::getMaxModulus()))
			compute_local_long (s, A, p, e);
		else
			compute_local_big (s, A, p, e);

		// normalize the answer
		for (std::vector<integer>::iterator p = s. begin(); p != s. end(); ++ p)
			*p = gcd (*p, m);
	}

	/* Compute the k-smooth part of the invariant factor, where k = 100.
	 * @param sev is the exponent part ...
	 * By local smith form and rank computation
	 * r >= 2;
	 */
	template <class Matrix>
	void SmithFormAdaptive::smithFormSmooth (std::vector<integer>& s, const Matrix& A, long r, const std::vector<long>& sev) {
		//....
		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation the k-smooth part of the invariant factors starts(via local and rank):" << std::endl;
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		linbox_check (s. size() >= (unsigned long)order);
		std::vector<long>::const_iterator sev_p; const long* prime_p; std::vector<integer>::iterator s_p;
		std::vector<integer> local(order); std::vector<integer>::iterator local_p;

		for (s_p = s. begin(); s_p != s. begin() + r; ++ s_p)
			*s_p = 1;
		for (; s_p != s. end(); ++ s_p)
			*s_p = 0;
		if (r == 0) return;

		for (sev_p = sev. begin(), prime_p = prime; sev_p != sev. begin() + NPrime; ++ sev_p, ++ prime_p) {
			int extra = 1;
			do {

				if ((*prime_p == 2) && (*sev_p < 32)) extra = 32 - *sev_p;
				integer m = 1;
				for (int i = 0; i < *sev_p + extra; ++ i) m *= * prime_p;
				report << "   Compute the local smith form mod " << *prime_p <<"^" << *sev_p + extra << std::endl;
				compute_local (local, A, *prime_p, *sev_p + extra);
				//check
				report << "   Check if it agrees with the rank: ";
				if ((local[r-1] % m != 0 ) && ((r == order) ||(local[r] % m == 0))) {report << "yes.\n"; break;}
				report << "no. \n";
				extra *= 2;
			} while (true);
			for (s_p = s. begin(), local_p = local. begin(); s_p != s. begin() + order; ++ s_p, ++ local_p) 
				*s_p *= *local_p;
		}
		report << "Computation of the smooth part is done.\n";
		
	}


#ifdef __LINBOX_HAVE_NTL			
	/* Compute the k-rough part of the invariant factor, where k = 100.
	 * By EGV+ algorithm or Iliopoulos' algorithm for Smith form.
	*/
	template <class Matrix>
	void SmithFormAdaptive::smithFormRough  (std::vector<integer>& s, const Matrix& A, integer m) {

		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Compuation of the k-rough part f the invariant factors starts(via EGV+ or Iliopolous):\n";
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		integer T; T = order; T <<= 20; T = pow (T, (int) sqrt((double)order));
		linbox_check ((s. size() >= (unsigned long)order) && (m > 0));
		if (m == 1) 
			report << "   Not rough part." << std::endl;
		else if ( m <=  PIRModular<int32>::getMaxModulus() ) {
			report << "    Elimination starts:\n";
			PIRModular<int32> R (m);
			DenseMatrix<PIRModular<int32> > A_ilio(R, A.rowdim(), A.coldim());
			MatrixHom::map (A_ilio, A, R);
			SmithFormIliopoulos::smithFormIn (A_ilio);
			int i; std::vector<integer>::iterator s_p;
			for (i = 0, s_p = s. begin(); s_p != s. begin() + order; ++ i, ++ s_p)
				R. convert(*s_p, A_ilio[i][i]);
			report << "    Elimination ends.\n";
		}
		// else if bisection possible
		else if (m > T)  {
			report << "   Big rough part, bisection starts:\n";
			typedef Modular<int> Field;
			typedef typename Matrix::Field Ring;
			typedef RationalSolverAdaptive Solver;
		   	typedef LastInvariantFactor<Ring, Solver> LIF;
			typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
			SmithFormBinary<Ring, OIF, MatrixRank<Ring, Field > > sf;;
			sf. setOIFThreshold (2);
			sf. setLIFThreshold (2);
			std::vector<long> primeL (prime, prime + NPrime);
			std::vector<typename Ring::Element> out (order);
			sf. smithForm (out, A, primeL);
			typename std::vector<typename Ring::Element>::iterator out_p;
			std::vector<integer>::iterator s_p;
			for (s_p = s. begin(), out_p = out. begin(); out_p != out. end(); ++ out_p, ++ s_p)
				A. field(). convert (*s_p, *out_p);
			report << "   Big rough part, bisection ends.\n";
		}
		else {
			report << "    Elimination start:\n";
			PIR_ntl_ZZ_p R (m);
			DenseMatrix<PIR_ntl_ZZ_p> A_ilio(R, A.rowdim(), A.coldim());
			MatrixHom::map (A_ilio, A, R);
			SmithFormIliopoulos::smithFormIn (A_ilio);
			int i; std::vector<integer>::iterator s_p;
			for (i = 0, s_p = s. begin(); s_p != s. begin() + order; ++ i, ++ s_p)
                            R. convert(*s_p, A_ilio[i][i]);
			report << "    Elimination ends.\n";
		}
		report << "Compuation of the k-rough part of the invariant factors finishes.\n";
	}
#endif

	/* Compute the Smith form via valence algorithms
	 * Compute the local Smtih form at each possible prime
	 * r >= 2;
	 */
	template <class Matrix>
	void SmithFormAdaptive::smithFormVal (std::vector<integer>&s, const Matrix& A, long r, const std::vector<long>& sev){
		//....
		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation the local smith form at each possible prime:\n";
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		linbox_check (s. size() >= (unsigned long)order);
		std::vector<long>::const_iterator sev_p; const long* prime_p; std::vector<integer>::iterator s_p;
		std::vector<integer> local(order); std::vector<integer>::iterator local_p;

		for (s_p = s. begin(); s_p != s. begin() + r; ++ s_p)
			*s_p = 1;
		for (; s_p != s. end(); ++ s_p)
			*s_p = 0;
		if (r == 0) return;

		for (sev_p = sev. begin(), prime_p = prime; sev_p != sev. begin() + NPrime; ++ sev_p, ++ prime_p) {
			if (*sev_p <= 0) continue;
			//only compute the local Smith form at each possible prime
			int extra = 2;
			if (*prime_p == 2) extra = 32;
			else {
				// cheating here, try to use the max word size modular
				extra = (int)(floor(log((double)PIRModular<int32>::getMaxModulus() - 1) / log (double(*prime_p))));
			}
			do {
				integer m = 1;
				for (int i = 0; i < extra; ++ i) m *= * prime_p;
				report << "   Compute the local smith form mod " << *prime_p <<"^" << extra << std::endl;
				compute_local (local, A, *prime_p, extra);
				//check
				report << "   Check if it agrees with the rank: ";
				if ((local[r-1] % m != 0 ) && ((r == order) ||(local[r] % m == 0))) {report << "yes.\n"; break;}
				report << "no. \n";
				extra *= 2;
			} while (true);
			for (s_p = s. begin(), local_p = local. begin(); s_p != s. begin() + order; ++ s_p, ++ local_p) 
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
	void SmithFormAdaptive::smithForm (std::vector<integer>& s, const Matrix& A) {
	 	//commentator.start ("Smith Form starts", "Smihtform");

		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation of the invariant factors starts (via an adaptive alg):" << std::endl;

		// compute the rank over a random prime field.
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		report << "Computation of the rank starts:\n";
		typedef typename Matrix::Field Ring;
		unsigned long r; 
		MatrixRank<Ring, Modular<int> > MR;
		r = MR. rank (A);
		report << "   Matrix rank over a random prime field: " << r << '\n';
		report << "Computation of the rank finished.\n";
		const long* prime_p;
		std::vector<long> e(NPrime); std::vector<long>::iterator e_p;

		report <<"   Compute the degree of min poly of AA^T: \n";
		typedef Modular<int> Field;
		integer Val; Field::Element v; unsigned long degree;
		RandomPrimeIterator rg ((int)(log( (double)(Field::getMaxModulus()) ) / M_LN2 - 2));
		Field F (*rg); 
		typename MatrixHomTrait<Matrix, Field>::value_type Ap(F, A.rowdim(), A.coldim());
		MatrixHom::map (Ap, A, F);
		Valence::one_valence (v, degree, Ap);
		report <<"   Degree of minimal polynomial of AA^T = " << degree << '\n';
		// if degree is small
		if (degree < sqrt(double(order))) {
			report << "   Computation of the valence starts:\n";
			Valence::valence (Val, degree, A);
			report << "      Valence = " << Val << std::endl;
			report << "   Computation of the valence of ends.\n";
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
		//Chnage A to DenseMatrix
		Ring R(A. field());
		DenseMatrix<Ring> DA(R,A.rowdim(),A.coldim()); 
                MatrixHom::map (DA, A, R);
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
		std::vector<integer> smooth (order), rough (order);
		smithFormRough (rough, DA, bonus);
		smithFormSmooth (smooth, A, r, e); 
		//fixed the rough largest invariant factor
		if (r > 0) rough[r-1] = r_mod;

		std::vector<integer>::iterator s_p, rough_p, smooth_p;

		/*
		report << "Smooth part\n";
		for (smooth_p = smooth. begin(); smooth_p != smooth. end(); ++ smooth_p)
			report<< *smooth_p << ' ';
		report<< '\n';
		report<<"Rough part\n";
		for (rough_p = rough. begin(); rough_p != rough. begin() + r; ++ rough_p)
			report<< *rough_p << ' ';
		report<< '\n';
		*/

		for (rough_p = rough. begin(); rough_p != rough. begin() + r; ++ rough_p) 
			if (* rough_p == 0) *rough_p = bonus;

		for (s_p = s. begin(), smooth_p = smooth. begin (), rough_p = rough. begin(); s_p != s. begin() + order; ++s_p, ++ smooth_p, ++ rough_p) 
			*s_p = *smooth_p * *rough_p;

		report << "Computation of the invariant factors ends." << std::endl;
		//commentator. stop ("done", NULL, "Smithform");
	}
	/* Compute the Smith form of a dense matrix
	 * By adaptive algorithm. 
	 * Compute the valence possible, or valence is rough, 
	 * Otherwise, compute the largest invariant factor,
	 * then based on that, compute the rough and smooth part, seperately.
	 */
	template <class IRing>
	void SmithFormAdaptive::smithForm (std::vector<integer>& s, const DenseMatrix<IRing>& A) {
	 	//commentator.start ("Smith Form starts", "Smithform");

		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation of the invariant factors starts (via an adaptive alg):" << std::endl;

		// compute the rank over a random prime field.
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		report << "Computation of the rank starts:\n";
		typedef typename DenseMatrix<IRing>::Field Ring;
		unsigned long r; 
		MatrixRank<Ring, Modular<int> > MR;
		r = MR. rank (A);
		report << "   Matrix rank over a random prime field: " << r << '\n';
		report << "Computation of the rank finished.\n";
		const long* prime_p;
		std::vector<long> e(NPrime); std::vector<long>::iterator e_p;

		report <<"   Compute the degree of min poly of AA^T: \n";
		typedef Modular<int> Field;
		integer Val; Field::Element v; unsigned long degree;
		RandomPrimeIterator rg ((int)(log( (double)(Field::getMaxModulus()) ) / M_LN2 - 2));
		Field F (*rg); 
		typename MatrixHomTrait<DenseMatrix<IRing>, Field>::value_type Ap(F,A.rowdim(),A.coldim());
                MatrixHom::map (Ap, A, F);
		Valence::one_valence (v, degree, Ap);
		report <<"   Degree of minial polynomial of AA^T = " << degree << '\n';
		// if degree is small
		if (degree < sqrt(double(order))) {
			report << "   Computation of the valence starts:\n";
			Valence::valence (Val, degree, A);
			report << "      Valence = " << Val << std::endl;
			report << "   Computation of the valence of ends.\n";
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
		OIF oif; oif. setThreshold  (10); oif.getLastInvariantFactor().setThreshold (6);
		typename Ring::Element _lif, _bonus; integer lif, bonus;
		//Chnage A to DenseMatrix
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
		std::vector<integer> smooth (order), rough (order);
		smithFormSmooth (smooth, A, r, e);
		smithFormRough (rough, A, bonus);
		//fixed the rough largest invariant factor
		if (r > 0) rough[r-1] = r_mod;

		std::vector<integer>::iterator s_p, rough_p, smooth_p;

		/*
		report << "Smooth part\n";
		for (smooth_p = smooth. begin(); smooth_p != smooth. end(); ++ smooth_p)
			report<< *smooth_p << ' ';
		report<< '\n';
		report<<"Rough part\n";
		for (rough_p = rough. begin(); rough_p != rough. begin() + r; ++ rough_p)
			report<< *rough_p << ' ';
		report<< '\n';
		*/

		for (rough_p = rough. begin(); rough_p != rough. begin() + r; ++ rough_p) 
			if (* rough_p == 0) *rough_p = bonus;

		for (s_p = s. begin(), smooth_p = smooth. begin (), rough_p = rough. begin(); s_p != s. begin() + order; ++s_p, ++ smooth_p, ++ rough_p) 
			*s_p = *smooth_p * *rough_p;

		report << "Computation of the invariant factors ends." << std::endl;
		//commentator. stop ("done", NULL, "Smithform");
	}
} 

#endif //__LINBOX_smith_form_adaptive_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
