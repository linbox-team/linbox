/* not sure filename
 * authors: bds and zw
 */
#include <vector>
#include <linbox/integer.h>
#include <linbox/util/debug.h>
#include <linbox/field/PIR-ntl-ZZ_p.h>
#include <linbox/field/PIR-modular-int32.h>
#include <linbox/field/local2_32.h>
#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/iliopoulos-elimination.h>
#include <linbox/algorithms/local-smith.h>
#include <linbox/algorithms/2local-smith.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/algorithms/last-invariant-factor.h>
#include <linbox/algorithms/one-invariant-factor.h>
#include <linbox/algorithms/matrix-rank.h>
#include <linbox/algorithms/matrix-mod.h>
#include <linbox/blackbox/random-matrix.h>
#include <linbox/blackbox/scompose.h>
#include <linbox/fflapack/fflapack.h>

namespace LinBox {

	/* Compute the local smith form at prime p, when modular (p^e) fits in long
	*/
	template <class Ring>
	void compute_local_long (std::vector <integer>& s, const DenseMatrix <Ring>& A, long p, long e) {
		
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		linbox_check ((s. size() >= (unsigned long)order) && (p > 0) && ( e >= 0));
		if (e == 0) return;

		if (p == 2) {
			 Local2_32 R;
			 DenseMatrix <Local2_32>* A_local; MatrixMod::mod (A_local, A, R);
			 std::list <Local2_32::Element> l;
			 LocalSmith<Local2_32> SF;
			 SF (l, *A_local, R);
			 std::list <Local2_32::Element>::iterator l_p;
			 std::vector <integer>::iterator s_p;
			 for (s_p = s. begin(), l_p = l. begin(); s_p != s. begin() + order; ++ s_p, ++ l_p) 
			 	*s_p = *l_p;
			delete A_local;
		}
		else if (e == 1) {
			Modular<double> F (p);
			int n = A. rowdim(); int m = A. coldim();
			Modular<double>::Element* A_local = new Modular<double>::Element [n * m];
			typename DenseMatrix <Ring>::ConstRawIterator raw_p;
			Modular<double>::Element* A_local_p;
			integer tmp;
			for (raw_p = A. rawBegin(), A_local_p = A_local; raw_p != A. rawEnd(); ++ raw_p, ++ A_local_p)
				//F. init (*A_local_p, *raw_p);
				{ A. field(). convert (tmp, *raw_p); F. init (*A_local_p, tmp);}
			unsigned int rank = FFLAPACK::Rank (F, n, m, A_local, m);
			delete[] A_local;
			std::vector <integer>::iterator s_p;
			for (s_p = s. begin(); s_p != s. begin() + (long) rank; ++ s_p)
				*s_p = 1;
			for (; s_p != s. begin() + order; ++ s_p)
				*s_p = 0;
		}
		else {
			long m = 1; int i = 0; for (i = 0; i < e; ++ i) m *= p;
			PIRModular <int32> R(m);
			DenseMatrix <PIRModular<int32> >* A_local; MatrixMod::mod (A_local, A, R);
			LocalSmith <PIRModular<int32> > SF;
			std::list <PIRModular<int32>::Element> l; SF (l, *A_local, R);
			std::list <PIRModular<int32>::Element>::iterator l_p;
			std::vector <integer>::iterator s_p;
			for (s_p = s. begin(), l_p = l. begin(); s_p != s. begin() + order; ++ s_p, ++ l_p)
				*s_p = *l_p;
			delete A_local;
		}
	
	}
	
	/* Compute the local smith form at prime p, when modular (p^e) doesnot fit in long
	*/
	template <class Ring>
	void compute_local_big (std::vector<integer>& s, const DenseMatrix<Ring>& A, long p, long e) {
		
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		linbox_check ((s. size() >= (unsigned long) order) && (p > 0) && ( e >= 0));

		std::cout << p << ' ' << e << std::endl;
		std::cerr << "Not implemented yet.\n";
		return;
	}


	/* Compute the local smith form at prime p
	*/
	template <class Ring>
	void compute_local (std::vector<integer>& s, const DenseMatrix<Ring>& A, long p, long e) {

		linbox_check ((p > 0) && ( e >= 0));
		integer m = 1; int i = 0; for ( i = 0; i < e; ++ i) m *= p;
		if (((p == 2) && (e <= 32)) || (m < PIRModular<int32>::getMaxModulus()))
			compute_local_long (s, A, p, e);
		else
			compute_local_big (s, A, p, e);
	}

	/* Compute the k-smooth part of the invariant factor, where k = 100.
	 * @param sev is the exponent part ...
	 * By local smith form and rank computation
	 * r >= 2;
	 */
	template <class Ring>
	void SmithFormSmooth (std::vector<integer>& s, const DenseMatrix<Ring>& A, long r, const std::vector<long>& sev) {
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
				if ((*prime_p == 2) && (*sev_p <= 32)) extra = 32 - *sev_p;
				integer m = 1;
				for (int i = 0; i < *sev_p + extra; ++ i) m *= * prime_p;
				report << "   Compute the local smith form mod " << *prime_p <<"^" << *sev_p + extra << std::endl;
				compute_local (local, A, *prime_p, *sev_p + extra);
				//check
				report << "   Check if it agrees with the rank: ";
				if ((local[r-1] % m != 0 ) && ((r == order) ||(local[r] % m == 0))) {report << "yes.\n"; break;}
				report << "no. \n";
				extra += 2;
			} while (true);
			for (s_p = s. begin(), local_p = local. begin(); s_p != s. begin() + order; ++ s_p, ++ local_p) 
				*s_p *= *local_p;
		}
		report << "Computation of the smooth part is done.\n";
		
	}
			
	/* Compute the k-rough part of the invariant factor, where k = 100.
	 * By EGV+ algorithm or Iliopoulos' algorithm for Smith form.
	*/
	template <class Ring>
	void SmithFormRough  (std::vector<integer>& s, const DenseMatrix<Ring>& A, integer m) {

		if (m == 1)  return;

		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Compuation of the k-rough part f the invariant factors starts(via EGV+ or Iliopolous):\n";
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		linbox_check ((s. size() >= (unsigned long)order) && (m > 0));
		if (m == 1) 
			report << "   Not rough part." << std::endl;
		else if ( m <  PIRModular<int32>::getMaxModulus() ) {
			PIRModular<int32> R (m);
			DenseMatrix<PIRModular<int32> >* A_ilio;
			MatrixMod::mod (A_ilio, A, R);
			IliopoulosElimination::smithIn (*A_ilio);
			int i; std::vector<integer>::iterator s_p;
			for (i = 0, s_p = s. begin(); s_p != s. begin() + order; ++ i, ++ s_p)
				R. convert(*s_p, (*A_ilio) [i][i]);
			delete A_ilio;
		}
		// else if bisection possible
		else {
			PIR_ntl_ZZ_p R (m);
			DenseMatrix<PIR_ntl_ZZ_p>* A_ilio;
			MatrixMod::mod (A_ilio, A, R);
			IliopoulosElimination::smithIn (*A_ilio);
			int i; std::vector<integer>::iterator s_p;
			for (i = 0, s_p = s. begin(); s_p != s. begin() + order; ++ i, ++ s_p)
				R. convert(*s_p, (*A_ilio) [i][i]);
			delete A_ilio;
		}
		report << "Compuation of the k-rough part of the invariant factors finishes.\n";
	}

	/* Compute the Smith form of a dense matrix
	 * By adaptive algorithm.
	 * Compute the largest invariant factor,
	 * then based on that, compute the rough and smooth part, seperately.
	 */
	template <class Ring>
	void SmithForm (std::vector<integer>& s, const DenseMatrix<Ring>& A) {

		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, PROGRESS_REPORT);
		report << "Computation of the invariant factors starts (via an adaptive alo):" << std::endl;

		// compute the rank over a random prime field.
		int order = A. rowdim() < A. coldim() ? A. rowdim() : A. coldim();
		unsigned long r; 
		MatrixRank<Ring, Modular<int> > MR;
		r = MR. rank (A);
		report << "Matrix rank over a random prime field: " << r << '\n';

		typedef Modular<int> Field;
		typedef RationalSolver<Ring, Field, RandomPrime> Solver;
	    typedef LastInvariantFactor<Ring, Solver> LIF;
		typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
		OIF oif; oif. setThreshold  (10); oif.getLastInvariantFactor().setThreshold (6);
		typename Ring::Element lif, bonus;
		//oif. oneInvariantFactor_bonus (lif, bonus, A, (int)r);
		oif. oneInvariantFactor (bonus, A, (int)r);
		report << "The largest invariant factor: " << bonus << '\n';
		// bonus = smooth * rough;
		long* prime_p;
		std::vector<long> e(NPrime); std::vector<long>::iterator e_p;
		integer rbonus; A. field(). convert (rbonus, bonus);
		for (prime_p = prime, e_p = e. begin(); e_p != e. end(); ++ prime_p, ++ e_p) {
			*e_p = 0;
			while (rbonus % *prime_p == 0) {
				++ *e_p;
				rbonus = rbonus / *prime_p;
			}
		}
		std::vector<integer> smooth (order), rough (order);
		SmithFormSmooth (smooth, A, r, e);
		SmithFormRough (rough, A, rbonus);

		std::vector<integer>::iterator s_p, rough_p, smooth_p;

		for (rough_p = rough. begin(); rough_p != rough. end(); ++ rough_p) 
			if (* rough_p == 0) *rough_p = rbonus;

		for (s_p = s. begin(), smooth_p = smooth. begin (), rough_p = rough. begin(); s_p != s. begin() + order; ++s_p, ++ smooth_p, ++ rough_p) 
			*s_p = *smooth_p * *rough_p;
	}
}

