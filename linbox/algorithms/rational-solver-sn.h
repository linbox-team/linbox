/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* sn-rational-solver.h */

/* Copyright (C) 2011 LinBox
 * Written Bryan Youse <>
 *
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_rational_solver_sn_H
#define __LINBOX_rational_solver_sn_H

#include <iostream>

#include "linbox/integer.h"
#include "linbox/field/param-fuzzy.h"
#include "linbox/solutions/methods.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/algorithms/dyadic-to-rational.h"
#include "linbox/blackbox/compose.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/vector-fraction.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/util/timer.h"
#include "linbox/field/PID-integer.h"

namespace LinBox {

	// bsd and mac problem
#undef _R

	/** \brief define the possible return status of the solver's computation.
	*/
	enum SNSolverReturnStatus {
		SNSS_OK, SNSS_FAILED, SNSS_SINGULAR, SNSS_INCONSISTENT
	};

	enum ShiftStatus {
		SHIFT_GROW, SHIFT_SHRINK, SHIFT_PEAK, SHIFT_SEARCH, SHIFT_MAX
	};

	/*
	 * A NumericSolver has
	 * init from a matrix A,
	 * solve(double* x, double* b)		// x = A^{-1}b
	 * apply(double* y, double * x);		// y = Ax
	 */

	template<class Ring, class NumericSolver>
	class RationalSolverSN {

	public:
		typedef typename Ring::Element Int;
		typedef std::vector<Int> IVector;
		// note: the type integer is also used.  For instance, we assume shift operator<< works on integer.
		typedef ParamFuzzy Field;
		typedef typename Field::Element Float;
		typedef std::vector<Float> FVector;
		typedef BlasBlackbox<Field> FMatrix;

	protected:
		Ring _R;
		VectorDomain<Ring> _VDR;
		Field _F;
		VectorDomain<Field> _VDF;
		NumericSolver _S;
		//inline static int check (int n, const double* M, integer* numx, integer& denx, double* b) ;
		//inline void update_r_xs (double* r, double* xs_int, double* xs_frac,
		//							int n, const double* M, double* x, int shift);
		//inline int rat_sol(IVector& numx, Int& denx, NumericSolver& _S, FVector& r, integer Bd);
		//inline void dyadicToRational(ZIVector& num, Int& den, vector<integer>& numx, integer& denx, integer Bd);
	private:
		size_t shift, shift_prev, shift_max, SHIFT_BOUND, HIT, MISS, iterations;
		ShiftStatus sstatus;
		bool searchPeak;
		double mnorm;
		bool exact_apply;
	public:

		RationalSolverSN(const Ring& R = Ring(), const NumericSolver& S = NumericSolver(),
				 bool ea=false) :
		       	_R(R), _VDR(R), _F(Field()), _VDF(Field()), _S(S), exact_apply(ea)
		{}

		/**
		 * IMatrix is matrix of integer type, eg. BlasBlackbox<PID-integer>
		 * IVector is linbox Vector of integer, eg. vector<PID-integer::Element>
		 * M is the matrix, b is rhs.
		 * num, den are the output  such that M*num = den*b (and den != 0 if successful).
		 */
		//  sparse matrix flag at the end, then avoid copying to DM as well ass
		//  new method to get hadamard bound and matrix norm!
		template <class IMatrix, class IVector>
		SNSolverReturnStatus solve(IVector& num, Int& den,
					   const IMatrix& M, const IVector& b)
		{
			Timer timer, solve_timer, rr_timer, tt;

			size_t n = b.size();
			// check basic feasiblility
			linbox_check((b.size() == M.rowdim()) && (num. size() == M.coldim()));

			// DM is M as matrix of doubles
			FMatrix DM(_F, n, n);
			//  Fix MatrixHom?
			//FMatrix* DMp = &DM;
			//MatrixHom::map<FMatrix, IMatrix, Field>(DMp, M, _F);

			if(n != M. rowdim() || n != M. coldim() || n != num.size()) {
				// std::cerr << "solve fail 1 - dimension mismatch" << std::endl;
				return SNSS_FAILED;
			}

			//  this is currently not used to check anything...
			integer entryBound = 1; entryBound <<= 49;  // nothing should exceed 2^50.
			SHIFT_BOUND = 52;

			//  why can't i put this in the for loop def???
			typename FMatrix::Iterator dm_p = DM.Begin();
			for (typename IMatrix::ConstIterator raw_p = M.Begin();
			     raw_p != M. End(); ++ raw_p, ++dm_p) {
				_F.init(*dm_p, *raw_p);
			}

			// build a numeric solver from new double matrix
			_S.init(DM);

			// r is b as vector of doubles.  (r is initial residual)
			FVector r(n);
			IVector bi(n);
			typename IVector::const_iterator b_p = b.begin();
			typename IVector::iterator bi_p = bi.begin();
			typename FVector::iterator r_p = r.begin();
			for (  ; b_p != b. begin() + n; ++b_p, ++r_p, ++bi_p) {
				*bi_p = *b_p;  //  copy original RHS
				_F.init(*r_p, *b_p);
			}

			//  denBound is the Hadamard bound, loopBound is roughly twice as much
			integer denBound, loopBound;
			zw_hbound (denBound, (int)n, (int)n, &*(DM.Begin()));
			loopBound = denBound*denBound;

			mnorm = zw_dOOnorm(&*(DM.Begin()), (int)n, (int)n);  //  infinity-norm of matrix
			//  set max shift to avoid exact applys
			size_t bits = 0;
			size_t mn2 = nextPower2((size_t)mnorm);
			for(;mn2;mn2>>=1, bits++);

			SHIFT_BOUND -= bits;
			//std::cerr << "BITS" << bits << "MAX" << SHIFT_BOUND << std::endl;

			loopBound *= (2*mnorm + zw_dmax((int)n, &*(r.begin()), 1));

			std::vector<integer> numx(n), tnum(n); // numerator of binary expansion
			integer denx = 1, tden; // denominator of binary expansion (denx is a power of 2).

			FVector x(n), xs_int(n), xs_frac(n);
			FVector lastr(n);
			IVector lastb(n);

			//set initial shift small.
			shift = 2;
			shift_prev = shift;
			shift_max = 0;
			searchPeak = false;
			sstatus = SHIFT_GROW;
			HIT = 0; MISS = 0;
			iterations = 0;
			integer ay, be;
			PID_integer Z;
			int ret;

			bool recon_success = false;
			int recon_status = 0;

			//timer.clear(); timer.start();
#ifdef SN_EARLY_TERM
			integer bound = denBound;
			//double it_cost = 0, rr_cost = 0;
#else
			integer bound = loopBound;
#endif
			//size_t rr_count = 0;
			//solve_timer.clear(); rr_timer.clear();
			do{
				//tt.clear(); tt.start();
				ret = rat_sol(numx, denx, xs_int, xs_frac, bi, lastb, r, lastr, x, bound, M);
				//tt.stop(); solve_timer += tt;

				if(ret == 1){
					// std::cerr << "numsym loop failed - likely lack of num accuracy" << std::endl;
					return SNSS_FAILED;
				}
				else if(ret == 2) denBound = denx; // zero residual

				// we're trying to early-term
				//std::cerr << bound << " " << loopBound << std::endl;
				if(bound < loopBound){
					//  update bound for next iteration (if applicable)
#if 0
				       	it_cost = solve_timer.realtime()/(double)iterations;
					   rr_cost = rr_timer.realtime()/(double)rr_count;
					   std::cerr << "iteration cost: " << it_cost << " v. rr cost: " << rr_cost << std::endl;
#endif
					Z.sqrt(bound, loopBound*bound);
					bound <<= 2;
					int rPos = rand()%(int)n;
					//std::cerr << "At iteration " << iterations << ", ";
					if(dyadicToRational(Z, ay, be, numx[rPos], denx, denBound) /*== 2*/){
						//std::cerr << "Random single worked!  ";
					}
					else{
						//std::cerr << "Random single failed." << std::endl;
						continue;
					}
				}

				//tt.clear(); tt.start();
				recon_status = dyadicToRational(Z, num, den, numx, denx, denBound);
				//tt.stop(); rr_timer += tt; ++rr_count;
				//std::cerr << "RRT: " << rr_timer << std::endl;

				recon_success = recon_status > 0;
				//if(!recon_success) std::cerr << "Full failed!" << std::endl;
				//else std::cerr << "Full worked!" << std::endl;
			} while((bound < loopBound) && !recon_success);

			//timer.stop(); std::cerr << "rat_sol time: " << solve_timer.realtime() << " rr time: " << rr_timer.realtime() << " Total: " << timer << std::endl;

#if 0
			writeVec(numx, "numx", 0, 10);
			std::cerr << denx << std::endl;
			writeVec(num, "num");
			std::cerr << "den: (large)" << std::endl;// << den << endl;
#endif

			if (recon_success) {
#if 0
				if(recon_status == 2) std::cerr << "reconstruction guaranteed" << std::endl;
				else std::cerr << "reconstruction speculative" << std::endl;

				std::cerr << "Solve success. Iterations: " << iterations << std::endl;
				std::cerr << HIT << " hits, " << MISS << " misses. (";
				fprintf(stderr,  "%.2f", (float)(HIT)/(float)(HIT+MISS)*100.0);
				std::cerr << "%) Maximum shift: " << shift_max << std::endl;
#endif
			}
			else{
				// std::cerr << "rat reconstruction asserts failure" << std::endl;
				// dumpData(M, b, numx, denx, denBound);
				return SNSS_FAILED;
			}

			if (_R.isZero(den)) {
				// std::cerr << "fail: zero denominator after rat-recons" << std::endl;
				return SNSS_FAILED;
			}

#if 0
			//  Answer checking
			IVector y(n), z(n);
			M.apply(y, num);
			_VDR.mul(z, b, den);
			if ( !_VDR.areEqual(y, z)) {
				std::cerr << "fail check: A*x != b exactly" << std::endl;
				dumpData(M, b, numx, denx, denBound);
				return SNSS_FAILED;
			}
#endif
			return SNSS_OK;

		} // solve

#include "rational-solver-sn.inl"

#if 0
		//embedded definitions now, so no declarations
		// functions used by solve()
		//protected:
		//print out a vector
		template <class Elt>
		inline static int printvec (const Elt* v, int n);
		/** Compute the OO-norm of a mtrix */
		inline static double zw_dOOnorm(const double* M, int m, int n);
		/** compute the maximam of absolute value of an array*/
		inline static double zw_dmax (const int N, const double* a, const int inc);
		/* apply  y <- Ax */
		inline static int zw_dapply (int m, int n, const double* A, const double* x, double* y);
		inline static int zw_mpzapply (int m, int n, const double* A, const integer* x, integer* y);
		//update the numerator; num = num * 2^shift + d;
		inline static int update_num (integer* num, int n, const double* d, int shift);
		//update r = r * shift - M d, where norm (r) < 2^32;
		inline static int update_r_int (double* r, int n, const double* M, const double* d, int shift);
		//update r = r * shift - M d, where 2^32 <= norm (r) < 2^53
		inline static int update_r_ll (double* r, int n, const double* M, const double* d, int shift);
		// compute  the hadamard bound
		inline static int zw_hbound (integer& b, int m, int n, const double* M);
		// compute the inverse of a general matrix
		inline static int zw_dgeinv(double* M, int n);
		/* solve Ax = b
		 * A, the integer matrix
		 * b, integer rhs
		 * Return value
		 * 0, ok.
		 * 1, the matrix is not invertible in floating point operations.
		 * 2, the matrix is not well conditioned.
		 * 3, incorrect answer, possible ill-conditioned.
		 */
		//inline int rsol (Ring& R, int n, const double* M, integer* numx, integer& denx, double* b);
#endif

	}; // class RationalSolverSN

} // namespace LinBox

#endif // __LINBOX_rational_solver_sn_H
