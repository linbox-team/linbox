/*
 * Written by Zhendong Wan  <wan@mail.eecis.udel.edu> 
 */

#ifndef __LINBOX_RATIONAL_SOLVER_ADAPTIVE_H
#define __LINBOX_RATIONAL_SOLVER_ADAPTIVE_H
#include <linbox/field/modular-int32.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/algorithms/rational-solver2.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/blackbox/dense.h>

namespace LinBox {

	class RationalSolverAdaptive {
	public:

		template<class IRing, class OutVector, class InVector>
		static SolverReturnStatus solveNonsingular(OutVector& num, typename IRing::Element& den, 
										const DenseMatrix<IRing>& M, const InVector& b) {

			linbox_check ((M. rowdim() == M. coldim()) && (b.size() == M.rowdim()) && (num. size() ==M.coldim()));
			std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
			report << "Rational solver start:\n";

			int n = M. rowdim();
			typedef typename IRing::Element IElement;
			typedef Modular<int32> Field;
			integer mnorm, bnorm; mnorm = 1; bnorm = 1;
			IRing R = M. field();
			typename InVector::const_iterator b_p; 
			typename OutVector::iterator num_p;
			IElement tmp_I; integer tmp;
			typename DenseMatrix<IRing>::ConstRawIterator raw_p;
			for (raw_p = M. rawBegin(); raw_p != M. rawEnd(); ++ raw_p) {
				R. convert (tmp, *raw_p);
				tmp = abs (tmp);
				if (tmp > mnorm) mnorm = tmp;
			}
			for (b_p = b. begin(); b_p != b.  end(); ++ b_p) {
				R. init (tmp_I, *b_p);
				R. convert (tmp, tmp_I);
				tmp = abs (tmp);
				if (tmp > bnorm) bnorm = tmp;
			}
				
			integer threshold; threshold = 1; threshold <<= 40;
			
			if ((mnorm < threshold) && (bnorm < threshold)) {

				report << "   Using numerical solver\n";
				double* DM = new double [n * n];
				double* Db = new double [n];
				double* DM_p, *Db_p;
				typename DenseMatrix<IRing>::ConstRawIterator raw_p;
				for (raw_p = M. rawBegin(), DM_p = DM; raw_p != M. rawEnd(); ++ raw_p, ++ DM_p) {
					R. convert (tmp, *raw_p);
					*DM_p = (double) tmp;
				}
				for (b_p = b. begin(), Db_p = Db; b_p != b. begin() + n; ++ b_p, ++ Db_p) {
					R. init (tmp_I, *b_p);
					R. convert (tmp, tmp_I);
					*Db_p = (double) tmp;
				}
				integer* numx = new integer[n];
				integer denx;
				int ret;
				ret = cblas_rsol (n, DM, numx, denx, Db);
				if (ret == 0) {
					report << "      Numerical solver success\n";
					R. init (den, denx);
					integer* numx_p;
					for (num_p = num. begin(), numx_p = numx; num_p != num. begin() + n; ++ num_p, ++ numx_p)
						R. init (*num_p, *numx_p);
				}
				else {
					report << "      Numerica method fails with error number " << ret <<"\n";
				}

				delete[] DM; delete[] Db; delete[] numx;
				if (ret == 0) {
					report << "Rational solver finished:\n";
					return SS_OK;
				}
			}
			report << "   Switch to Dixon lifting\n";

			RationalSolver<IRing, Field, RandomPrime> solver;
			SolverReturnStatus ret;
			ret = solver. solveNonsingular (num, den, M, b);
			report << "Rational solver finished:\n";
			return ret;
		}
	};
}

#endif
