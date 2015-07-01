/*
 * Written by Zhendong Wan  <wan@mail.eecis.udel.edu> 
 */

#ifndef __LINBOX_RATIONAL_SOLVER_ADAPTIVE_H
#define __LINBOX_RATIONAL_SOLVER_ADAPTIVE_H
#include <linbox/field/modular-int32.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/blackbox/dense.h>

namespace LinBox {

	class RationalSolverAdaptive {
	public:

		template<class IRing, class OutVector, class InVector>
		static SolverReturnStatus solveNonsingular(OutVector& num, typename IRing::Element& den, 
										const DenseMatrix<IRing>& M, const InVector& b) {

			linbox_check ((M. rowdim() == M. coldim()) && (b.size() == M.rowdim()) && (num. size() ==M.coldim()));
			typedef Modular<int32> Field;
			RationalSolver<IRing, Field, RandomPrime, NumericalTraits> numerical_solver;
			SolverReturnStatus ret;
			ret = numerical_solver. solve(num, den, M, b);
			if (ret != SS_OK) {
				RationalSolver<IRing, Field, RandomPrime> solver;
				ret = numerical_solver. solve(num, den, M, b);
			}

			return ret;
		}
	};
}

#endif
