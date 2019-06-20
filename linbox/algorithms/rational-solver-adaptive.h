/* Copyright (C) 2007 LinBox
 *
 * Written by Zhendong Wan  <wan@mail.eecis.udel.edu>
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


#ifndef __LINBOX_rational_solver_adaptive_H
#define __LINBOX_rational_solver_adaptive_H

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"


#include "linbox/ring/modular.h"
#include "linbox/algorithms/rational-solver.h"
#include "linbox/randiter/random-prime.h"

namespace LinBox
{

	// Generic non-numerical solver requires conversion of the vector
	template<class IRing, class OutVector, class InVector>
	struct RationalSolverAdaptiveClass {
		static SolverReturnStatus solveNonsingular(OutVector& num, typename IRing::Element& den, const BlasMatrix<IRing>& M, const InVector& b)
		{
			linbox_check ((M. rowdim() == M. coldim()) && (b.size() == M.rowdim()) && (num. size() ==M.coldim()));
			typedef Givaro::Modular<int32_t> Field;
			// typedef Givaro::Modular<double> Field;
			DixonSolver<IRing, Field, PrimeIterator<IteratorCategories::HeuristicTag>, Method::SymbolicNumericNorm> numerical_solver;
			//DixonSolver<IRing, Field, PrimeIterator<IteratorCategories::HeuristicTag>, Method::SymbolicNumericOverlap> numerical_solver;
			SolverReturnStatus ret;
			ret = numerical_solver. solve(num, den, M, b);

			if (ret != SS_OK) {
				DixonSolver<IRing, Field, PrimeIterator<IteratorCategories::HeuristicTag>> solver;
				BlasVector<IRing> Ib(M.field()); Ib.reserve(b.size());
				typename IRing::Element tmp;
				for(typename InVector::const_iterator biter = b.begin();
				    biter != b.end();
				    ++biter) {
					Ib.push_back( M.field().init(tmp, *biter) );
				}
				ret = solver. solve(num, den, M, Ib);
			}

			return ret;
		}
	};


	// Specialization when the vector is already over the ring
	template<class IRing, class OutVector, template<typename T> class Container>
	struct RationalSolverAdaptiveClass<IRing, OutVector, Container<typename IRing::Element> > {
		static SolverReturnStatus solveNonsingular(OutVector& num, typename IRing::Element& den, const BlasMatrix<IRing>& M, const Container<typename IRing::Element> & b) {
			linbox_check ((M. rowdim() == M. coldim()) && (b.size() == M.rowdim()) && (num. size() ==M.coldim()));
			typedef Givaro::Modular<int32_t> Field;
			// typedef Givaro::Modular<double> Field;
			DixonSolver<IRing, Field, PrimeIterator<IteratorCategories::HeuristicTag>, Method::SymbolicNumericOverlap> numerical_solver;
			SolverReturnStatus ret;
			ret = numerical_solver. solve(num, den, M, b);

			if (ret != SS_OK) {
				DixonSolver<IRing, Field, PrimeIterator<IteratorCategories::HeuristicTag> > solver;
				ret = solver. solve(num, den, M, b);
			}

			return ret;
		}
	};


	class RationalSolverAdaptive {
	public:
		template<class IRing, class OutVector, class InVector>
		static SolverReturnStatus solveNonsingular(OutVector& num, typename IRing::Element& den, const BlasMatrix<IRing>& M, const InVector& b) {
			return RationalSolverAdaptiveClass<IRing,OutVector,InVector>::solveNonsingular(num, den, M, b);
		}
	};

}

#endif //__LINBOX_rational_solver_adaptive_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
